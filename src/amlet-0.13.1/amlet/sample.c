/**
 * sample.c
 *
 * OPHELIA Library
 * Copyright (C) 2009-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/oratio.h>

#include "amlet_driver.h"
#include "lattice.h"

#include <string.h>

/**
 * Allocate memory for (quasi-)random draws.
 *
 * @param pop is the population size.
 * @param n is the number of random vectors to draw per individual.
 * @param s is the size of each vector.
 *
 * @return an allocated tabular of the required size to store the
 *   draws, with each line corresponding to one individual.
 */
double **amlet_sample_alloc(int pop, int n, int s)
{
  int i;
  double **S;

  S = _mallocs(pop, double *);
  for (i = 0; i < pop; i++) S[i] = _mallocs(n*s, double);

  return S;
}

/**
 * Deallocate memory associated to (quasi-)random draws.
 *
 * @param S is the (allocated) array of draws.
 * @param pop is the population size.
 */
void amlet_sample_free(double **S, int pop)
{
  int i;

  for (i = 0; i < pop; i++) free(S[i]);
  free(S);
}

/**
 * Sample generation for use in (quasi-)Monte Carlo simulation
 *
 * @param f is the file where the U(0,1) draws will be saved, if not NULL.
 * @param ran is an initialized random numbers generator.
 * @param method is the sampling method.
 * @param n is the number of vectors to draw.
 * @param s is the size of draws vectors.
 * @param sample contains, on output, the generated sample.
 */
void amlet_sample_generate_01(FILE *f, Random *ran, int np, int method,
			      int n, int s, double **sample)
{
  int i, j;

  switch (method) {
  case MOCA_HALTON: {
    MOCA_Halton *h;
    double shift[s];

    h = moca_halton_new(s, MOCA_HALTON_NO_SCRAMBLING);
    moca_halton_sample(h, n, sample[0], s);

    for (i = 1; i < np; i++) memcpy(sample[i], sample[0], n*s*sizeof(double));
    for (i = 0; i < np; i++) {
      /* Apply the random shift. */
      moca_shift_init32(ran, s, shift);
      for (j = 0; j < n; j++)
	moca_shift_apply(s, &(sample[i][j*s]), shift);
    }

    moca_halton_free(h);

    break;
  }
  case MOCA_RANK1_LATTICE: {
    MOCA_Rank1Lattice *lattice;
    int a[s];
    double shift[s];
    int inc, d;
    int *l;

    l = amlet_get_default_lattice();

    inc = amlet_get_default_lattice_max_dimension()+1;
    d = ds_find_inter(l, amlet_get_default_lattice_number(),
		      inc, sizeof(int), &n, ds_search_int);
    if (s <= 30)
      for (i = 0; i < s; i++) a[i] = l[d*inc + (i+1)];
    else {
      Random *r = ran_random();
      /* Switch to a Monte-Carlo sample. */
      amlet_sample_generate_01(NULL, ran, np, MOCA_MONTE_CARLO, n, s, sample);
      ran_random_free(r);
      break;
    }

    for (i = 0; i < s; i++) printf("%d ", a[i]);
    printf("\n");

    lattice = moca_rank1_lattice_new(n, s, a);

    moca_rank1_lattice_sample(lattice, n, sample[0], s);
    for (i = 1; i < np; i++) memcpy(sample[i], sample[0], n*s*sizeof(double));
    for (i = 0; i < np; i++) {
      /* Apply the random shift. */
      moca_shift_init32(ran, s, shift);
      for (j = 0; j < n; j++) {
	moca_shift_apply(s, &(sample[i][j*s]), shift);
	/* Apply Baker transformation. */
	moca_baker_point(s, &(sample[i][j*s]), 1);
      }
    }
    moca_rank1_lattice_free(lattice);
    break;
  }
  case MOCA_SOBOL: {
    MOCA_Sobol *sbl = moca_sobol_new(0, s);

    for (i = 0; i < np; i++) moca_sobol_sample(sbl, n, s, sample[i], s);

    moca_sobol_free(sbl);
    break;
  }
  default:
    /* Generate a Monte Carlo sample */
    for (i = 0; i < np; i++)
      for (j = 0; j < n*s; j++)
	sample[i][j] = ran_random_get_val(ran);
    break;
  }

  if (f)
    for (i = 0; i < np; i++) ran_export_sample(f, n, s, sample[i], s);
}

/**
 * @param f is the file where the U(0,1) draws will be saved, if not NULL.
 * @param ran is an initialized random numbers generator.
 * @param n_random_coeff is the number of random coefficients in the model.
 */
void amlet_sample_generate(FILE *f, Random *ran, AMLET_Driver *dr, int new,
			   double **sample, int *types, int n_random_coeff)
{
  int np = dr->var[NP].i;
  int i, j, k;

  /* Generate the draws, QMCREP times in order to allow error
     estimation with Quasi-Monte Carlo draws. */
  if (new) {
    for (i = 0; i < dr->var[QMCREP].i; i++)
      amlet_sample_generate_01(f, ran, np, dr->var[SAMPLING].i,
			       dr->var[NREP].i, n_random_coeff,
			       sample+i*np);
  }
  /* Convert [0-1] draws to the required types. */
  for (i = 0; i < np*dr->var[QMCREP].i; i++)
    for (k = 0; k < dr->var[NREP].i; k++)
      for (j = 0; j < n_random_coeff; j++)
	sample[i][k*n_random_coeff+j] =
	  moca_convert_random_value(types[j], sample[i][k*n_random_coeff+j]);
}

