/**
 * test_uniform.c
 *
 * ORATIO Library
 * Copyright (C) 1999-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "test_ran.h"

#include <oratio/stat.h>

#define SIZE 20000

/**
 * @param f is a write-open file.
 * @param r is an initialized random numbers generator.
 */
void uniform_test(FILE *f, Random *r, int n)
{
  int *obs;
  int df, nclasses;

  double v;
  double chi, chi_th;

  int i;

  nclasses = 4096;
  obs = _mallocs(nclasses, int);

  for (i = 0; i < nclasses; i++) obs[i] = 0;

  for (i = 0; i < SIZE; i++) {
    v = (ran_random_get_val(r))*nclasses;
    df = v;
    obs[df]++;
  }

  chi = 0;
  v = ((double)SIZE)/nclasses;
  for (i = 0; i < nclasses; i++) chi+= ((obs[i]-v)*(obs[i]-v))/v;

  chi_th=1-st_chi2_cdf(chi, nclasses-1);
  /*  chi_th=1-st_chi2_cdf(4211, nclasses-1); */

  /* Graphical representation on the unit square. */
  if (f) for (i = 0; i < n*n; i++)
      fprintf(f, "%g %g\n", ran_random_get_val(r), ran_random_get_val(r));
}

/**
 * Performs some checkings with the streams-version of generator MRG32K3a.
 *
 * @param f is a write-open file.
 * @param r is an initialized MRG32K3a random numbers generator.
 */
void mrg32k3a_check(FILE *f, Random *r)
{
  printf("Tests on the MRG32K3a generator.\n");
  uniform_test(f, r, 64);
}

/**
 * Perform some tests on uniform U(0,1) distributions.
 *
 * @param f is a write-open file.
 * @param r is an initialized random numbers generator.
 */
void uniforms_check(FILE *f, Random *r)
{
  unsigned long state;
  int id;

  for (id = 1; id < END_RANDOM_GENERATORS; id++) {
    state = 0;
    ran_random_set_generator(r, id, &state, 1);

    printf("%d. Tests on the generator %s.\n",
	   id, ran_random_generator_name(r));

    f = fopen(ran_random_generator_name(r), "w");
    uniform_test(f, r, 64);
    fclose(f);
  }

  ran_random_set_generator(r, RAN_MRG32K3A, &state, 1);
  mrg32k3a_check(NULL, r);
  ran_random_resetnextsubstream(r);
  mrg32k3a_check(NULL, r);
  ran_random_resetstartsubstream(r);
  mrg32k3a_check(NULL, r);
  ran_random_resetstartstream(r);
  mrg32k3a_check(NULL, r);
}
