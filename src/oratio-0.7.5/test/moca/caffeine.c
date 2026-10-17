/**
 * caffeine.c
 *
 * ORATIO Library
 * Copyright (C) 2005-2007 <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <stdio.h>

#include <oratio/ran.h>
#include <oratio/moca.h>

#include "caffeine.h"

#define MAX_SIZE 100000

/**
 * Faure tests.
 *
 * @param ran is an initialized random numbers generator.
 */
int moca_faure_check(double **sample, Random *ran)
{
  FILE *f;

  double quasi[MAX_SIZE];
  int dimen, atmost;

  int i;

  MOCAFaure *mcf;

  f = fopen("FaureDraws", "w");

  mcf = moca_faure_new(0);

  dimen = 3;
  atmost = 100;
  if (!moca_faure_init(mcf, dimen, atmost)) ;

  for (i = 0; i < 10; i++) {
    moca_faure_next(mcf, quasi);
    fprintf(f, "%lg %lg %lg\n", quasi[0], quasi[1], quasi[2]);
  }

  fclose(f);

  moca_faure_free(0);

  return 0;
}

/**
 * Halton tests
 *
 * @param ran is an initialized random numbers generator.
 */
int moca_halton_check(double **sample, Random *ran)
{
  FILE *f;
  MOCA_Halton *h;
  double x[3];

  int i, j;

  f = fopen("HaltonDraws", "w");

  for (i = 1; i <= 100; i++) fprintf(f, "%lg\n", moca_halton(2, 1, i, 0, 0));

  h = moca_halton_new(3, MOCA_HALTON_BW_SCRAMBLING);

  for (i = 0; i < 10; i++) {
    moca_halton_next(h, x);
    for (j = 0; j < 3; j++) printf("%lg ", x[j]);
    printf("\n");
  }

  moca_halton_free(h);

  fclose(f);

  return 0;
}

int main()
{
  int err = 0;

  Random *r;

  double **sample;

  sample = nt_matrix_new(10,1000);

  r = ran_random_new(0, NULL, 0);

  /* QMC tests */

  /*  err = moca_faure_check(sample, r);

      err = moca_halton_check(sample, r); */

  err = moca_sobol_check(r);

  //  err = moca_lattice_check(r);

  /* RQMC test */

  //  err = moca_test_rqmc(r);

  /* Distribution tests */

  //  moca_test_distributions(sample);

  ran_random_free(r);

  nt_matrix_free(sample);

  return err;
}

