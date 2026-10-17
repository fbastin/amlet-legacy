/**
 * test_sobol.c
 *
 * ORATIO Library
 * Copyright (C) 2005-2007 <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "caffeine.h"

/**
 * Sobol tests
 *
 * @param ran is an initialized random numbers generator.
 *
 * @return 0 if the test has been successful, an error code otherwise.
 */
int moca_sobol_check(Random *ran)
{
  FILE *f;

  int n;

  double *quasi;
  int dimen = 2;

  int i;

  MOCA_Sobol *sbl;

  n = 100000;
  quasi = _mallocs(dimen, double);

  sbl = moca_sobol_new(0, dimen);

  f = fopen("SobolDraws", "w");

  moca_sobol_init(sbl);
  for (i = 0; i < n; i++) {
    moca_sobol_next(sbl, quasi);
    ran_write_vector(f, quasi, dimen);
  }

  if (f) fclose(f);

  free(quasi);
  moca_sobol_free(sbl);

  return 0;
}
