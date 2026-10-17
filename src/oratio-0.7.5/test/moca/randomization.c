/**
 * randomization.c
 *
 * Monte-Carlo and Quasi-Monte Carlo Library
 * Copyright (C) 2008 <fabian.bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "caffeine.h"

#include <oratio/oratio.h>

/**
 * Some basic tests on randomized quasi-Monte Carlo sequences.
 *
 * @param r is an initialized random numbers generator.
 *
 * @return 0 if everything has gone well, some error code otherwise.
 */
int moca_test_rqmc(Random *r)
{
  MOCA_Sobol *sbl;

  int i, j, nrep;
  double *x, *shift;
  int s, n;

  nrep = 10;

  s = 3;
  n = 10;

  sbl = moca_sobol_new(0, s);
  x = _mallocs(s, double);
  shift = _mallocs(s, double);

  moca_sobol_init(sbl);
  for (i = 0; i < nrep; i++) {
    moca_shift_init(r, s, shift);
    for (j = 0; j < n; j++) {
      moca_sobol_next(sbl, x);
      moca_shift_apply(s, x, shift);
      ran_write_vector(stdout, x, s);
    }
    moca_sobol_reset(sbl);
  }

  free(x);
  free(shift);

  moca_sobol_free(sbl);

  return 0;
}
