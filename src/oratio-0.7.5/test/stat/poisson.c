/**
 * poisson.c
 *
 * Stat Library
 * Copyright (C) 2006-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/stat.h"

#include <stdio.h>

/**
 * Tests on the Poisson distribution.
 */
int st_poisson_test()
{
  printf("Tests on the Poisson distribution\n");
  printf("---------------------------------\n");

  /* PMF */
  printf("Some values of Poisson PMF.\n");
  printf("x=%d, lambda=%g: %g\n", 0, 2.5, st_poisson_pmf(0, 2.5));
  printf("x=%d, lambda=%g: %g\n", 1, 2.5, st_poisson_pmf(1, 2.5));
  printf("x=%d, lambda=%g: %g\n", 2, 2.5, st_poisson_pmf(2, 2.5));
  printf("x=%d, lambda=%g: %g\n", 3, 0.5, st_poisson_pmf(3, 0.5));

  /* CDF */
  printf("CDF\n");
  printf("x=%d, lambda=%g: %g\n", 0, 2.5, st_poisson_cdf(0, 2.5));
  printf("x=%d, lambda=%g: %g\n", 1, 2.5, st_poisson_cdf(1, 2.5));

  return 0;
}
