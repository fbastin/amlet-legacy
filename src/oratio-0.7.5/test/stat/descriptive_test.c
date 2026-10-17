/**
 * descriptive_test.c
 *
 * Stat Library
 * Copyright (C) 2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/oratio.h"

#include <stdio.h>

/**
 * Tests on the descriptive statistics functions.
 */
int st_descriptive_test()
{
  double q[3];
  double obs[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  int n;

  printf("Test on the descriptive statistics functions\n");
  printf("--------------------------------------------\n");

  n = sizeof(obs)/sizeof(obs[0]);

  printf("%d observations: ", n);
  nt_vect_print(stdout, "obs", obs, n);

  printf("Median: %g\n", st_median(obs, 1, 0, n));
  printf("Median restricted to observed values: %g\n",
	 st_median_observation(obs, 1, 0, n));

  st_quantiles(n, obs, 1, 3, q);
  nt_vect_print(stdout, "Quantiles", q, 3);
  printf("Quantile 0.5: %g\n", st_quantile(n, obs, 1, 0.5));

  return 0;
}
