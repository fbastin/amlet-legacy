/**
 * lognormal.c
 *
 * Stat Library
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <oratio/stat.h>

#include <stdio.h>

int st_lognormal_test()
{
  double mu, sigma;

  printf("Test on the lognormal distribution\n");
  printf("-----------------------------------\n");

  mu = 0.0;
  sigma = 0.5;
  printf("X lognormally distributed, with parameters %g and %g\n", mu, sigma);
  printf("Mean: %g, Standard Deviation: %g, Variance: %g\n",
	 st_lognormal_mean(mu, sigma), st_lognormal_stdv(mu, sigma),
	 st_lognormal_variance(mu, sigma));

  return 0;
}
