/**
 * logistic.c
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <math.h>

/**
 * Compute the cumulative distribution function of a Rayleigh variate.
 *
 * @param x is the point where cdf has to be evaluated.
 * @param beta is the Rayleigh parameter.
 *
 * @return the cdf of a Rayleigh variate.
 */
double st_rayleigh_cdf(int x, double beta)
{
  return 1-exp(-x*x/beta);
}

/**
 * Compute the inverse cumulative distribution function of a Rayleigh
 * variate.
 *
 * @param u is the point where inverse cdf has to be evaluated.
 * @param beta is the Rayleigh parameter.
 *
 * @return the inverse cdf of a Rayleigh variate.
 */
double st_rayleigh_icdf(double u, double beta)
{
  return sqrt(-beta*log(1-u));
}
