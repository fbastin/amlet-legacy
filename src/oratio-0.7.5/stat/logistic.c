/**
 * logistic.c
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <math.h>

/**
 * Compute the cumulative distribution function of a logistic variate.
 *
 * @param x is the point where cdf has to be evaluated.
 * @param mu is the location parameter, also the mean.
 * @param s is the scale parameter.
 *
 * @return the cdf of a logistic variate.
 */
double st_logistic_cdf(double x, double mu, double s)
{
  return 1/(1+exp(-(x-mu)/s));
}

/**
 * Compute the inverse cumulative distribution function of a logistic
 * variate.
 *
 * @param u is the point where inverse cdf has to be evaluated.
 * @param mu is the location parameter, also the mean.
 * @param s is the scale parameter.
 *
 * @return the inverse cdf of a logistic variate.
 */
double st_logistic_icdf(double u, double mu, double s)
{
  return mu - s*log(1/u-1);
}

