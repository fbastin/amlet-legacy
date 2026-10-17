/**
 * weibull.c
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <math.h>

/**
 * Compute the cumulative distribution function of a Weibull variate.
 *
 * @param x is the point where cdf has to be evaluated.
 * @param alpha is the shape parameter; alpha > 0.
 * @param beta is the scale parameter; beta > 0.
 *
 * @return the cdf of a Weibull of parameters alpha and beta.
 */
double st_weibull_cdf(double x, double alpha, double beta)
{
  return 1-exp(-pow(x/beta, alpha));
}

/**
 * Compute the inverse cumulative distribution function of a Weibull variate.
 *
 * @param u is the point where inverse cdf has to be evaluated.
 * @param alpha is the shape parameter; alpha > 0.
 * @param beta is the scale parameter; beta > 0.
 *
 * @return the inverse cdf of a Weibull of parameters alpha and beta.
 */
double st_weibull_icdf(double u, double alpha, double beta)
{
  return beta*pow(log(1/(1-u)), 1.0/alpha);
}

