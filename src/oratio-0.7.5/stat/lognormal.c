/**
 * lognormal.c
 *
 * ORATIO Library - statistical module
 * Copyright (C) 1999-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

/**
 * Lognormal distribution
 */

#include "oratio/ntools.h"
#include "oratio/stat.h"

#include <math.h>
#include <float.h>

#define _min(a, b) (((a) < (b)) ? a : b)
#define _max(a, b) (((a) < (b)) ? b : a)

/**
 * Evaluate the probability density function of the lognormal variate
 * L of parameters mu and sigma, the mean and standard deviation of log L.
 *
 * @param x is the point where the probability density function is
 *   evaluated.
 *
 * @param mu is the mean of log L.
 * @param sigma is the standard deviation of log L.
 *
 * @return the PDF value of L at x.
 */
double st_lognormal_pdf(double x, double mu, double sigma)
{
  double temp = log(x)-mu;

  if (x > 0.0)
    return 1/(x*sigma*sqrt(2*M_PI))*exp(-(temp*temp)/(2*sigma*sigma));
  else
    return 0.0;
}

/**
 * Evaluate the cumulative density function of the lognormal L of
 * parameters mu and sigma, the mean and standard deviation of log L.
 *
 * @param x is the point where the cumulative distribution function is
 *   evaluated. x must be positive.
 * @param mu is the mean of log L.
 * @param sigma is the standard deviation of log L.
 * @param ccum contains on output the complement of result, i.e. ccum
 *   = 1.0 - result.
 *
 * @return the value of the CDF of L at x.
 *
 * Reference:
 *  Evans, Merran, Hastings, Nicholas and Peacock, Brian, "Statistical
 *  Distributions, Second Edition", Wiley, 1993, pp 102-105
 */
double st_lognormal_cdf(double x, double mu, double sigma)
{
  /* Return NaN if x < 0 of sigma <= 0.0 */
  if (x < 0.0 || sigma <= 0.0) return nt_nan();

  if (x == 0.0)
    return 0.0;
  else
    /* we ensure that the returned number does not exceed 1.0 (for
       instance due to rounding errors). */
    return _min(0.5*(1+erf((log(x)-mu)/(sigma*M_SQRT2))),1.0);
}

double st_lognormal_icdf(double q, double mu, double sigma)
{
  return exp(mu + sigma*st_gaussian_icdf(q));
}

double st_lognormal_mean(double mu, double sigma)
{
  /**
   * Ref.: http://en.wikipedia.org/wiki/Log_normal_distribution.
   */

  return exp(mu+0.5*sigma*sigma);
}

double st_lognormal_stdv(double mu, double sigma)
{
  double sigma2 = sigma*sigma;
  double mean = st_lognormal_mean(mu, sigma);

  return sqrt(exp(sigma2)-1)*mean;
}

double st_lognormal_variance(double mu, double sigma)
{
  double sigma2 = sigma*sigma;
  double mean = st_lognormal_mean(mu, sigma);

  return (exp(sigma2)-1)*mean*mean;
}
