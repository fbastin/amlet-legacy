/**
 * chi2.c
 *
 * STAT Library
 * Copyright (C) 1999-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/st_utils.h"

#include <math.h>

/**
 * Return the value of a chi-square probabibility density function,
 * at x.
 *
 * @param x is the point where the probability density is evaluated.
 * @param k is the number of degrees of freedom.
 *
 * @return the value of the chi 2 pdf at x.
 */
double st_chi2_pdf(double x, double k)
{
  double d = k/2.0;

  return (pow(x, d-1.0)*exp(-x/2))/(pow(2, d)*exp(lgamma(d)));
}

/**
 * Compute the chi-square cumulative distribution function.
 *
 * @param x is the point where the cumulative density is evaluated.
 * @param k is the number of degrees of freedom.
 */
double st_chi2_cdf(double x, double k)
{
  return st_gamma_p(k/2, x/2);
}

/**
 * Return the mean of a chi-square distribution with k degrees of freedom.
 *
 * @param k is the number of degrees of freedom.
 */
double st_chi2_mean(double k)
{
  return k;
}

/**
 * Return the variance of a chi-square distribution with k degrees
 * of freedom.
 *
 * @param k is the number of degrees of freedom.
 */
double st_chi2_variance(double k)
{
  return 2*k;
}

/**
 * Return the skewness of a chi-square distribution with k degrees of
 * freedom.
 *
 * @param k is the number of degrees of freedom.
 *
 * @return the skewness of the distribution.
 */
double st_chi2_skewness(double k)
{
  return sqrt(8.0/k);
}
