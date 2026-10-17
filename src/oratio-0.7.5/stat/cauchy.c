/**
 * cauchy.c
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2008-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

/**
 * Cauchy distribution
 */

#include "oratio/ntools.h"

#include <math.h>

/**
 * Compute the value of a Cauchy probabibility density function at x.
 *
 * @param x is the point where the Cauchy pdf has to be evaluated.
 * @param a is the location parameter (the median).
 * @param b is the scale parameter.
 *
 * @return the value of the Cauchy pdf at x.
 */
double st_cauchy_pdf(double x, double a, double b)
{
  double q = (x-a)/b;
  return 1/(b*M_PI*(1+q*q));
}

/**
 * Compute the value of a Cauchy cumulative distribution function at x.
 *
 * @param x is the point where the Cauchy cdf has to be evaluated.
 * @param a is the location parameter (the median).
 * @param b is the scale parameter.
 *
 * @return the value of the Cauchy cdf at x.
 */
double st_cauchy_cdf(double x, double a, double b)
{
  return 0.5 + M_1_PI*atan((x-a)/b);
}

/**
 * Inverse cumulative distribution of a Cauchy random variable.
 *
 * @param u is the value corresponding to the desired quantile.
 * @param a is the location parameter (the median).
 * @param b is the scale parameter.
 *
 * @return the value of the Cauchy icdf at u.
 */
double st_cauchy_icdf(double u, double a, double b)
{
  /* Produce infinite values if u is set to 0.0 or 1.0 */
  if (u == 1.0) return nt_infinity();
  else if (u == 0.0) return -nt_infinity();

  return a+b*tan(M_PI*(u-0.5));
}

/**
 * Compute the median of a Cauchy distribution.
 *
 * @param a is the location parameter (the median).
 * @param b is the scale parameter.
 *
 * @return the median of the Cauchy distribution.
 */
double st_cauchy_median(double a, double b)
{
  return a;
}

/**
 * Compute the mode of a Cauchy distribution.
 *
 * @param a is the location parameter (the median).
 * @param b is the scale parameter.
 *
 * @return the mode of the Cauchy distribution.
 */
double st_cauchy_mode(double a, double b)
{
  return a;
}
