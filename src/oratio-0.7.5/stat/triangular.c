/**
 * triangular.c
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <math.h>

/**
 * Compute the density of a Triangular distribution at x.
 *
 * @param x is the point where the pdf has to be evaluated.
 * @param a is the lower limit.
 * @param b is the upper limit.
 * @param c is the mode of the distribution.
 *
 * @return the density of the Triangular distribution evaluated at x.
 */
double st_triangular_pdf(double x, double a, double b, double c)
{
  if (x <= a || x >= b) return 0.0;
  else if (x <= c) return 2*(x-a)/((b-a)*(c-a));
  else return 2*(b-x)/((b-a)*(b-c));
}

/**
 * Compute the cumulative distribution function of a Triangular
 * distribution at x.
 *
 * @param x is the point where the cdf has to be evaluated.
 * @param a is the lower limit.
 * @param b is the upper limit.
 * @param c is the mode of the distribution.
 *
 * @return the cumulative distribution function of the Triangular
 *   distribution evaluated at x.
 */
double st_triangular_cdf(double x, double a, double b, double c)
{
  if (x <= a) return 0.0;
  else if (x >= b) return 1.0;
  else if (x <= c) return ((x-a)*(x-a))/((b-a)*(c-a));
  else return 1.0 - ((b-x)*(b-x))/((b-a)*(b-c));
}

/**
 * Compute the inverse cumulative distribution function of a Triangular
 * distribution.
 *
 * @param p is the value (0 <= p <= 1) with respect we want to compute
 *   the associated quantile.
 * @param a is the lower limit.
 * @param b is the upper limit.
 * @param c is the mode of the distribution.
 *
 * @return the inverse cumulative distribution function of the Triangular
 *   distribution evaluated at p.
 */
double st_triangular_icdf(double u, double a, double b, double c)
{
  if (u <= (c-a)/(b-a))
    return a+sqrt((b-a)*(c-a)*u);
  else
    return b-sqrt((b-a)*(b-c)*(1.0-u));
}

/**
 * Compute the mode of a Triangular distribution.
 *
 * @param a is the lower limit.
 * @param b is the upper limit.
 * @param c is the mode of the distribution.
 *
 * @return the mode of the Triangular distribution.
 */
double st_triangular_mean(double a, double b, double c)
{
  return (a+b+c)/3.0;
}

/**
 * Compute the variance of a Triangular distribution.
 *
 * @param a is the lower limit.
 * @param b is the upper limit.
 * @param c is the mode of the distribution.
 *
 * @return the variance of the Triangular distribution.
 */
double st_triangular_variance(double a, double b, double c)
{
  return (a*(a-b-c)+b*(b-c)+c*c)/18.0;
}

/**
 * Compute the standard deviation of a Triangular distribution.
 *
 * @param a is the lower limit.
 * @param b is the upper limit.
 * @param c is the mode of the distribution.
 *
 * @return the standard deviation of the Triangular distribution.
 */
double st_triangular_stdv(double a, double b, double c)
{
  return sqrt(st_triangular_variance(a, b, c));
}

/**
 * Compute the mode of a Triangular distribution.
 *
 * @param a is the lower limit.
 * @param b is the upper limit.
 * @param c is the mode of the distribution.
 *
 * @return the mode of the Triangular distribution.
 */
double st_triangular_mode(double a, double b, double c)
{
  return c;
}

/**
 * Compute the median of a Triangular distribution.
 *
 * @param a is the lower limit.
 * @param b is the upper limit.
 * @param c is the mode of the distribution.
 *
 * @return the median of the Triangular distribution.
 */
double st_triangular_median(double a, double b, double c)
{
  if (c > (b-a)/2.0)
    return a+sqrt(((b-a)*(c-a))/2.0);
  else
    return b-sqrt(((b-a)*(b-c))/2.0);
}
