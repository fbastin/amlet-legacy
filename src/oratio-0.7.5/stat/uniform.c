/**
 * uniform.c
 *
 * STAT Library
 * Copyright (C) 1999-2006 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Compute the probability distribution function of a variate X
 * uniformly distributed on ]a,b[.
 *
 * @param x is a real where the function has to be evaluated.
 * @param a is the interval lower bound.
 * @param b is the interval upper bound.
 *
 * @return the probability distribution function at x.
 */ 
double st_uniform_pdf(double x, double a, double b)
{
  if (x > a && x < b)
    return 1/(b-a);
  else
    return 0.0;
}

/**
 * Mean of a variate X uniformly distributed on ]a,b[.
 *
 * @param a is the interval lower bound.
 * @param b is the interval upper bound.
 *
 * @return the mean of an uniform distribution on ]a,b[.
 */
double st_uniform_mean(double a, double b)
{
  return (b-a)/2;
}

/**
 * Variance of a variate X uniformly distributed on ]a,b[.
 *
 * @param a is the interval lower bound.
 * @param b is the interval upper bound.
 *
 * @return the variance of an uniform distribution on ]a,b[.
 */
double st_uniform_variance(double a, double b)
{
  return ((b-a)*(b-a))/12;
}

/**
 * Median of a variate X uniformly distributed on ]a,b[.
 *
 * @param a is the interval lower bound.
 * @param b is the interval upper bound.
 *
 * @return the median of an uniform distribution on ]a,b[.
 */
double st_unforme_median(double a, double b)
{
  return (b-a)/2;
}
