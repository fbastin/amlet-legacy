/**
 * exponential.c
 *
 * ORATIO Library - statistical module
 * Copyright (C) 1999-2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the MIT License
 */

#include "oratio/constants.h"
#include "oratio/ntools.h"
#include "oratio/st_descriptive.h"

#include <math.h>

/**
 * Compute the exponential probability distribution function.
 *
 * @param x is the point where to evaluate the probability distribution.
 * @param lambda is the exponential distribution parameter.
 *
 * @return the probability distribution of an exponential random
 *   variable, evaluated at x.
 */
double st_exponential_pdf(double x, double lambda)
{
  return (x >= 0) ? lambda*exp(-lambda*x) : 0;
}

/**
 * Compute the first derivative of exponential probability distribution
 * function with respect to the parameter distribution.
 *
 * @param x is the point where to evaluate the derivative of the probability
 *   distribution.
 * @param lambda is the exponential distribution parameter.
 *
 * @return the first derivative of the probability distribution of an
 *   exponential random variable, evaluated at x.
 */
double st_exponential_pdf_llderivative(double x, double lambda)
{
  return (x >= 0) ? exp(-lambda*x)*(1-lambda*x) : 0;  
}

/**
 * Compute the exponential cumulative distribution function.
 *
 * @param x is the point where to evaluate the distribution.
 * @param lambda is the exponential distribution parameter.
 *
 * @return the cumulative distribution of an exponential random
 *   variable, evaluated at x.
 */
double st_exponential_cdf(double x, double lambda)
{
  return (x >= 0) ? 1-exp(-lambda*x) : 0;
}

/**
 * Compute the inverse distribution function of an exponential random
 * variable.
 *
 * @param u is the point where to evaluated the inverse distribution
 *   function.
 * @param lambda is the exponential distribution parameter.
 *
 * @return the inverse distribution function of an exponential
 *   distribution, evaluated at u.
 */
double st_exponential_icdf(double u, double lambda)
{
  return (u >= 1.0) ? nt_infinity() : -log(1-u)/lambda;
}

/**
 * Compute the mean of an exponential random variable.
 *
 * @param lambda is the exponential distribution parameter.
 *
 * @return the mean of an exponential random variable of parameter lambda.
 */
double st_exponential_mean(double lambda)
{
  return 1.0/lambda;
}

/**
 * Compute the standard deviation of an exponential random variable.
 *
 * @param lambda is the exponential distribution parameter.
 *
 * @return the standard deviation of an exponential random variable of
 *   parameter lambda.
 */
double st_exponential_stdv(double lambda)
{
  return 1.0/lambda;
}

/**
 * Compute the variance of an exponential random variable.
 *
 * @param lambda is the exponential distribution parameter.
 *
 * @return the variance of an exponential random variable of parameter lambda.
 */
double st_exponential_variance(double lambda)
{
  return 1.0/(lambda*lambda);
}

/**
 * Compute the median of an exponential random variable.
 *
 * @param lambda is the exponential distribution parameter.
 *
 * @return the median of an exponential random variable of parameter lambda.
 */
double st_exponential_median(double lambda)
{
  return M_LN2/lambda;
}

/**
 * Compute the maximum likelihood estimator of the exponential parameter.
 *
 * @param n is the number of observations.
 */
double st_exponential_ml_estimator(int n, double *x, int incx, double *lambda)
{
  *lambda = 1.0/st_mean(x, incx, 0, n);

  return *lambda;
}
