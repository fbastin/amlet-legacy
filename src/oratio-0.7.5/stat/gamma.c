/**
 * gamma.c
 *
 * STAT Library
 * Copyright (C) 1999-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/stat.h"

#include <math.h>

#define PLIM 1.0e5

/**
 * Gamma probability distribution function.
 *
 * @param x is the point where the probability distribution is evaluated.
 * @param alpha > 0 is the shape parameter of the gamma distribution.
 * @param lambda > 0 is the scale parameter of the gamma distribution.
 *
 * @return the value at x of the density of a gamma variate of shape
 *   parameter k and scale parameter theta.
 */
double st_gamma_pdf(double x, double alpha, double lambda)
{
  if (alpha == 1.0) return st_exponential_pdf(x, lambda);

  return (x <= 0) ? 0.0 :
    exp(alpha*log(lambda*x)-lambda*x-lgamma(alpha))/x;
}

/* inline double st_gamma_cdf(double x, double alpha, double lambda)
{
  if (x <= 0.0)
    return 0.0;
  if (alpha == 1.0)
    return st_exponential_cdf(x, 1.0);

  if (alpha >= PLIM) {
    double d2, S, z;

    d2 = x + 1.0/3.0 - alpha - 0.02/alpha;
    S = alpha - 1.0/2.0;
    z = d2 * sqrt((1 + belog(S/x))/x);

    return st_gaussian_cdf(z);
  }

  if (x <= 1.0 || x < alpha) {
    double factor, gin, rn, term;
    double acu = EPSARRAY[d];

    factor = exp(alpha*(x) - x - st_gamma_ln(alpha));

    gin = 1.0;
    term = 1.0;
    rn = alpha;
    do {
      rn += 1.0;
      term *= x/rn;
      gin += term;
    } while (term >= acu * gin);
    return gin*factor/alpha;
  }
  else
    return 1.0 - barF (alpha, d, x);
    } */

 /* inline double st_gamma_icdf(double u, double alpha, double lambda)
{
} */

/**
 * Mean of a Gamma distribution function of parameters alpha and lambda.
 *
 * @param alpha > 0 is the shape parameter of the gamma distribution.
 * @param lambda > 0 is the scale parameter of the gamma distribution.
 *
 * @return the mean of a Gamma distribution function of parameters
 *   alpha and lambda.
 */
double st_gamma_mean(double alpha, double lambda)
{
  return alpha/lambda;
}

/**
 * Variance of a Gamma distribution function of parameters alpha and lambda.
 *
 * @param alpha > 0 is the shape parameter of the gamma distribution.
 * @param lambda > 0 is the scale parameter of the gamma distribution.
 *
 * @return the variance of a Gamma distribution function of parameters
 *   alpha and lambda.
 */
double st_gamma_variance(double alpha, double lambda)
{
  return alpha/(lambda*lambda);
}

/**
 * Standard of a Gamma distribution function of parameters alpha and
 * lambda.
 *
 * @param alpha > 0 is the shape parameter of the gamma distribution.
 * @param lambda > 0 is the scale parameter of the gamma distribution.
 *
 * @return the standard deviation of a Gamma distribution function of
 *   parameters alpha and lambda.
 */
double st_gamma_stdv(double alpha, double lambda)
{
  return sqrt(alpha)/lambda;
}

/**
 * Mode of a Gamma distribution function of parameters alpha and lambda.
 *
 * @param alpha >= 1.0 is the shape parameter of the gamma distribution.
 * @param lambda > 0 is the scale parameter of the gamma distribution.
 *
 * @return the mode of a Gamma distribution function of parameters
 *   alpha and lambda.
 */
double st_gamma_mode(double alpha, double lambda)
{
  return (alpha-1.0)/lambda;
}
