/**
 * gamma.c
 *
 * ORATIO Library
 * Copyright (C) 2005-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/ntools.h>

#include <float.h>
#include <math.h>

/**
 * Ratio of two gamma functions.
 *
 * @param x is a positive real; the numerator is a gamma function evaluated
 *   x+1/2, while the denominator is a gamma function evaluated at x.
 *
 * @return the value of the ratio.
 */
double nt_gamma_ratio_half(double x)
{
  if (x <= 0.0) return nt_nan();

  if (x <= 10.0) {
    double y = lgamma(x + 0.5) - lgamma (x);
    return exp(y);
  }

  double sum;
  if (x <= 300.0) {
    /* The sum converges very slowly for small x, but faster as x increases */
    double term = 1.0;
    int i = 1;

    sum = 1.0;
    while (term > DBL_EPSILON*sum) {
      term *= (i - 1.5)*(i - 1.5) /(i*(x + i - 1.5));
      sum += term;
      i++;
    }
    return sqrt((x - 0.5)*sum);
  }

  /* Asymptotic series for Gamma(x + 0.5) / (Gamma(x) * Sqrt(x)) */
  double y = 1.0 / (8.0*x);
  sum = 1.0 + y*(-1.0 + y*(0.5 + y*(2.5 - y*(2.625 + 49.875*y))));
  return sum*sqrt(x);
}
