/**
 * poisson.c
 *
 * ORATIO Library - statistical module
 * Copyright (C) 1999-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/ntools.h"
#include "oratio/st_utils.h"

#include <math.h>

double factorial_ln(int n)
{
  static float a[101];  /* a statically allocated array is always
			   initialized to 0. */

  if (n <= 1) return 0.0;
  if (n <= 100)
    return a[n] ? a[n] : (a[n] = lgamma(n+1.0));
  else
    return lgamma(n+1.0);
}

/**
 * Compute the probability function of a Poisson variate of parameter lambda.
 *
 * @param x is a positive integer where the function has to be evaluated.
 * @param lambda is the mean of the Poisson variable.
 *
 * @return the probability function of a Poisson of parameter lambda
 *   at x, i.e. P[X = x].
 */
double st_poisson_pmf(int x, double lambda)
{
  double p;

  p = exp(log(lambda)*x-factorial_ln(x)-lambda);

  return p;
}

/**
 * Compute the inverse cumulative distribution function of a Poisson variate
 * of parameter lambda.
 *
 * @param u is a real number in (0,1).
 * @param lambda is the mean of the Poisson variable.
 *
 * @return the smallest integer x such that F(x) >= u.
 */
int st_poisson_icdf(double u, double lambda)
{
  int x = 0;
  /* p is the probability mass; f is the cdf */
  double p, f;

  p = f = exp(lambda);

  while (f < u) {
    x++;
    p = lambda/x*p;
    f = f+p;
  }

  return x;
}

/**
 * Compute the cumulative distribution function of a Poisson variate
 * of parameter lambda.
 *
 * @param x is a positive integer where the function has to be evaluated.
 * @param lambda is the mean of the Poisson variable.
 *
 * @return the cdf of a Poisson of parameter lambda at x, i.e. P[X <= x].
 */
double st_poisson_cdf(int x, double lambda)
{
  return st_gamma_q(x+1.0, lambda);
}
