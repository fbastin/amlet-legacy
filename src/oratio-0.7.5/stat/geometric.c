/**
 * geometric.c
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Geometric distribution
 */

#include "oratio/ntools.h"

/**
 * Probability mass function of a geometric random variable.
 *
 * @param k is a natural where the probability function has to be evaluated.
 * @param p is the Bernouilli parameter of the random distributiun; 0 < p <= 1.
 *
 * @return the probability function at k of a geometric random
 *   variable of parameter p.
 */
double st_geometric_pmf(int k, double p)
{
  return p*nt_ipow(1-p, k);
}

/**
 * Cumulative distribution function of a geometric random variable.
 *
 * @param k is a natural where the distribution function has to be evaluated.
 * @param p is the Bernouilli parameter of the random distributiun; 0 < p <= 1.
 *
 * @return the distribution function at k of a geometric random
 *   variable of parameter p.
 */
double st_geometric_cdf(int k, double p)
{
  return 1-nt_ipow(1-p, k+1);
}

/**
 * Inverse cumulative distribution function of a geometric random variable.
 *
 * @param u is a real number in (0,1).
 * @param p is the Bernouilli parameter of the random distributiun; 0 < p <= 1.
 *
 * @return the u-quantile of the geometric variate of parameter p.
 */
int st_geometric_icdf(double u, double p)
{
  if (p == 0 || p == 1) return p;
  return floor(log(1-u)/log(1-p));
}

/**
 * Mean of a geometric random variable.
 *
 * @param p is the Bernouilli parameter of the random distributiun; 0 < p <= 1.
 *
 * @return the mean of a geometric random variable of parameter p.
 */
double st_geometric_mean(double p)
{
  return (1-p)/p;
}

/**
 * Mean of a geometric random variable.
 *
 * @param p is the Bernouilli parameter of the random distributiun; 0 < p <= 1.
 *
 * @return the mean of a geometric random variable of parameter p.
 */
double st_geometric_variance(double p)
{
  return (1-p)/(p*p);
}

/**
 * Mean of a geometric random variable.
 *
 * @param p is the Bernouilli parameter of the random distributiun; 0 < p <= 1.
 *
 * @return the mean of a geometric random variable of parameter p.
 */
double st_geometric_mode(double p)
{
  return 0;
}
