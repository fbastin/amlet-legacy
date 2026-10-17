/**
 * poisson.h
 *
 * ORATIO Library
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_STAT_POISSON_H__
#define __ORATIO_STAT_POISSON_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Poisson distribution.
   */

  /**
   * Probability mass function of a Poisson random variable
   *
   * @param x is a positive integer where the function has to be evaluated.
   * @param lambda is the mean of the Poisson variable.
   *
   * @return the probability function of a Poisson of parameter lambda
   *   at x, i.e. P[X = x].
   */
  double st_poisson_pmf(int x, double lambda);

  /**
   * Compute the cumulative distribution function of a Poisson variate
   * of parameter lambda.
   *
   * @param x is a positive integer where the function has to be evaluated.
   * @param lambda is the mean of the Poisson variable.
   *
   * @return the cdf of a Poisson of parameter lambda at x, i.e. P[X <= x].
   */
  double st_poisson_cdf(int x, double lambda);

  /**
   * Compute the inverse cumulative distribution function of a Poisson variate
   * of parameter lambda.
   *
   * @param u is a real number in (0,1).
   * @param lambda is the mean of the Poisson variable.
   *
   * @return the smallest integer x such that F(x) >= u.
   */
  int st_poisson_icdf(double u, double lambda);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_STAT_POISSON_H__ */
