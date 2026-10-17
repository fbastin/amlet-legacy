/**
 * st_geometric.h
 *
 * ORATIO Library
 * Copyright (C) 2005-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_STAT_GEOMETRIC_H__
#define __ORATIO_STAT_GEOMETRIC_H__

/**
 * Cauchy distribution.
 */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Geometric distribution.
   */

  /**
   * Probability mass function of a geometric random variable.
   *
   * @param k is a natural where the probability function has to be evaluated.
   * @param p is the bernouilli parameter of the random distributiun;
   *  0 < p <= 1.
   *
   * @return the probability function at k of a geometric random
   *   variable of parameter p.
   */
  double st_geometric_pmf(int k, double p);

  /**
   * Cumulative distribution function of a geometric random variable.
   *
   * @param k is a natural where the distribution function has to be evaluated.
   * @param p is the bernouilli parameter of the random distributiun;
   *   0 < p <= 1.
   *
   * @return the distribution function at k of a geometric random
   *   variable of parameter p.
   */
  double st_geometric_cdf(int k, double p);

  /**
   * Inverse cumulative distribution function of a geometric random variable.
   *
   * @param u is a real number in (0,1).
   * @param p is the bernouilli parameter of the random distributiun;
   *  0 < p <= 1.
   *
   * @return the u-quantile of the geometric variate of parameter p.
   */
  int st_geometric_icdf(double u, double p);

  /**
   * Mean of a geometric random variable.
   *
   * @param p is the bernouilli parameter of the random distributiun;
   *  0 < p <= 1.
   *
   * @return the mean of a geometric random variable of parameter p.
   */
  double st_geometric_mean(double p);

  /**
   * Mean of a geometric random variable.
   *
   * @param p is the bernouilli parameter of the random distributiun;
   *   0 < p <= 1.
   *
   * @return the mean of a geometric random variable of parameter p.
   */
  double st_geometric_variance(double p);

  /**
   * Mean of a geometric random variable.
   *
   * @param p is the bernouilli parameter of the random distributiun;
   *  0 < p <= 1.
   *
   * @return the mean of a geometric random variable of parameter p.
   */
  double st_geometric_mode(double p);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_STAT_GEOMETRIC_H__ */
