/**
 * pareto.h
 *
 * ORATIO Library
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_STAT_PARETO_H__
#define __ORATIO_STAT_PARETO_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Compute the cumulative distribution function of a Pareto variate.
   *
   * @param x is the point where cdf has to be evaluated.
   * @param alpha is the Pareto index.
   * @param beta is the mode of the distribution.
   *
   * @return the cdf of a Pareto variate of mode beta and index alpha.
   */
  double st_pareto_cdf(double x, double alpha, double beta);

  /**
   * Compute the inverse cumulative distribution function of a Pareto variate.
   *
   * @param u is the point where inverse cdf has to be evaluated.
   * @param alpha is the Pareto index.
   * @param beta is the mode of the distribution.
   *
   * @return the inverse cdf of a Pareto variate of mode beta and index alpha.
   */
  double st_pareto_icdf(double u, double alpha, double beta);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_STAT_PARETO_H__ */
