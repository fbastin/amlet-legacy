/**
 * st_weibull.h
 *
 * ORATIO Library
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_STAT_WEIBULL_H__
#define __ORATIO_STAT_WEIBULL_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Compute the cumulative distribution function of a Weibull variate.
   *
   * @param x is the point where cdf has to be evaluated.
   * @param alpha is the shape parameter; alpha > 0.
   * @param beta is the scale parameter; beta > 0.
   *
   * @return the cdf of a Weibull of parameters alpha and beta.
   */
  double st_weibull_cdf(double x, double alpha, double beta);

  /**
   * Compute the inverse cumulative distribution function of a Weibull variate.
   *
   * @param u is the point where inverse cdf has to be evaluated.
   * @param alpha is the shape parameter; alpha > 0.
   * @param beta is the scale parameter; beta > 0.
   *
   * @return the inverse cdf of a Weibull of parameters alpha and beta.
   */
  double st_weibull_icdf(double u, double alpha, double beta);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_STAT_WEIBULL_H__ */
