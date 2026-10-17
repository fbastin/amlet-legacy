/**
 * logistic.h
 *
 * ORATIO Library
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_STAT_LOGISTIC_H__
#define __ORATIO_STAT_LOGISTIC_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Compute the cumulative distribution function of a logistic variate.
   *
   * @param x is the point where cdf has to be evaluated.
   * @param mu is the location parameter, also the mean.
   * @param s is the scale parameter.
   *
   * @return the cdf of a logistic variate.
   */
  double st_logistic_cdf(int x, double mu, double s);

  /**
   * Compute the inverse cumulative distribution function of a logistic
   * variate.
   *
   * @param u is the point where inverse cdf has to be evaluated.
   * @param mu is the location parameter, also the mean.
   * @param s is the scale parameter.
   *
   * @return the inverse cdf of a logistic variate.
   */
  double st_logistic_icdf(double u, double mu, double s);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_STAT_LOGISTIC_H__ */
