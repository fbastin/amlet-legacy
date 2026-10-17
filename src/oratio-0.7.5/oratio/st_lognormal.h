/**
 * st_lognormal.h
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2002-2014 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

/**
 * Lognormal distribution
 */

#ifndef __ORATIO_STAT_LOGNORMAL_H__
#define __ORATIO_STAT_LOGNORMAL_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Evaluate the probability density function of the lognormal variate
   * L of parameters mu and sigma, the mean and standard deviation of log L.
   *
   * @param x is the point where the probability density function is
   *   evaluated.
   *
   * @param mu is the mean of log L.
   * @param sigma is the standard deviation of log L.
   *
   * @return the PDF value of L at x.
   */
  double st_lognormal_pdf(double x, double mu, double sigma);

  /**
   * Evaluate the cumulative density function of the lognormal L of
   * parameters mu and sigma, the mean and standard deviation of log L.
   *
   * @param x is the point where the cumulative distribution function is
   *   evaluated. x must be positive.
   * @param mu is the mean of log L.
   * @param sigma is the standard deviation of log L.
   * @param ccum contains on output the complement of result, i.e. ccum
   *   = 1.0 - result.
   *
   * @return the value of the CDF of L at x.
   */
  double st_lognormal_cdf(double x, double mu, double sigma);

  double st_lognormal_icdf(double q, double mu, double sigma);

  double st_lognormal_mean(double mu, double sigma);

  double st_lognormal_stdv(double mu, double sigma);

  double st_lognormal_variance(double mu, double sigma);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_STAT_LOGNORMAL_H__ */
