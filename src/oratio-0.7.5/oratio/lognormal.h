/**
 * lognormal.h
 *
 * STAT Library
 * Copyright (C) 1999-2000 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the MIT License
 */

/**
 * ROUTINES CONCERNING THE LOGNORMAL DISTRIBUTION
 */

#ifndef STAT_LOGNORMAL_H
#define STAT_LOGNORMAL_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

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
#endif /* __cplusplus */

#endif  /* __STAT_LOGNORMAL_H__ */
