/**
 * st_normal.h
 *
 * ORATIO Library
 * Copyright (C) 1999-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the BSD License
 */

#ifndef __ORATIO_ST_NORMAL_H__
#define __ORATIO_ST_NORMAL_H__

#ifdef __cplusplus
extern "C" {
#endif


  /**
   * Normal distribution.
   */

  /**
   * Returns the value of a Gaussian probability distribution function
   * evaluated at x.
   *
   * @param x is the point where the probability distribution is evaluated.
   *
   * @return the value of the Gaussian pdf evaluated at x.
   */
  double st_gaussian_pdf(double x);

  /**
   * Evaluate the density function of a normal function of mean mu and
   * standard deviation sigma, at x.
   *
   * @param x is the point where the probability distribution is evaluated.
   * @param mu is the mean of the distribution.
   * @param sigma is the standard deviation of the distribution.
   *
   * @return the value of the density fonction at x.
   */
  double st_normal_pdf(double x, double mu, double sigma);

  /**
   * Returns the value of a Gaussian cumulative distribution function
   * evaluated at x, using Cody's algorithm.
   *
   * @param x is the upper limit of integration
   *
   * @return the value of a Gaussian distribution at x.
   */
  double st_gaussian_cody_cdf(double x);

  /**
   * Returns the value of a Gaussian cumulative distribution function
   * evaluated at x.
   *
   * @param x: upper limit of integration
   *
   * @return the Gaussian cdf at x.
   */
  double st_gaussian_cdf(double x);

  /**
   * Returns the value at x of a Gaussian complementary cumulative distribution
   * (ccdf) function. Marsaglia "fast" implementation.
   *
   * @param x is the point where the ccdf has to be evaluated.
   *
   * @return the Gaussian ccdf value at x.
 */
  double st_gaussian_ccdf(double x);

  /**
   * Compute an approximation of the cumulative function of a normal
   * distribution, of mean mu and standard deviation sigma, at x.
   *
   * @param x is the point where the cdf has to be evaluated.
   * @param mu is the mean of the distribution.
   * @param sigma is the standard deviation of the distribution.
   */
  double st_normal_cdf(double x, double mu, double sigma);

  /**
   * Compute the Gaussian cumulative distribution function.
   *
   * @param p is the value (0 <= p <= 1) with respect we want to compute
   *   the associated quantile.
   *
   * @return the value of the inverse cumulative function of a
   *   gaussian distribution evaluated at p.
   */
  double st_gaussian_icdf(const double p);

  /**
   * Compute the Normal inverse cumulative distribution function.
   *
   * @param p is the value (0 <= p <= 1) with respect we want to compute
   *   the associated quantile.
   * @param mu is the mean of the normal distribution.
   * @param sigma is the standard deviation of the normal distribution.
   *
   * @return the value of the inverse cumulative function of a
   *   normal distribution of mean mu and standard deviation sigma.
   */
  double st_normal_icdf(const double q, double mu, double sigma);

  /**
   * Return the mean of a normal distribution.
   *
   * @param mu is the mean of the normal distribution.
   * @param sigma is the standard deviation of the normal distribution.
   *
   * @return the mean of the normal distribution.
   */
  double st_normal_mean(double mu, double sigma);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_ST_NORMAL_H__ */
