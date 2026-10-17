/**
 * st_laplace.h
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_ST_LAPLACE_H__
#define __ORATIO_ST_LAPLACE_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Compute the Laplace probability distribution function.
   *
   * @param x is the point where to evaluate the probability distribution.
   * @param mu is the location parameter.
   * @param b is the scale parameter; b > 0.
   *
   * @return the probability distribution of a Laplace random variable,
   *   evaluated at x.
   */
  double st_laplace_pdf(double x, double mu, double b);

  /**
   * Compute the Laplace cumulative distribution function.
   *
   * @param x is the point where to evaluate the cdf.
   * @param mu is the location parameter.
   * @param b is the scale parameter; b > 0.
   *
   * @return the cdf of a Laplace random variable, evaluated at x.
   */
  double st_laplace_cdf(double x, double mu, double b);

  /**
   * Compute the inverse distribution function of a Laplace random
   * variable.
   *
   * @param u is the point where to evaluated the inverse distribution
   *   function.
   * @param mu is the location parameter.
   * @param b is the scale parameter; b > 0.
   *
   * @return the inverse distribution function of an exponential
   *   distribution, evaluated at u.
   */
  double st_laplace_icdf(double u, double mu, double b);

  /**
   * Compute the mean of an Laplace random variable.
   *
   * @param mu is the location parameter.
   * @param b is the scale parameter; b > 0.
   *
   * @return the mean of a Laplace random variable of parameters mu et b.
   */
  double st_laplace_mean(double mu, double b);

  /**
   * Compute the mean of an Laplace random variable.
   *
   * @param mu is the location parameter.
   * @param b is the scale parameter; b > 0.
   *
   * @return the mean of a Laplace random variable of parameters mu et b.
   */
  double st_laplace_median(double mu, double b);

  /**
   * Compute the mean of an Laplace random variable.
   *
   * @param mu is the location parameter.
   * @param b is the scale parameter; b > 0.
   *
   * @return the mode of a Laplace random variable of parameters mu et b.
   */
  double st_laplace_mode(double mu, double b);

  /**
   * Compute the standard deviation of an exponential random variable.
   *
   * @param mu is the location parameter.
   * @param b is the scale parameter; b > 0.
   *
   * @return the standard deviation of a Laplace random variable of parameters
   *   mu et b.
   */
  double st_laplace_stdv(double mu, double b);

  /**
   * Compute the variance of an Laplace random variable.
   *
   * @param mu is the location parameter.
   * @param b is the scale parameter; b > 0.
   *
   * @return the variance of a Laplace random variable of parameters mu et b.
   */
  double st_laplace_variance(double mu, double b);

#ifdef __cplusplus
}
#endif

#endif  /* _ORATIO_ST_LAPLACE_H__ */
