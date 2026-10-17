/**
 * st_triangular.h
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2002-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Triangular distribution
 */

#ifndef __ORATIO_STAT_TRIANGULAR_H__
#define __ORATIO_STAT_TRIANGULAR_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Triangular distribution.
   */

  /**
   * Compute the density of a Triangular distribution at x.
   *
   * @param x is the point where the pdf has to be evaluated.
   * @param a is the lower limit.
   * @param b is the upper limit.
   * @param c is the mode of the distribution.
   *
   * @return the density of the Triangular distribution evaluated at x.
   */
  double st_triangular_pdf(double x, double a, double b, double c);

  /**
   * Compute the cumulative distribution function of a Triangular
   * distribution at x.
   *
   * @param x is the point where the cdf has to be evaluated.
   * @param a is the lower limit.
   * @param b is the upper limit.
   * @param c is the mode of the distribution.
   *
   * @return the cumulative distribution function of the Triangular
   *   distribution evaluated at x.
   */
  double st_triangular_cdf(double x, double a, double b, double c);

  /**
   * Compute the inverse cumulative distribution function of a Triangular
   * distribution.
   *
   * @param p is the value (0 <= p <= 1) with respect we want to compute
   *   the associated quantile.
   * @param a is the lower limit.
   * @param b is the upper limit.
   * @param c is the mode of the distribution.
   *
   * @return the inverse cumulative distribution function of the Triangular
   *   distribution evaluated at p.
   */
  double st_triangular_icdf(double u, double a, double b, double c);

  /**
   * Compute the mode of a Triangular distribution.
   *
   * @param a is the lower limit.
   * @param b is the upper limit.
   * @param c is the mode of the distribution.
   *
   * @return the mode of the Triangular distribution.
   */
  double st_triangular_mean(double a, double b, double c);

  /**
   * Compute the variance of a Triangular distribution.
   *
   * @param a is the lower limit.
   * @param b is the upper limit.
   * @param c is the mode of the distribution.
   *
   * @return the variance of the Triangular distribution.
   */
  double st_triangular_variance(double a, double b, double c);

  /**
   * Compute the standard deviation of a Triangular distribution.
   *
   * @param a is the lower limit.
   * @param b is the upper limit.
   * @param c is the mode of the distribution.
   *
   * @return the standard deviation of the Triangular distribution.
   */
  double st_triangular_stdv(double a, double b, double c);

  /**
   * Compute the mode of a Triangular distribution.
   *
   * @param a is the lower limit.
   * @param b is the upper limit.
   * @param c is the mode of the distribution.
   *
   * @return the mode of the Triangular distribution.
   */
  double st_triangular_mode(double a, double b, double c);

  /**
   * Compute the median of a Triangular distribution.
   *
   * @param a is the lower limit.
   * @param b is the upper limit.
   * @param c is the mode of the distribution.
   *
   * @return the median of the Triangular distribution.
   */
  double st_triangular_median(double a, double b, double c);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_STAT_TRIANGULAR_H__ */
