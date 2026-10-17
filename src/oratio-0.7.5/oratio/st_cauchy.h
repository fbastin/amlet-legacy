/**
 * st_cauchy.h
 *
 * ORATIO Library
 * Copyright (C) 2008-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_STAT_CAUCHY_H__
#define __ORATIO_STAT_CAUCHY_H__

/**
 * Cauchy distribution.
 */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Compute the value of a Cauchy probabibility density function at x.
   *
   * @param x is the point where the Cauchy pdf has to be evaluated.
   * @param a is the location parameter (the median).
   * @param b is the scale parameter.
   *
   * @return the value of the Cauchy pdf at x.
   */
  double st_cauchy_pdf(double x, double a, double b);

  /**
   * Compute the value of a Cauchy cumulative distribution function at x.
   *
   * @param x is the point where the Cauchy cdf has to be evaluated.
   * @param a is the location parameter (the median).
   * @param b is the scale parameter.
   *
   * @return the value of the Cauchy cdf at x.
   */
  double st_cauchy_cdf(double x, double a, double b);

  /**
   * Inverse cumulative distribution of a Cauchy random variable.
   *
   * @param u is the value corresponding to the desired quantile.
   * @param a is the location parameter (the median).
   * @param b is the scale parameter.
   *
   * @return the value of the Cauchy icdf at u.
   */
  double st_cauchy_icdf(double u, double a, double b);
  
  /**
   * Compute the median of a Cauchy distribution.
   *
   * @param a is the location parameter (the median).
   * @param b is the scale parameter.
   *
   * @return the median of the Cauchy distribution.
   */
  double st_cauchy_median(double a, double b);

  /**
   * Compute the mode of a Cauchy distribution.
   *
   * @param a is the location parameter (the median).
   * @param b is the scale parameter.
   *
   * @return the mode of the Cauchy distribution.
   */
  double st_cauchy_mode(double a, double b);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_STAT_CAUCHY_H__ */
