/**
 * st_gumbel.h
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2002-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Gumbel distribution
 */

#ifndef __ORATIO_STAT_GUMBEL_H__
#define __ORATIO_STAT_GUMBEL_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Cumulative distribution function of a Gumbel variate.
   *
   * @param x is the point where one wants to evaluate the cumulative
   *   distribution function.
   * @param a is the location parameter (the mode).
   * @param b is the scale parameter. b must be strictly positive.
   *
   * @return the cumulative distribution function of a Gumbel variable
   *   of mode a and scale parameter b, evaluated at x.
   */
  double st_gumbel_cdf(double x, double a, double b);
	
  /**
   * Inverse cumulative distribution function of a Gumbel variate.
   *
   * @param u is the value of the cumulative distribution function.
   * @param a is the location parameter (the mode).
   * @param b is the scale parameter. b must be strictly positive.
   *
   * @return the inverse cumulative distribution function of a Gumbel
   *   variate of mode a and scale parameter b, evaluated at u.
   */
  double st_gumbel_icdf(double u, double a, double b);
	
  /**
   * Probability distribution function of a Gumbel variate.
   *
   * @param a is the location parameter (the mode).
   * @param b is the scale parameter. b must be strictly positive.
   *
   * @return the probability distribution function of a Gumbel variable
   *   of mode a and scale parameter b, evaluated at x.
   */
   double st_gumbel_pdf(double x, double a, double b);

  /**
   * Compute the mean of a Gumbel variate
   *
   * @param a is the location parameter (the mode).
   * @param b is the scale parameter. b must be strictly positive.
   *
   * @return the mean of a Gumbel variate of parameters (a,b).
   */
  double st_gumbel_mean(double a, double b);

  /**
   * Compute the mean of a Gumbel variate
   *
   * @param a is the location parameter (the mode).
   * @param b is the scale parameter. b must be strictly positive.
   *
   * @return the mean of a Gumbel variate of parameters (a,b).
   */
  double st_gumbel_variance(double a, double b);

  /**
   * Compute the mode of a Gumbel variate.
   *
   * @param a is the location parameter (the mode).
   * @param b is the scale parameter. b must be strictly positive.
   *
   * Return the mode of a Gumbel variate of parameters (a,b).
   */
  double st_gumbel_mode(double a, double b);

  /**
   * Compute the median of a Gumbel variate.
   *
   * @param a is the location parameter (the mode).
   * @param b is the scale parameter. b must be strictly positive.
   *
   * @return the median of a Gumbel variate of parameters (a,b).
   */
  double st_gumbel_median(double a, double b);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_STAT_GUMBEL_H__ */
