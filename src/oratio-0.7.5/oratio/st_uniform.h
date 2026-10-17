/**
 * st_uniform.h
 *
 * ORATIO Library
 * Statistical routines
 * Copyright (C) 1999-2015 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_STAT_UNIFORM_H__
#define __ORATIO_STAT_UNIFORM_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * UNIFORM DISTRIBUTION
   */

  /**
   * Compute the probability distribution function of a variate X
   * uniformly distributed on ]a,b[.
   *
   * @param x is a real where the function has to be evaluated.
   * @param a is the interval lower bound.
   * @param b is the interval upper bound.
   *
   * @return the probability distribution function at x.
   */ 
 double st_uniform_pdf(double x, double a, double b);

  /**
   * Mean of a variate X uniformly distributed on ]a,b[.
   *
   * @param a is the interval lower bound.
   * @param b is the interval upper bound.
   *
   * @return the mean of an uniform distribution on ]a,b[.
   */
  double st_uniform_mean(double a, double b);

  /**
   * Variance of a variate X uniformly distributed on ]a,b[.
   *
   * @param a is the interval lower bound.
   * @param b is the interval upper bound.
   *
   * @return the variance of an uniform distribution on ]a,b[.
   */
  double st_uniform_variance(double a, double b);

  /**
   * Median of a variate X uniformly distributed on ]a,b[.
   *
   * @param a is the interval lower bound.
   * @param b is the interval upper bound.
   *
   * @return the median of an uniform distribution on ]a,b[.
   */
  double st_unforme_median(double a, double b);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_STAT_UNIFORM_H__ */
