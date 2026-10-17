/**
 * chi2.h
 *
 * ORATIO Library
 * Copyright (C) 1999-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_STAT_CHI2_H__
#define __ORATIO_STAT_CHI2_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Chi 2 distribution.
   */

  /**
   * Return the value of a chi-square probabibility density function,
   * at x.
   *
   * @param x is the point where the probability density is evaluated.
   * @param k is the number of degrees of freedom.
   *
   * @return the value of the chi 2 pdf at x.
   */
  double st_chi2_pdf(double x, double k);

  /**
   * Compute the chi-square cumulative distribution function.
   *
   * @param k is the number of degrees of freedom.
   */
  double st_chi2_cdf(double x, double k);

  /**
   * Return the mean of a chi-square distribution with k degrees of freedom.
   *
   * @param x is the point where the cumulative density is evaluated.
   * @param k is the number of degrees of freedom.
   */
  double st_chi2_mean(double k);

  /**
   * Return the variance of a chi-square distribution with k degrees
   * of freedom.
   *
   * @param k is the number of degrees of freedom.
   */
  double st_chi2_variance(double k);

  /**
   * Return the skewness of a chi-square distribution with k degrees of
   * freedom.
   *
   * @param k is the number of degrees of freedom.
   *
   * @return the skewness of the distribution.
   */
  double st_chi2_skewness(double k);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_STAT_CHI2_H__ */
