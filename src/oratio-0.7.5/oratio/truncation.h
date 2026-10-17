/**
 * truncation.h
 *
 * STAT Library
 * Copyright (C) 2006 Fabian Bastin <bastin@cerfacs.fr>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __STAT_TRUNCATION_H__
#define __STAT_TRUNCATION_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Evaluation of the probability distribution of a truncated
   * distribution.
   *
   * @param Fmin is the value of the CDF at the lower bound of the
   *   interval.
   * @param Fmax is the value of the CDF at the upper bound of the
   *   interval.
   *
   * @return the probability function value of the truncated distribution.
   */
  double st_truncated_pdf(double fx, double Fmin, double Fmax);

  /**
   * Evaluation of the cumulative distribution function of a truncated
   * distribution.
   *
   * @param Fmin is the value of the CDF at the lower bound of the
   *   interval.
   * @param Fmax is the value of the CDF at the upper bound of the
   *   interval.
   *
   * @return the CDF value of the truncated distribution.
   */
  double st_truncated_cdf(double Fx, double Fmin, double Fmax);

#ifdef __cplusplus
}
#endif

#endif  /* __STAT_TRUNCATION_H__ */
