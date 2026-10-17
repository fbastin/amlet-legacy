/**
 * rayleigh.h
 *
 * ORATIO Library
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_STAT_RAYLEIGH_H__
#define __ORATIO_STAT_RAYLEIGH_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Compute the cumulative distribution function of a Rayleigh variate.
   *
   * @param x is the point where cdf has to be evaluated.
   * @param beta is the Rayleigh parameter.
   *
   * @return the cdf of a Rayleigh variate.
   */
  double st_rayleigh_cdf(int x, double beta);

  /**
   * Compute the inverse cumulative distribution function of a Rayleigh
   * variate.
   *
   * @param u is the point where inverse cdf has to be evaluated.
   * @param beta is the Rayleigh parameter.
   *
   * @return the inverse cdf of a Rayleigh variate.
   */
  double st_rayleigh_icdf(double u, double beta);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_STAT_RAYLEIGH_H__ */
