/**
 * gamma.h
 *
 * ORATIO Library
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_NTOOLS_GAMMA_H__
#define __ORATIO_NTOOLS_GAMMA_H__

#ifdef __cplusplus
extern "C" {
#endif  /* __cplusplus */

  /**
   * Ratio of two gamma functions.
   *
   * @param x is a positive real; the numerator is a gamma function evaluated
   *   x+1/2, while the denominator is a gamma function evaluated at x.
   *
   * @return the value of the ratio.
   */
  double nt_gamma_ratio_half(double x);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __ORATIO_NTOOLS_GAMMA_H__ */
