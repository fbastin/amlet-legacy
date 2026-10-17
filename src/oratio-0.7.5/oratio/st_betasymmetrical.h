/**
 * st_betasymmetrical.h
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_ST_BETASYMMETRICAL_H__
#define __ORATIO_ST_BETASYMMETRICAL_H__

/**
 * Symmetrical Beta distribution
 */

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Compute the mean of a symmetrical beta distribution of parameter alpha.
   *
   * @param alpha is the (shape) parameter of the symmetrical beta
   *   distribution; alpha must be positive.
   *
   * @return the mean of the distribution.
   */
  double st_betasymmetrical_mean(double alpha);

  /**
   * Compute the inverse cumulative distribution function of a symmetrical
   * beta distribution of parameter alpha.
   *
   * @param u is the point where the icdf has to be evaluated; u must be in
   *   [0,1].
   * @param alpha is the (shape) parameter of the symmetrical beta
   *   distribution; alpha must be positive.
   *
   * @return the value of the icdf at u.
   */
  double st_betasymmetrical_icdf(double u, double alpha);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_ST_BETASYMMETRICAL_H__ */
