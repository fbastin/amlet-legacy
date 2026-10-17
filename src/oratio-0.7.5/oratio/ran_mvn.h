/**
 * ran_mvn.h
 *
 * ORATIO Library - Random Numbers Generation module
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Random numbers generator.
 */

#ifndef __ORATIO_RAN_MVN_H__
#define __ORATIO_RAN_MVN_H__

/**
 * Multivariate normal generation.
 */

#include "ran_rng.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Draw from a multivariate Gaussian distribution.
   *
   * @param ran is an initialized random numbers generator.
   * @param n is the dimension of the random vector.
   * @param x is on output a realization of the multivariate Gaussian
   *   random vector.
   * @param incx is the increment of x.
   */
  void ran_mvn_gaussian(Random *ran, int n, double *x, int incx);

  /**
   * Draw from a multivariate Gaussian distribution (Box-Muller variant).
   *
   * @param ran is an initialized random numbers generator.
   * @param n is the dimension of the random vector.
   * @param x is on output a realization of the multivariate Gaussian
   *   random vector.
   * @param incx is the increment of x.
   */
  void ran_mvn_gaussian_box_muller(Random *ran, int n,
				   double *x, int incx);

  /**
   * Draw from a multivariate normal distribution.
   *
   * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
   * @param ran is an initialized random numbers generator.
   * @param n is the dimension of the random vector.
   * @param L is a lower triangular matrix such that LL' gives the
   *   variance-covariance matrix.
   * @param ldL is the leading dimension of L.
   * @param mu is the mean vector. If mu is set to NULL, a zero mean vector is
   *   assumed.
   * @param incmu is the increment of mu.
   * @param x is on output a realization of the multivariate normal
   *   random vector.
   * @param incx is the increment of x.
   */
  void ran_mvn(Random *ran,
	       const enum CBLAS_ORDER Order, int n, double *L, int ldL,
	       double *mu, int incmu, double *x, int incx);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_RAN_MVN_H__ */
