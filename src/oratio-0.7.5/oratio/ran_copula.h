/**
 * ran_copula.h
 *
 * ORATIO Library - Random Number Generation Module
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the BSD License.
 */

#ifndef __ORATIO_RAN_COPULA_H__
#define __ORATIO_RAN_COPULA_H__

#include "cblas.h"
#include "ran_rng.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Generate a Gaussian copula with uniform marginals.
   *
   * @param ran is an initialized random numbers generator.
   * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
   * @param n is the dimension of the random vector.
   * @param L is a lower triangular matrix such that LL' gives the
   *   correlation  matrix.
   * @param ldL is the leading dimension of L.
   * @param x is on output a realization of the multivariate normal
   *   random vector.
   * @param incx is the increment of x.
   */
  void ran_copula_gaussian_unif(Random *ran,
				const enum CBLAS_ORDER Order,
				int n, double *L, int ldL,
				double *x, int incx);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_RAN_COPULA_H__ */
