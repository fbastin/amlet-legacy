/**
 * st_mvn.h
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

/**
 * Multivariate normal distribution
 */

#ifndef __ORATIO_STAT_MVN_H__
#define __ORATIO_STAT_MVN_H__

#include "cblas.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Convert a standard multivariate normal observation into an observation
   * of a multivariate normal with given mean and variance-covariance matrix.
   *
   * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
   * @param n is the dimension of the random vector.
   * @param L is a lower triangular matrix such that LL' gives the
   *   variance-covariance matrix.
   * @param ldL is the leading dimension of L.
   * @param mu is the mean vector. If mu is set to NULL, a zero mena vector is
   *   assumed.
   * @param incmu is the increment of mu.
   * @param x is the observation from a standard multivariate normal
   *   distribution on input, and an observation from the desired multivariate
   *   normal on output.
   * @param incx is the increment of x.
   */
  void st_mvn_convert(const enum CBLAS_ORDER Order, int n,
		      double *L, int ldL,
		      double *mu, int incmu,
		      double *x, int incx);

  /**
   * Convert an uniform vector into a vector from a multivariate normal with
   * given mean and variance-covariance matrix.
   *
   * @param n is the dimension of the random vector.
   * @param u is a random vector of dimension n whose components are
   *   observations taken from independent U(0,1).
   * @param x is a on output a vector of independant gaussian variates;
   *   x and u can point to the same memory address on input.
   * @param incx is the increment of x.
   */
  void st_mvn_standard_icdf(int n, double *u, int incu,
			    double *x, int incx);

  /**
   * Density function of a multivariate standard normal.
   *
   * @param s is the dimension of the random vector.
   * @param x is the point where the density function has to be evaluated.
   * @param incx is the increment of x.
   *
   * @return the value of the density function of a multivariate standard
   *   normal, evaluated at x.
   */
  double st_mvn_standard_pdf(int s, double *x, int incx);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_STAT_MVN_H__ */
