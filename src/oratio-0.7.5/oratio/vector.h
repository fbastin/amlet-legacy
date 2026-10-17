/**
 * vector.h
 *
 * NTools Library
 * Copyright (C) 2002-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_NTOOLS_VECTOR_H__
#define __ORATIO_NTOOLS_VECTOR_H__

#include "norm.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Vector operations.
   *
   * The following functions aim to provide some convenient tools for
   * the manipulation of vectors, but for classical linear algebra
   * functions, please use a BLAS library.
   */

  /**
   * Set all vector components to zero.
   *
   * @param n is the vector dimension.
   * @param x is a real vector. On output, all components of x are equal to 0.
   * @param incx is the vector increment.
   */
  void nt_vect_zeros(int n, double *x, int incx);

  /**
   * Set all vector components to one.
   *
   * @param n is the vector dimension.
   * @param x is a real vector. On output, all components of x are equal to 1.
   * @param incx is the vector increment.
   */
  void nt_vect_ones(int n, double *x, int incx);

  int nt_vect_insert_vect(double x[], int nx, double y[], int ny, int a);

  /**
   * Fill selected elements of the vector x with the content of y.
   */
  void nt_vect_select(double x[], double y[], int *selector, int n);

  /**
   * Artithmetic operations
   */

  /**
   * Compute the distance between vectors x and y
   */
  double nt_vect_dist(double x[], double y[], int n);

  /**
   * Hadamard product
   *
   * y = alpha*x.*y
   */
  void nt_vect_dhad(int n, double alpha,
		    double x[], int incx, double y[], int incy);

  /**
   * Inverse Hadamard product
   *
   * y = alpha*x^(-1).*y
   */
  void nt_vect_dhad_inv(int n, double alpha,
			double x[], int incx, double y[], int incy);

  /**
   * Comparison functions
   */
  int nt_vect_is_zero(double x[], int n);

  int nt_vect_is_one(double x[], int n);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_NTOOLS_VECTOR_H__ */
