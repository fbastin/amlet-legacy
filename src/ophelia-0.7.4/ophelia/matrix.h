/**
 * matrix.h
 *
 * OPHELIA Library
 * Copyright (C) 2005-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __OPHELIA_MATRIX_H__
#define __OPHELIA_MATRIX_H__

#include <oratio/f2c.h>
#include <oratio/cblas.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Copy the off-diagonal lower part of a square matrix in the off-diagonal
   * upper-part.
   *
   * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
   * @param n is the matrix order.
   * @param A is a real square matrix nxn.
   * @param ldA is the leading dimension of A.
   */
  void op_matrix_copy_lower(const enum CBLAS_ORDER Order,
				   int n, double *A, int ldA);

  /**
   * Copy the off-diagonal lower part of a square matrix in the off-diagonal
   * upper-part.
   *
   * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
   * @param n is the matrix order.
   * @param A is a real square matrix nxn.
   * @param ldA is the leading dimension of A.
   */
  void op_matrix_copy_upper(const enum CBLAS_ORDER Order,
				   int n, double *A, int ldA);

  /**
   * Compute the inverse of a symmetric indefinite square matrix.
   *
   * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
   * @param Uplo designs which half-matrix is filled (CblasUpper or
   *   CblasLower).
   * @param I is the matrix to invert. The storage on input and output
   *   is given by Order and Uplo.
   * @param npar is the matrix main dimension.
   *
   * @return the information value delivered but Lapack routines used to
   *   compute the inverse.
   */
  int op_matrix_syin_inverse(const enum CBLAS_ORDER Order,
			     const enum CBLAS_UPLO Uplo,
			     double *I, int npar);

  /**
   * Compute the inverse of a symmetric positive definite square matrix.
   *
   * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
   * @param Uplo designs which half-matrix is filled (CblasUpper or
   *   CblasLower).
   * @param I is the matrix to invert. The storage on input and output
   *   is given by Order and Uplo.
   * @param npar is the matrix main dimension.
   *
   * @return the information value delivered but Lapack routines used to
   *   compute the inverse.
   */
  int op_matrix_sydp_inverse(const enum CBLAS_ORDER Order,
			     const enum CBLAS_UPLO Uplo,
			     double *I, int npar);

  /**
   * Inverse a real symmetric matrix.
   *
   * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
   * @param Uplo designs which half-matrix is filled (CblasUpper or
   *   CblasLower).
   * @param I is the matrix to invert. The storage on input and output
   *   is given by Order and Uplo.
   * @param n is the matrix main dimension.
   *
   * @return the information value delivered but Lapack routines used to
   *   compute the inverse.
   */
  int op_matrix_sy_inverse(const enum CBLAS_ORDER Order,
			   const enum CBLAS_UPLO Uplo,
			   double *I, int n);

#ifdef __cplusplus
}
#endif

#endif  /* __OPHELIA_MATRIX_H__ */
