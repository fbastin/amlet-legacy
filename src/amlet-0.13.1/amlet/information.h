/**
 * information.h
 *
 * AMLET (Another Mixed Logit Estimation Tool) Library
 * Copyright (C) 2002-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __AMLET_INFORMATION_H__
#define __AMLET_INFORMATION_H__

#include <oratio/cblas.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * t-test
   *
   * @param n is the number of parameters.
   * @param theta if the vector of parameters.
   * @param I is a matrix having parameters variances on its diagonal
   *   and zero as other entries.
   *        alpha: level used for the test
   * @param t is a vector of size n on input, and contains the
   *        t-statistics on output
   *
   * @return the percentile used for the test.
   */
  double amlet_t_statistics(int n, double *theta, double *hypothetical,
			    double **I, double alpha, double *t);

  /**
   * Compute the inverse of the likehood Hessian.
   *
   * @param m is an estimated model, including the maximum likelihood Hessian.
   */
  void amlet_likelihood_Hinv(AMLET_Model *m);

  /**
   * Compute the correlation matrix.
   *
   * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
   * @param Uplo designs which half-matrix is filled (CblasUpper or
   *   CblasLower).
   * @param V is the variance/covariance matrix.
   *
   * #POST: The correlation values are stocked in triangular off-diagonal of
   *        of V that is not used for covariance values.
   */
  void amlet_correlation(const enum CBLAS_ORDER Order,
			 const enum CBLAS_UPLO Uplo,
			 double **V, int n);

  void amlet_weak_covariance_matrix(AMLET_Model *m, int n,
				    double scale);

  /**
   * Compute the information and correlation matrices.
   */
  void amlet_robust_covariance_matrix(AMLET_Model *m, int n,
				      double scale);

  /**
   * Compute the information and correlation matrices.
   *
   * @param n is the number of parameters.
   * @param H is the Hessian matrix (or an approximation) of the
   *   objective function; it can be set to NULL if the estimation is
   *   no robust.
   * @param B is the outer-product based matrix.
   * @param robust should be set to 1 for a robust matrix computation, 0
   *   else.
   * @param scale is a scale factor, usually set as 1.0/nind, where nind
   *   is the number of individuals/
   */
  int amlet_covariance_matrix(AMLET_Model *m,
			       int n, double **B, int robust,
			       double scale);

  void amlet_dump_information_matrices(AMLET_Model *m);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_INFORMATION_H__ */
