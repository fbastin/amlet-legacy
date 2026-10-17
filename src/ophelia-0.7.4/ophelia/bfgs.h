/**
 * bfgs.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2002-2003 Laurent Martin
 *                         Emilie Wanufelle <ewa@math.fundp.ac.be>,
 *                         Fabian Bastin    <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __NLP2_BFGS_H__
#define __NLP2_BFGS_H__

#include "bfgs_struct.h"
#include "io.h"

#include <oratio/ntools.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * BFGS method initialization.
   *
   * @param bfgs is a memory allocated BFGS structure.
   * @param n is the main problem dimension.
   * @param m is the secondary problem dimension, usually set to 0.
   */
  void nlp_bfgs_init(BFGS *bf, int n, int m);

  /**
   * Positive definite secant update (BFGS), unfactored form.
   *
   * Ref.: J.E. Dennis, Jr. and Robert B. Schnabel, "Numerical Methods
   * for Unconstrained Optimization and Nonlinear Equations", SIAM,
   * 1983, Algorithm 9.4.1, page 355.
   *
   * The implementation in Algorithm 9.4.1 is more conservative than the
   * current, due to the use of an additional parameter eta, described
   * in the guidelines 2, in particular section "fdigits".
   *
   * @param Order designs the orientation of the matrix Bk. The
   *   possible values are CblasRowMajor and CblasColMajor.
   * @param Uplo designs the triangular part that has to be updated. The
   *   possible value CblasUpper and CblasLower.
   * @param ldbk is the leading dimension of matrix Bk.
   * @param n is the problem dimension.
   * @param yk is the gradient difference (g+-gk).
   * @param incyk is the increment in vector yk.
   * @param sk is the current step: sk = x+-xk.
   * @param incsk is the increment in vector sk.
   * @param work is a working vector of size n.
   *
   * @return 1 is the updated has been performed, 0 else.
   */
  int bfgs_update(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
		  const int n, double *Bk, const int ldbk,
		  const double *yk, const int incyk,
		  const double *sk, const int incsk,
		  double *work);

  /**
   * Stopping test for BFGS linesearch method.
   *
   * @param bf is an initialized BFGS structure.
   */
  int nlp_bfgs_stopping_test(BFGS *bf);

  /**
   * BFGS method. Expert routine.
   * Function computing the minimum of a function of class C1
   *
   * @param bf is an initialized BFGS structure.
   * @param Order designs the orientation of the matrix Bk. The
   *   possible values are CblasRowMajor and CblasColMajor.
   * @param Uplo designs the triangular part that has to be updated. The
   *   possible value CblasUpper and CblasLower.
   * @param f is the function to minimize.
   * @param n is the primary problem dimension.
   * @param x is on input the starting point. It is the solution on output.
   * @param output is the output file for results and iteration informations
   * @param Bk is the initial secant matrix.
   * @param data represents the data needed by function f (NULL if no
   *   data is necessary).
   * @param is a working vector of size at least 3*n.
   *
   * @return the optimal value.
   */
  double nlp_bfgsx(BFGS *bf,
		   const enum CBLAS_ORDER Order,
		   const enum CBLAS_UPLO Uplo,
		   C_GENERIC f, int f_type,
		   double (*fg)(double *, int, void *, double *),
		   void *data, FILE *output,
		   double **Bk, int n, double *x,
		   double *work);

  /**
   * BFGS secant update.
   *
   * @param Order designs the orientation of the matrix Bk. The
   *   possible values are CblasRowMajor and CblasColMajor.
   * @param Uplo designs the triangular part that has to be updated. The
   *   possible value CblasUpper and CblasLower.
   * @param n is the primary problem dimension.
   * @param m is the secondary problem dimension.
   */
  int nlp_bfgs_update(const enum CBLAS_ORDER Order,
		      const enum CBLAS_UPLO Uplo,
		      int n, int m, double *H,
		      double *values, double *valueplus,
		      double *grad, double *gradplus, double *step,
		      int *info, double *work, void *data);

  /**
   * BFGS method. Simple routine.
   * Function computing the minimum of a function of class C1
   *
   * @param work is a double precision working vector of size at least 3*n.
   */
  double nlp_bfgs(BFGS *bf, C_GENERIC f, C_GENERIC g, void *data,
		  FILE *output, double **hessian, double *work);

#ifdef __cplusplus
}
#endif

#endif
