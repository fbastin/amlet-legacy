/**
 * sr1.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2003 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

/**
 * SR-1 approximation
 */

#include <sr1.h>

#include <math.h>

#include <oratio/ntools.h>

#define R_SR1 1e-8  /* tolerance used to reject too small denominators
		       Proposition of Nocedal and Wright: 1e-8 */

/**
 * SR1 Hessian Update
 *
 * @param Order designs the orientation of the matrix Bk. The
 *   possible values are CblasRowMajor and CblasColMajor.
 * @param Uplo designs the triangular part that has to be updated. The
 *   possible value CblasUpper and CblasLower.
 * @param n is the number of variables.
 * @param yk is the gradients difference.
 * @param sk is the step.
 *
 * @return 1 if the approximation Bk has been updated, 0 else.
 */
int sr1_update(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	       const int n, double *Bk, const int ldbk,
	       const double *yk, const int incyk,
	       const double *sk, const int incsk,
	       double *work)
{
  /**
   * SR1 Hessian Update
   * Ref.: Nocedal and Wright, Numerical Optimization
   *
   * @return 1 if the approximation Bk has been updated, 0 else.
   */

  /* Would be useful to compare this implementation to other ones */

  double dot;

  cblas_dcopy(n, yk, incyk, work, 1);

  /* Compute y_k-B_ks_k */
  cblas_dsymv(Order, Uplo, n, -1.0, Bk, ldbk, sk, incsk, 1.0, work, 1);

  /* Compute s_k^T(y_k-B_ks_k) */
  dot = cblas_ddot(n, sk, incsk, work, 1);

  /* Check the denominator */
  if (fabs(dot) <= R_SR1*cblas_dnrm2(n, sk, incsk)*cblas_dnrm2(n, work, 1))
    return 0;

  /* Compute B_{k+1} */
  cblas_dsyr(Order, Uplo, n, 1/dot, work, 1, Bk, ldbk);

  return 1;
}

int nlp_sr1_update(const enum CBLAS_ORDER Order,
		   const enum CBLAS_UPLO Uplo,
		   int n, int m, double *H,
		   double *values, double *valuseplus,
		   double *grad, double *gradplus, double *step,
		   int *info, double *work, void *data)
{
  cblas_dcopy(n, gradplus, 1, work, 1);
  cblas_daxpy(n, -1.0, grad, 1, work, 1);

  return sr1_update(Order, Uplo, n, H, n, work, 1, step, 1, work+n);
}
