/**
 * nlp_least_squares.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006 Fabian Bastin <bastin@cerfacs.fr>
 *
 * Licensed under the MIT License
 */

#include "nlp_least_squares.h"

#include <oratio/ntools.h>

#include <stdio.h>
#include <string.h>

/**
 * Compute the least-squares gradient from the Jacobian components.
 *
 * @param m is the number of terms in the least-squares sum.
 * @param n is the problem dimension.
 * @param r is the residuals vector.
 * @param transpose is equal to 1 if the gradient has to be transposed, 0 else.
 */
void nlp_ls_grad(const int m, const int n, double *J, int ldJ,
			double *r, int incr, int transpose)
{
  int i;

  if (transpose)
    for (i = 0; i < n; i++) J[i] = cblas_ddot(m, r, incr, J+n+i*ldJ, 1);
  else
    for (i = 0; i < n; i++) J[i] = cblas_ddot(m, r, incr, J+ldJ+i, ldJ);
}

int nlp_ls_gn(const enum CBLAS_ORDER Order,
	      const enum CBLAS_UPLO Uplo,
	      int n, int m, double *hessian,
	      double *values, double *valueplus,
	      double *grad, double *gradplus,
	      double *step,
	      int *info, double *work, void *data)
{
  int i;
  double *g;

  g = (gradplus) ? gradplus : grad;

  nt_ma_zeros(Order, n, n, hessian, n);

  /* How to take account of the transposition??? */
  for (i = 1; i <= m; i++)
    cblas_dsyr(Order, Uplo, n, 1.0, g+i*n, 1, hessian, n);
  /*
    USE info>
  if (transpose) {
    for (i = 0; i < m; i++)
      cblas_dsyr(Order, Uplo, n, 1.0, g+n+i, m, hessian, n);
  }
  */

  return 0;
}

int nlp_ls_gn_transpose(const enum CBLAS_ORDER Order,
			const enum CBLAS_UPLO Uplo,
			int n, int m, double *hessian,
			double *values, double *valueplus,
			double *grad, double *gradplus,
			double *step,
			int *info, double *work, void *data)
{
  int i;
  double *g;

  g = (gradplus) ? gradplus : grad;

  nt_ma_zeros(Order, n, n, hessian, n);

  /* How to take account of the transposition??? */
  for (i = 0; i < m; i++)
    cblas_dsyr(Order, Uplo, n, 1.0, g+n+i, m, hessian, n);

  return 0;
}

/**
 * Least-squares secant update.
 *
 * @param Order designs the orientation of the matrix Bk. The
 *   possible values are CblasRowMajor and CblasColMajor.
 * @param Uplo designs the triangular part that has to be updated. The
 *   possible value CblasUpper and CblasLower.
 * @param n is the number of variables.
 *
 * @return 1 if the approximation Bk has been updated, 0 else.
 */
int nlp_ls_secant(const enum CBLAS_ORDER Order,
		  const enum CBLAS_UPLO Uplo,
		  int n, int m, double *Ak,
		  double *values, double *valuesplus,
		  double *grad, double *gradplus,
		  double *step,
		  int *info, double *work, void *data)
{
  /**
   * Ref. Dennis J.E. and Schnabel R.B., "Numerical methods for unconstrained
   * optimizationand nonlinear equations", SIAM, 1983, Chapter 10.
   */

  double snorm, ynorm, ys, scale, tmp;

  /* y, stored in work, represent the gradients difference */
  cblas_dcopy(n, gradplus, 1, work, 1);
  cblas_daxpy(n, -1.0, grad, 1, work, 1);  

  ys = cblas_ddot(n, step, 1, work, 1);
  snorm = cblas_dnrm2(n, step, 1);
  ynorm = cblas_dnrm2(n, work, 1);

  if (ys > sqrt(DBL_EPSILON)*snorm*ynorm) {
    /* Scaling factor */
    cblas_dsymv(Order, Uplo, n, 1.0, Ak, n, step, 1, 0.0, work+2*n, 1);
    scale = nt_min(1.0, ys/cblas_ddot(n, step, 1, work+2*n, 1));

    ys = 1.0/ys;

    /* Compute y#, stored in work+n */
    cblas_dcopy(n, gradplus, 1, work+n, 1);
    cblas_dgemv(CblasRowMajor, CblasTrans, m, n, -1.0, grad+n, n,
		valuesplus, 1, 1.0, work+n, 1);
    /* Compute y#-Ac*sc, stored in work+n */
    cblas_daxpy(n, -scale, work+2*n, 1, work+n, 1);

    /* To revise */
    if (scale < 1.0) nt_ma_scale(Order, n, n, scale, Ak, n);

    /* Ac := Ac + 1/yc'sc*((y#-Ac*sc)*yc'+yc*(y#-Ac*sc)') */
    cblas_dsyr2(Order, Uplo, n, ys, work, 1, work+n, 1, Ak, n);

    /* Ac := Ac -<y#-Ac*sc,sc>/(yc'sc)^2*yc*yc'*/
    tmp = -cblas_ddot(n, work+n, 1, step, 1)*ys*ys;
    cblas_dsyr(Order, Uplo, n, tmp, work, 1, Ak, n);
  }

  return 0;
}

int nlp_ls_secant_init(const enum CBLAS_ORDER Order,
		       const enum CBLAS_UPLO Uplo,
		       int n, int m, double *hessian,
		       double *values, double *valuesplus,
		       double *grad, double *gradplus,
		       double *step,
		       int *info, double *work, void *data)
{
  nt_ma_zeros(Order, n, n, hessian, n);

  return 0;
}

double **nlp_ls_secant_get_hessian(NlpSecant *ns, const enum CBLAS_ORDER Order,
				   const enum CBLAS_UPLO Uplo, int n, int m)
{
  int i;

  // memcpy(ns->data, ns->ns->H, n*n*sizeof(double));
  memcpy(*((double **)(ns->data)), *(ns->ns->H), n*n*sizeof(double));

  /* Add ns->H to ns->data and store the result in ns->data. */
  /* We perform the addition row-by-row */
  for (i = 0; i < n; i++) {
    /*    cblas_daxpy(n, 1.0, ns->H[i], 1,
	  ((double **)(ns->data))[i], 1); */
    cblas_daxpy(n, 1.0, ns->H[i], 1,
		((double **)(ns->data))[i], 1);
  }

  return (double **)(ns->data);
}

NlpSecant *nlp_ls_secant_new(int n, int id)
{
  NlpSecant *ns;

  ns = nlp_secant_new_default(n, id);
  ns->init = nlp_ls_secant_init;
  ns->ns = NULL;  /* ns should be set to the Gauss-newton update later */
  // ns->data = _mallocs(n*n, double);
  ns->data = nt_matrix_new(n, n);

  return ns;
}

void nlp_ls_secant_free(NlpSecant *ns)
{
  if (ns) {
    nt_matrix_free(ns->H);
    nt_matrix_free((double **)(ns->data));
    free(ns);
  }
}
