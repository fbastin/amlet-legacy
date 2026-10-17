/**
 * nlp_least_squares.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/cblas.h>

#include "secant.h"

#ifndef __NLP2_LEAST_SQUARES_H__
#define __NLP2_LEAST_SQUARES_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Compute the least-squares gradient from the Jacobian components.
   *
   * @param m is the number of terms in the least-squares sum.
   * @param n is the problem dimension.
   * @param J is the Jacobian matrix.
   * @param r is the residuals vector.
   * @param transpose is equal to 1 if the gradient has to be transposed,
   *   0 else.
   */
  void nlp_ls_grad(const int m, const int n, double *J, int ldJ,
		   double *r, int incr, int transpose);

  double **nlp_ls_secant_get_hessian(NlpSecant *ns,
				     const enum CBLAS_ORDER Order,
				     const enum CBLAS_UPLO Uplo, int n, int m);

  /**
   * Gauss-Newton Hessian approximation.
   *
   * @param Order designs the orientation of the matrix Bk. The
   *   possible values are CblasRowMajor and CblasColMajor.
   * @param Uplo designs the triangular part that has to be updated. The
   *   possible value CblasUpper and CblasLower.
   * @param n is the number of variables.
   */
  int nlp_ls_gn(const enum CBLAS_ORDER Order,
		const enum CBLAS_UPLO Uplo,
		int n, int m, double *hessian,
		double *values, double *valuesplus,
		double *grad, double *gradplus,
		double *step,
		int *info, double *work, void *data);

  int nlp_ls_gn_transpose(const enum CBLAS_ORDER Order,
			  const enum CBLAS_UPLO Uplo,
			  int n, int m, double *hessian,
			  double *values, double *valueplus,
			  double *grad, double *gradplus,
			  double *step,
			  int *info, double *work, void *data);

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
		    int *info, double *work, void *data);


  /**
   *
   * @param Order designs the orientation of the matrix Bk. The
   *   possible values are CblasRowMajor and CblasColMajor.
   * @param Uplo designs the triangular part that has to be updated. The
   *   possible value CblasUpper and CblasLower.
   * @param n is the number of variables.
   */
  int nlp_ls_secant_init(const enum CBLAS_ORDER Order,
			 const enum CBLAS_UPLO Uplo,
			 int n, int m, double *hessian,
			 double *values, double *valuesplus,
			 double *grad, double *gradplus,
			 double *step,
			 int *info, double *work, void *data);

  /**
   * Define a new secant structure, specialized to nonlinear
   * least-squares problems.
   *
   * @param n is the number of variables.
   */
  NlpSecant *nlp_ls_secant_new(int n, int id);

  void nlp_ls_secant_free(NlpSecant *ns);

#ifdef __cplusplus
}
#endif

#endif  /* __NLP2_LEAST_SQUARES_H__ */
