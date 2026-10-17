/**
 * wrapper.h
 *
 * OPHELIA Library
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ATHOS_WRAPPER_H__
#define __ATHOS_WRAPPER_H__

#include <oratio/ntools.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef double (*nlp_wrapper_fct)(C_GENERIC f, int f_type,
				    double *x, int n, int m, void *data,
				    double *values, double *grad, double **H,
				    void *prob_info, int type_prob, int *info);

  typedef void (*nlp_wrapper_grad)(C_GENERIC f, int f_type,
				   C_GENERIC g, int g_type,
				   double *x , int n, void *data,
				   double val, double *grad,
				   double *scaling, void *prob_info,
				   int prob_type, int *info,
				   double *h, double *work);

  typedef void (*nlp_wrapper_hess)(const enum CBLAS_ORDER Order,
				   const enum CBLAS_UPLO Uplo,
				   C_GENERIC f, int f_type,
				   C_GENERIC g, int g_type,
				   C_GENERIC H, int h_type,
				   double *x , int n, void *data,
				   double *grad, double **hessian,
				   double *scaling, void *prob_info,
				   int prob_type, int *info, double *work);
  /**
   * BTR wrapper for the function evaluation.
   *
   * @param f is the objective function to evaluate; f is assumed to be
   *   twice-continuously derivate in a neighbourhood of x.
   * @param x is the point where the function has to be evaluated.
   * @param n is the primary problem dimension.
   */
  double nlp_wrapper_fct_dflt(C_GENERIC f, int f_type,
			      double *x, int n, int m, void *data,
			      double *values, double *grad, double **H,
			      void *prob_info, int type_prob, int *info);

  /**
   * BTR wrapper for the gradient evaluation.
   *
   * @param f is the objective function to evaluate; f is assumed to be
   *   twice-continuously derivate in a neighbourhood of x.
   * @param n is the primary problem dimension.
   */
  void nlp_wrapper_grad_dflt(C_GENERIC f, int f_type,
			     C_GENERIC g, int g_type,
			     double *x , int n, void *data, double val,
			     double *grad, double *scaling,
			     void *prob_info, int prob_type, int *info,
			     double *h, double *work);
  /**
   * BTR wrapper for the Hessian evaluation.
   *
   * @param Order designs the orientation of the matrix Bk. The
   *   possible values are CblasRowMajor and CblasColMajor.
   * @param Uplo designs the triangular part that has to be updated. The
   *   possible value CblasUpper and CblasLower.
   * @param f is the objective function to evaluate; f is assumed to be
   *   twice-continuously derivate in a neighbourhood of x.
   * @param x is the point where the Hessian is evaluated.
   * @param n is the primary problem dimension.
   */
  void nlp_wrapper_hess_dflt(const enum CBLAS_ORDER Order,
			     const enum CBLAS_UPLO Uplo,
			     C_GENERIC f, int f_type,
			     C_GENERIC g, int g_type,
			     C_GENERIC H, int h_type,
			     double *x , int n, void *data,
			     double *grad, double **hessian, double *scaling,
			     void *prob_info, int prob_type, int *info,
			     double *work);

#ifdef __cplusplus
}
#endif

#endif  /* __ATHOS_WRAPPER_H__ */
