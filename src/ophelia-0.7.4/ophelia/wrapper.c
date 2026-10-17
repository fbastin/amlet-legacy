/**
 * wrapper.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006 Fabian Bastin <bastin@cerfacs.fr>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "wrapper.h"

#include <float.h>

#include "btr.h"
#include "stopping.h"

/**
 * BTR wrapper for the function evaluation.
 *
 * @param f is the objective function to evaluate; f is assumed to be
 *   twice-continuously derivate in a neighbourhood of x.
 * @param x is the point where the function has to be evaluated.
 */
double nlp_wrapper_fct_dflt(C_GENERIC f, int f_type,
				   double *x, int n, int m, void *data,
				   double *values, double *grad, double **H,
				   void *prob_info, int type_prob, int *info)
{
  double val = DBL_MAX;
  /* This function should be modulated */
  /* The type should indicated if we evaluate the function AND the
     gradient or only the function */

  switch (f_type) {
  case NlpF:
    val = ((NT_f)(f))(x, n, data);
    break;
  case NlpFval:
    ((NT_fval)(f))(x, n, data, &val);
    break;
  case NlpFG:
    val = ((NT_fg)(f))(x, n, data, grad);
    break;
  case NlpFGH:
    break;
  case NlpMFG:
    val = ((NT_mfJ)(f))(x, n, data, values, m, grad);
    break;
  default:
    *info = NLP_STOPPING_PARAMETERS;
  }

  return val;
}

/**
 * BTR wrapper for the gradient evaluation.
 *
 * @param f is the objective function to evaluate; f is assumed to be
 *   twice-continuously derivate in a neighbourhood of x.
 */
void nlp_wrapper_grad_dflt(C_GENERIC f, int f_type,
				  C_GENERIC g, int g_type,
				  double *x , int n, void *data, double val,
				  double *grad, double *scaling,
				  void *prob_info, int prob_type, int *info,
				  double *h, double *work)
{
  switch (g_type) {
  case NlpFval:
    ((NT_fval)(g))(x, n, data, grad);
    break;
  case NlpFG:
    val = ((NT_fg)(g))(x, n, data, grad);
    break;
  default:
    nt_derive_grad_fd((NT_f)(f), x, grad, h, 0.0,
		      scaling, n, NULL, work, data);
    break;
  }
}

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
 */
void nlp_wrapper_hess_dflt(const enum CBLAS_ORDER Order,
				  const enum CBLAS_UPLO Uplo,
				  C_GENERIC f, int f_type,
				  C_GENERIC g, int g_type,
				  C_GENERIC H, int h_type,
				  double *x , int n, void *data,
				  double *grad, double **hessian,
				  double *scaling, void *prob_info,
				  int prob_type, int *info, double *work)
{
  if (h_type == NlpExact)
    ((NT_fval)(H))(x, n, data, *hessian);
  else {
    switch (f_type) {
    case NlpFG:
      nt_derive_hess_fd((NT_fg)(f), x, grad, hessian,
			work, 0.0, scaling, n, NULL, NULL, data);
      break;
    default:
      nt_derive_hess_cd((NT_f)(f), x, hessian, work,
			0.0, scaling, n, NULL, NULL, data);
      break;
    }
  }
}
