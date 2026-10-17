/**
 * nlp_secant.c
 *
 * OPHELIA Library
 * Copyright (C) 2005-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "secant.h"

#include "bfgs.h"
#include "sr1.h"
#include "nlp_least_squares.h"

#include <oratio/ctools.h>
#include <oratio/ntools.h>

/**
 * General framework for secant approximations.
 */

void *nlp_set_hessian_approximation(int id, int info)
{
  switch (id) {
  case NLP_SECANT_BFGS:
    return (C_GENERIC)bfgs_update;
    break;
  case NLP_SECANT_SR1:
    return (C_GENERIC)sr1_update;
    break;
  case NLP_GAUSS_NEWTON:
    if (info) return (C_GENERIC)nlp_ls_gn_transpose;
    return (C_GENERIC)nlp_ls_gn;
    break;
  case NLP_SECANT_LS:
    return (C_GENERIC)nlp_ls_secant;
    break;
  default:
    return NULL;
  }
}

NlpSecant *nlp_secant_new_default(int n, int id)
{
  NlpSecant *ns;

  ns = _malloc(NlpSecant);
  ns->H = nt_matrix_new(n, n);
  ns->n = n;
  ns->id = id;
  ns->init = nlp_secant_init_default;
  ns->initialized = 0;
  ns->ns = NULL;

  return ns;
}

NlpSecant *nlp_secant_new(int n, int id)
{
  NlpSecant *ns;

  switch(id) {
  case NLP_SECANT_LS:
    ns = nlp_ls_secant_new(n, id);
    break;
  default:
    ns = nlp_secant_new_default(n, id);
    break;
  }

  return ns;
}

void nlp_secant_free(NlpSecant *ns)
{
  if (ns) {
    switch(ns->id) {
    case NLP_SECANT_LS:
      nlp_ls_secant_free(ns);
      break;
    default:
      nt_matrix_free(ns->H);
      free(ns);
    }
  }
}

/**
 * Set the update function.
 */
int nlp_secant_set_fct(NlpSecant *ns)
{
  switch (ns->id) {
  case NLP_SECANT_BFGS:
    ns->update = nlp_bfgs_update;
    break;
  case NLP_SECANT_SR1:
    ns->update = nlp_sr1_update;
    break;
  case NLP_GAUSS_NEWTON:
    ns->update = nlp_ls_gn;
    ns->init = nlp_ls_gn;
    break;
  case NLP_SECANT_LS:
    ns->init = nlp_ls_secant_init;
    ns->update = nlp_ls_secant;
    break;
  default:
    ns->update = NULL;
    return 0;
  }

  return 1;
}

/**
 * Default initialization for secant approximations.
 */
int nlp_secant_init_default(const enum CBLAS_ORDER Order,
			    const enum CBLAS_UPLO Uplo,
			    int n, int m, double *H,
			    double *values, double *valuesplus,
			    double *grad, double *gradplus,
			    double *step, int *info, double *work, void *data)
{
  double ng = cblas_dnrm2(n, grad, 1);

  if (fabs(ng) < DBL_EPSILON) ng = 1.0;

  nt_ma_diag(Order, n, n, H, n, 1.0/ng);

  return 0;
}

void nlp_secant_init(NlpSecant *ns, const enum CBLAS_ORDER Order,
		     const enum CBLAS_UPLO Uplo,
		     int n, int m, double *grad, double *gradplus,
		     double *step, int *info, double *work, void *data)
{
  if (!ns->initialized) {
    ns->initialized = 1;
    nlp_secant_init_default(Order, Uplo, n, m, *(ns->H),
			    NULL, NULL, grad, gradplus,
			    step, info, work, data);
  }
}

/**
 * Define the Hessian matrix associated to the particular secant method.
 *
 * @return a pointer to the Hessian matrix corresponding to the
 *   secant method defined in ns.
   */
double **nlp_secant_get_hessian(NlpSecant *ns, const enum CBLAS_ORDER Order,
				const enum CBLAS_UPLO Uplo, int n, int m)
{
  switch (ns->id) {
  case NLP_SECANT_LS:
    return nlp_ls_secant_get_hessian(ns, Order, Uplo, n, m);
    break;
  default:
    return ns->H;
  }
}
