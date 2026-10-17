/**
 * step.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 1999-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "step.h"

#include "bounds.h"
#include "model.h"
#include "pcg.h"
#include "projection.h"

#include <math.h>
#include <string.h>

#include <oratio/ctools.h>
#include <oratio/ntools.h>

/**
 * Brief descriptions of the various step techniques.
 */
c_symbol nlp_step[] = {
  {"", PERSONAL},
  {"Cauchy step", BTR_STEP_CAUCHY},
  {"Steihaug-Toint step", BTR_STEP_STEIHAUG_TOINT},
  {"Restricted Steihaug-Toint step", BTR_STEP_RESTRICTED_STEIHAUG_TOINT},
  {"projected Cauchy step", BTR_STEP_PROJ_CAUCHY},
  {"simple projected conjugated gradient", BTR_STEP_SPCG_ST},
};

/**
 * Return a short description of the step method used in the
 * trust-region algorithm.
 *
 * @param btr is a initialized BTR structure.
 *
 * @return the name of the step technique used in btr.
 */
char *btr_get_step_name(BTR *btr)
{
  return c_symbol_get_name(btr->step_type, nlp_step, _sizeenum(nlp_step));
}

/**
 * Truncated conjugated gradient algorithm implementation
 *
 * Ref.: Andrew R. Conn, Nicholas I. M. Gould, Philippe L. Toint,
 *   "Trust-Region Methods", SIAM, pp 202-207, 2000.
 *
 * @param btr is an initialized BTR structure.
 */
void btr_step_steihaug_toint(BTR *btr,
			     const enum CBLAS_ORDER Order,
			     const enum CBLAS_UPLO Uplo,
			     double **hessian, double *work,
			     void *data)
{

  /* We currently only consider the 2-norm */
  double *g, *v, *p, *sum, *Hp;
  double alpha, beta, kappa, root, sigma, normg0;

  double sqrnrm, tmp;

  int k;   /* iteration counter */
  size_t size;

  int n;

  n = btr->n;

  /* allocation de memoire */
  g = work;
  v = work+n;
  p = v+n;
  sum = p+n;
  Hp = _mallocs(n, double);

  /* initialisation des variables */
  size = n*sizeof(double);
  memcpy(g, btr->grad, size);
  memcpy(v, btr->grad, size);
  memcpy(p, btr->grad, size);

  for (k = 0; k < n; k++) {
    btr->step[k] = 0;
    p[k] = -p[k];
  }
  normg0 = cblas_dnrm2(n, g, 1);
  k = 1;

  while (cblas_dnrm2(n, g, 1) > normg0*nt_min(0.01, sqrt(normg0)) && k <= n) {
    cblas_dsymv(Order, Uplo, n, 1.0, *hessian, n, p, 1, 0.0, Hp, 1);
    kappa = cblas_ddot(n, p, 1, Hp, 1);
    if (kappa <= 0) {
      sigma = cblas_ddot(n, btr->step, 1, p, 1);
      sqrnrm = cblas_ddot(n, p, 1, p, 1);
      root = sqrt(sigma*sigma+sqrnrm*
		  (btr->delta*btr->delta - 
		   cblas_ddot(n, btr->step, 1, btr->step, 1)));
      sigma = (root-sigma)/sqrnrm;
      cblas_daxpy(n, sigma, p, 1, btr->step, 1);
      k = n+1;
    }
    else {
      alpha = cblas_ddot(n, g, 1, v, 1)/kappa;
      cblas_dcopy(n, btr->step, 1, sum, 1);
      cblas_daxpy(n, alpha, p, 1, sum, 1);
      if (cblas_ddot(n, sum, 1, sum, 1) >= btr->delta*btr->delta) {
	/* we are leaving the trust region */
	sigma = cblas_ddot(n, btr->step, 1, p, 1);
	sqrnrm = cblas_ddot(n, p, 1, p, 1);
	root = sqrt(sigma*sigma+sqrnrm*
		    (btr->delta*btr->delta - 
		     cblas_ddot(n, btr->step, 1, btr->step, 1)));
	sigma = (root-sigma)/sqrnrm;
	cblas_daxpy(n, sigma, p, 1, btr->step, 1);
	k = n+1; /* go out */
      }
      else {
	memcpy(btr->step, sum, size);
	tmp = cblas_ddot(n, g, 1, v, 1);
	/* Compute g_k+1 */
	/* IT SHOULD BE POSSIBLE TO SIMPLIFY THIS LINE */
	//	cblas_dsymv(Order, Uplo, n, alpha, *hessian, n, p, 1, 1.0, g, 1);
	cblas_daxpy(n, alpha, Hp, 1, g, 1);
	memcpy(v, g, size);
	beta = cblas_ddot(n, g, 1, v, 1) / tmp;
	cblas_dscal(n, beta, p, 1);
	cblas_daxpy(n, -1, v, 1, p, 1);
	/*	printf("Gradient norm: %g\n", cblas_dnrm2(n, g, 1)); */
      }
    }
    k++;
  }

  free(Hp);
}

/**
 * Bounds constraints
 */

double nlp_bounds_maximal_step_size(double *bounds)
{
  return 1.0;
}

/**
 * Allow to attach a predefined step computation function.
 *
 * @param btr is a initialized BTR structure.
 * @param step_type is a step identifier constant.
 */
void btr_set_step(BTR *btr, const int step_type)
{
  btr->step_type = step_type;

  switch (step_type) {
  case BTR_STEP_CAUCHY:
    btr->step_function = btr_step_cauchy;
    break;
  case BTR_STEP_STEIHAUG_TOINT:
    btr->step_function = btr_step_steihaug_toint;
    break;
  case BTR_STEP_RESTRICTED_STEIHAUG_TOINT:
    c_error("%s not yet implemented",
	    nlp_step[BTR_STEP_RESTRICTED_STEIHAUG_TOINT].name);
    break;
  case BTR_STEP_PROJ_CAUCHY:
    btr->step_function = btr_proj_cauchy;
    break;
  case BTR_STEP_SPCG_ST:
    btr->step_function = nlp_spcg_st;
    break;
  }
}
