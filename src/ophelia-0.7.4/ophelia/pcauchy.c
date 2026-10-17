/**
 * pcauchy.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "pcauchy.h"

#include "bounds.h"
#include "model.h"

#include <oratio/ntools.h>

/**
 * Projected Cauchy step
 */

double btr_proj_cauchy_step(double t, double *step, double *grad,
				   double **hessian, BTR *btr,
				   const enum CBLAS_ORDER Order,
				   const enum CBLAS_UPLO Uplo,
				   void *data)
{
  int n;     /* problem dimension */
  double m = 0.0;  /* model value at the projected point: NOT COMPUTED */

  n = btr->n;

  /* compute the point to project: x-t*grad */
  cblas_dcopy(n, grad, 1, step, 1);
  cblas_dscal(n, -t, step, 1);
  cblas_daxpy(n, 1.0, btr->x, 1, step, 1);

  /* compute point on the projected gradient path */
  btr->projection_wrapper(btr, step);

  /* step computation: step := step-x */
  cblas_daxpy(n, -1.0, btr->x, 1, step, 1);

  return -m;
}

void nlp_proj_cauchy_direction(double t, double *step, double *grad,
				      double **hessian, BTR *btr,
				      const enum CBLAS_ORDER Order,
				      const enum CBLAS_UPLO Uplo,
				      void *data)
{
  int n;     /* problem dimension */

  /* Difference with nlp_proj_cauchy_step? */

  n = btr->n;

  /* compute the point to project */
  cblas_dcopy(n, grad, 1, step, 1);
  cblas_dscal(n, -t, step, 1);
  cblas_daxpy(n, 1.0, btr->x, 1, step, 1);

  /* compute point on the projected gradient path */
  btr->projection_wrapper(btr, step);

  /* step computation: step := step-x */
  cblas_daxpy(n, -1.0, btr->x, 1, step, 1);
}

/**
 * Ref.: Lin and More, "Newton's Method For Large Bound-Constrained
 * Optimization Problems", Siam Journal on Optimization 9(4), 1999.
 */
void nlp_proj_cauchy_bounds_st(BTR *b, double *alpha,
		     const enum CBLAS_ORDER Order,
		     const enum CBLAS_UPLO Uplo,
		     double **hessian, void *data)
{
  int n;             /* problem dimension */
  double *p;
  /* minimal and maximal breakpoints */
  double brptmin, brptmax, gts, m;
  double deltasqr;   /* square of the trust-region radius */
  int interp;        /* 1 for interpolation, 0 for extrapolation */
  int search;        /* 1 while search has to be conducted, 0 when it
			is completed. */

  NlpBounds *bounds = (NlpBounds *)b->projection_data;

  /* This is somewhat crude, isn't it? */
  double mu0 = 0.01;
  double interpf = 0.1;
  double extrapf = 10;
  double alphas;

  n = b->n;
  p = b->work;
  deltasqr = b->delta*b->delta;

  /* Compute the initial point to project */
  cblas_dcopy(n, b->grad, 1, p, 1);
  cblas_dscal(n, -1.0, p, 1);

  /* Find the minimal and maximal break-point on x-alpha*g */
  nlp_bounds_breakpt_st(n, b->x, bounds->xl, bounds->xu, p,
			&brptmin, &brptmax);

  /* Evaluate the initial alpha amd decide if the algorithm must
     interpolate or extrapolate. */
  nlp_bounds_gpstep_st(n, b->x, bounds->xl, bounds->xu,
		       -*alpha, b->grad, b->step);

  if (cblas_ddot(n, b->step, 1, b->step, 1) > deltasqr)
    /* We are outside the trust-region. */
    interp = 1;
  else {
    /* Compute the quadratic model. */
    cblas_dsymv(Order, Uplo, n, 1.0, *hessian, n, b->step, 1, 0.0, p, 1);
    gts = cblas_ddot(n, b->grad, 1, b->step, 1);
    m = 0.5*cblas_ddot(n, b->step, 1, p, 1)+gts;
    interp = (m >= mu0*gts);
  }

  /* Either interpolate or extrapolate to find a successful step. */

  if (interp) {
    /* Reduce alpha until a successful step is found. */
    search = 1;

    while (search) {
      /* Tron comment: this is a crude interpolation procedure, to
	 improve. */

      *alpha *= interpf;
      nlp_bounds_gpstep_st(n, b->x, bounds->xl, bounds->xu,
			   -*alpha, b->grad, b->step);
      if (cblas_ddot(n, b->step, 1, b->step, 1) <= deltasqr) {
	cblas_dsymv(Order, Uplo, n, 1.0, *hessian, n, b->step, 1, 0.0, p, 1);
	gts = cblas_ddot(n, b->grad, 1, b->step, 1);
	m = 0.5*cblas_ddot(n, b->step, 1, p, 1)+gts;
	search = (m > mu0*gts);
      }
    }
  }
  else {
    /* Increase alpha until a successful step is found. */
    search = 1;
    alphas = *alpha;

    while (search) {
      /* Tron comment: this is a crude extrapolation procedure, to
	 improve. */
      *alpha *= extrapf;
      nlp_bounds_gpstep_st(n, b->x, bounds->xl, bounds->xu,
			   -*alpha, b->grad, b->step);
      if (cblas_ddot(n, b->step, 1, b->step, 1) > deltasqr) {
	cblas_dsymv(Order, Uplo, n, 1.0, *hessian, n, b->step, 1, 0.0, p, 1);
	gts = cblas_ddot(n, b->grad, 1, b->step, 1);
	m = 0.5*cblas_ddot(n, b->step, 1, p, 1)+gts;
	if (m < mu0*gts) {
	  search = 1;
	  alphas = *alpha;
	}
      }
      else
	search = 0;
    }

    /* Recover the last successful step.*/
    *alpha = alphas;
    nlp_bounds_gpstep_st(n, b->x, bounds->xl, bounds->xu,
			 -*alpha, b->grad, b->step);
  }
}

/**
 * Projected Cauchy step
 *
 * @param btr is am initialized BTR structure.
 * @param Order designs the orientation of the matrix Bk. The
 *   possible values are CblasRowMajor and CblasColMajor.
 * @param Uplo designs the triangular part that has to be updated. The
 *   possible value CblasUpper and CblasLower.
 */
void btr_proj_cauchy(BTR *btr,
		     const enum CBLAS_ORDER Order,
		     const enum CBLAS_UPLO Uplo,
		     double **hessian, double *work,
		     void *data)
{
  double k_ubs, k_lbs, k_frd, k_frd2;
  double k_tfr, k_lbt; /* constantes sur t */
  double t, l, u;
  double *s_l, *s_u;
  int s_0_l, s_0_u;
  double m, delta, delta2;
  double step_norm2;

  double pr;  /* used in dot products */

  /* arret_proj_cauchy will be set to 1 if the stopping criterium is
     satisfied, 0 else */
  int arret_proj_cauchy, ok;

  int n; /* problem dimension */
  int i;

  /* Initialisation */
  arret_proj_cauchy = 0;

  n = btr->n;
  delta = btr->delta;
  delta2 = delta*delta;

  k_lbs = k_frd = 0.85;
  k_tfr = 0.95;
  k_lbt = 0.99;
  k_ubs = 0.25;
  k_frd2 = k_frd*k_frd;
  l = 0;
  u = t = delta;

  /* Pas initial */
  s_l = _mallocs(n, double);
  s_u = _mallocs(n, double);
  s_0_l = s_0_u = 1;

  for (i = 0; i < n; i++) s_l[i] = 0.0;

  cblas_dcopy(n, btr->grad, 1, s_u, 1);
  m = cblas_dnrm2(n, btr->grad, 1);
  for (i = 0; i < n; i++) s_u[i] *= -(delta+1)/m;

  do {
    do {
      nlp_proj_cauchy_direction(t, btr->step, btr->grad, hessian, btr,
				Order, Uplo, data);
      /* Evaluate the model at the projected point */
      m = -model(Order, Uplo, btr->step, btr->grad, *hessian, n, data);
      
      step_norm2 = cblas_ddot(n, btr->step, 1, btr->step, 1);
      if (step_norm2 > delta2) {
	u = t;
	t /= 2;
      }

    } while (step_norm2 > delta2);

    /* Stopping conditions */
    pr = cblas_ddot(n, btr->grad, 1, btr->step, 1);
    if (m > k_ubs*pr) {
      u = t;
      s_0_u = 0;
    }
    else if ( (m < k_lbs*pr) &&
	      (t < nt_min(k_tfr*delta,k_lbt)) &&
	      (step_norm2 < k_frd*k_frd*delta2) ) {
      l = t;
      s_0_l = 0;
    }
    else
      arret_proj_cauchy = 1;

    /* Compute a new value for the parameter */
    if (!arret_proj_cauchy) {
      /* Cas 1 */
      t = 0.5*(l+u);
      if (s_0_l && s_0_u) {
	/* Cas 2 */
	ok = 0;
	do {
	  t = 0.5*(l+u);
	  nlp_proj_cauchy_direction(t, btr->step, btr->grad, hessian,
				    btr, Order, Uplo, data);
	  /* Evaluate the model at the projected point */
	  m = model(Order, Uplo, btr->step, btr->grad, *hessian, n, data);
	  pr = cblas_ddot(n, btr->grad, 1, btr->step, 1);
	  if (m > k_ubs*pr)
	    u = t;
	  else if ( (m < k_lbs*pr) &&
		    (cblas_ddot(n, btr->step, 1, btr->step, 1) <
		     k_frd*k_frd*delta*delta) &&
		    (t < nt_min(k_tfr*delta,k_lbt)) )
	    l = t;
	  else
	    ok = 1;
	}
	while (!ok);
      }
    }

  } while (!arret_proj_cauchy);

  free(s_l);
  free(s_u);
}

