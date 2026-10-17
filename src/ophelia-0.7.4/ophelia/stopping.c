/**
 * stopping.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 1999-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "stopping.h"

#include "btr.h"

#include <math.h>

#include <oratio/ntools.h>

static const char *nlp_stopping_text[] = {
  "stopping criterion not fulfilled",
  /* NLP_STOPPING_GRADIENT */
  "gradient norm",
  /* NLP_STOPPING_PROJECTED_GRADIENT_PATH */
  "projected gradient norm",
  /* NLP_STOPPING_RELATIVE_GRADIENT */
  "relative gradient norm",
  /* NLP_STOPPING_PROJECTED_RELATIVE_GRADIENT */
  "projected relative gradient norm",
  /* NLP_STOPPING_DECREASE */
  "too small decrease",
  /* NLP_STOPPING_STEP */
  "step length",
  /* NLP_STOPPING_DIST */
  "distance to the feasible set",
  /* NLP_STOPPING_UNDEFINED_ERROR */
  "undefined error",
  /* NLP_STOPPING_DELTA */
  "trust-region radius",
  /* NLP_STOPPING_ITERATIONS */
  "number of iterations",
  /* NLP_STOPPING_FAILURE_STEP */
  "failure in computation of a suitable step",
  /* NLP_STOPPING_INFEASIBILITY */
  "infeasible point",
  /* NLP_STOPPING_ROUNDING_ERRORS */
  "rounding errors",
  /* NLP_STOPPING_PARAMETERS */
  "invalide parameters",
  /* NLP_STOPPING_USER */
  "user origin"
};

int btr_successful(BTR *btr, double rho)
{
  return (rho >= btr->eta1) ? 1 : 0;
}

int btr_very_successful(BTR *btr, double rho)
{
  return (rho >= btr->eta2) ? 1 : 0;
}

const char *nlp_stopping_msg(int id)
{
  return nlp_stopping_text[id];
}

/**
 * Convergence indicator.
 *
 * @return 1 is the convergence (usually to a first-order critical
 * point) has been observed, 0 else.
 */
int nlp_has_converged(int id)
{
  if (id >= NLP_STOPPING_GRADIENT &&
      id <= NLP_STOPPING_PROJECTED_RELATIVE_GRADIENT)
    /* Convergence has been observed */
    return 1;
  else /* if (id >= NLP_STOPPING_DECREASE && NLP_STOPPING_PARAMETERS) */
    return 0;
}

/**
 * Compute the relative gradient
 *
 * Ref: Dennis, J. E. and Schnabel, R. B. (1983) Numerical Methods
 * for Unconstrained Optimization and Nonlinear Equations,
 * Prentice-Hall, Englewood Cliffs, NJ, page 160.
 *
 * @param val designs the objective function value
 * @param x is the point where the relative gradient has to be evaluated
 * @param grad is the gradient at x
 * @param typx is a typical value of x, or NULL if it is not available
 * @param n is the dimension of x and grad
 * @param typf is a typical value of the objective function
 *
 * @return the relative gradient value
 */
double nlp_relative_gradient(double val, double *x, double *grad,
				    double *typx, int n, double typf)
{
  int i;
  double gmax = 0.0;
  double typxi;

  for (i = 0; i < n; i++) {
    typxi = (typx) ? typx[i] : 1.0;
    gmax = nt_max(gmax,
		fabs(grad[i]*
		     nt_max(fabs(x[i]), typxi))/nt_max(fabs(val), typf));
  }

  return gmax;
}

/**
 * Compute the relative change (infinite norm) in x
 *
 * Ref: Dennis, J. E. and Schnabel, R. B. (1983) Numerical Methods
 * for Unconstrained Optimization and Nonlinear Equations,
 * Prentice-Hall, Englewood Cliffs, NJ, page 160.
 *
 * @param x is the point where the relative gradient has to be evaluated
 * @param typx is a typical value of x, or NULL if it is not available
 * @param n is the dimension of x and grad
 */
double nlp_relative_change(double *x, double *xplus,
				  double *typx, int n)
{
  /* THIS FUNCTION SHOULD BE CHECKED! */

  int i;
  double relx = 0.0;
  double typxi;

  for (i = 0; i < n; i++) {
    typxi = (typx) ? typx[i] : 1.0;
    relx = nt_max(relx, fabs(xplus[i]-x[i])/nt_max(fabs(xplus[i]), typxi));
  }

  return relx;
}

/**
 * Generic tests
 */

int nlp_stopping_gradient_2norm_test(double *grad,
					       double n, double tol)
{
  double norm;

  norm = cblas_ddot(n, grad, 1, grad, 1);

  return (norm < tol*tol) ? NLP_STOPPING_GRADIENT : NLP_STOPPING_NOT;
}

int nlp_stopping_radius_test(double radius, double tol)
{
  return (radius < tol) ? NLP_STOPPING_DELTA : NLP_STOPPING_NOT;
}

/**
 * BFGS tests
 */

/**
 * Stopping test for BFGS linesearch method.
 *
 * @param bf is an initialized BFGS structure.
 */
int nlp_bfgs_stopping_test(BFGS *bf)
{
  double relative_gradient;

  if (bf->info > 1) {
    switch (bf->info) {
    case 2:
      printf("Stop due to too small interval for the linesearch.\n");
      break;
    case 3:
      printf("Stop due to too large step.\n");
      break;
    case 4:
      printf("Stop due to too small step.\n");
      break;
    case 6:
      printf("Rounding errors prevent further progress.\n");
      break;
    default:
      printf("Failure of the linesearch (error %d\n).\n", bf->info);
      break;
    }
    return NLP_STOPPING_FAILURE_STEP;
  }

  if (!(bf->info))
    c_warning("Wolfe conditions not fullfilled by the found step\n");

  relative_gradient = 
    nlp_relative_gradient(bf->value, bf->x, bf->grad, NULL, bf->n, 1.0);

  /*  if (bf->output && bf->verbose)
      fprintf(bf->output, "Relative gradient: %g\n", relative_gradient); */

  if (relative_gradient < bf->tol) return NLP_STOPPING_GRADIENT;

  return NLP_STOPPING_NOT;
}

/**
 * BTR tests
 */
double nlp_btr_relative_step(double rho, double *x, double *step, int n)
{
  int i;
  double relx = 0.0;
  /*  double typxi; */

  for (i = 0; i < n; i++) {
    /*    typxi = (typx) ? typx[i] : 1.0; */
    relx = nt_max(relx, fabs(step[i])/nt_max(fabs(x[i]), 1.0));
  }

  return relx;
}

int btr_constrained_stopping_test(BTR *btr, double rho)
{
  if (nlp_relative_gradient(btr->value, btr->x, btr->grad, NULL, btr->n, 1.0)
      < btr->tol)
    return NLP_STOPPING_GRADIENT;

  if (btr->delta < btr->tol) return NLP_STOPPING_DELTA;
  if (cblas_dnrm2(btr->n, btr->step, 1) < 0.1*btr->tol)
    return NLP_STOPPING_STEP;

  return NLP_STOPPING_NOT;
}

/**
 * Standard stopping test for constrained optimization relying on projections.
 *
 * Ref. Conn, Gould and Toint, "Trust-Region Methods", SIAM, 2000,
 * Section 12.1, especially page 450.
 */
int btr_constrained_projection_stopping_test(BTR *btr, double rho)
{
  double g, relative_gradient;
  
  double *y;

  y = _mallocs(btr->n, double);

  /* projection(x-grad(x)) - x = 0 */
  cblas_dcopy(btr->n, btr->x, 1, y, 1);
  cblas_daxpy(btr->n, -1.0, btr->grad, 1, y, 1);
  // nt_vect_print(stdout, "y", y, btr->n);
  btr->projection_wrapper(btr, y);
  cblas_daxpy(btr->n, -1.0, btr->x, 1, y, 1);
  // nt_vect_print(stdout, "y", y, btr->n);

  /*  nt_vect_print(stdout, "g", btr->grad, btr->n);
  nt_vect_print(stdout, "x", btr->x, btr->n);
  nt_vect_print(stdout, "y", y, btr->n); */

  g = cblas_dnrm2(btr->n, y, 1);
  //  if (btr->verbose) {
    printf("Projected gradient path norm: %g\n", g);
    printf("Step norm: %g\n", cblas_dnrm2(btr->n, btr->step, 1));
    //  }

  free(y);

  if (g < btr->tol) return NLP_STOPPING_PROJECTED_GRADIENT_PATH;  
  if (btr->k == btr->max_iter) return NLP_STOPPING_ITERATIONS;
  if (cblas_dnrm2(btr->n, btr->step, 1) < 0.1*btr->tol)
    return NLP_STOPPING_STEP;

  relative_gradient =
    nlp_relative_gradient(btr->value, btr->x, btr->grad, NULL, btr->n, 1.0);
  if (relative_gradient < btr->tol)
    return NLP_STOPPING_RELATIVE_GRADIENT;
 
 return NLP_STOPPING_NOT;
}

/**
 * Ref. Conn, Gould and Toint, "Trust-Region Methods", SIAM, 2000,
 * Section 12.1, especially page 450.
 */
int btr_stop_projection_relative(BTR *btr, double rho,
				 double val, double *typx, double typf)
{
  double relative_gradient;
  double gmax = 0.0;
  double sqrtol = btr->tol*btr->tol;
  int i, n;
  double typxi;
  
  double *y;

  n = btr->n;
  y = _mallocs(n, double);

  cblas_dcopy(n, btr->x, 1, y, 1);
  cblas_daxpy(n, -1.0, btr->grad, 1, y, 1);
  btr->projection_wrapper(btr, y);
  cblas_daxpy(n, -1.0, btr->x, 1, y, 1);

  for (i = 0; i < n; i++) {
    typxi = (typx) ? typx[i] : 1.0;
    gmax = nt_max(gmax,
		  fabs(y[i]*nt_max(fabs(btr->x[i]), typxi))/
		  nt_max(fabs(val), typf));
  }

  if (gmax < sqrtol) return NLP_STOPPING_PROJECTED_RELATIVE_GRADIENT;  

  free(y);

  if (btr->k == btr->max_iter) return NLP_STOPPING_ITERATIONS;

  relative_gradient =
    nlp_relative_gradient(btr->value, btr->x, btr->grad, NULL, btr->n, 1.0);
  if (relative_gradient < btr->tol)
    return NLP_STOPPING_RELATIVE_GRADIENT;
 
 return NLP_STOPPING_NOT;
}

int nlp_btr_stopping_test(BTR *btr, double rho)
{
  double relative_gradient, relative_step;

  relative_gradient =
    nlp_relative_gradient(btr->value, btr->x, btr->grad, NULL, btr->n, 1.0);

  if (btr->output && btr->verbose)
    fprintf(btr->output, "Relative gradient: %g\n", relative_gradient);

  if (relative_gradient < btr->tol)
    return NLP_STOPPING_RELATIVE_GRADIENT;

  if (btr->delta < btr->tol) return NLP_STOPPING_DELTA;

  /* step... */
  if (btr_very_successful(btr, rho)) {
    relative_step = nlp_btr_relative_step(rho, btr->x, btr->step, btr->n);
    if (relative_step < btr->tol*btr->tol) return NLP_STOPPING_STEP;
  }

  if (btr->k >= btr->max_iter) return NLP_STOPPING_ITERATIONS;

  if (btr->delta < btr->tol) return NLP_STOPPING_DELTA;

  return NLP_STOPPING_NOT;
}

int nlp_btr_primal_stopping_test(BTR *btr, double rho)
{
  double relative_gradient, relative_step;

  relative_gradient =
    nlp_relative_gradient(btr->value, btr->x, btr->grad, NULL, btr->n, 1.0);
  if (relative_gradient < btr->tol) {
    printf("NLP_STOPPING_RGRAD");
    return NLP_STOPPING_RELATIVE_GRADIENT;
  }

  if (btr->delta < btr->tol) {
    printf("NLP_STOPPING_DELTA");
    return NLP_STOPPING_DELTA;
  }

  /* step... */
  relative_step = nlp_btr_relative_step(rho, btr->x, btr->step, btr->n);
  if (btr_successful(btr, rho) && relative_step < btr->tol) {
    printf("NLP_STOPPING_STEP");
    return NLP_STOPPING_STEP;
  }

  /* distance to the feasible set */
  /*  if (btr->dist)
    if (btr->dist(btr->x, btr->n, btr->data) < btr->tol) {
      printf("NLP_STOPPING_DIST");
      return NLP_STOPPING_DIST;
      } */

  /*  if (btr->output && btr->verbose)
      fprintf(btr->output, "Relative gradient: %g\n", relative_gradient); */

  return 0;
}
