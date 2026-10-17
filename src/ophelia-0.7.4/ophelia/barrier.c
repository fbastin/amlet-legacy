/**
 * barrier.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2005 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "barrier.h"

#include <float.h>
#include <math.h>

#include <oratio/ntools.h>

/**
 * Update the barrier parameter.
 *
 * Ref: Andrew R. Conn, Nicholas M. Gould, Dominique Orban, and
 * Philippe L. Toint, "A Primal-Dual Trust-Region Algorithm for
 * Non-convex Nonlinear Programming", Mathematical Programming, 87(2):
 * 215-149, 2000.
 *
 * @param mu, the barrier parameter
 */
double nlp_update_barrier_parameter(double mu)
{
  return nt_min(0.1*mu, pow(mu, 1.5));
}

/**
 * Compute the barrier function associated to bounds constraints.
 *
 * @param mu, the barrier parameter
 */
double nlp_log_barrier_bounds(double mu,
			      double *x, void *bounds, int n)
{
  int i;
  double barrier = 0.0;
  double *b = (double *)bounds;

  for (i = 0; i < n; i++) {
    if (b[2*i] != -DBL_MAX) barrier += log(nt_max(0, x[i]-b[2*i]));
    if (b[2*i+1] != DBL_MAX) barrier += log(nt_max(0, b[2*i+1]-x[i]));
  }

  /*  if (isinf(barrier)) {
    printf("Zut\n");
    nt_vect_print(stdout, "x", x, n);
    } */

  return -mu*barrier;
}

/**
 * Add to grad the log-barrier gradient, for bounds constraints.
 */
void nlp_log_barrier_bounds_d1(double mu,
			       double *x, void *bounds, int n,
			       double *grad)
{
  int i;
  double *b = (double *)bounds;

  for (i = 0; i < n; i++) {
    if (b[2*i] != -DBL_MAX) grad[i] -= mu/(x[i]-b[2*i]);
    if (b[2*i+1] != DBL_MAX) grad[i] -= mu/(b[2*i+1]-x[i]);
  }
}

/**
 * Add to hessian the log-barrier hessian, for bounds constraints.
 */
void nlp_log_barrier_bounds_d2(double mu,
			       double *x, void *bounds, int n,
			       double **hessian)
{
  int i;
  double *b = (double *)bounds;
  double t;

  for (i = 0; i < n; i++) {
    if (b[2*i] != -DBL_MAX) {
      t = 1/(x[i]-b[2*i]);
      hessian[i][i] += mu*t*t;
    }
    if (b[2*i+1] != DBL_MAX) {
      t = 1/(b[2*i+1]-x[i]);
      hessian[i][i] += mu*t*t;
    }
  }
}

/**
 * Compute the barrier function associated to positivity constraints.
 *
 * @param mu, the barrier parameter
 */
double nlp_log_barrier_positivity(double mu,
				  double *x, double n)
{
  int i;
  double barrier = 0.0;

  for (i = 0; i < n; i++) barrier += log(x[i]);

  return -mu*barrier;
}
