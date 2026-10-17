/**
 * barrier.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2005 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __NLP2_BARRIER_H__
#define __NLP2_BARRIER_H__

#ifdef __cplusplus
extern "C" {
#endif

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
  double nlp_update_barrier_parameter(double mu);

  /**
   * Compute the barrier function associated to bounds constraints.
   *
   * @param mu, the barrier parameter
   * @param n, the problem dimension
   */
  double nlp_log_barrier_bounds(double mu,
				double *x, void *bounds, int n);

  /**
   * Add to grad the log-barrier gradient, for bounds constraints.
   */
  void nlp_log_barrier_bounds_d1(double mu,
				 double *x, void *bounds, int n,
				 double *grad);

  /**
   * Add to hessian the log-barrier hessian, for bounds constraints.
   */
  void nlp_log_barrier_bounds_d2(double mu,
				 double *x, void *bounds, int n,
				 double **hessian);

  /**
   * Compute the barrier function associated to positivity constraints.
   *
   * @param mu, the barrier parameter
   */
  double nlp_log_barrier_positivity(double mu,
				    double *x, double n);

#ifdef __cplusplus
}
#endif

#endif  /* __NLP2_BARRIER_H__ */
