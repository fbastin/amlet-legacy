/**
 * stopping.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 1999-2004 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the MIT License
 */

#ifndef __NLP2_STOPPING_H__
#define __NLP2_STOPPING_H__

#include "btr_struct.h"
#include "bfgs_struct.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Constants regarding the stopping test that is fulfilled
   */
  enum {
    NLP_STOPPING_NOT = 0,
    /* gradient 2-norm */
    NLP_STOPPING_GRADIENT,
    /* projected gradient path 2-norm */
    NLP_STOPPING_PROJECTED_GRADIENT_PATH,
    /* relative gradient: see Dennis and Schnabel (1983) */
    NLP_STOPPING_RELATIVE_GRADIENT,
    NLP_STOPPING_PROJECTED_RELATIVE_GRADIENT,
    /* too small decrease of the objective */
    NLP_STOPPING_DECREASE,
    /* too small step */
    NLP_STOPPING_STEP,
    /* distance to the feasible set */
    NLP_STOPPING_DIST,
    /* unidentified error */
    NLP_STOPPING_UNDEFINED_ERROR,
    /* too small trust-region radius */
    NLP_STOPPING_DELTA,
    /* number of iterations */
    NLP_STOPPING_ITERATIONS,
    /* no step found (linesearch) */
    NLP_STOPPING_FAILURE_STEP,
    /* infeasibility (feasible approaches) */
    NLP_STOPPING_INFEASIBILITY,
    /* no possible progress due to rounding errors */
    NLP_STOPPING_ROUNDING_ERRORS,
    /* incorrect parameters; the optimization cannot be performed */
    NLP_STOPPING_PARAMETERS,
    /* user stopping */
    NLP_STOPPING_USER
  };

  /* definir message variable */
  const char *nlp_stopping_msg(int id);

  /**
   * Convergence indicator.
   *
   * @return 1 is the convergence (usually to a first-order critical
   * point) has been observed, 0 else.
   */
  int nlp_has_converged(int id);

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
			       double *typx, int n, double typf);

  /**
   * Generic tests
   */

  /**
   * BFGS tests
   */

  /**
   * Stopping test for BFGS linesearch method.
   */
  int nlp_bfgs_stopping_test(BFGS *bf);

  /**
   * BTR specific stopping tests.
   */
  int btr_constrained_projection_stopping_test(BTR *btr, double rho);

  int nlp_btr_primal_stopping_test(BTR *btr, double rho);

  int btr_constrained_stopping_test(BTR *btr, double rho);

  int nlp_btr_stopping_test(BTR *btr, double rho);

#ifdef __cplusplus
}
#endif

#endif /* __NLP2_STOPPING_H__ */
