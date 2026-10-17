/**
 * dcl.h
 *
 * AMLET Library
 * Copyright (C) 2010-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __AMLET_DCL_H__
#define __AMLET_DCL_H__

#include <ophelia/ophelia.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * STRUCTURE DEFINITIONS
   */

  /**
   * Constants for optimization techniques.
   */
  enum {
    NO_OPTIMIZATION,
    MCBTR,
    MCRTR,
    MCRTRDA,
    MCBTRDA,
    MCBFGS,
    MCBFGSDA,
    MCBTRC
  };

  /**
   * Constants for Hessian approximations.
   */
  enum {
    HESSIAN_BFGS,
    HESSIAN_BHHH,
    HESSIAN_SR1
  };

  union _optimethod {
    BTR  *b;   /* used for BTR/RTR algorithm, in the optimization process */
    BFGS *bf;  /* used for BTR algorithm, in the optimization process */
  };

  typedef union _optimethod DclOptiMethod;

  /**
   * Apply a safeguard to the utilities, such that the logit formula is
   * mathematically equivalent to the original one, but is numerically
   * more stable.
   *
   * @param utilities[i] is the utility value of alternative i; and on
   *   return utilities[i] is the exponential of the original value.
   * @param den is the sum of the exponential of utilities values on
   *   return; its value on input in undefined.
   * @param availability[i] is equal to 1 if the alternative i is
   *   available, 0 else.
   * @param nalt is the number of alternatives.
   * @param choice is the alternative choice if it superior or equal to
   *   0; a negative value can be use to indicate that the choice value is
   *   unknown.
   *
   * @return the probability choice if choice is positive (or null),
   *   -1.0 else.
   */
  double dcl_utilities_safeguard(double *utilities, double *den,
				 int *availability, int nalt, int choice);

  /**
   * @param n is the number of parameters to estimate.
   */
  void dcl_ml_init_solver(NTLog *mclog, int n,
			  int method, DclOptiMethod *opti,
			  double tol, int maxiter, int hessian_method);

  /**
   * Solve the maximum likelihood problem.
   *
   * @param n is the number of parameters to estimate.
   * @param method is the optimisation methods that has to be used.
   * @param n_random_coeff is the number of random coefficients that are
   *   present in the model.
   * @param verbose is set to 1 if the function is verbose, displaying
   *   various comments on screen, 0 else.
   * @param log is a ntools log file.
   */
  double dcl_ml_solve(NTLog *log, int n,
		      int method, DclOptiMethod *opti,
		      int n_random_coeff,
		      double **H, void *data,
		      int verbose);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_DCL_H__ */
