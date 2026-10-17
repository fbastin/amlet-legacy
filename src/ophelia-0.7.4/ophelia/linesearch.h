/**
 * linesearch.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2002-2006 Laurent Martin
 *                         Emilie Wanufelle <ewa@math.fundp.ac.be>,
 *                         Fabian Bastin    <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __NLP2_LINESEARCH_H__
#define __NLP2_LINESEARCH_H__

#include "common_struct.h"
#include "projection.h"
#include "wrapper.h"

#include <oratio/ctools.h>

#ifdef __cplusplus
extern "C" {
#endif

#define LINESEARCH_STEP_MAX 1.e3
#define LINESEARCH_STEP_MIN 1.e-9

  enum {
    NLP_NEWTON_STEP = 0,
    NLP_BACKSTRAP_STEP
  };

  enum {
    NLP_LINESEARCH_WOLFE_CONDITIONS = 1,
    NLP_LINESEARCH_SUFFICIENT_DECREASE,
    NLP_LINESEARCH_SMALL_INTERVAL,
    NLP_LINESEARCH_LARGE_STEP,
    NLP_LINESEARCH_SMALL_STEP,
    NLP_LINESEARCH_ROUNDING_ERRORS,
    NLP_LINESEARCH_DESCENT_DIRECTION,
    NLP_LINESEARCH_IMPROPER_INPUT,
    NLP_LINESEARCH_FAILURE,
  };

  /**
   * Find a step which satisfies a sufficient decrease condition and a
   * curvature condition.
   *
   *   at each stage the subroutine updates an interval of 
   *   uncertainty with endpoints stx and sty. the interval of 
   *   uncertainty is initially chosen so that it contains a 
   *   minimizer of the modified function 
   *        f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s). 
   *
   *   if a step is obtained for which the modified function 
   *   has a nonpositive function value and nonnegative derivative, 
   *   then the interval of uncertainty is chosen so that it 
   *   contains a minimizer of f(x+stp*s). 
   *   the algorithm is designed to find a step which satisfies 
   *   the sufficient decrease condition 
   *         f(x+stp*s) .le. f(x) + ftol*stp*(gradf(x)'s), 
   *   and the curvature condition 
   *         abs(gradf(x+stp*s)'s)) .le. gtol*abs(gradf(x)'s). 
   *   if ftol is less than gtol and if, for example, the function 
   *   is bounded below, then there is always a step which satisfies 
   *   both conditions. if no step can be found which satisfies both 
   *   conditions, then the algorithm usually stops when rounding 
   *   errors prevent further progress. in this case stp only 
   *   satisfies the sufficient decrease condition. 
   *
   *   Parameters
   *     n is a positive int input variable set to the number 
   *       of variables. 
   *     x is an array of length n. on input it must contain the 
   *       base point for the line search. on output it contains 
   *       x + stp*s. 
   *     f is a variable. on input it must contain the value of f 
   *       at x. on output it contains the value of f at x + stp*s. 
   *     g is an array of length n. on input it must contain the 
   *       gradient of f at x. on output it contains the gradient 
   *       of f at x + stp*s. 
   *     s is an input array of length n which specifies the 
   *       search direction. 
   *     stp is a nonnegative variable. on input stp contains an 
   *       initial estimate of a satisfactory step. on output 
   *       stp contains the final estimate. 
   *     ftol and gtol are nonnegative input variables. 
   *       termination occurs when the sufficient decrease 
   *       condition and the directional derivative condition are 
   *       satisfied. 
   *       ftol should be smaller than 5.e-1
   *       suggested value = 1.e-4 for newton methods
   *                       = 1.e-1 for more exact line searches
   *       Default Value = 1.e-4
   *
   *       gtol should be greater than 1.e-4 
   *       Default Value = 0.9
   *     xtol is a nonnegative input variable. termination occurs 
   *       when the relative width of the interval of uncertainty 
   *       is at most xtol. 
   *       Default Value = 2.2e-16
   *     stpmin and stpmax are nonnegative input variables which 
   *       specify lower and upper bounds for the step. (in this reverse 
   *       communication implementatin they are defined in a common 
   *       statement). 
   *       stpmin Default Value = 1.e-9
   *       stpmax Default Value = 1.e3
   *     info is an integer output variable set as follows: 
   *       info =-1 improper input parameters. 
   *       info = 1  the sufficient decrease condition and the 
   *                 directional derivative condition hold. 
   *       info = 2  relative width of the interval of uncertainty 
   *                 is at most xtol. 
   *       info = 4  the step is at the lower bound stpmin. 
   *       info = 5  the step is at the upper bound stpmax. 
   *       info = 6  rounding errors prevent further progress. 
   *                 there may not be a step which satisfies the 
   *                 sufficient decrease and curvature conditions. 
   *                 tolerances may be too small.
   *     work is a working vector of minimum size n.
   * 
   *   functions called 
   *     dcstep 
   *   Initial version:
   *   Argonne National Laboratory. Minpack project. June 1983 
   *   Jorge J. More', David J. Thuente 
   */
  int dcsrch(int n, double *grad, double *fvalue, double *x,
	     nlp_wrapper_fct wf, C_GENERIC f, int f_type,
	     double (*fg)(double *, int, void *, double *),
	     void *data, deriv1 d1,
	     double *scaling, double *s,
	     double *stp, double ftol, double xtol, 
	     double gtol, double stpmax, double stpmin,
	     double *work);

  /**
   * Piecewise linesearch along a (projected) direction p), using a
   * quadratic model.
   *
   * @param Order designs the orientation of the matrix Bk. The
   *   possible values are CblasRowMajor and CblasColMajor.
   * @param Uplo designs the triangular part that has to be updated. The
   *   possible value CblasUpper and CblasLower.
   * @param n is the problem dimension.
   * @param x is the current iterate.
   * @param g is the gradient at x.
   * @param H is the hessian at x.
   * @param p is the direction search.
   * @param maxiter is the allowed maximum number of iterations.
   * @param data is a pointer to problem-specific data.
   *
   * @return the step length in direction p.
   */
  double nlp_projected_search(const enum CBLAS_ORDER Order,
			      const enum CBLAS_UPLO Uplo,
			      NlpProjection proj,
			      int n, double *x, double alpha,
			      double *g, double *H,
			      double *p,int maxiter,
			      double *s,
			      double *vplus, double *splus, double *gplus,
			      double *work, void *data, double alup,
			      void *data_projection);

#ifdef __cplusplus
}
#endif

#endif  /* __NLP2_LINESEARCH_H__ */
