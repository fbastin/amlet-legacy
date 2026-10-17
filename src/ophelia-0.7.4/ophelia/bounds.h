/**
 * bounds.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 1999-2004 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __NLP2_BOUNDS_H__
#define __NLP2_BOUNDS_H__

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
    double *xl;  /* lower bounds */
    double *xu;  /* upper bounds */
  } NlpBounds;

  /**
   * Reset lower and upper bounds
   */
  void nlp_bounds_reset(double *bounds, int n);

  double *nlp_bounds_init(double **bounds, int n);

  /**
   * Test the feasibility of point x for bounds constraints.
   *
   * @return 1 is x is feasible, 0 else.
   */
  int nlp_bounds_feasibility(double *x, int n, void *data);

  int nlp_bounds_strict_feasibility(double *x, int n, void *data);

  /**
   * Compute the distance between a point x and the bounds described in data.
   *
   * We assume that x is feasible.
   */
  double nlp_bounds_dist(double *x, int n, void *data);

  void nlp_bounds_set(double *bounds, int i, double lower, double upper);

  /**
   * Project the point x onto the bounds defined by data.
   *
   * @param n is the problem dimension.
   * @param x is the point to project; on output, x is unchanged,
   *   unless x and x_proj design the same pointer.
   * @param w is a weighting vector (unused here).
   * @param data is a pointer to a double array, containing the bounds.
   * @param x_proj is a initialized vector of dimension n; it may be
   *   the same as x.
   * @param on output x_proj[i] is in [ data[2*i], data[2*i+1] ].
   *
   * @return 0 (no error code).
   */
  int nlp_bounds_projection(int n, double *x, double *w, double *x_proj,
			    void *data);

  /**
   * Bounds management when two double arrays are used.
   */

  /**
   * This function computes the number of break-points, and the minimal
   * and maximal break-points if the projection of x+alpha*w on the
   * n-dimensional interval [xl, xu].
   *
   * @param n is the problem dimension.
   * @param x is the vector x, assumed to be feasible; it is unchanged
   *   on output.
   * @param xl is the vector of lower bounds.
   * @param xu is the vector of upper bounds.
   * @param w is a direction vector.
   * @param *brptmin does not need to be specified on input, and is the
   *   minimal breakpoint on output.
   * @param *brptmin does not need to be specified on input, and is the
   *   maximal breakpoint on output.
   *
   * @return the number of break points.
   */
  int nlp_bounds_breakpt_st(int n, double *x, double *xl, double *xu,
			    double *w, double *brptmin, double *brptmax);

  /**
   * Compute the gradient projection step, for bounds constraints:
   * s = P[x+alpha*w]-x.
   *
   * @param n is the problem dimension.
   * @param x is the vector x, assumed to be feasible; it is unchanged
   *   on output.
   * @param xl is the vector of lower bounds.
   * @param xu is the vector of upper bounds.
   * @param w is a direction vector.
   * @param s does not need to be specified on input; it is the step
   * vector on output.
   */
  void nlp_bounds_gpstep_st(int n, double *x, double *xl, double *xu,
			    double alpha, double *w, double *s);

  /**
   * Project the point x onto the bounds defined by bounds, represented
   * by a structure nlp_bounds.
   *
   * @param n is the problem dimension.
   * @param x is the point to project; on output, x is unchanged,
   *   unless x and x_proj design the same pointer.
   * @param x_proj is a initialized vector of dimension n; it may be
   *   same as x.
   *
   * @return 0 (no error code).
   */
  int nlp_bounds_projection_st(int n, double *x, double *w, double *x_proj,
			       void *bounds);

#ifdef __cplusplus
}
#endif

#endif  /* __NLP2_BOUNDS__ */
