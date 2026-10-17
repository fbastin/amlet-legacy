/**
 * bounds.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 1999-2002 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "bounds.h"

#include "projection.h"

#include <float.h>
#include <math.h>

#include <oratio/ctools.h>
#include <oratio/ntools.h>

/**
 * Reset lower and upper bounds
 */
void nlp_bounds_reset(double *bounds, int n)
{
  int i;

  for (i = 0; i < n; i++) {
    bounds[2*i] = -DBL_MAX;
    bounds[2*i+1] = DBL_MAX;
  }
}

double *nlp_bounds_init(double **bounds, int n)
{
  *bounds = _mallocs(2*n, double);
  nlp_bounds_reset(*bounds, n);

  return *bounds;
}

/**
 * Test the feasibility of point x for bounds constraints.
 *
 * @return 1 is x is feasible, 0 else.
 */
int nlp_bounds_feasibility(double *x, int n, void *data)
{
  double *bounds = (double *)data;
  int i;

  for (i = 0; i < n; i++)
    if ((x[i] < bounds[2*i]) || (x[i] > bounds[2*i+1])) return 0;

  return 1;
}

int nlp_bounds_strict_feasibility(double *x, int n, void *data)
{
  double eps = 1e-8;
  double *bounds = (double *)data;
  int i;

  for (i = 0; i < n; i++)
    if ((x[i] < bounds[2*i] + eps) || (x[i] > bounds[2*i+1] + eps)) return 0;

  return 1;
}

/**
 * Compute the distance between a point x and the bounds described in data.
 *
 * We assume that x is feasible.
 */
double nlp_bounds_dist(double *x, int n, void *data)
{
  double *bounds = (double *)data;
  double dist, t;
  int i;

  dist = DBL_MAX;

  for (i = 0; i < n; i++) {
    t = nt_min(fabs(bounds[2*i]-x[i]), fabs(bounds[2*i+1]-x[i]));
    dist = nt_min(t, dist);
  }

  return dist;
}

void nlp_bounds_set(double *bounds, int i, double lower, double upper)
{
  bounds[2*i] = lower;
  bounds[2*i+1] = upper;
}

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
				 void *data)
{
  int i;
  double *bounds = (double *)data;

  for (i = 0; i < n; i++)
    x_proj[i] = nt_max(bounds[2*i], nt_min(bounds[2*i+1], x[i]));

  return 0;
}

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
			  double *w, double *brptmin, double *brptmax)
{
  /**
   * Ref.: Lin and More, "Newton's Method For Large Bound-Constrained
   * Optimization Problems", Siam Journal on Optimization 9(4), 1999.
   * TRON v1.2.
   */

  int i;
  int nbrpt = 0;
  double brpt;

  *brptmin = -DBL_MAX;
  *brptmax = DBL_MAX;

  for (i = 0; i < n; i++) {
    if (x[i] < xu[i] && w[i] > 0.0) {
      nbrpt++;
      brpt = (x[i]-xu[i])/w[i];
      *brptmin = nt_min(brpt, *brptmin);
      *brptmax = nt_min(brpt, *brptmax);
    }
    else if (x[i] > xl[i] && w[i] < 0.0) {
      nbrpt++;
      brpt = (xl[i]-x[i])/w[i];
      *brptmin = nt_min(brpt, *brptmin);
      *brptmax = nt_min(brpt, *brptmax);
    }
  }

  /* Handle the exceptional case */
  if (nbrpt == 0) *brptmin = *brptmax = 0.0;

  return nbrpt;
}

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
			  double alpha, double *w, double *s)
{
  /**
   * Ref.: Lin and More, "Newton's Method For Large Bound-Constrained
   * Optimization Problems", Siam Journal on Optimization 9(4), 1999.
   * TRON v1.2.
   */

  int i;
  double tmp;

  for (i = 0; i < n; i++) {
    tmp = alpha*w[i];
    if (x[i] + tmp < xl[i])
      s[i] = xl[i] - x[i];
    else if (x[i] + tmp > xu[i])
      s[i] = xu[i] - x[i];
    else
      s[i] = tmp;
  }
}

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
int nlp_bounds_projection_st(int n, double *x, double *w,
				    double *x_proj, void *bounds)
{
  int i;

  for (i = 0; i < n; i++)
    x_proj[i] = nt_max(((NlpBounds *)bounds)->xl[i],
		       nt_min(((NlpBounds *)bounds)->xu[i], x[i]));

  return 0;
}
