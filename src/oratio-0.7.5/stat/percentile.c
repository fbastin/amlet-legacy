/**
 * percentile.c
 *
 * STAT Library
 * Copyright (C) 1999-2005 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "percentile.h"

#include <float.h>
#include <math.h>

#include <oratio/ctools.h>

/**
 * Bissection procedure.
 */
double st_percentile_bissect(double q, STCDF f, double *pars,
			     double lo, double up, double tol)
{
  double mid, alpha;

  alpha = q+2*q*tol;
  mid = 0.5*(lo+up);

  while (fabs(q-alpha) > tol*q) {
    printf("Coucou 4 %g %g %g\n", lo, up,  alpha);
    mid = 0.5*(lo+up);
    alpha = f(mid, pars);
    if (alpha < q)
      lo = mid;
    else
      up = mid;
  }

  return mid;
}

/**
 * Computation of percentiles of arbitrary continuous cumulative
 * distribution functions.
 *
 * @param alpha is the target qunatile.
 * @param pars is the vector of parameters associated to the random variate.
 * @param tol is the used tolerance; if tol <= 0.0, a default
 *   accuracy will be set.
 */
double st_percentile(double alpha, STCDF f, double *pars,
		     double median, double sigma, double tol)
{
  double q;
  double lo, up;   /* lower and upper bounds of the interalpha */
  double k;

  /**
   * Strategy: we first aim to define some interalpha contained the
   * desired point. Once this is done we proceed to a standard
   * bissection.
   */

  printf("Coucou 1\n");

  if (alpha <= 0.0 || alpha >= 1.0) {
    c_warning("invalid value %g for routine st_percentile_func");
    return 0.0;
  }

  if (tol <= 0.0) tol = 100*DBL_EPSILON; /* quite arbitrary */
  lo = median;
  k = sigma;

  q = f(median, pars);

  printf("3. Coucou %g %g\n", median, q);

  if (fabs(alpha-q) <= tol*alpha) return median;

  if (alpha > q) {
    /* compute the upper bound */
    up = median+sigma;
    while ((q = f(up, pars)) < alpha) {
      k *= 2;
      up += k;
    }
  }
  else {
    /* compute the lower bound */
    up = median;
    lo = median-sigma;
    while (f(lo, pars) > alpha) {
      k *= 2;
      lo -= k;
    }
  }

  printf("Coucou 3\n");

  /* Bissection */
  return st_percentile_bissect(alpha, f, pars, lo, up, tol);
}
