/**
 * activesets.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <math.h>
#include <float.h>
#include <stdio.h>

#include "activesets.h"
#include "bounds.h"

/**
 * Bounds
 */

/**
 * Free variables initialization. By default, all variables are free.
 *
 * @param n is the problem dimension.
 */
void nlp_as_free_init_dflt(int n, int *free, int incfree)
{
  int i;

  for (i = 0; i < n; i += incfree) free[i] = 1.0;
}

/**
 * Determine the active set for bounds constraints.
 *
 * @param n is the problem dimension.
 * @param x is the point where the active set is determined.
 * @param incx is the increment of x.
 */
int nlp_as_bounds(void *as, int n, double *x, int incx,
		  double *g, int incg, int *free, int incfree)
{
  NlpBounds *b = (NlpBounds *)as;

  int nfree, i;

  nfree = n;

  for (i = 0; i < n; i++) {
    if ((fabs(b->xl[i]-x[i*incx]) <= DBL_EPSILON && g[i*incg] > 0) ||
	(fabs(x[i*incx]-b->xu[i]) <= DBL_EPSILON && g[i*incg] < 0)) {
      free[i*incfree] = 0;
      nfree--;
    }
  }

  return nfree;
}

/**
 * Determine the active set for bounds constraints and step validation.
 *
 * @param n is the problem dimension.
 * @param x is the point where the active set is determined.
 * @param incx is the increment of x.
 */
int nlp_as_bounds_step(void *as, int n, double *x, int incx,
		       double *step, int incs, double *p, int incp,
		       int *free, int incfree, int *nfree)
{
  NlpBounds *b = (NlpBounds *)as;

  int i;

  /* Update the indexes of free variables */
  for (i = 0; i < n; i++) {
    if (free[i*incfree]) {
      if (p[i*incp] > 0 && fabs(x[i*incx]+step[i*incs]-b->xu[i])
	  <= DBL_EPSILON) {
	step[i*incs] = b->xu[i]-x[i*incx];
	free[i*incfree] = 0;
	(*nfree)--;
      }
      else if (p[i*incp] < 0 && fabs(x[i*incx]+step[i*incs]-b->xl[i])
	       <= DBL_EPSILON ) {
	step[i*incs] = b->xl[i]-x[i*incx];
	free[i*incfree] = 0;
	(*nfree)--;
      }
    }
  }

  return *nfree;
}

/**
 * Compute the restricted gradient under bounds constraints.
 *
 * @param n is the problem dimension.
 */
void nlp_as_rg_bounds(int n, double *p, int incp,
			     double *g, int incg, int *xfree, int incf,
			     int nfree)
{
  int i;

  /* Set the gradient to 0 at the active bounds. */
  for (i = 0; i < n; i++) p[i*incp] = g[i*incg]*xfree[i*incf];
}

/**
 * Monotonicity
 */

/**
 * Determine the active set for monotonicity constraints and step validation.
 *
 * @param n is the problem dimension.
 * @param x is the point where the active set is determined.
 * @param incx is the increment of x.
 */
int nlp_as_monotonicity_step(NLP_Constraint *c, int n, double *x, int incx,
			     double *step, int incs, double *p, int incp,
			     int *nsets, int *sets, int incsets)
{
  int i, j, k, l;

  //  double mean;

  j = 0;
  l = 0;

  for (i = 1; i < *nsets; i++) {
    k = sets[i*incsets];
    if (j >= c->cidx[l] && k < c->cidx[l+1]) {
      if ((x[k*incx]+step[k*incs]) - (x[j*incx]+step[j*incs]) <= DBL_EPSILON) {
	if ((fabs(p[k*incp]) <= DBL_EPSILON && fabs(p[j*incp]
						    <= DBL_EPSILON))
	    || (p[j*incp] >= p[k*incp]))
	  /* The variables will become linked. */
	  sets[i*incsets] = -sets[i*incsets];
      }
    }
    else if (j >= c->cidx[l]) l += 2;
    j = k;
  }

  /* Update the active sets. */
  j = 1;
  for (i = 0; i < *nsets-1; i++) {
    if (sets[(i+1)*incsets] > 0) {
      sets[(j++)*incsets] = sets[(i+1)*incsets];
    }
  }
  sets[j*incsets] = n;
  *nsets = j;

  /* We normalize the steps. */
  /* CORRECT??? */
  /*  for (i = 0; i < *nsets; i++) {
    mean = 0.0;
    if (sets[(i+1)*incsets]-sets[i*incsets]) {
      for (j = sets[i*incsets]; j < sets[(i+1)*incsets]; j++)
	mean += (x[j*incx]+step[j*incs]-mean)/(j+1-sets[i*incsets]);
      for (j = sets[i*incsets]; j < sets[(i+1)*incsets]; j++)
	step[j*incs] = mean - x[j*incx];
	}
	} */

  return *nsets;
}

/**
 * Determine the active set under monotonicity constraints.
 *
 * @param c is a monotonicity constraint.
 * @param n is the problem dimension.
 * @param x is the point where the active set is determined.
 * @param incx is the increment of x.
 * @param p is the search direction.
 * @param sets[i] designs the start of the i-th active set.
 * @param incsets is the increment of sets.
 */
int nlp_as_monotonicity(NLP_Constraint *c, int n, double *x, int incx,
			double *p, int incp, int *sets, int incsets)
{
  int free, nsets;
  int i, j;

  if (!c) {
    for (i = 0; i <= n; i++) sets[i*incsets] = i;
    return 1;
  }

  nsets = 0;
  for (i = 0; i <= c->cidx[0]; i++) sets[(nsets++)*incsets] = i;
  
  for (j = 0; j < c->nidx; j++) {

    /* Start a new active set. */
    if (j) {
      for (i = c->cidx[2*j-1]; i <= c->cidx[2*j]; i++)
	sets[(nsets++)*incsets] = i;
    }

    /* Get the extent of the active set. */
    for (i = c->cidx[2*j]+1; i < c->cidx[2*j+1]; i++) {
      free = 1;
      if (fabs(x[(i-1)*incx] - x[i*incx]) <= DBL_EPSILON) {
	if (fabs(p[(i-1)*incp]) <= DBL_EPSILON &&
	    fabs(p[i*incp] <= DBL_EPSILON))
	  free = 0;
	else if (p[(i-1)*incp] >= p[i]) free = 0;
      }
      if (free) sets[(nsets++)*incsets] = i;
    }
  }

  for (i = c->cidx[2*c->nidx-1]; i < n; i++) sets[(nsets++)*incsets] = i;
  sets[nsets*incsets] = n;

  return nsets;
}

/**
 * Compute the restricted gradient under monotonicity constraints.
 *
 * @param n is the problem dimension.
 * @param incg is the increment in g.
 */
void nlp_as_rg_monotonicity(int n, double *p, int incp,
			    double *g, int incg, int *sets, int incs,
			    int nsets)
{
  int i, j;

  for (i = 0; i < nsets; i++) {
    p[sets[i]] = g[sets[i]];
    for (j = sets[i*incs]+1; j < sets[(i+1)*incs]; j++) p[sets[i*incs]] += g[j];
    for (j = sets[i*incs]+1; j < sets[(i+1)*incs]; j++) p[j] = p[sets[i*incs]];
  }
}

/**
 * Compute the 2-norm of the reduced gradient.
 */
double nlp_as_rg_nrm2(int n, double *g, int incg, int nsets, int *sets, int incs)
{
  double nrm = 0.0, tmp;

  int i;

  for (i = 0; i < nsets; i++) {
    tmp = g[sets[i*incs]*incg];
    nrm += tmp*tmp;
  }

  return sqrt(nrm);
}

double nlp_as_rdot(int n, double *x, int incx, double *y, int incy,
		   int nsets, int *sets, int incs)
{
  double tmp = 0.0;

  int i;

  for (i = 0; i < nsets; i++) {
    tmp = x[sets[i*incs]*incx]*y[i*sets[incs]*incy];
  }

  return tmp;
}

void nlp_as_monotonicity_print(int nset, int *sets)
{
  int i;

  for (i = 0; i < nset; i++) printf("%d ", sets[i]);
  printf("%d\n", sets[nset]);
}
