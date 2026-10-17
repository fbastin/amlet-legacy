/**
 * bsplines.c
 *
 * ORATIO Library
 * Copyright (C) 2005-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <oratio/ctools.h>
#include <oratio/bsplines.h>

#include <string.h>
#include <stdio.h>

/**
 * Compute the B-splines space dimension.
 *
 * @param p is the degree of the spline basis functions (order p+1).
 * @param m is the dimension of the open knot vector, of the form (2.13)
 *   in Piegl and Tiller, p. 66, that is, the first and last knots have
 *   multiplicity p+1, while interior knots have multiplicity 1.
 *
 * @return the B-splines space dimension.
 */
int nt_bsplines_space_dimension(int p, int m)
{
  /**
   * Reference: P2.8, p. 70, in Piegl and Tiller.
   */

  return m-p-1;
}

/**
 * Create a new b-spline object.
 *
 * @param p is the degree of the spline basis functions (order p+1).
 * @param nk is the dimension of the knots vector (with all
 *   multiplicities set to 1).
 * @param knots is the vectors of knots points (with all
 *   multiplicities set to 1); it is unchanged on output.
 * @param inck is the increment of knots (typically 1).
 */
NT_BSpline *nt_bsplines_new(int p, int nk,
			    double *knots, int inck)
{
  NT_BSpline *spline;

  spline = _malloc(NT_BSpline);

  spline->p = p;
  spline->n = nt_bsplines_knots_create(p, nk, knots, &(spline->n),
				       &(spline->U));
  spline->m = nt_bsplines_space_dimension(p, spline->n);

  spline->N = _mallocs(p+1, double);

  return spline;
}

/**
 * Free the memory associated to a B-spline.
 *
 * @param bspline is the B-spline object to free.
 */
void nt_bsplines_free(NT_BSpline *bspline)
{
  if (bspline) {
    free(bspline->U);
    free(bspline->N); 
    free(bspline);
  }
}

/**
 * Create a knot span.
 *
 * @param p is the degree of the spline basis functions (order p+1).
 * @param nk is the dimension of the vector x.
 * @param k is the knot vector with all multiplicities set to one,
 *   and increasing components.
 * @param nU is, on output, the dimension of the vector U.
 * @param U is on output a newly memory-allocated knot vector,
 *   with increasing components, of the form (2.13) in Piegl and
 *   Tiller, p. 66, that is, the first and last knots have multiplicity
 *   p+1, while interior knots have multiplicity 1.
 *
 * @return the value of nU, i.e. the knots vector dimension.
 */
int nt_bsplines_knots_create(int p, int nx, double *x,
			     int *nU, double **U)
{
  int i;

  *nU = nx+2*p;

  *U = _mallocs(*nU, double);

  memcpy(*U+p, x, nx*sizeof(double));

  for (i = 0; i < p; i++) (*U)[i] = (*U)[p];
  for (i = p; i > 0; i--) (*U)[*nU-i] = (*U)[*nU-p-1];

  return *nU;
}

/**
 * Compute the knot vector dimension, according to the form (2.13) in
 * Piegl and Tiller, p. 66
 *
 * @param p is the degree of the spline basis functions (order p+1).
 * @param nk is the dimension of the knots vector (with all
 *   multiplicities set to 1).
 *
 * @return the B-splines knot vector dimension.
 */
int nt_bsplines_knots_dimension(int p, int nk)
{
  return nk+2*p;
}

/**
 * Allocate the memory needed for the knots vector.
 *
 * @param p is the degree of the spline basis functions (order p+1).
 * @param nk is the dimension of the knots vector, with all
 *   multiplicities set to 1.
 *
 * @return a newly allocated vector with the right dimension.
 */
double *nt_bsplines_knots_allocate(int p, int nk)
{
  return _mallocs(nt_bsplines_knots_dimension(p, nk), double);
}

/**
 * Initialize a knot span.
 *
 * @param p is the degree of the spline basis functions (order p+1).
 * @param n is, on input the dimension of the vector *U, where each
 *   component has a multiplicity equal to 1, and, on output, the
 *   dimension of the vector *U with repetitions.
 * @param U is on output aa knot vector, with increasing components,
 *   of the form (2.13) in Piegl and Tiller, p. 66, that is, the first and
 *   last knots have multiplicity p+1, while interior knots have
 *   multiplicity 1.
 * @param realloc is equal on input to 1 if a reallocation of the
 *   vector *U has to be performed, 0 else (in which case it is the user
 *   responsability to check that *U has a sufficient size).
 *
 * @return the value of n.
 */
int nt_bsplines_knots_initialize(int p, int *n, double **U,
				 int realloc)
{
  int i;

  if (realloc) *U = _reallocs(*U, *n+2*p, double);
  memmove(*U+p, *U, (*n)*sizeof(double));
  *n += 2*p;

  for (i = 0; i < p; i++) (*U)[i] = (*U)[p];
  for (i = p; i > 0; i--) (*U)[*n-i] = (*U)[*n-p-1];

  return *n;
}

/**
 * Determine the knot span index.
 *
 * @param n is the dimension of the knot vector.
 * @param p is the degree of the spline basis functions (order p+1).
 * @param u is some real belonging to the domain of interest.
 * @param U is a knot vector, with increasing components, of the form
 *   (2.13) in Piegl and Tiller, p. 66, that is, the first and last
 *   knots have multiplicity p+1, while interior knots have
 *   multiplicity 1.
 *
 * @return the span index, i.e. the index such that u belongs to
 *   [U[i], U[i+1]), expect if u is equal to U[n+1], in which case we
 *   set the index to n.
 */
int nt_bsplines_find_span(int n, int p, double u, double U[])
{
  /**
   * Algorithm A2.1
   */

  int low, mid, high;

  /* Special case */
  if (u == U[n-p]) return n-p-1;

  /* Do binary search */
  low = p;
  high = n-p;
  mid = (low+high)/2;

  while (u < U[mid] || u >= U[mid+1]) {
    if (u < U[mid])
      high = mid;
    else
      low = mid;
    mid = (low+high)/2;
  }

  return mid;
}

/**
 * Computes the nonvanishing basis functions and stores them in the
 * array N[0],...N[p].
 *
 * @param i is the span index for the real u.
 * @param u is the point where the basis functions have to be evaluated.`
 * @param p is the degree of the spline basis functions (order p+1).
 * @param U is a knot vector, with increasing components, of the form
 *   (2.13) in Piegl and Tiller, p. 66, that is, the first and last
 *   knots have multiplicity p+1, while interior knots have
 *   multiplicity 1.
 * @param N must be of dimension p+1.
 */
void nt_bsplines_basis_func(int i, double u, int p, double U[], double N[])
{
  /**
   * Reference: Algorithm A2.2, p. 70, in Piegl and Tiller.
   */

  double *left, *right;
  double saved, temp;

  int j, r;

  left = _mallocs(p+1, double);
  right = _mallocs(p+1, double);

  N[0] = 1.0; /* the others elements can be assumed to be set to 0 */

  for (j = 1; j <= p; j++) {
    left[j] = u - U[i+1-j];
    right[j] = U[i+j] - u;
    saved = 0.0;
    for (r = 0; r < j; r++) {
      temp = N[r]/(right[r+1]+left[j-r]);
      N[r] = saved+right[r+1]*temp;
      saved = left[j-r]*temp;
    }
    N[j] = saved;
  }

  free(left);
  free(right);
}

/**
 * Evaluate the basis function Nip at u.
 *
 * @param p is the degree of the spline basis functions (order p+1).
 * @param m is the high index of U (m+1 knots).
 * @param i is the index of basis function.
 * @param u is the real where the basis function is evaluated.
 * @param U is a knot vector, with increasing components, of the form
 *   (2.13) in Piegl and Tiller, p. 66, that is, the first and last
 *   knots have multiplicity p+1, while interior knots have
 *   multiplicity 1.
 *
 * @return Nip(u).
 */
double nt_bsplines_one_basis_func(int p, int m,
				      double U[], int i, double u)
{
  /**
   * Reference: Algorithm A2.4, pp. 74-75, in Piegl and Tiller.
   */

  double *N;
  double Uleft, Uright;
  double saved, temp;

  int j, k;

  /* Special cases: we face a bound */
  if ((i == 0 && u == U[0]) || (i == m-p-1 && u == U[m])) return 1.0;

  /* Local property (see p. 55) */
  if (u < U[i] || u >= U[i+p+1]) return 0.0;

  N = _mallocs(p+1, double);

  /* Initialize zeroth-degree functions */
  for (j = 0; j <= p; j++) {
    if (u >= U[i+j] && u < U[i+j+1])
      N[j] = 1.0;
    else
      N[j] = 0.0;
  }

  /* Compute triangular table */
  for (k = 1; k <= p; k++) {
    saved = (N[0] == 0.0) ? 0.0 : ((u-U[i])*N[0])/(U[i+k]-U[i]);
    for (j = 0; j < p-k+1; j++) {
      Uleft = U[i+j+1];
      Uright = U[i+j+k+1];
      if (N[j+1] == 0.0) {
	N[j] = saved;
	saved = 0.0;
      }
      else {
	temp = N[j+1]/(Uright-Uleft);
	N[j] = saved+(Uright-u)*temp;
	saved = (u-Uleft)*temp;
      }
    }
  }

  temp = N[0];
  free(N);

  return temp;
}

/**
 * Evaluate a spline curve at some point.
 *
 * @param n is the dimension of the knot vector (under the form (2.3)).
 * @param p is the degree of the spline basis functions (order p+1).
 * @param u is the point where to evaluate the spline.
 * @param U is an open knot vector, with increasing components, of
 *   dimension n, following the form (2.13) in Piegl and Tiller, p. 66.
 * @param P is the vector of control points, of dimension n-p.
 * @param N is a double vector of dimension at least equal to p+1.
 *
 * @return the spline curve value.
 */
double nt_bsplines_curve_point(int n, int p, double u, double U[],
			       double P[], double N[])
{
  /**
   * Reference: Algorithm A3.1, p. 82, in Piegl and Tiller.
   */
  double val = 0.0;

  int span, i;

  span = nt_bsplines_find_span(n, p, u, U);

  nt_bsplines_basis_func(span, u, p, U, N);

  for (i = 0; i <= p; i++) val += N[i]*P[span-p+i];

  return val;
}

/**
 * Evaluation a B-Spline at a particular point.
 *
 * @param b is a space of B-splines.
 * @param P is the vector of control points defining a particular B-spline.
 * @param u is the point where to evaluate the B_spline.
 */
double nt_bsplines_evaluate(NT_BSpline *b, double *P, double u)
{
  return nt_bsplines_curve_point(b->n, b->p, u, b->U, P, b->N);
}
