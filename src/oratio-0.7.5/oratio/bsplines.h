/**
 * bsplines.h
 *
 * ORATIO Library
 * Copyright (C) 2006-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_BSPLINES_H__
#define __ORATIO_BSPLINES_H__

#ifdef __cplusplus
extern "C" {
#endif

/**
 * B-Spline estimation and interpolation
 *
 * The following code is an adaptation of functions presented in Piegl
 * and Tiller (1996).
 *
 * Reference:
 *   The NURBS Book
 *   Les Piegl and Wayne Tiller
 *   Springer-Verlag, New York, NY
 *   Second Edition, 1996
 */

  typedef struct {
    int p; /* degree of the spline basis functions (order p+1) */
    int n; /* dimension of the open knot vector */
    int m; /* present for compatibility issue with amlet - deprecated */
    double *U; /* U is the open knot vector associated to the spline */
    double *N;
  } NT_BSpline;

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
  NT_BSpline *nt_bsplines_new(int p, int nk, double *knots, int inck);

  /**
   * Free the memory associated to a B-spline.
   *
   * @param bspline is the B-spline object to free.
   */
  void nt_bsplines_free(NT_BSpline *bspline);

  /**
   * Evaluation a B-Spline at a particular point.
   *
   * @param b is a space of B-splines.
   * @param P is the vector of control points defining a particular B-spline.
   * @param u is the point where to evaluate the B_spline.
   */
  double nt_bsplines_evaluate(NT_BSpline *b, double *P,
			      double u);

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
  int nt_bsplines_space_dimension(int p, int n);

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
  int nt_bsplines_knots_dimension(int p, int n);

  /**
   * Allocate the memory needed for the knots vector.
   *
   * @param p is the degree of the spline basis functions (order p+1).
   * @param nk is the dimension of the knots vector, with all
   *   multiplicities set to 1.
   *
   * @return a newly allocated vector with the right dimension.
   */
  double *nt_bsplines_knots_allocate(int p, int n);

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
			       int *nU, double **U);

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
				   int realloc);

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
  int nt_bsplines_find_span(int n, int p, double u, double *U);

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
  void nt_bsplines_basis_func(int i, double u, int p,
			      double *U, double *N);

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
  double nt_bsplines_one_basis_func(int p, int m, double *U,
					int i, double u);

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
  double nt_bsplines_curve_point(int n, int p, double u, double *U,
				 double *P, double *N);

  /**
   * The following routines are defined for the moca library but could/should
   * be part of the ntools library.
   */

  /**
   * Initialize the data for a single spline.
   *
   * @param nU is a pointer to the size of the vector x
   (on input and on output).
   * @param p is the degree of the spline basis functions (order p+1).
   * @param U is the knot vector, without repetitions on input if
   *   mulknots is set to 0, with repetitions on output.
   * @param mulknots is a pointer to an integer indicating is x is a
   *   knot vector without repetitions if it is equal to 0, or with
   *   repetitions if it is equal to 1. On ouput the pointer content is
   *   set to 1.
   *
   * @return the b-splines space dimension.
   */
  int moca_spline_init(int *nU, int p, double *U, int *mulknots);

#define MOCA_DEFAULT_SPLINE_DEGREE 3

#ifdef __cplusplus
}
#endif

#endif /* __ORATIO_BSPLINES_H__ */
