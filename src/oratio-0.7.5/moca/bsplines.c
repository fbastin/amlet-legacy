/**
 * bsplines.c
 *
 * MOCA (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <oratio/ntools.h>

/*
  TOTO: combine spline pointwise and gradient evaluation.
 */

/**
 * Evaluate the spline gradient, with respect to the vector of control points.
 *
 * @param m is the dimension of P, i.e. the dimension of the splines space.
 * @param nU is the dimension of the knot vector.
 * @param p is the degree of the spline basis functions (order p+1)
 * @param u is the point where the gradient is evaluated.
 * @param U is a knot vector, with increasing components, of the form
 *   (2.13) in Piegl and Tiller, p. 66, that is, the first and last
 *   knots have multiplicity p+1, while interior knots have
 *   multiplicity 1.
 * @param P is the control vector.
 * @param N is a double vector of dimension at least equal to p+1.
 * @param p is the degree of the spline basis functions (order p+1)
 */
void moca_bsplines_curve_point_grad(int m,
				    int nU, int p, double u, double U[],
				    double P[], double N[],
					   double g[])
{
  int span, i;

  span = nt_bsplines_find_span(nU, p, u, U);

  nt_bsplines_basis_func(span, u, p, U, N);

  nt_vect_zeros(m, g, 1);
  for (i = 0; i <= p; i++) g[span-p+1] = N[i];
}

/**
 * Initialize the data for a single spline.
 *
 * @param nU is a pointer to the size of the vector x (on input and on output).
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
int moca_spline_init(int *nU, int p, double *U, int *mulknots)
{
  int m;

  nt_bsplines_knots_initialize(p, nU, &U, 1);
  m = nt_bsplines_space_dimension(p, *nU);

  return m;
}
