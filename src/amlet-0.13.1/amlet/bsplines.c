/**
 * bsplines.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/ntools.h>

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
 * @param N is a double vector of dimension at least equal to p+1.
 * @param p is the degree of the spline basis functions (order p+1)
 */
void amlet_bsplines_curve_point_grad(int m, int nU, int p,
				     double u, double U[],
				     double N[], double g[], double factor)
{
  int span, i;

  span = nt_bsplines_find_span(nU, p, u, U);

  nt_bsplines_basis_func(span, u, p, U, N);

  for (i = 0; i <= p; i++) g[span-p+i] += N[i]*factor;
}

