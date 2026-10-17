/**
 * bsplines.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2006 Fabian Bastin <bastin@cerfacs.fr>
 *
 * Licensed under the MIT License
 */

#ifndef __AMLET_BSPLINES_H__
#define __AMLET_BSPLINES_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Evaluate the spline gradient, with respect to the vector of
   * control points.
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
				       double N[], double g[], double factor);

#ifdef __cplusplus
}
#endif

#endif /* __AMLET_BSPLINES_H__ */
