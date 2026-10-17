/**
 * gss.h
 *
 * OPHELIA Library
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the BSD License
 */

#ifndef __OPHELIA_GSS_H__
#define __OPHELIA_GSS_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Golden Section Search method.
   *
   * @parameter f is a pointer to the function to minimize with respect to the
   *   first argument x. data is a generic pointer that can be used to pass
   *   any data required for the computation of f. One assumes that f(x) is
   *   continuous over [a,b] and f(x) is unimodal over [a,b].
   * @parameter x is on output the solution of the problem.
   * @paramater data is a pointer to the data to be used when computing f;
   *   data can be set to NULL if no data are required.
   * @parameter a is the interval left side.
   * @parameter b is the interval right side.
   * @parameter tol is a tolerance parameter: the function stops is the
   *   interval width is less than tol. If tol is less or equal to zero, the
   *   function continues until the the maximum number of allowed iterations is
   *   reached, or if the width is set to 0 (due to rounding errors).
   * @parameter imax is the maximum numbers of allowed iterations. If imax < 0,
   *   this criteria is ignored.
   *
   * @return the optimal value, computed at the optimal solution *x.
   */
  double op_gss(double (*f)(double x, void *data),
		double *x, void *data,
		double a, double b, double tol, int imax);

#ifdef __cplusplus
}
#endif

#endif  /* __OPHELIA_GSS_H__ */
