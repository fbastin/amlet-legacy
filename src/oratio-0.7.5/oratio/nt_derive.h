/**
 * nt_derive.h
 *
 * ORATIO Library
 * Copyright (C) 2002-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_NTOOLS_DERIVE_H__
#define __ORATIO_NTOOLS_DERIVE_H__

#include <oratio/ctools.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Function type definition for step computation when approximating
   * a derivative.
   * The first parameter corresponds to the point when the derivative
   * will be evaluated, while the second one corresponds to the tolerance.
   */
  typedef double (*NT_H)(double, double);

  /**
   * Derivation in dimension one
   */

  /* First Derivate */

  /**
   * Approximate the first derivate by finite difference.
   *
   * Ref.: Conn, Gould and Toint, "Trust-Region Methods", SIAM, 2000,
   *       pp 297-302.
   *
   * @param f is the function to differentiate.
   * @param x is the point at which the derivate is evaluated.
   * @paran *h is set on output to the value of the used step.
   * @param tol is the numerical precision in the function evaluation;
   *   if tol is set to 0, the machine precision is used.
   * @param scaling is the scaling factor; if no scaling is used in the
   *   function evaluation, the scaling factor should be set to 1.0
   *        dh: routine used to compute the step in one dimension
   *            (NULL if no routine is provided)
   *        data: data needed by function f (NULL if no data is necessary)
   *
   * @return the derivate approximation
   */
  double nt_derive_fd(double (*f)(double, void *),
		      double x, double *h, double tol, double scaling,
		      NT_H dh, void *data);

  /**
   * Approximate the first derivate by finite difference.
   *
   * Ref.: Conn, Gould and Toint, "Trust-Region Methods", SIAM, 2000,
   *       pp 297-302.
   *
   * @param f is the function to differentiate.
   * @param x is the point at which the derivate is evaluated.
   * @param tol is the numerical precision in the function evaluation;
   *   if tol is set to 0, the machine precision is used.
   * @param scaling is the scaling factor; if no scaling is used in the
   *   function evaluation, the scaling factor should be set to 1.0
   *        scaling: scaling factor; if no scaling is used in the
   *                 function evaluation, the scaling factor should be
   *                 set to 1.0
   *        dh: routine used to compute the step in one dimension
   *            (NULL if no routine is provided)
   *        data: data needed by function f (NULL if no data is necessary)
   *
   * #POST: *h: value of the step used
   *        return the derivate approximation
   */
  double nt_derive_fd_with_val(double (*f)(double, void *), double x,
			       double val, double *h, double tol,
			       double scaling, NT_H dh, void *data);

  /**
   * Approximate the first derivate by central difference.
   *
   * Ref.: Conn, Gould and Toint, "Trust-Region Methods", SIAM, 2000,
   *       pp 297-302.
   *
   * @param f is the function to differentiate.
   * @param x is the point where to differentiate the function.
   *        tol: numerical precision in the function evaluation; if tol
   *            is set to 0, the machine precision is used.
   *        scaling: scaling factor; if no scaling is used in the
   *                 function evaluation, the scaling factor should be
   *                 set to 1.0
   *        dh: routine used to compute the step in one dimension
   *            (NULL if no routine is provided)
   *        data: data needed by function f (NULL if no data is necessary)
   *
   * #POST: *h: value of the step used
   *        return the derivate approximation
   */
  double nt_derive_cd(double (*f)(double, void *),
		      double x, double *h, double tol, double scaling,
		      NT_H dh, void *data);

  /**
   * Approximate the first derivate using the Ridder's algorithm.
   *
   * Ref.:
   * 1) Press, Teukolsky, Vetterlin and Flannery, "Numerical
   * Recipes in C (second edition), Cambridge University Press, 1992,
   * pp 186-189
   * 2) Conn, Gould and Toint, "Trust-Region Methods", SIAM, 2000,
   * pp 297-302.
   *
   * @param f is the function to differentiate.
   * @param x is the point where to differentiate the function.
   *        tol: numerical precision in the function evaluation; if tol
   *            is set to 0, the machine precision is used.
   *        scaling: scaling factor; if no scaling is used in the
   *                 function evaluation, the scaling factor should be
   *                 set to 1.0
   *        dh: routine used to compute the step in one dimension
   *            (NULL if no routine is provided)
   *        data: data needed by function f (NULL if no data is necessary)
   *
   * #POST: *h: value of the step used
   *        return the derivative approximation
   */
  double nt_derive_ridder(double (*f)(double, void *),
			  double x, double *h, double *err, double tol,
			  double scaling, NT_H dh, void *data);

  /* Second derivative */

  /**
   * Approximate the second derivate by central difference.
   *
   * Ref.:
   * 1) Press, Teukolsky, Vetterlin and Flannery, "Numerical
   * Recipes in C (second edition), Cambridge University Press, 1992,
   * pp 186-189
   * 2) Conn, Gould and Toint, "Trust-Region Methods", SIAM, 2000,
   * pp 297-302.
   *
   * @param f is the function to differentiate.
   * @param x is the point where to differentiate the function.
   *        tol: numerical precision in the function evaluation; if tol
   *            is set to 0, the machine precision is used.
   *        scaling: scaling factor; if no scaling is used in the
   *                 function evaluation, the scaling factor should be
   *                 set to 1.0
   *        dh: routine used to compute the step in one dimension
   *            (NULL if no routine is provided)
   *        data: data needed by function f (NULL if no data is necessary)
   *
   * #POST: *h: value of the step used
   *        return the derivate approximation
   */
  double nt_derive_2nd_cd(double (*f)(double, void *),
			  double x, double *h, double tol, double scaling,
			  NT_H dh, void *data);

  /**
   * Multidimensional differentiation
   */

  /* First derivative */

  /**
   * Approximate the gradient by finite differences.
   *
   * Ref.:
   * Conn, Gould and Toint, "Trust-Region Methods", SIAM, 2000,
   * pp 297-302.
   *
   * @param f is the function to differentiate.
   * @param x is the point where to differentiate the function.
   *        tol: numerical precision in the function evaluation; if tol
   *             is set to 0, the machine precision is used
   *        scaling: scaling factors; if no scaling is used in the
   *                 function evaluation, the scaling factors should be
   *                 set to 1.0
   *        n: dimension of the problem
   *        dh: routine used to compute the step in one dimension
   *            (NULL if no routine is provided)
   *        data: data needed by function f (NULL if no data is necessary)
   *
   * #POST: h[i] containts the step used in the computation of the ith
   *        partial derivate
   *        gradient points to an array containing the gradient of f
   *        at x
   */
  void nt_derive_grad_fd(double (*f)(double *, int, void *),
			 double *x, double *gradient, double *h,
			 double tol, double *scaling, int n,
			 NT_H dh, double *work, void *data);

  /**
   * Approximate the gradient by finite differences.
   *
   * Ref.:
   * Conn, Gould and Toint, "Trust-Region Methods", SIAM, 2000,
   * pp 297-302.
   *
   * @param f is the function to differentiate.
   * @param x is the point where to differentiate the function.
   *        val: value of f at x
   *        tol: numerical precision in the function evaluation; if tol
   *             is set to 0, the machine precision is used
   *        scaling: scaling factors; if no scaling is used in the
   *                 function evaluation, the scaling factors should be
   *                 set to 1.0
   *        n: dimension of the problem
   *        dh: routine used to compute the step in one dimension
   *            (NULL if no routine is provided)
   *        data: data needed by function f (NULL if no data is necessary)
   *
   * #POST: h[i] containts the step used in the computation of the ith
   *        partial derivate
   *        gradient points to an array containing the gradient of f
   *        at x
   */
  void nt_derive_grad_fd_with_val(double (*f)(double *, int, void *),
				  double *x, double val,
				  double *gradient, double *h,
				  double tol, double *scaling, int n,
				  NT_H dh, double *work, void *data);

  /**
   * Approximate the gradient by central differences.
   *
   * Ref.:
   * Conn, Gould and Toint, "Trust-Region Methods", SIAM, 2000,
   * pp 297-302.
   *
   * @param f is the function to differentiate.
   * @param x is the point where to differentiate the function.
   *        tol: numerical precision in the function evaluation; if tol
   *             is set to 0, the machine precision is used
   *        scaling: scaling factors; if no scaling is used in the
   *                 function evaluation, the scaling factors should be
   *                 set to 1.0
   *        n: dimension of the problem
   *        dh: routine used to compute the step in one dimension
   *            (NULL if no routine is provided)
   * @param work is a working vector of size at least equal to 2*n.
   *        data: data needed by function f (NULL if no data is necessary)
   *
   * #POST: h[i] containts the step used in the computation of the ith
   *        partial derivate
   *        gradient points to an array containing the gradient of f
   *        at x
   */
  void nt_derive_grad_cd(double (*f)(double *, int, void *),
			 double *x, double *gradient, double *h,
			 double tol, double *scaling, int n,
			 NT_H dh, double *work, void *data);

  /**
   * Approximate the gradient with the Ridder's algorithm
   *
   * High cost in function evaluations!
   *
   * Ref.:
   * 1) Press, Teukolsky, Vetterlin and Flannery, "Numerical
   * Recipes in C (second edition), Cambridge University Press, 1992,
   * pp 186-189
   * 2) Conn, Gould and Toint, "Trust-Region Methods", SIAM, 2000,
   * pp 297-302.
   *
   * @param f is the function to differentiate.
   * @param x is the point where to differentiate the function.
   * @param work is a working vector of size at least equal to 2n if not NULL,
   *   otherwise is has to be set to NULL.
   *        tol: numerical precision in the function evaluation; if tol
   *             is set to 0, the machine precision is used
   *        scaling: scaling factors; if no scaling is used in the
   *                 function evaluation, the scaling factors should be
   *                 set to 1.0
   *        n: dimension of the problem
   *        dh: routine used to compute the step in one dimension
   *            (NULL if no routine is provided)
   *        data: data needed by function f (NULL if no data is necessary)
   *
   * #POST: h[i] containts the step used in the computation of the ith
   *        partial derivate
   *        gradient points to an array containing the gradient of f
   *        at x
   */
  void nt_derive_grad_ridder(double (*f)(double *, int, void *),
			     double *x, double *gradient, double *h,
			     double *err, double tol, double *scaling,
			     int n, NT_H dh, double *work, void *data);

  /* Second derivative */

  /**
   * Approximate the Hessian by finite differences on the gradient.
   *
   * Ref.:
   * 1) Press, Teukolsky, Vetterlin and Flannery, "Numerical
   * Recipes in C (second edition), Cambridge University Press, 1992,
   * pp 186-189
   * 2) Conn, Gould and Toint, "Trust-Region Methods", SIAM, 2000,
   * pp 297-302.
   *
   * @param f is the function to differentiate.
   * @param x is the point where to differentiate the function.
   * @param h is a pointer to an integer array that contains on output
   *   the used steps.
   * @param tol is the numerical precision in the function evaluation;
   *   if tol is set to 0, the machine precision is used.
   * @param scaling is a scaling factor; if no scaling is used in the
   *   function evaluation, the scaling factor should be set to 1.0.
   * @param dh is the routine used to compute the step in one
   *   dimension (NULL if no routine is provided).
   * @param work is a working vector of size at least equal to n.
   * @param data corresponds to the data needed by function f (NULL if
   *   no data is necessary).
   *
   * @return the derivative approximation.
   */
  void nt_derive_hess_fd(double (*g)(double *, int, void *,
				   double *grad),
			 double *x, double *grad, double **hessian,
			 double *h, double tol,
			 double *scaling, int n, NT_H dh,
			 double *work, void *data);

  /**
   * @param g is the function to differentiate.
   * @param x is the point where to differentiate the function.
   * @param h is a pointer to an integer array that contains on output
   *   the used steps, if h is not initially set to NULL.
   */
  void nt_derive_hess_fd_gradval(double (*g)(double *, int, void *,
					     double *grad),
				 double *x, double *grad, double **hessian,
				 double *h, double tol,
				 double *scaling, int n, NT_H dh,
				 double *work, void *data);

  /**
   * Approximate the Hessian by central differences.
   *
   * Ref.:
   * Conn, Gould and Toint, "Trust-Region Methods", SIAM, 2000,
   * pp 297-302.
   *
   * @param f is the function to differentiate.
   * @param x is the point where to differentiate the function.
   * @param h is a pointer to an integer array that contains on output
   *   the used steps, if h is not initially set to NULL.
   * @param tol is the numerical precision in the function evaluation; if tol
   *   is set to 0, the machine precision is used
   *        scaling: scaling factors; if no scaling is used in the
   *                 function evaluation, the scaling factors should be
   *                 set to 1.0
   * @param n is the problem dimension.
   *        dh: routine used to compute the step in one dimension
   *            (NULL if no routine is provided)
   * @param data: data needed by function f (NULL if no data is necessary)
   * @param work is a working vector, of size at least equal to 4*n,
   *   otherwise it has to be set to NULL.
   */
  void nt_derive_hess_cd(double (*f)(double *, int, void *),
			 double *x, double **hessian, double *h,
			 double tol, double *scaling, int n,
			 NT_H dh, double *work, void *data);

#ifdef __cplusplus
}
#endif

#endif  /* __NTOOLS_DERIVE_H__ */
