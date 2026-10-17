/**
 * derive.c
 *
 * ORATIO Library
 * Copyright (C) 2002-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/nt_derive.h"
#include "oratio/nt_utilities.h"
#include "float.h"

#include "cblas.h"

#include <stdio.h>
#include <string.h>

static double eps = 0.0;  /* still useful ? */

/* Constants for the Ridder's algorithm */
#define NTAB   5
#define REDUCE 1.2
#define ERR    1e100
#define SAFE   2.0

/**
 * Derivation in dimension 1
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
		    NT_H dh, void *data)
{
  return nt_derive_fd_with_val(f, x, (*f)(x, data), h, tol, scaling, dh, data);
}

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
			     double scaling, NT_H dh, void *data)
{
  if (!tol) {
    if (!eps)
      eps = DBL_EPSILON;
    tol = eps;
  }

  if (scaling <= 0.0) scaling = 1.0;

  if (dh)
    *h = dh(x, tol);
  else {
    tol = sqrt(tol);
    *h = nt_sign(x)*tol*nt_max(fabs(x), scaling);
  }
  *h = (x+(*h))-x;

  return ((*f)(x+(*h), data)-val)/(*h);
}

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
double nt_derive_cd(double (*f)(double, void *), double x,
		    double *h, double tol,
		    double scaling, NT_H dh, void *data)
{
  if (!tol) {
    if (!eps)
      eps = DBL_EPSILON;
    tol = eps;
  }

  if (dh)
    *h = dh(x, tol);
  else {
    tol = pow(tol, 1.0/3.0);
    *h = nt_sign(x)*tol*nt_max(fabs(x), 1);
  }
  *h = (x+(*h))-x;

  return ((*f)(x+(*h), data)-(*f)(x-(*h), data))/(2*(*h));
}

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
			double scaling, NT_H dh, void *data)
{
  double array[NTAB][NTAB];
  double mul, hh, fact, reducesqr, err_temp, df = 0;
  int i, j;

  if (!tol) {
    if (!eps)
      eps = DBL_EPSILON;
    tol = eps;
  }

  mul = pow(REDUCE, NTAB/1.5);
  *err = ERR;
  reducesqr = REDUCE*REDUCE;

  if (scaling <= 0.0) scaling = 1.0;

  if (dh)
    *h = dh(x, tol);
  else {
    tol = pow(tol, 1.0/3.0);
    *h = nt_sign(x)*tol*nt_max(fabs(x), scaling);
  }
  *h *= mul;
  *h = (x+(*h))-x;
  hh = *h;

  array[0][0] = ((*f)(x+hh, data)-(*f)(x-hh, data))/(2.0*(hh));

  for (i = 1; i < NTAB; i++) {
    hh /= REDUCE;
    /* Initialise la première ligne à l'aide de différences centrales */
    array[0][i] = ((*f)(x+hh, data)-(*f)(x-hh, data))/(2.0*hh);
    fact = reducesqr;
    /* Les colonnes successives dans le tableau prennent des pas de plus en
     * plus petits, avec des ordres d'extrapolation plus grands. */
    for (j = 1; j <= i; j++) {
      array[j][i] = (array[j-1][i]*fact-array[j-1][i-1])/(fact-1.0);
      fact *= reducesqr;
      err_temp = nt_max(array[j][i]-array[j-1][i], array[j][i]-array[j-1][i-1]);
      if (err_temp < *err) {
	*err = err_temp;
	df = array[j][i];
	*h = hh;
      }
    }
    /* Si l'ordre supérieur est pire par un facteur significatif SAFE,
     * alors sort plus rapidement. */
    if (fabs(array[i][i]-array[i-1][i-1]) >= SAFE*(*err)) break;
  }

  return df;
}

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
			NT_H dh, void *data)
{
  if (!tol) {
    if (!eps)
      eps = DBL_EPSILON;
    tol = eps;
  }

  if (scaling <= 0.0) scaling = 1.0;

  if (dh)
    *h = dh(x, tol);
  else {
    tol = pow(tol, 0.25);
    *h = nt_sign(x)*tol*nt_max(fabs(x), scaling);
  }
  *h = (x+(*h))-x;

  return ((*f)(x+(*h), data)-2*(*f)(x, data)+(*f)(x-(*h), data))/((*h)*(*h));
}

/**
 * Multidimensional differentiation.
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
		       NT_H dh, double *work, void *data)
{
  double *xplus, *h0;
  double hh;

  double val;

  int i;

  xplus = (work) ? work : _mallocs(n, double);
  h0 = (h) ? h : _mallocs(n, double);

  if (!tol) {
    if (!eps)
      eps = DBL_EPSILON;
    tol = eps;
  }

  if (dh)
    for (i = 0; i < n; i++) h0[i] = dh(x[i], tol);
  else {
    tol = sqrt(tol);
    for (i = 0; i < n; i++) {
      if (scaling [i] <= 0.0) scaling[i] = 1.0;
      h0[i] = nt_sign(x[i])*tol*nt_max(fabs(x[i]), scaling[i]);
    }
  }

  /* slightly improve the precision */
  for (i = 0; i < n; i++) h0[i] = (x[i]+h0[i])-x[i];

  memcpy(xplus, x, n*sizeof(double));

  val = (*f)(x, n, data);

  for (i = 0; i < n; i++) {
    hh = h0[i];
    xplus[i] += hh;
    gradient[i] = ((*f)(xplus, n, data) - val)/hh;
    xplus[i] = x[i];
  }

  if (!work) free(xplus);
  if (!h) free(h0);
}

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
				NT_H dh, double *work, void *data)
{
  double *xplus, *h0;
  double hh;

  int i;

  xplus = (work) ? work : _mallocs(n, double);
  h0 = (h) ? h : _mallocs(n, double);

  if (!tol) {
    if (!eps)
      eps = DBL_EPSILON;
    tol = eps;
  }

  if (dh)
    for (i = 0; i < n; i++) h0[i] = dh(x[i], tol);
  else {
    tol = sqrt(tol);
    for (i = 0; i < n; i++) {
      if (scaling [i] <= 0.0) scaling[i] = 1.0;
      h0[i] = nt_sign(x[i])*tol*nt_max(fabs(x[i]), scaling[i]);
    }
  }

  /* slightly improve the precision */
  for (i = 0; i < n; i++) h0[i] = (x[i]+h0[i])-x[i];

  memcpy(xplus, x, n*sizeof(double));

  for (i = 0; i < n; i++) {
    hh = h0[i];
    xplus[i] += hh;
    gradient[i] = ((*f)(xplus, n, data) - val)/hh;
    xplus[i] = x[i];
  }

  if (!work) free(xplus);
  if (!h) free(h0);
}

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
		       NT_H dh, double *work, void *data)
{
  double *xplus, *xmin, *h0;
  double hh;
  int i;

  h0 = (h) ? h : _mallocs(n, double);
  xplus = (work) ? work : _mallocs(2*n, double);
  xmin = xplus+n;

  if (!tol) {
    if (!eps)
      eps = DBL_EPSILON;
    tol = eps;
  }

  if (dh)
    for (i = 0; i < n; i++) h0[i] = dh(x[i], tol);
  else {
    tol = pow(tol, 1.0/3.0);
    for (i = 0; i < n; i++) {
      if (scaling [i] <= 0.0) scaling[i] = 1.0;
      h0[i] = nt_sign(x[i])*tol*nt_max(fabs(x[i]), scaling[i]);
    }
  }

  /* slightly improve the precision */
  for (i = 0; i < n; i++) h0[i] = (x[i]+h0[i])-x[i];

  memcpy(xmin, x, n*sizeof(double));
  memcpy(xplus, x, n*sizeof(double));

  for (i = 0; i < n; i++) {
    hh = h0[i];
    xmin[i] -= hh;
    xplus[i] += hh;

    gradient[i] = ((*f)(xplus, n, data) - (*f)(xmin, n, data))/(2.0*hh);

    xmin[i] = xplus[i] = x[i];
  }

  if (!work) free(xplus);
  if (!h) free(h0);
}

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
			   int n, NT_H dh, double *work, void *data)
{
  double array[NTAB][NTAB], *xplus, *xmin, *h0;
  double mul, hh, fact, reducesqr, err_temp;
  int i, j, dim;

  xplus = (work) ? work : _mallocs(2*n, double);
  h0 = (h) ? h : _mallocs(n, double);
  xmin = xplus+n;

  if (!tol) {
    if (!eps)
      eps = DBL_EPSILON;
    tol = eps;
  }

  mul = pow(REDUCE, NTAB/1.5);
  reducesqr = REDUCE*REDUCE;

  if (dh)
    for (i = 0; i < n; i++) h0[i] = dh(x[i], tol);
  else {
    tol = pow(tol, 1.0/3.0);
    for (i = 0; i < n; i++) h0[i] = nt_sign(x[i])*tol*nt_max(fabs(x[i]), 1);
  }

  for (i = 0; i < n; i++) {
    h0[i] *= mul;
    h0[i] = (x[i]+h0[i])-x[i];
  }

  for (dim = 0; dim < n; dim++) {
    *err = ERR;
    hh = h0[dim];
    memcpy(xplus, x, n*sizeof(double));
    memcpy(xmin, x, n*sizeof(double));
    xplus[dim] += hh;
    xmin[dim] -= hh;
    array[0][0] = ((*f)(xplus, n, data)-(*f)(xmin, n, data))/(2.0*(hh));

    for (i = 1; i < NTAB; i++) {
      xplus[dim] -= hh;
      xmin[dim] += hh;
      hh /= REDUCE;
      xplus[dim] += hh;
      xmin[dim] -= hh;
      /* Initialise la première ligne à l'aide de différences centrales */
      array[0][i] = ((*f)(xplus, n, data)-(*f)(xmin, n, data))/(2.0*hh);
      fact = reducesqr;
      /* Les colonnes successives dans le tableau prennent des pas de plus en
       * plus petits, avec des ordres d'extrapolation plus grands. */
      for (j = 1; j <= i; j++) {
	array[j][i] = (array[j-1][i]*fact-array[j-1][i-1])/(fact-1.0);
	fact *= reducesqr;
	err_temp = nt_max(array[j][i]-array[j-1][i], array[j][i]-array[j-1][i-1]);
	if (err_temp < *err) {
	  *err = err_temp;
	  gradient[dim] = array[j][i];
	  h0[dim] = hh;
	}
      }
      /* Si l'ordre supérieur est pire par un facteur significatif SAFE,
       * alors sort plus rapidement. */
      if (fabs(array[i][i]-array[i-1][i-1]) >= SAFE*(*err)) break;
    }
  }

  if (!work) free(xplus);
  if (!h) free(h0);
}

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
		       double *h, double tol, double *scaling,
		       int n, NT_H dh, double *work, void *data)
{
  int i, j;
  double *u, *h0;

  u = (work) ? work : _mallocs(n, double);
  h0 = (h) ? h : _mallocs(n, double);

  if (!tol) {
    if (!eps)
      eps = DBL_EPSILON;
    tol = eps;
  }

  if (dh)
    for (i = 0; i < n; i++) h0[i] = dh(x[i], tol);
  else {
    tol = sqrt(tol);
    for (i = 0; i < n; i++) {
      if (scaling [i] <= 0.0) scaling[i] = 1.0;
      h0[i] = nt_sign(x[i])*tol*nt_max(fabs(x[i]), scaling[i]);
    }
  }

  /* slightly improve the precision */
  for (i = 0; i < n; i++) h0[i] = (x[i]+h0[i])-x[i];

  memcpy(u, x, n*sizeof(double));

  /* compute the gradient at x */
  g(x, n, data, grad);

  for (i = 0; i < n; i++) {
    u[i] += h0[i];
    g(u, n, data, hessian[i]);

    cblas_daxpy(n, -1.0, grad, 1, hessian[i], 1);
    cblas_dscal(n, 1/h0[i], hessian[i], 1);

    u[i] -= h0[i];
  }

  /* regularization of the result */
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      hessian[j][i] = hessian[i][j] = (hessian[i][j]+hessian[j][i])/2.0;

  if (!work) free(u);
  if (!h) free(h0);
}

/**
 * @param g is the function to differentiate.
 * @param x is the point where to differentiate the function.
 * @param h is a pointer to an integer array that contains on output
 *   the used steps, if h is not initially set to NULL.
 */
void nt_derive_hess_fd_gradval(double (*g)(double *, int, void *,
					   double *grad),
			       double *x, double *grad, double **hessian,
			       double *h, double tol, double *scaling, int n,
			       NT_H dh, double *work, void *data)
{
  int i, j;
  double *u, *h0;

  u = (work) ? work : _mallocs(n, double);
  h0 = (h) ? h : _mallocs(n, double);

  if (!tol) {
    if (!eps)
      eps = DBL_EPSILON;
    tol = eps;
  }

  if (dh)
    for (i = 0; i < n; i++) h0[i] = dh(x[i], tol);
  else {
    tol = sqrt(tol);

    for (i = 0; i < n; i++) {
      if (scaling[i] <= 0.0) scaling[i] = 1.0;
      h0[i] = nt_sign(x[i])*tol*nt_max(fabs(x[i]), scaling[i]);
    }
  }

  /* slightly improve the precision */
  for (i = 0; i < n; i++) h0[i] = (x[i]+h0[i])-x[i];

  memcpy(u, x, n*sizeof(double));

  for (i = 0; i < n; i++) {
    u[i] += h0[i];
    g(u, n, data, hessian[i]);
    cblas_daxpy(n, -1.0, grad, 1, hessian[i], 1);
    cblas_dscal(n, 1/h0[i], hessian[i], 1);
    u[i] -= h0[i];
  }

  /* regularization of the result */
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      hessian[j][i] = hessian[i][j] = (hessian[i][j]+hessian[j][i])/2.0;

  if (!work) free(u);
  if (!h) free(h0);
}

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
		       NT_H dh, double *work, void *data)
{
  /* TODO: implement a better stockage of the Hessian in order to take
     into account of the symmetry */

  int i, j;
  double *l, *u, *v, *w, *h0 = NULL;

  double val;

  h0 = (h) ? h : _mallocs(n, double);
  l = (work) ? work : _mallocs(4*n, double);
  u = l+n;
  v = l+2*n;
  w = l+3*n;

  if (!tol) {
    if (!eps)
      eps = DBL_EPSILON;
    tol = eps;
  }

  if (dh)
    for (i = 0; i < n; i++) h0[i] = dh(x[i], tol);
  else {
    tol = pow(tol, 0.25);
    for (i = 0; i < n; i++)
      if (scaling[i] <= 0.0) scaling[i] = 1.0;
    for (i = 0; i < n; i++)
      h0[i] = nt_sign(x[i])*tol*nt_max(fabs(x[i]), scaling[i]);
  }
  for (i = 0; i < n; i++) h0[i] = (x[i]+h0[i])-x[i];

  /* initialisation des vecteurs l, u, v, w */
  memcpy(l, x, n*sizeof(double));
  memcpy(u, x, n*sizeof(double));
  memcpy(v, x, n*sizeof(double));
  memcpy(w, x, n*sizeof(double));

  val = (*f)(u, n, data);

  for (i = 0; i < n; i++) {

    /* Treat separately the case where i = j, in order to save one
       function evaluation */
    l[i] += h0[i];
    w[i] -= h0[i];

    hessian[i][i] = ((*f)(l, n, data) - 2*val +
		     (*f)(w, n, data))/(h0[i]*h0[i]);

    u[i] += h0[i];
    v[i] -= h0[i];

    for (j = i+1; j < n; j++) {
      l[j] += h0[j];
      u[j] -= h0[j];
      v[j] += h0[j];
      w[j] -= h0[j];

      hessian[i][j] = ((*f)(l, n, data) - (*f)(u, n, data)-
		       (*f)(v, n, data) + (*f)(w, n, data))/(4*h0[i]*h0[j]);
      if (i != j) hessian[j][i] = hessian[i][j];

      l[j] -= h0[j];
      u[j] += h0[j];
      v[j] -= h0[j];
      w[j] += h0[j];
    }
    l[i] -= h0[i];
    u[i] -= h0[i];
    v[i] += h0[i];
    w[i] += h0[i];

  }

  if (!work) free(l);
  if (!h) free(h0);
}
