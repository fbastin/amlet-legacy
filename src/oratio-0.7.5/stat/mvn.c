/**
 * mvn.c
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2009-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * This file is licensed under BSD License.
 */

#include <math.h>

#include "oratio/cblas.h"
#include "oratio/stat.h"

/**
 * Convert a standard multivariate normal observation into an observation
 * of a multivariate normal with given mean and variance-covariance matrix.
 *
 * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
 * @param n is the dimension of the random vector.
 * @param L is a lower triangular matrix such that LL' gives the
 *   variance-covariance matrix.
 * @param ldL is the leading dimension of L.
 * @param mu is the mean vector. If mu is set to NULL, a zero mean vector is
 *   assumed.
 * @param incmu is the increment of mu.
 * @param x is the observation from a standard multivariate normal
 *   distribution on input, and an observation from the desired multivariate
 *   normal on output.
 * @param incx is the increment of x.
 */
void st_mvn_convert(const enum CBLAS_ORDER Order, int n,
			   double *L, int ldL,
			   double *mu, int incmu,
			   double *x, int incx)
{
  /* Apply L. */
  cblas_dtrmv(Order, CblasLower, CblasNoTrans, CblasNonUnit,
	      n, L, ldL, x, incx);

  /* Do the translation. */
  if (mu) cblas_daxpy(n, 1.0, mu, incmu, x, incx);
}

/**
 * Convert an uniform vector into a vector from a multivariate normal with
 * given mean and variance-covariance matrix.
 *
 * @param n is the dimension of the random vector.
 * @param u is a random vector of dimension n whose components are
 *   observations taken from independent U(0,1).
 * @param x is a on output a vector of independant gaussian variates;
 *   x and u can point to the same memory address on input.
 * @param incx is the increment of x.
 */
void st_mvn_standard_icdf(int n, double *u, int incu,
				 double *x, int incx)
{
  int i;

  for (i = 0; i < n; i++) x[i*incx] = st_gaussian_icdf(u[i*incu]);
}


/**
 * Density function of a multivariate standard normal.
 *
 * @param s is the dimension of the random vector.
 * @param x is the point where the density function has to be evaluated.
 * @param incx is the increment of x.
 *
 * @return the value of the density function of a multivariate standard
 *   normal, evaluated at x.
 */
double st_mvn_standard_pdf(int s, double *x, int incx)
{
  int d = s/2;
  int i;

  double den = 1.0;

  for (i = 0; i < d; i++) den *= 2*M_PI;
  if (s % 2) den *= sqrt(2*M_PI);

  return exp(-0.5*cblas_dnrm2(s, x, incx))/den;
}
