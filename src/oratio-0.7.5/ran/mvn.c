/**
 * mvn.c
 *
 * ORATIO Library
 * Copyright (C) 2009-2012 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/cblas.h"
#include "oratio/ran.h"
#include "oratio/stat.h"

/**
 * Multivariate normal generation.
 */

/**
 * Draw from a multivariate Gaussian distribution.
 *
 * @param ran is an initialized random numbers generator.
 * @param n is the dimension of the random vector.
 * @param x is on output a realization of the multivariate Gaussian
 *   random vector.
 * @param incx is the increment of x.
 */
void ran_mvn_gaussian(Random *ran, int n, double *x, int incx)
{
  int i;

  for (i = 0; i < n; i++) x[i*incx] = ran_gaussian_icdf(ran);
}

/**
 * Draw from a multivariate Gaussian distribution (Box-Muller variant).
 *
 * @param ran is an initialized random numbers generator.
 * @param n is the dimension of the random vector.
 * @param x is on output a realization of the multivariate Gaussian
 *   random vector.
 * @param incx is the increment of x.
 */
void ran_mvn_gaussian_box_muller(Random *ran, int n,
				 double *x, int incx)
{
  int i;

  for (i = 0; i < n; i++) x[i*incx] = ran_gaussian_box_muller(ran);
}

/**
 * Draw from a multivariate normal distribution.
 *
 * @param ran is an initialized random numbers generator.
 * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
 * @param n is the dimension of the random vector.
 * @param L is a lower triangular matrix such that LL' gives the
 *   variance-covariance matrix.
 * @param ldL is the leading dimension of L.
 * @param mu is the mean vector. If mu is set to NULL, a zero mena vector is
 *   assumed.
 * @param incmu is the increment of mu.
 * @param x is on output a realization of the multivariate normal
 *   random vector.
 * @param incx is the increment of x.
 */
void ran_mvn(Random *ran,
	     const enum CBLAS_ORDER Order, int n, double *L, int ldL,
	     double *mu, int incmu, double *x, int incx)
{
  ran_mvn_gaussian(ran, n, x, incx);
  st_mvn_convert(Order, n, L, ldL, mu, incmu, x, incx);
}
