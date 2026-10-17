/**
 * copula.c
 *
 * ORATIO - Random Numbers Module
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/ran_copula.h"
#include "oratio/ran_mvn.h"
#include "oratio/stat.h"

/**
 * Generate a Gaussian copula with uniform marginals.
 *
 * @param ran is an initialized random numbers generator.
 * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
 * @param n is the dimension of the random vector.
 * @param L is a lower triangular matrix such that LL' gives the
 *   correlation  matrix.
 * @param ldL is the leading dimension of L.
 * @param x is on output a realization of the multivariate normal
 *   random vector.
 * @param incx is the increment of x.
 */
void ran_copula_gaussian_unif(Random *ran,
			      const enum CBLAS_ORDER Order,
			      int n, double *L, int ldL,
			      double *x, int incx)
{
  int i;

  ran_mvn(ran, Order, n, L, n, NULL, 1, x, 1);
  for (i = 0; i < n; i++) x[i] = st_gaussian_cdf(x[i]);
}
