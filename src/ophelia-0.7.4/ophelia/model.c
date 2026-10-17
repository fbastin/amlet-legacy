/**
 * model.c
 *
 * ATHOS Library
 * Copyright (C) 1999-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "model.h"

#include <oratio/ctools.h>
#include <oratio/ntools.h>

/**
 * Computation of the quadratic model
 */
double model(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	     double *step, double *gradient, double *hessian, int n,
	     void *data)
{
  double value;
  double *h;

  /* Is it necessary to create h? */
  h = _mallocs(n, double);

  /* Quadratic model */

  /* h = hessian*step */
  cblas_dsymv(Order, Uplo, n, 1.0, hessian, n, step, 1, 0.0, h, 1);
  /*  cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, hessian, n,
      step, 1, 0.0, h, 1);  */
  value = cblas_ddot(n, gradient, 1, step, 1) +
    0.5*cblas_ddot(n, step, 1, h,  1);

  free(h);

  return -value;
}

/**
 * Computation of the quadratic model.
 *
 * @param n is the problem dimension.
 * @param work is a double working vector of size at least equal to n.
 * @param hessian is the Hessian, or an approximation of it, of the
 *   objective function.
 */
double btr_quadratic_model(const enum CBLAS_ORDER Order,
				  const enum CBLAS_UPLO Uplo,
				  int n, double *step,
				  double *gradient, double *hessian,
				  void *data, double *work)
{
  double value;

  /* Quadratic model */

  /* work = hessian*step */
  cblas_dsymv(Order, Uplo, n, 1.0, hessian, n, step, 1, 0.0, work, 1);
  value = cblas_ddot(n, gradient, 1, step, 1) +
    0.5*cblas_ddot(n, step, 1, work,  1);

  return -value;
}

/**
 * Computation of the quadratic retrospective model.
 *
 * @param n is the problem dimension.
 * @param work is a double working vector of size at least equal to n.
 * @param hessian is the Hessian, or an approximation of it, of the
 *   objective function.
 */
double btr_quadratic_retro_model(const enum CBLAS_ORDER Order,
					const enum CBLAS_UPLO Uplo,
					int n, double *step,
					double *gradient, double *hessian,
					void *data, double *work)
{
  double value;

  /* Quadratic model */

  /* work = hessian*step */
  cblas_dsymv(Order, Uplo, n, 1.0, hessian, n, step, 1, 0.0, work, 1);
  value = cblas_ddot(n, gradient, 1, step, 1) +
    0.5*cblas_ddot(n, step, 1, work,  1);

  return -value;
}
