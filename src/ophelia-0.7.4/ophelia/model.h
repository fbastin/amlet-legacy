/**
 * model.h
 *
 * ATHOR Library
 * Copyright (C) 1999-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __NLP2_MODEL_H__
#define __NLP2_MODEL_H__

#include <oratio/cblas.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Computation of the quadratic model
   */
  double model(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	       double *step, double *gradient, double *hessian, int n,
	       void *data);

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
			     void *data, double *work);

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
				   void *data, double *work);

#ifdef __cplusplus
}
#endif

#endif  /* __ATHOR_MODEL_H__ */
