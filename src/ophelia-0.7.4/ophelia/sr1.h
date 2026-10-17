/**
 * sr1.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2003 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

/**
 * SR-1 Hessian approximation
 */

#ifndef __NLP2_SR1_H__
#define __NLP2_SR1_H__

#include <oratio/cblas.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * SR1 Hessian Update
   *
   * @param Order designs the orientation of the matrix Bk. The
   *   possible values are CblasRowMajor and CblasColMajor.
   * @param Uplo designs the triangular part that has to be updated. The
   *   possible value CblasUpper and CblasLower.
   * @param n is the number of variables.
   * @param yk is the gradients difference.
   * @param sk is the step.
   *
   * @return 1 if the approximation Bk has been updated, 0 else.
   */
  int sr1_update(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
		 const int n, double *Bk, const int ldbk,
		 const double *yk, const int incyk,
		 const double *sk, const int incsk,
		 double *work);


  int nlp_sr1_update(const enum CBLAS_ORDER Order,
		     const enum CBLAS_UPLO Uplo,
		     int n, int m, double *H,
		     double *values, double *valuesplus,
		     double *grad, double *gradplus, double *step,
		     int *info, double *work, void *data);

#ifdef __cplusplus
}
#endif

#endif  /* __NLP2_SR1_H__ */
