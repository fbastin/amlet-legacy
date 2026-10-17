/**
 * pcauchy.h
 *
 * OPHELIA Library
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __OPHELIA_PCAUCHY_H__
#define __OPHELIA_PCAUCHY_H__

#include "btr_struct.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Projected Cauchy step
   *
   * @param btr is am initialized BTR structure.
   * @param Order designs the orientation of the matrix Bk. The
   *   possible values are CblasRowMajor and CblasColMajor.
   * @param Uplo designs the triangular part that has to be updated. The
   *   possible value CblasUpper and CblasLower.
   */
  void btr_proj_cauchy(BTR *btr,
		       const enum CBLAS_ORDER Order,
		       const enum CBLAS_UPLO Uplo,
		       double **hessian, double *work,
		       void *data);

#ifdef __cplusplus
}
#endif

#endif  /* __OPHELIA_PCAUCHY_H__ */
