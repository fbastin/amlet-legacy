/**
 * pcg.h
 *
 * ATHOS Library
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ATHOS_PCG_H__
#define __ATHOS_PCG_H__

#include "btr_struct.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * A very simple (bound-)projected truncated conjugate-gradient algorithm.
   *
   * @param Order designs the orientation of the matrix Bk. The
   *   possible values are CblasRowMajor and CblasColMajor.
   * @param Uplo designs the triangular part that has to be updated. The
   *   possible value CblasUpper and CblasLower.
   * @param b is am initialized BTR structure.
   */
  void nlp_spcg_st(BTR *b,
		   const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
		   double **H, double *work, void *data);

#ifdef __cplusplus
}
#endif

#endif  /* __ATHOS_PCG_H__ */
