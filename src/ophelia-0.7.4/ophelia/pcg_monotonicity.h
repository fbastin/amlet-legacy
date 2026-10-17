/**
 * pcg_monotonicity.h
 *
 * ATHOS Library
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ATHOS_PCGM_H__
#define __ATHOS_PCGM_H__

#include "btr_struct.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * A very simple (monotonicity-)projected truncated
   * conjugate-gradient algorithm.
   *
   * @param b is am initialized BTR structure.
   * @param Order designs the orientation of the matrix Bk. The
   *   possible values are CblasRowMajor and CblasColMajor.
   * @param Uplo designs the triangular part that has to be updated. The
   *   possible value CblasUpper and CblasLower.
   */
  void nlp_spcgm_st(BTR *b,
		    const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
		    double **H, double *work, void *data);

#ifdef __cplusplus
}
#endif

#endif  /* __ATHOS_PCGM_H__ */
