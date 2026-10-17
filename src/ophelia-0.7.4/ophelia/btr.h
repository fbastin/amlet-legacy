/**
 * btr.h
 *
 * OPHELIA Library
 * Copyright (C) 1999-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __OPHELIA_BTR_H__
#define __OPHELIA_BTR_H__

#include "btr_struct.h"

#ifdef __cplusplus
extern "C" {
#endif

  enum {
    NlpNONE = 0,
    NlpF = 1,
    NlpFval,
    NlpFG,
    NlpFGH,
    NlpMFG
  };

  enum {
    NlpNone = 0,
    NlpComputed,
    NlpExact,
    NlpApprox
  };

  /**
   * BTR implementation. Expert routine
   *
   * @param f is the objective function; if fg is not null, f is not used.
   * @param f_type is the objective function type.
   * @param g is the gradient evaluation function.
   * @param g_type is the gradient function type.
   * @param H is the hessian evaluation function.
   * @param h_type is the hessian function type.
   * @param data: problem data.
   * @param tol: tolerance used in the function evaluation; if tol is
   *   set to 0.0, the machine precision is used.
   * @param output: output file for results and iteration informations
   * @param hessian: preallocated array of size btr->n*btr-btr->n, if
   *   the hessian is needed
   * @param Order: CblasRowMajor or CblasColMajor
   * @param Uplo: CblasUpper or CblasLower
   * @param work: working vector of size btr->n
   *
   * @return the function value at x, who should satisfy some
   *   criticality criteria (depending on the step function)
   */
  double nlp_btrx(BTR *btr, C_GENERIC f, int f_type,
		  C_GENERIC g, int g_type,
		  C_GENERIC H, int h_type,
		  void *data, double tol,
		  FILE *output, double **hessian,
		  const enum CBLAS_ORDER Order,
		  const enum CBLAS_UPLO Uplo,
		  double *work);

  /**
   * BTR implementation. Simple routine.
   */
  double nlp_btr(BTR *btr, C_GENERIC f, C_GENERIC g, void *data,
		 FILE *output, double **hessian, double *work);

#ifdef __cplusplus
}
#endif

#endif   /* __OPHELIA_BTR_H__ */
