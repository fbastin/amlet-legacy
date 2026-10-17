/**
 * btr_step.h
 *
 * OPHELIA Library
 * Copyright (C) 1999-2004 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __OPHELIA_BTR_STEP_H__
#define __OPHELIA_BTR_STEP_H__

#include "btr_struct.h"

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Constants for step computation techniques.
   */
  enum {
    PERSONAL = 0,
    /* No constraints */
    BTR_STEP_CAUCHY,
    BTR_STEP_STEIHAUG_TOINT,
    /* Bounds constraints */
    BTR_STEP_RESTRICTED_STEIHAUG_TOINT,
    /* With constraints, projected cauchy step */
    BTR_STEP_PROJ_CAUCHY,
    /* Simple projected conjugated gradient step */
    BTR_STEP_SPCG_ST
  };

  enum {
    INTERIOR,
    BOUNDARY_TR,
    ACTIVE_CSTR
  };

  /**
   * Step functions
   */

  /**
   * Allow to attach a predefined step computation function.
   *
   * @param btr is a initialized BTR structure.
   * @param step_type is a step identifier constant.
   */
  void btr_set_step(BTR *btr, int step_type);

  /**
   * Return a short description of the step method used in the
   * trust-region algorithm.
   *
   * @param btr is a initialized BTR structure.
   *
   * @return the name of the step technique used in btr.
   */
  char *btr_get_step_name(BTR *btr);

  /**
   * Cauchy step implementation
   *
   * Ref.: Jorge Nocedal and Stephen J. Wright, "Numerical Optimization",
   *   Springer, pp 69-70, 1999
   */
  void btr_step_cauchy(BTR *b,
		       const enum CBLAS_ORDER Order,
		       const enum CBLAS_UPLO Uplo,
		       double **hessian, double *work,
		       void *data);

  /**
   * Truncated conjugated gradient algorithm implementation
   *
   * Ref.: Andrew R. Conn, Nicholas I. M. Gould, Philippe L. Toint,
   *   "Trust-Region Methods", SIAM, pp 202-207, 2000.
   */
  void btr_step_steihaug_toint(BTR *btr,
			       const enum CBLAS_ORDER Order,
			       const enum CBLAS_UPLO Uplo,
			       double **hessian, double *work,
			       void *data);

  /**
   * Projected Cauchy step
   */
  void btr_proj_cauchy(BTR *btr,
		       const enum CBLAS_ORDER Order,
		       const enum CBLAS_UPLO Uplo,
		       double **hessian, double *work,
		       void *data);

  void btr_step_steihaug_toint_bounds(BTR *btr,
				      const enum CBLAS_ORDER Order,
				      const enum CBLAS_UPLO Uplo,
				      double **hessian, double *work,
				      void *data);

#ifdef __cplusplus
}
#endif

#endif  /* __OPHELIA_BTR_STEP_H__ */
