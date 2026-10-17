/**
 * common_struct.h
 *
 * OPHELIA Library
 * Copyright (C) 1999-2003 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __OPHELIA_COMMON_STRUCT_H__
#define __OPHELIA_COMMON_STRUCT_H__

#include <stdio.h>
#include <sys/types.h>

#include <oratio/cblas.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Derivation functions types
   */
  typedef double (*deriv1)(double *x, int n, void *data, double *grad);
  typedef void (*deriv2)(double *x, int n, void *data, double **hessian);
  typedef int (*deriv2_approx)(const enum CBLAS_ORDER Order,
			       const enum CBLAS_UPLO Uplo,
			       const int n, double *Hk, const int ldhk,
			       const double *yk, const int incyk,
			       const double *sk, const int incsk,
			       double *work);

#ifdef __cplusplus
}
#endif

#endif /* __OPHELIA_COMMON_STRUCT_H__ */
