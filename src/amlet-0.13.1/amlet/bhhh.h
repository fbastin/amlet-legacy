/**
 * bhhh.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2004 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __AMLET_BHHH_H__
#define __AMLET_BHHH_H__

#include <oratio/cblas.h>

#ifdef __cplusplus
extern "C" {
#endif

  int bhhh_update(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
		  const int n, double *Bk, const int ldbk,
		  const double *yk, const int incyk,
		  const double *sk, const int incsk,
		  double *work);

#ifdef __cplusplus
}
#endif

#endif

