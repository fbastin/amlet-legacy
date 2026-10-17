/**
 * psb.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2003 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the MIT License
 */

/**
 * Powell-symmetric-Broyden approximation
 */

/* NOT COMPLETED!!! */

#include <oratio/cblas.h>

#include <math.h>

int psb_update(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	       const int n, double *Bk, const int ldbk,
	       const double *yk, const int incyk,
	       const double *sk, const int incsk,
	       double *work)
{
  double dot;
  double dot2;

  cblas_dcopy(n, yk, incyk, work, 1);

  /* Compute y_k-B_ks_k */
  cblas_dsymv(Order, Uplo, n, -1.0, Bk, ldbk, sk, incsk, 1.0, work, 1);

  /* Compute s_k^T(y_k-B_ks_k) */
  dot = cblas_ddot(n, sk, incsk, work, 1);

  /* Compute s_k^Ts_k */
  dot2 = cblas_ddot(n, sk, incsk, sk, incsk);

  /* Compute B_k+dot/dot2^2*sk*s_k^T */
  cblas_dsyr(Order, Uplo, n, dot/(dot2*dot2), sk, incsk, Bk, ldbk);

  /*  cblas_dsyr2(Order, Uplo, n, 1.0/dot2, sk, incsk, Bk, ldbk); */

  return 0;
}
