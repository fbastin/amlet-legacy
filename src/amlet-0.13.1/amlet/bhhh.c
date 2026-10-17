/**
 * bhhh.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2004 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "bhhh.h"

#include <amlet/amlet.h>
#include "amlet_config.h"

/* Define the minimum probability acceptable for md->sprob[i] */
#define MinSProb 1e-100

int bhhh_update(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
		const int n, double *BHHH, const int ldbk,
		const double *yk, const int incyk,
		const double *sk, const int incsk,
		double *work)
{
  int i, j, np;

  /* Initialize the BHHH matrix */
  for (i = 0; i < n*n; i++) BHHH[i] = 0.0;

  np = 0;

  /* Update the BHHH matrix */
  for (i = 0; i < data.np; i++) {
    if (!data.n_random_coeff) {
      /* Multinomial logit case */
      cblas_dsyr(Order, Uplo, n, 1.0,
		 data.sgrad[i], 1, BHHH, n);
      np++;
    }
    else {
      /*      nt_vect_print(stdout, "g", data.sgrad[i], n);
	      printf("Prob %g %g\n", data.sprob[i], MinSProb); */
      if (data.sprob[i] >= 0.0001/* MinSProb */) {
	/* Mixed logit case */
	cblas_dsyr(Order, Uplo, n, 1.0/(data.sprob[i]*data.sprob[i]),
		   data.sgrad[i], 1, BHHH, n);
	np++;
      }
    }
  }

  for (i = 1; i < n; i++)
    for (j = 0; j < i; j++) BHHH[i*n+j] = BHHH[j*n+i];

  /* Divide the elements of BHHH matrix */
  for (i = 0; i < n*n; i++) BHHH[i] /= np;

  return 1;
}
