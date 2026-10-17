/**
 * amlet_stopping.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2005-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "amlet_stopping.h"

#include <math.h>
#include <oratio/ntools.h>

#define RATIO_GRADIENT  0.1

int amlet_btrda_stopping_test(BTR *btr, double *step, int nlast,
			      dcl_data *md, double accuracy,
			      double rfct, double rho)
{
  double tmp = 0.0;
  double sqrtol = btr->tol*btr->tol;
  double relgr;

  if (md->size_sampling == md->max_size || accuracy*accuracy <= sqrtol) {
    /* gradient norm */
    relgr = nlp_relative_gradient(btr->value, btr->x, btr->grad,
				  NULL, btr->n, 1.0);
    //    if (verbose) printf("Relative gradient: %g\n", relgr);
    /*    g = cblas_ddot(btr->n, btr->grad, 1, btr->grad, 1); */
    tmp = RATIO_GRADIENT*accuracy;
    if (cblas_dnrm2(btr->n, btr->grad, 1) < tmp) {
      //      if (verbose) printf("Stop due to gradient\n");
      return NLP_STOPPING_GRADIENT;
    }
    if (relgr < tmp) {
      //      if (verbose) printf("Stop due to relative gradient\n");
      return NLP_STOPPING_RELATIVE_GRADIENT;
    }
    if (btr->delta < btr->tol) {
      /*      if (verbose)
	printf("Stop due to too small trust-region radius: %g\n",
	btr->delta); */
      return NLP_STOPPING_DELTA;
    }
    if (rho > btr->eta1 && cblas_dnrm2(btr->n, btr->step, 1) < 0.01*btr->tol) {
      /*      if (verbose)
	printf("Stop due to too small step: %g %g\n",
	cblas_dnrm2(btr->n, btr->step, 1), btr->tol); */
      return 1;
    }
  }

  return 0;
}
