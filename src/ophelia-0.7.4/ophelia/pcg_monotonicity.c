/**
 * pcg_monotonicity.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "pcg.h"

#include "activesets.h"
#include "bounds.h"
#include "linesearch.h"
#include "model.h"
#include "step.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

#include <oratio/ctools.h>
#include <oratio/ntools.h>

/**
 * A very simple (monotonicity-)projected truncated
 * conjugate-gradient algorithm.
 *
 * @param Order designs the orientation of the matrix Bk. The
 *   possible values are CblasRowMajor and CblasColMajor.
 * @param Uplo designs the triangular part that has to be updated. The
 *   possible value CblasUpper and CblasLower.
 * @param b is am initialized BTR structure.
 */
void nlp_spcgm_st(BTR *b,
		  const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
		  double **H, double *work, void *data)
{
  int n;           /* problem dimension */
  int nset, nset0; /* number of free variables */
  int *sets;
  int j, k, l;

  double *p, *Hp;
  double *g, *rg, *gplus, *splus, vplus = 0.0;

  double eps_term, ng;
  double epsi = 1e-5, epsi_max = 1e-15;
  double kappa, beta;
  double alpha, aguess;

  double omega = 0.01;

  int i, its, iface;
  /* outside will be set to INTERIOR if the iterate stays inside the
     feasible region, BOUNDARY_TR if it goes on the trust-region
     boundary, or ACTIVE_CSTR if it actives one or more additional
     constraints. */
  int outside;

  int cgits, cgitstart, cgitmax;

  double minsig, bound;
  double value, ngfree, ngfree0, root;

  NLP_Constraint *c = (NLP_Constraint *)(b->constraints->data);

  /* Initialization */

  n = b->n;
  cgitstart = 0;
  cgitmax = n;
  outside = 0;
  aguess = 0;
  value = 0.0;
  alpha = 0.0;

  nset0 = nset = n;

  /* We should use a working vector */

  p = _mallocs(n, double);
  Hp = _mallocs(n, double);
  g = work;
  rg = _mallocs(n, double);
  sets = _mallocs(n+1, int);
  gplus = _mallocs(n, double);
  splus = _mallocs(n, double);
  
  nt_vect_zeros(n, b->step, 1);
  cblas_dcopy(n, b->grad, 1, g, 1);
  cblas_dcopy(n, b->grad, 1, rg, 1);
  cblas_dscal(n, -1.0, rg, 1);

  /* Define the active set. */
  nset = nlp_as_monotonicity(c, n, b->x, 1, rg, 1, sets, 1);

  /*  nlp_as_monotonicity_print(nset, sets); */

  /* Compute the restricted negative gradient. */
  nlp_as_rg_monotonicity(n, p, 1, rg, 1, sets, 1, nset);
  /* 2-norm of the restricted gradient. */
  ng = ngfree0 = nlp_as_rg_nrm2(n, p, 1, nset, sets, 1);  
  printf("Restricted gradient: %g\n", ng);

  /* VALUE OF EPI_TERM? */
  eps_term = nt_max(nt_min(epsi, sqrt(ng))*ng, 0.95*epsi_max);

  cgits = 0;

  /* Loop on successive sets of CG iterations (in different faces). */
  for (iface = 0; iface < cgitmax; iface++) {
#ifdef DEBUG
    printf("IFACE %d\n", iface);
#endif
    /*    nlp_as_monotonicity_print(nset, sets); */
    nset0 = nset;
    /* Reinitialize the first direction to the reduced negative gradient
       in the current face. */
    nlp_as_rg_monotonicity(n, p, 1, g, 1, sets, 1, nset);
    cblas_dscal(n, -1.0, p, 1);

    /* Nothing more to do in this face if the projected gradient (in
       p) is small enough. */
    if (iface > 0) ng = nlp_as_rg_nrm2(n, p, 1, nset, sets, 1);
    if (ng < eps_term) {
      cblas_dcopy(n, g, 1, gplus, 1);
      goto END;
    }

    /* CG loop */
    for (its = cgitstart; its < cgitmax; its++) {
      cgits++;

      cblas_dsymv(Order, Uplo, n, 1.0, *H, n, p, 1, 0.0, Hp, 1);
      kappa = cblas_ddot(n, p, 1, Hp, 1);

      minsig = 1.0e99;
      bound = 0;

      cblas_dcopy(n, b->x, 1, work+2*n, 1);
      cblas_daxpy(n, 1.0, b->step, 1, work+2*n, 1);

      l = 0;

      // PASSAGE A VERIFIER CORRECT? CA A L'AIR !!!!
      for (i = 0; i < nset; i++) {
	j = sets[i];
	if (j >= c->cidx[l]) {
	  k = sets[i+1];
	  if (k < c->cidx[l+1]) {
	    if (p[j]-p[k] > DBL_EPSILON) {
	      bound = ((work+2*n)[k]-(work+2*n)[j])/(p[j]-p[k]);
	      minsig = nt_min(minsig, bound);
	    }
	  }
	  else {
	    l += 2;
	    if (l == 2*c->nidx) break;
	  }
	}
      }

      /* Check if negative curvature is encountered. */
      if (kappa > 0) {
	outside = 0;
	alpha = aguess = -cblas_ddot(n, g, 1, p, 1)/kappa;
	if (minsig < aguess) {
	  outside = 2;
	  aguess = minsig;
	}
	else {
	/* Check the trust-region boundary. */
	  cblas_dcopy(n, b->step, 1, work+2*n, 1);
	  cblas_daxpy(n, alpha, p, 1, work+2*n, 1);
	  if (cblas_ddot(n, work+2*n, 1, work+2*n, 1) >= b->delta*b->delta)
	    outside = 1;
	}
      }
      else {
	/* kappa <= 0 */
	printf("NEGATIVE CURVATURE\n");
	outside = 2;  /* ??? */
	aguess = minsig;
      }

      if (outside) {
	/* Compute the alpha in order to reach the boundary. */
	beta = cblas_ddot(n, b->step, 1, p, 1);
	bound = cblas_ddot(n, p, 1, p, 1);
	root = sqrt(beta*beta+bound*
		    (b->delta*b->delta-cblas_ddot(n, b->step, 1, b->step, 1)));
	beta = (root-beta)/bound;
	if (beta < aguess) {
	  aguess = beta;
	  outside = 1;
	}

	/* Piecewise linesearch loop along the current (projected) CG
	   direction. */
	cblas_daxpy(n, aguess, p, 1, b->step, 1);
	cblas_daxpy(n, aguess, Hp, 1, g, 1); 
	value += vplus;
	/* Update the indexes of free variables */
	if (outside == 2)
	  nlp_as_monotonicity_step(c, n, b->x, 1, b->step, 1, p, 1,
				   &nset, sets, 1);	  
	break;  /* we terminate the CG */
      }
      else {
#ifdef DEBUG
	printf("INSIDE\n");
#endif
	/* Compute the next iterate */
	cblas_daxpy(n, alpha, p, 1, b->step, 1);
	cblas_daxpy(n, alpha, Hp, 1, g, 1);
	value += alpha*cblas_ddot(n, g, 1, p, 1) + 0.5*alpha*alpha*kappa;

	/* Compute the new free gradient and its norm. */
	nlp_as_rg_monotonicity(n, rg, 1, g, 1, sets, 1, nset);
	ngfree = nlp_as_rg_nrm2(n, rg, 1, nset, sets, 1);
#ifdef DEBUG
	printf("Residual norm: %g\n", ngfree);
#endif

	if (cgits >= cgitmax) {
	  /* Too many CG iterations */
	  cblas_dcopy(n, g, 1, gplus, 1);
	  goto END;
	}

	if (ngfree < ngfree0*nt_min(omega, sqrt(ngfree0))) {
	  cblas_dcopy(n, g, 1, gplus, 1);
	  goto END;
	}

	/* CG body: compute new search direction. */
	/* The beta computation has to be modified if a
	   preconditioning is used */
	beta = ngfree/ng;
	beta *= beta;
	cblas_dscal(n, beta, p, 1);
	cblas_daxpy(n, -1.0, rg, 1, p, 1);
	ng = ngfree;
      }

    }

    if (nset == 0 || nset == nset0) break;

    cgitstart = cgits;
  }

  cblas_dcopy(n, g, 1, gplus, 1);

 END:

  /* Memory liberation */
  free(p);
  free(Hp);
  free(gplus);
  free(splus);
  free(sets);
}
