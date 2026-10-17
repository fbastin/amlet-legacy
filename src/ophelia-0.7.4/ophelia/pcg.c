/**
 * pcg.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "pcg.h"

#include "activesets.h"
#include "bounds.h"
#include "linesearch.h"
#include "model.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

#include <oratio/ctools.h>
#include <oratio/ntools.h>

/**
 * A very simple (bound-)projected truncated conjugate-gradient algorithm.
 *
 * @param b is am initialized BTR structure.
 * @param Order designs the orientation of the matrix Bk. The
 *   possible values are CblasRowMajor and CblasColMajor.
 * @param Uplo designs the triangular part that has to be updated. The
 *   possible value CblasUpper and CblasLower.
 */
void nlp_spcg_st(BTR *b,
		 const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
		 double **H, double *work, void *data)
{
  int n;             /* problem dimension */
  int nfree, nfree0; /* number of free variables */
  int *xfree;        /* 0-1 vector indicating which variables are
		        free variables */
  double *p, *pr, *sum;
  double *g, *gplus, *splus, vplus = 0.0;

  double eps_term, ng;
  double epsi = 1e-5, epsi_max = 1e-15;
  double kappa, beta;
  double alpha, aguess;

  NlpBounds *bounds = (NlpBounds *)(b->projection_data);

  int i, its, iface;
  int interior, outside;

  int cgits, cgitstart, cgitmax;

  double minsig, maxsig, bound;
  double value, ngfree, root;

  /* Memory allocation */

  n = b->n;
  cgitstart = 0;
  cgitmax = n;
  outside = 0;
  aguess = 0;
  value = 0.0;
  alpha = 0.0;

  nfree0 = nfree = n;

  /* We should use a working vector */

  xfree = _mallocs(n, int);
  p = _mallocs(n, double);
  pr = _mallocs(n, double);
  g = work;
  gplus = _mallocs(n, double);
  splus = _mallocs(n, double);
  sum = _mallocs(n, double);

  /* We set all elements of xfree to 1 */;
  nlp_as_free_init_dflt(n, xfree, 1);

  nt_vect_zeros(n, b->step, 1);
  cblas_dcopy(n, b->grad, 1, g, 1);

  /* Define the active set. */
  nfree = nlp_as_bounds(b->projection_data, n, b->x, 1, g, 1, xfree, 1);

  /* Compute the restricted gradient. */
  nlp_as_rg_bounds(n, p, 1, g, 1, xfree, 1, nfree);
  /* 2-norm of the restricted gradient. */
  ng = cblas_dnrm2(n, p, 1);  

  eps_term = nt_max(nt_min(epsi, sqrt(ng))*ng, 0.95*epsi_max);

  interior = 1;
  cgits = 0;

  /* Loop on successive sets of CG iterations (in different faces). */
  for (iface = 0; iface < cgitmax; iface++) {
#ifdef DEBUG
    printf("IFACE %d\n", iface);
#endif
    nfree0 = nfree;
    cblas_dcopy(n, b->grad, 1, g, 1);
    /* Reinitialize the first direction to the free negative gradient
       in the current face. */
    nlp_as_rg_bounds(n, p, 1, g, 1, xfree, 1, nfree);
    cblas_dscal(n, -1.0, p, 1);

    /* Nothing more to do in this face if the projected gradient (in
       p) is small enough. */
    if (iface > 0) ng = cblas_dnrm2(n, p, 1);
    if (ng < eps_term) {
      cblas_dcopy(n, g, 1, gplus, 1);
      goto END;
    }

    /* CG loop */
    for (its = cgitstart; its < cgitmax; its++) {
      cgits++;

      cblas_dsymv(Order, Uplo, n, 1.0, *H, n, p, 1, 0.0, pr, 1);
      kappa = cblas_ddot(n, p, 1, pr, 1);
      /*      nt_vect_print(stdout, "p", p, n);
	      nt_vect_print(stdout, "pr", pr, n); */

      /* Find the step to the first bound - similar to breakpoints */
      maxsig = 0;
      minsig = 1.0e99;

      for (i = 0; i < n; i++) {
	if (xfree[i]) {
	  bound = 0;  /* should not remain equal to 0 */
	  if (p[i] > 0)
	    bound = (bounds->xu[i]-b->x[i]-b->step[i])/p[i];
	  else if (p[i] < 0)
	    bound = (bounds->xl[i]-b->x[i]-b->step[i])/p[i];
	  minsig = nt_min(minsig, bound);
	  maxsig = nt_max(bound, maxsig);
	}
      }

      if (kappa > 0) {
	alpha = aguess = -cblas_ddot(n, g, 1, p, 1)/kappa;
	outside = 0;
	/* Check the bounds. */
	if (minsig < alpha) {
	  outside = 1;
	  aguess = minsig;
	}
	/* Check the trust-region boundary. */
	cblas_dcopy(n, b->step, 1, sum, 1);
	cblas_daxpy(n, alpha, p, 1, sum, 1);
	if (cblas_ddot(n, sum, 1, sum, 1) >= b->delta*b->delta)
	  outside = 1;
      }
      else {
	/* kappa <= 0 */
	outside = 1;
	aguess = maxsig;
      }

      if (outside) {
	/* Compute the alpha in order to reach the boundary. */
	beta = cblas_ddot(n, b->step, 1, p, 1);
	bound = cblas_ddot(n, p, 1, p, 1);
	root = sqrt(beta*beta+bound*
		    (b->delta*b->delta-cblas_ddot(n, b->step, 1, b->step, 1)));
	beta = (root-beta)/bound;
	aguess = nt_min(beta, aguess);
	/* Piecewise linesearch loop along the current (projected) CG
	   direction. */
	alpha = nlp_projected_search(Order, Uplo, nlp_bounds_projection_st,
				     n, b->x, aguess, g,
				     *H, p, n,
				     b->step, &vplus,
				     splus, gplus, work+2*n,  data, beta,
				     b->projection_data);
	cblas_dcopy(n, splus, 1, b->step, 1);
	value += vplus;
	/* Update the indexes of free variables */
	nlp_as_bounds_step(b->projection_data, n, b->x, 1,
			   b->step, 1, p, 1, xfree, 1, &nfree);
	break;  /* we terminate the CG */
      }
      else {
	/* Compute the next iterate */
	cblas_daxpy(n, alpha, p, 1, b->step, 1);
	cblas_daxpy(n, alpha, pr, 1, g, 1);
	value += alpha*cblas_ddot(n, g, 1, p, 1) + 0.5*alpha*alpha*kappa;

	/* Compute the new free gradient and its norm. */
	nlp_as_rg_bounds(n, g, 1, g, 1, xfree, 1, nfree);
	ngfree = cblas_dnrm2(n, g, 1);

	/* CG termination? */
	if (ngfree < eps_term) {
	  /* Desired accuracy achieved. */
	  goto END;
	}
	if (cgits >= cgitmax) {
	  /* Too many CG iterations */
	  cblas_dcopy(n, gplus, 1, g, 1);
	  goto END;	  
	}

	/* CG body: compute new search direction. */
	beta = ngfree/ng;
	beta *= beta;
	cblas_dscal(n, beta, p, 1);
	cblas_daxpy(n, -1.0, g, 1, p, 1);
	ng = ngfree;
      }
    }

    if (nfree == 0 || nfree == nfree0) break;

    cgitstart = cgits;
  }

  cblas_dcopy(n, gplus, 1, g, 1);

 END:

  /* Memory liberation */
  free(p);
  free(pr);
  free(xfree);
  free(gplus);
  free(splus);
  free(sum);
}
