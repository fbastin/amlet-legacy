/**
 * linesearch.c
 *
 * OPHELIA Library
 *
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "linesearch.h"

#include <oratio/ntools.h>

#include "bounds.h"

/**
 * Piecewise linesearch along a (projected) direction p, using a
 * quadratic model, for a bounds constrained problem.
 *
 * @param Order designs the orientation of the matrix Bk. The
 *   possible values are CblasRowMajor and CblasColMajor.
 * @param Uplo designs the triangular part that has to be updated. The
 *   possible value CblasUpper and CblasLower.
 * @param n is the problem dimension.
 * @param x is the current iterate.
 * @param g is the gradient at x.
 * @param H is the hessian at x.
 * @param p is the direction search, and must corresponds to a descent
 *   direction.
 * @param maxiter is the allowed maximum number of iterations.
 * @param work is a working vector of size at least equal to 2*n.
 * @param data is a pointer to problem-specific data.
 * @param alup is an upper bound on alpha (typically the maximum value
 *   to reach the boundary of a trust-region).
 *
 * @return the step length in direction p.
 */
double nlp_projected_search(const enum CBLAS_ORDER Order,
			    const enum CBLAS_UPLO Uplo,
			    NlpProjection proj,
			    int n, double *x, double alpha,
			    double *g, double *H,
			    double *p,int maxiter,
			    double *s,
			    double *vplus, double *splus, double *gplus,
			    double *work, void *data, double alup,
			    void *data_projection)
{
  double abest;              /* best found alpha */
  double allw;               /* lower bound on alpha */
  double sls;
  double vbest;
  double tmp;

  double minlg = 1.0e-10;    /* maximum length of interval [allw, alup] */

  int above = 0;

  int j;

  allw = 0;
  abest = vbest = *vplus = 0;

  if (!maxiter) maxiter = n;

  for (j = 0; j < maxiter; j++) {
    /* Compute the trial point, the associated gradient and model_value */
    /* We ensure that alpha is not greater than alup. */

    /* We project the step. */
    /* work = x+s+alpha*p */
    cblas_dcopy(n, x, 1, work, 1);
    cblas_daxpy(n, 1.0, s, 1, work, 1);
    cblas_daxpy(n, alpha, p, 1, work, 1);    
    proj(n, work, NULL, work,
	 (NlpBounds *)data_projection);
    /* splus */
    cblas_daxpy(n, -1.0, x, 1, work, 1);
    cblas_dcopy(n, work, 1, splus, 1);
    /* We compute the step from s to splus and store the result in work. */
    cblas_daxpy(n, -1.0, s, 1, work, 1);

    tmp = cblas_dnrm2(n, splus, 1);
    if (tmp < DBL_EPSILON) break;

    sls = cblas_ddot(n, g, 1, work, 1);

    if (sls < 0) {
      /* we are able to lower the model at the first order */
      /* gplus = g+Hs */
      cblas_dcopy(n, g, 1, gplus, 1);
      cblas_dsymv(Order, Uplo, n, 1.0, H, n, work, 1, 1.0, gplus, 1);
      /* vplus = sls + 0.5s'(gplus-g) */
      cblas_dcopy(n, gplus, 1, work+n, 1);
      cblas_daxpy(n, -1.0, g, 1, work+n, 1);
      *vplus = sls + 0.5*cblas_ddot(n, work, 1, work+n, 1);
      if (*vplus < vbest) {
	/* We have found a better solution for the model. */
	vbest = *vplus;
	abest = alpha;
      }
      if (*vplus > 0.01*sls) {
	/* Poor reduction; we have an upper bound on alpha. */
	alup = alpha;
	above = 1;
	/* We proceed by bissection. */
	alpha = 0.5*(allw+alup);
      }
      else if (!above && *vplus < 0.9*sls) {
	/*	tmp = 0.0;
	for (i = 0; i < n; i++)
	  tmp = nt_max(tmp,
	  (((((NlpBounds *)data_projection)->xu)[i]-x[i])-splus[i])*
	  (splus[i]-(((NlpBounds *)data_projection)->xl)[i]-x[i]));
	if (tmp == 0.0)
	  break;
	  else { */
	  /* Define a lower bound on alpha. */
	  allw = alpha;
	  alpha = /*(upset) ? */
	    /* We proceed by bissection */
	    0.5*(allw+alup); /* : */
	    /* Extrapolation */
	    /*	    2*alpha; */
	  /*	} */
      }
      else
	break;
    }
    else {
      /* The model is locally increasing, so we have an upper bound on
	 alpha. */
      alup = alpha;
      //      upset = 1;
      /* We proceed by bissection. */
      alpha = 0.5*(allw+alup);
    }
    if (/*upset && */vbest < 0) {
      /* Check the interval length. */
      if (alup-allw < minlg) break;
    }
  }

  if (vbest < *vplus) {
    cblas_dcopy(n, x, 1, work, 1);
    cblas_daxpy(n, 1.0, s, 1, work, 1);
    cblas_daxpy(n, abest, p, 1, work, 1);    
    proj(n, work, NULL, work,
	 (NlpBounds *)data_projection);
    /* splus */
    cblas_daxpy(n, -1.0, x, 1, work, 1);
    cblas_dcopy(n, work, 1, splus, 1);
    /* gplus = g+H*splus */
    cblas_dcopy(n, g, 1, gplus, 1);
    cblas_dsymv(Order, Uplo, n, 1.0, H, n, splus, 1, 1.0, gplus, 1);
    *vplus = vbest;
  }

  return abest;
}
