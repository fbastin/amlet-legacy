/**
 * cauchy.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 1999-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "step.h"

#include <oratio/ctools.h>
#include <oratio/ntools.h>

/**
 * Cauchy step implementation
 *
 * Ref.: Jorge Nocedal and Stephen J. Wright, "Numerical Optimization",
 *   Springer, pp 69-70, 1999
 */
void btr_step_cauchy(BTR *btr,
		     const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
		     double **hessian, double *work, void *data)
{
  double curvature, tau, norm;
  double *pr;

  int n;

  n = btr->n;
  pr = _mallocs(n, double);

  cblas_dsymv(Order, Uplo, n, 1.0, *hessian, n, btr->grad, 1, 0.0, pr, 1);
  curvature = cblas_ddot(n, btr->grad, 1, pr, 1);
  norm = cblas_dnrm2(n, btr->grad, 1);

  if (curvature <= 0)
    tau = 1.0;
  else
    tau = nt_min(norm*norm*norm/(btr->delta*curvature),1.0);

  cblas_dcopy(n, btr->grad, 1, btr->step, 1);
  cblas_dscal(n, -tau*btr->delta/norm, btr->step, 1);

  free(pr);
}
