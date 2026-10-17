/**
 * radius.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 1999-2003 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "radius.h"

#include <float.h>

#include <oratio/ntools.h>

/**
 * Accept or reject the iterate candidate.
 *
 * @param b is an initialized BTR structure.
 * @param rho is the ratio between actuel and predicted reductions.
 *
 * @return 1 if the iterate candidate is accepted, 0 else.
 */
int nlp_btr_accept_candidate(BTR *b, double rho)
{
  return (rho >= b->eta1) ? 1 : 0;
}

/**
 * Radius initialization procedures
 */

/**
 * Compute the initial trust-region radius on basis of the gradient value.
 *
 * Ref.: Lin and More, "Newton's Method For Large Bound-Constrained
 * Optimization Problems", Siam Journal on Optimization 9(4), 1999.
 *
 * @param b is an initialized BTR structure.
 */
double btr_radius_init_grad(BTR *b, void *data)
{
  /* This is inspired form Lancelot and reflects the fact that a too
     large value could lease to further problems. */
  const double factor = 0.1;

  return factor*cblas_dnrm2(b->n, b->grad, 1);
}

/**
 * Radius update procedures
 */

void btr_radius_update_basic(BTR *b, double rho)
{
  if (rho >= b->eta1) {
    if (rho >= b->eta2)
      /* the trust-region radius is doubled */
      b->delta *= 2;
    else
      b->delta *= (b->gamma2+1)/2;
  }
  else /* trust-region radius reduction */
    b->delta *= (b->gamma1+b->gamma2)/2;

  /* This line protects us against infinity in IEEE artithmetic mode */
  b->delta = nt_min(DBL_MAX, b->delta);
}

void btr_radius_update_usual(BTR *b, double rho)
{
  double stepnrm;

  if (rho >= b->eta1) {
    stepnrm = cblas_dnrm2(b->n, b->step, 1);

    if (b->k == 0) b->delta = nt_min(stepnrm, b->delta);

    if (rho >= b->eta2)
      b->delta = nt_min(10e20, nt_max(4*stepnrm, b->delta));
    else
      b->delta *= b->gamma2;
  }
  else /* trust-region radius reduction */
    b->delta *= b->gamma1;
}
