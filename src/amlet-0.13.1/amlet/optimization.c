/**
 * optimization.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2004-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "amlet/optimization.h"

#include <amlet/dcl.h>

#include "amlet_btr.h"
#include "amlet_stopping.h"
#include "bhhh.h"
#include "amlet_driver.h"
#include "info.h"
#include "io.h"
#include "log.h"

void amlet_btr_set_hessian_approx(BTR *b, int hessian_id)
{
  switch (hessian_id) {
  case HESSIAN_BHHH:
    /*    btr_set_hessian_approximation(b, NLP_SECANT_SR1); */
    b->d2approx = (C_GENERIC)bhhh_update;
    break;
  case HESSIAN_SR1:
    btr_set_hessian_approximation(b, NLP_SECANT_SR1, 0);
    break;
  case HESSIAN_BFGS:
  default:
    btr_set_hessian_approximation(b, NLP_SECANT_BFGS, 0);
    break;
  }
}

void amlet_init_btr(BTR *b, AMLET_Driver *dr)
{
  b->printer = amlet_print_iteration_btr;
  b->tol = dr->var[EPS].d;

  /* Simulation experiments tends to show that the radius reduction is
     too slow */
  //  btr_set_constants(b, 0.01, 0.9, 0.2, 0.5);
  //  btr_set_constants(b, 0.01, 0.75, 0.5, 0.5);
  btr_set_constants(b, 0.01, 0.9, 0.3, 0.5);

  b->max_iter = dr->var[NITER].i;

  /* TODO: define a strategy to select the initial radius */
  b->delta = AMLET_DEFAULT_TRUST_REGION_RADIUS;
  b->d1 = likelihood_gradient;

  b->data = 0;

  b->output = mclog->f;

  /* TODO: better initialize the approximations */
  /* Hessian approximation */
  amlet_btr_set_hessian_approx(b, dr->var[HESSIAN_UPDATE].i);
}

void amlet_init_bfgs(BFGS *bf, AMLET_Driver *dr)
{
  bf->printer = amlet_print_iteration_bfgs;
  bf->tol = dr->var[EPS].d;
  bf->max_iter = dr->var[NITER].i;

  bf->d1 = likelihood_gradient;
}
