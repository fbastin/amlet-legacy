/**
 * btr_struct.c
 *
 * OPHELIA Library
 * Copyright (C) 1999-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <stdlib.h>

#include "bounds.h"
#include "btr_defaults.h"
#include "defaults.h"
#include "io.h"
#include "radius.h"
#include "secant.h"
#include "step.h"
#include "stopping.h"
#include "utilities.h"
#include "wrapper.h"

#include "model.h"
#include "pcg.h"
#include "projection.h"

#include <oratio/ctools.h>
#include <oratio/ntools.h>

/**
 * Set the default constants of the BTR algorithm.
 *
 * @param btr is an initialized BTR structure.
 */
void btr_set_constants(BTR *btr, double eta1, double eta2,
			      double gamma1, double gamma2)
{
  btr->eta1 = eta1;
  btr->eta2 = eta2;
  btr->gamma1 = gamma1;
  btr->gamma2 = gamma2;
}

/**
 * Set the trust-region radius.
 *
 * @param btr is an initialized BTR structure.
 * @param delta is the desired value for the trust-region radius.
 */
void btr_set_delta(BTR *btr, double delta)
{
  btr->delta = delta;
}

/**
 * Set the tolerance used to test convergence.
 *
 * @param btr is an initialized BTR structure.
 */
void btr_set_tol(BTR *btr, double tol)
{
  btr->tol = tol;
}

/**
 * Set the projection information.
 *
 * @param btr is an initialized BTR structure.
 */
void btr_set_projection_data(BTR *btr, void *data)
{
  btr->projection_data = data;
}

/**
 * Set the maximum iterations number.
 *
 * @param btr is an initialized BTR structure.
 */
void btr_set_max_iter(BTR *btr, int max_iter)
{
  btr->max_iter = max_iter;
}

/**
 *
 * @param btr is an initialized BTR structure.
 */
void btr_set_hessian_approx(BTR* btr, C_GENERIC fct)
{
  btr->d2approx = fct;
}

/**
 *
 * @param btr is an initialized BTR structure.
 */
void nlp_btr_wrappers_default(BTR *btr)
{
  btr->wrapper_fct = nlp_wrapper_fct_dflt;
  btr->wrapper_grad = nlp_wrapper_grad_dflt;
  btr->wrapper_hess = nlp_wrapper_hess_dflt;
}

/**
 * Constraints management.
 */

/**
 * Initialize the constraints element of a BTR structure.
 *
 * @param btr is an initialized BTR structure.
 */
void nlp_btr_constraints_init(BTR *btr)
{
  btr->constraints = ds_list_new();
  btr->nconstraints = 0;
}

/**
 * Register a constraints block.
 *
 * @param btr is an initialized BTR structure.
 *
 * @return 1 if the registration has been successful 0 else.
 */
int nlp_btr_constraints_register(BTR *btr, NLP_Constraint *c)
{
  DS_List *l = btr->constraints;
  NLP_Constraint *c2;

  while (l) {
    l = l->next;
    c2 = (NLP_Constraint *)(l->data);
    if (c2->id == c->id) return 0;
  }

  btr->constraints = ds_list_prepend(btr->constraints, c);
  btr->nconstraints++;

  return 1;
}

/**
 * Deregister a constraint block.
 *
 * @param btr is an initialized BTR structure.
 * @param id is the constraint identifier returned by nlp_btr_constraint
 *
 * @return 1 if the deregistration has been successful 0 else.
 */
int nlp_btr_constraints_deregister(BTR *btr, int id)
{
  DS_List *l = btr->constraints;
  NLP_Constraint *c;

  /* Search the constraints blocs. */
  while (l) {
    c = (NLP_Constraint *)(l->data);
    if (c->id == id) {
      ds_list_remove_node(l, l);
      btr->nconstraints--;
      return 1;
    }
    l = l->next;
  }

  return 0;
}

/**
 * @param btr is an initialized BTR structure.
 */
void nlp_btr_constraints_deregister_all(BTR *btr)
{
  ds_list_free(btr->constraints, NULL);
  btr->nconstraints = 0;
}

/**
 * Secant updates management.
 *
 * @param btr is an initialized BTR structure.
 */
double **nlp_btr_secant_register(BTR *btr, int secant_id)
{
  NlpSecant *ns;
  btr->secant_updates = _reallocs(btr->secant_updates,
				  ++(btr->n_secant_updates),
				  struct _NlpSecant *);
  ns = nlp_secant_new(btr->n, secant_id);

  btr->secant_updates[btr->n_secant_updates-1] = ns;

  if (secant_id == NLP_SECANT_LS) {
    int i;

    /* Search NLP_GAUSS_NEWTON */
    for (i = 0; i < btr->n_secant_updates-1; i++) {
      if (btr->secant_updates[i]->id == NLP_GAUSS_NEWTON) {
	ns->ns = btr->secant_updates[i];
	break;
      }
    }
  }
  return ns->H;
}

/**
 *
 * @param btr is an initialized BTR structure.
 */
void nlp_btr_secant_deregister(BTR *btr, int secant_id)
{
  int i;

  for (i = 0; i < btr->n_secant_updates; i++) {
    if (btr->secant_updates[i]->id == secant_id) {
      nlp_secant_free(btr->secant_updates[i]);
      btr->secant_updates[i] = btr->secant_updates[btr->n_secant_updates-1];
      break;
    };
  }

  if (--(btr->n_secant_updates))
    btr->secant_updates = _reallocs(btr->secant_updates, btr->n_secant_updates,
				    NlpSecant *);
}

/*
 *
 * @param btr is an initialized BTR structure.
 */
double **btr_set_hessian_approximation(BTR *b, int id, int info)
{
  double **H;
  b->d2approx = nlp_set_hessian_approximation(id, info);

  H = nlp_btr_secant_register(b, id);
  nlp_secant_set_fct(b->secant_updates[b->n_secant_updates-1]);

  return H;
}

/**
 * Initialization functions
 */

/**
 * Initialize a BTR structure.
 *
 * @param btr is a memory allocated BTR structure.
 * @param n is the problem dimension.
 * @param m is the secondary, problem dependent, dimension, usually set to 0.
 */
void op_btr_init(BTR *btr, int n, int m)
{
  int i;

  btr_set_constants(btr, BTR_ETA_1, BTR_ETA_2, BTR_GAMMA_1, BTR_GAMMA_2);
  btr->retro = 0;

  /* Problem dimension */
  btr->n = n;
  btr->m = m;

  /* Iterate */
  btr->x = _reallocs(btr->x, n, double);
  btr->scaling = _reallocs(btr->scaling, n, double);
  btr->grad = _reallocs(btr->grad, n*(m+1), double);
  btr->step = _reallocs(btr->step, n, double);

  /* Some initialization */
  for (i = 0; i < n; i++) btr->x[i] = 0.0;
  for (i = 0; i < n; i++) btr->scaling[i] = 1.0;
  for (i = 0; i < n*(m+1); i++) btr->grad[i] = 0.0;
  for (i = 0; i < n; i++) btr->step[i] = 0.0;

  btr->value = DBL_MAX;

  /* Initial trust region radius */
  btr_set_delta(btr, BTR_INITIAL_DELTA);
  btr->update_radius = btr_radius_update_usual;

  btr_set_tol(btr, NLP_TOL);
  btr_set_max_iter(btr, NLP_MAX_ITER);

  btr->k = 0;  /* iterations counter */

  /* Constraints sets. */
  nlp_btr_constraints_init(btr);

  while (btr->n_secant_updates)
    nlp_btr_secant_deregister(btr, btr->secant_updates[0]->id);

  /* derivation routines */
  btr->d2approx         = NULL;
  btr->n_secant_updates = 0;
  btr->secant_updates   = NULL;
  nlp_btr_wrappers_default(btr);

  nlp_btr_constraints_deregister_all(btr);
  nlp_btr_constraints_init(btr);

  /* step */
  btr_set_step(btr, BTR_STEP_STEIHAUG_TOINT);

  /* We assume by default that no projection has to be performed. */
  btr->projection_wrapper      = NULL;
  btr->projection_norm_wrapper = NULL;
  btr->projection_data         = NULL;

  btr->feasible = nlp_feasible_default;

  /* stopping test */
  btr->stop = nlp_btr_stopping_test;

  /* trust-region radius initialization */
  btr->radius_init = btr_radius_init_grad;

  btr->output = stdout;
  btr->printer = NULL;

  btr->info = 0;

  btr->work = btr->data = NULL;

  btr->dist = NULL;
}

/**
 * Reset a btr structure.
 *
 * @param btr is an initialized BTR structure.
 */
void op_btr_reset(BTR *btr)
{
  int i;

  btr->value = DBL_MAX;
  for (i = 0; i < btr->n*(btr->m+1); i++) btr->grad[i] = DBL_MAX;
  for (i = 0; i < btr->n; i++) btr->step[i] = 0.0;

  btr->k = 0;  /* iterations counter */

  /* Initial trust region radius */
  btr_set_delta(btr, BTR_INITIAL_DELTA);

  btr->info = 0;
}

/**
 * Create a new BTR structure.
 */
BTR *op_btr_new()
{
  BTR *btr = _malloc(BTR);

  btr->n = btr->m = 0;

  btr_set_constants(btr, BTR_ETA_1, BTR_ETA_2, BTR_GAMMA_1, BTR_GAMMA_2);
  btr->retro = 0;

  /* Iterate */
  btr->x = NULL;
  btr->x_type = NULL;
  btr->scaling = NULL;
  btr->grad = NULL;
  btr->step = NULL;

  btr->value = DBL_MAX;

  /* Initial trust region radius */
  btr_set_delta(btr, BTR_INITIAL_DELTA);
  btr->update_radius = btr_radius_update_usual;

  btr_set_tol(btr, NLP_TOL);
  btr_set_max_iter(btr, NLP_MAX_ITER);

  btr->k = 0;  /* iterations counter */

  /* Constraints sets. */
  nlp_btr_constraints_init(btr);

  /* derivation routines */
  btr->d2approx         = NULL;
  btr->n_secant_updates = 0;
  btr->secant_updates   = NULL;
  nlp_btr_wrappers_default(btr);

  /* step */
  btr_set_step(btr, BTR_STEP_STEIHAUG_TOINT);

  /* We assume by default that no projection has to be performed. */
  btr->projection_wrapper      = NULL;
  btr->projection_norm_wrapper = NULL;
  btr->projection_data         = NULL;

  btr->feasible = nlp_feasible_default;

  /* stopping test */
  btr->stop = nlp_btr_stopping_test;

  /* trust-region radius initialization */
  btr->radius_init = btr_radius_init_grad;

  btr->output = stdout;
  btr->printer = NULL;

  btr->info = 0;

  btr->work = btr->data = NULL;

  btr->dist = NULL;

  btr->verbose = SILENT;

  return btr;
}

/**
 * Free memory associated to a BTR object.
 *
 * @param btr is an initialized BTR structure.
 */
void op_btr_free(BTR *btr)
{
  free(btr->x);
  free(btr->scaling);
  free(btr->step);
  free(btr->grad);

  while (btr->n_secant_updates) {
    nlp_btr_secant_deregister(btr, btr->secant_updates[0]->id);
  }

  nlp_btr_constraints_deregister_all(btr);
}
