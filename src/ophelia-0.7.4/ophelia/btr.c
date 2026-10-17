/**
 * btr.c
 *
 * OPHELIA Library
 * Copyright (C) 2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Faire référence au papier Mathematical Programming.
 */

#include "btr.h"

#include <string.h>

#include "btr_struct.h"
#include "io.h"
#include "radius.h"
#include "step.h"
#include "stopping.h"
#include "utilities.h"
#include "wrapper.h"

#include "defaults.h"
#include "model.h"

#include <oratio/oratio.h>

/**
 * Define if derivatives are computed along with the function or not.
 */
void nlp_set_derivatives_type(BTR *b, int f_type, int *d1, int *d2)
{
  *d1 = *d2 = NlpNone;

  switch (f_type) {
  case 0:
    break;
  case NlpFG:
  case NlpMFG:
    *d1 = NlpComputed;
    break;
  case NlpFGH:
    *d1 = NlpComputed;
    *d2 = NlpComputed;
    break;
  default:
    *d1 = NlpApprox;
    break;
  }

  if (b->n_secant_updates) *d2 = NlpApprox;
}

/**
 * BTR implementation. Expert routine
 *
 * @param f is the objective function.
 * @param f_type is the objective function type.
 * @param g is the gradient evaluation function.
 * @param g_type is the gradient function type.
 * @param H is the hessian evaluation function.
 * @param h_type is the hessian function type.
 * @param data: problem data.
 * @param tol: tolerance used in the function evaluation; if tol is
 *   set to 0.0, the machine precision is used.
 * @param output: output file for results and iteration informations
 * @param hessian: preallocated array of size btr->n*btr-btr->n, if
 *   the hessian is needed
 * @param Order: CblasRowMajor or CblasColMajor
 * @param Uplo: CblasUpper or CblasLower
 * @param work: working vector of size 4*btr->n.
 *
 * @return the function value at x, who should satisfy some
 *   criticality criteria (depending on the step function)
 */
double nlp_btrx(BTR *btr, C_GENERIC f, int f_type,
		C_GENERIC g, int g_type,
		C_GENERIC H, int h_type,
		void *data, double tol,
		FILE *output, double **hessian,
		const enum CBLAS_ORDER Order,
		const enum CBLAS_UPLO Uplo,
		double *work)
{
  double new_value, rho = 0.0;
  double *x_cand = NULL;  /* candidate */
  double tmp, val, best;
  double *values = NULL, *valuesplus = NULL;

  double *newgrad = NULL;
  int n;
  int i, newwork = 0;

  int has_newgrad, has_newHessian;
  int id_secant = -1, new = -1;

  int d1, d2;

  double **Horig = hessian;

  /* Initialization */
  n = btr->n;
  btr->k = 0; /* set the iteration counter */

  /* The 2-norm is currently the only allowed norm */

  nlp_set_derivatives_type(btr, f_type, &d1, &d2);

  /* Parameters checking.  */
  if (!(f && (d1 || g)) || !(hessian || btr->n_secant_updates)) {
    btr->info = NLP_STOPPING_PARAMETERS;
    goto End;
  }

  if (!work) {
    work = _mallocs(4*n, double);
    newwork = 1;
  }
  btr->work = work;

  /* Memory allocation */
  x_cand = _mallocs(n, double);
  /* if (btr->d2approx || fg) */
  newgrad = _mallocs(n*(btr->m+1), double);
  if (btr->m) {
    values = _mallocs(btr->m, double);
    valuesplus = _mallocs(btr->m, double);
  }

  for (i = 0; i < n; i++) btr->step[i] = DBL_MAX;

  //  printf("BTR: %d\n", btr->n_secant_updates);

  /* Compute the function value and its gradient. */
  btr->value = btr->wrapper_fct(f, f_type, btr->x, n, btr->m, data,
				values, btr->grad,
				hessian, NULL, 0, &(btr->info));
  if (btr->verbose == HIGHLY_VERBOSE)
    printf("Value at starting point: %g\n", btr->value);

  if (d1 != NlpComputed) {
    btr->wrapper_grad(f, f_type, g, g_type, btr->x, n, data, btr->value,
		      btr->grad, btr->scaling,
		      NULL, 0, &(btr->info), x_cand, work);
  }

  if (d2 != NlpApprox)
    btr->wrapper_hess(Order, Uplo, f, f_type, g, g_type, H, h_type,
		      btr->x, n, data, work, hessian, btr->scaling,
		      NULL, 0, &(btr->info), x_cand);
  else
    id_secant = 0;

  for (i = 0; i < btr->n_secant_updates; i++) {
    btr->secant_updates[i]->init(Order, Uplo, btr->n,
				 btr->m, *(btr->secant_updates[i]->H),
				 values, valuesplus,
				 btr->grad, NULL, NULL,
				 &(btr->info), work, data);
  }

  /* Trust-region radius initialization */
  if (btr->radius_init) btr->delta = btr->radius_init(btr, data);

  while (!(btr->info = btr->stop(btr, rho)) && !btr->info) {

    if (btr->verbose && btr->n_secant_updates > 1) {
      printf("Secant method: %d %d\n", id_secant,
	     btr->secant_updates[id_secant]->id);
    }
    if (btr->printer && output) btr->printer(output, btr, data);

    has_newgrad = has_newHessian = 0;

    if (id_secant >= 0)
      hessian = nlp_secant_get_hessian(btr->secant_updates[id_secant],
				       Order, Uplo, n, btr->m);
    btr->step_function(btr, Order, Uplo, hessian, work, data);

    cblas_dcopy(n, btr->x, 1, x_cand, 1);
    cblas_daxpy(n, 1.0, btr->step, 1, x_cand, 1);

    if (!btr->feasible(x_cand, btr->n, btr->data)) {
      /* It would be better to refine the step computation */
      rho = 0.0;
      new_value = DBL_MAX;
      if (btr->verbose) c_warning("unfeasible point encountered");
    }
    else {
      new_value = btr->wrapper_fct(f, f_type, x_cand, n, btr->m,
				   data, valuesplus, newgrad,
				   hessian, NULL, 0, &(btr->info));

      if (d1 == NlpComputed) has_newgrad = 1;
      if (d2 == NlpComputed) has_newHessian = 1;

      rho = 0.0;
      tmp = 0.0;
      if (finite(new_value)) {
	tmp = btr->value-new_value;
	best = rho = tmp /
	  btr_quadratic_model(Order, Uplo, n, btr->step,
			      btr->grad, *hessian, data, work);
	btr->rho = rho;
	if (id_secant >= 0) {
	  new = id_secant;
	  for (i = 0; i < btr->n_secant_updates; i++) {
	    if (i != id_secant) {
	      hessian = nlp_secant_get_hessian(btr->secant_updates[i],
					       Order, Uplo, n, btr->m);
	      val = tmp /
		btr_quadratic_model(Order, Uplo, n, btr->step, btr->grad,
				    *hessian, data, work);
	      if (val > best) {
		best = val;
		new = i;
	      }
	    }
	  }
	  id_secant = new;
	}
      }

      if (d2 == NlpApprox && finite(new_value)) {
	if (d1 != NlpComputed)
	  btr->wrapper_grad(f, f_type, g, g_type, x_cand, n, data, new_value,
			    newgrad, btr->scaling, NULL, 0,
			    &(btr->info), work, work+n);

	for (i = 0; i < btr->n_secant_updates; i++) {
	  btr->secant_updates[i]->update(Order, Uplo, btr->n,
					 btr->m, *(btr->secant_updates[i]->H),
					 values, valuesplus,
					 btr->grad, newgrad,
					 btr->step, &(btr->info),
	  				 work, data);
	}
	has_newgrad = 1;
	has_newHessian = 1;
      }
    }

    /* Update the iterate */
    if (nlp_btr_accept_candidate(btr, rho)) {
      /* Set the new iterate to the candidate point. */
      cblas_dcopy(btr->n, x_cand, 1, btr->x, 1);
      btr->value = new_value;
      if (values) memcpy(values, valuesplus, (btr->m)*sizeof(double));

      if (has_newgrad)
	cblas_dcopy(n*(btr->m+1), newgrad, 1, btr->grad, 1);
      else
	btr->wrapper_grad(f, f_type, g, g_type, btr->x, n, data, btr->value,
			  btr->grad, btr->scaling,
			  NULL, 0, &(btr->info), x_cand, work);
      if (!has_newHessian)
	btr->wrapper_hess(Order, Uplo, f, f_type, g, g_type, H, h_type,
			  btr->x, n, data, work, hessian, btr->scaling,
			  NULL, 0, &(btr->info), x_cand);
    }

    if (btr->retro && (rho >= btr->eta1)) {
      rho = tmp /
	btr_quadratic_retro_model(Order, Uplo, n, btr->step,
				  btr->grad, *hessian, data, work);
    }

    /* Update the trust region radius */
    btr->update_radius(btr, rho);

    if ((++(btr->k)) == btr->max_iter) btr->info = NLP_STOPPING_ITERATIONS;

  }

 End:

  if (id_secant >= 0) {
    hessian = nlp_secant_get_hessian(btr->secant_updates[id_secant],
				     Order, Uplo, n, btr->m);
    if (Horig) nt_matrix_insert_matrix(Horig, n, n, hessian, n, n, 0, 0);
  }

  /* Print results */
  if (output) btr_print_results(output, btr, data);

  /* memory deallocation */

  free(x_cand);
  if (newgrad) free(newgrad);

  if (newwork) free(work);
  if (values) {
    free(values);
    free(valuesplus);
  }

  return btr->value;
}

/**
 * BTR algorithm. Simple routine.
 *
 * @param g if not NULL, used to evaluate the function too.
 */
double nlp_btr(BTR *btr, C_GENERIC f, C_GENERIC g, void *data,
		      FILE *output, double **hessian, double *work)
{
  if (g)
    return nlp_btrx(btr, f, NlpF, g, NlpFG, NULL, NlpApprox, data, 0.0,
		    output, hessian, CblasRowMajor, CblasUpper, work);
  else
    return nlp_btrx(btr, f, NlpF, NULL, 0, NULL, 0, data, 0.0,
		    output, hessian, CblasRowMajor, CblasUpper, work);
}
