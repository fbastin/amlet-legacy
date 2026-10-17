/**
 * outer.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2005 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "outer.h"

#include "btr_struct.h"

#include "barrier.h"
#include "io.h"

#include <stdio.h>

#include <oratio/ctools.h>
#include <oratio/ntools.h>

double nlp_outer_objective(double *x, int n, void *data)
{
  double val = 1.0;
  nlp_outer *no = (nlp_outer *)data;

  /* Compute the regular objective */
  val = no->pbl->objective(x, n, no->data);

  /* Add the penalty term */
  if (no->penalty) val += no->penalty(no->mu, x, no->pbl->bounds, n);

  return val;
}

void nlp_outer_d1(double *x, int n, void *data, double *grad)
{
  double *work;
  nlp_outer *no = (nlp_outer *)data;

  work = _mallocs(n, double);

  /*  nt_vect_print(stdout, "x", x, n); */

  /* Compute the regular gradient */
  /*  val = no->pbl->objective(x, n, no->data); */
  nt_derive_grad_cd(no->pbl->objective, x, grad, work, 0.0,
		    no->b->scaling, n, NULL, NULL, no->data);
  /*  nt_vect_print(stdout, "g1", grad, n); */

  /*  if (btr->d1)
    (btr->d1)(btr->x, n, data, btr->grad);
    else */

  /* Add the penalty term */
  if (no->penalty_d1) {
    no->penalty_d1(no->mu, x, no->pbl->bounds, n, grad);
  }
  else {
    double *tmp;

    tmp = _mallocs(n, double);

    nt_derive_grad_cd(no->pbl->objective, x, tmp, work, 0.0,
		      no->b->scaling, n, NULL, NULL, no->data);

    free(tmp);
  }
  /*  nt_vect_print(stdout, "g2", grad, n); */

  /*  nt_derive_grad_cd(nlp_outer_objective, x, grad, work, 0.0,
      no->b->scaling, n, NULL, data);
      nt_vect_print(stdout, "g1", grad, n); */

/*  else {
  if (no->penalty) val += no->penalty(no->mu, x, no->pbl->bounds, n);
  nt_derive_grad_cd(nlp_outer_objective, x, grad, work, 0.0,
  no->b->scaling, n, NULL, data);
  }*/
  /*    if (btr->d1)
      (btr->d1)(btr->x, n, data, btr->grad);
    else
      nt_derive_grad_cd(f, btr->x, btr->grad, x_cand, 0.0,
      btr->scaling, n, NULL, data); */

  free(work);
}

void nlp_outer_d2(double *x, int n, void *data, double **hessian)
{
  double *work;
  nlp_outer *no = (nlp_outer *)data;

  work = _mallocs(n, double);

  /*  if (btr->d2)
    (btr->d2)(btr->x, n, data, hessian);
    else if (!btr->d2approx) { */
  /*  nt_derive_hess_cd(nlp_outer_objective, x, hessian, work,
      0.0, no->b->scaling, n, NULL, data); */
  /*  nt_derive_grad_cd(no->pbl->objective, x, grad, work, 0.0,
      no->b->scaling, n, NULL, no->data); */
  nt_derive_hess_cd(no->pbl->objective, x, hessian, work,
		    0.0, no->b->scaling, n, NULL, NULL, no->data);
  /* } */

  /* Add the penalty term */
  if (no->penalty_d2) {
    no->penalty_d2(no->mu, x, no->pbl->bounds, n, hessian);
  }
  else {

  }

  free(work);
}

nlp_outer *nlp_outer_new(nlp_problem *pbl, BTR *b)
{
  nlp_outer *no;

  no = _malloc(nlp_outer);

  no->pbl = pbl;
  no->b = b;
  no->k = no->t = 0;
  no->mu = 10;

  no->penalty = NULL;
  no->penalty_d1 = NULL;
  no->penalty_d2 = NULL;

  no->output = NULL;
  no->printer = NULL;

  no->verbose = 1;

  return no;
}

void nlp_outer_set_verbose(nlp_outer *no)
{
  if (!(no->output)) no->output = stdout;
  if (!(no->printer)) no->printer = nlp_outer_print_iteration;

  no->verbose = 1;
}

void nlp_outer_init_barrier_parameter(nlp_outer *no, double mu)
{
  no->mu = mu;
}

void nlp_outer_free(nlp_outer *no)
{
  if (no) free(no);
}

int stopping_outer(double mu)
{
  double tol = 1e-8;

  if (mu < tol) return 1;

  return 0;
}

/**
 * Ref.: Jorge Nocedal and Steven Wright, "Numerical Optimization",
 * Springer, 19999, p. 505.
 */
void nlp_outer_iteration(nlp_outer *no, double **H, double *work)
{
  int stop = 0;

  /*  double *primal = NULL;
      double *dual = NULL; */

#ifdef DEBUG
  /* We record the barrier/penalty parameter evolution */
  FILE *barrier;

  barrier = fopen("barrier.txt", "w");
#endif

  /*  no->b->d1 = nlp_outer_d1;
      no->b->d2 = nlp_outer_d2; */

  do {
    if (no->printer) nlp_outer_print_iteration(no->output, no, no->data);

    op_btr_reset(no->b);
    nt_matrix_identity(no->b->n, no->b->n, *H, no->b->n);

    /* approximately minimize the augmented function */
    /*    nlp_btr(no->b, (C_GENERIC)nlp_outer_objective, NULL, no, stdout,
     *H,  work); */

    /* Stopping test */
    if (stopping_outer(no->mu))
      stop = 1;
    else {
      /* update zeta_k */
      no->k++;
      no->t += no->b->k;

      /* update mu */
      no->mu = nlp_update_barrier_parameter(no->mu);

      /* update the starting point */
    }
  } while (!stop);

  /* Print results */
  printf("Total number of inner iterations: %d\n", no->t);

#ifdef DEBUG
  fclose(barrier);
#endif
}
