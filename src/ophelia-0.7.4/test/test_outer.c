/**
 * test_btr_outer.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006 Fabian Bastin <bastin@cerfacs.fr>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "test_outer.h" 

#include <stdio.h>

#include <oratio/ntools.h>

#include "bounds.h"

/**
 * Test the barrier parameter evolution.
 */
int nlp_test_barrier_parameter()
{
  int k = 0;
  double mu = 1000000;

  FILE *f = fopen("test_barrier.txt", "w");
  if (!f) return 1;

  while (mu > 1e-8) {
    fprintf(f, "%d %g\n", k++, mu);
    mu = nlp_update_barrier_parameter(mu);
  }
  fprintf(f, "%d %g\n", k++, mu);

  fclose(f);

  return 0;
}

int nlp_outer_tests(BTR *b)
{
  nlp_problem *pbl;
  nlp_outer *no;

  double work[2];
  int n = 2;
  double **H;
  double *bounds;

  nlp_bounds_init(&bounds, n);

  pbl = nlp_problem_new(n);

  n = 2;
  op_btr_init(b, n, 0);

  /* CGT-1 */
  /*  pbl->objective = cgt_1;
  pbl->bounds = bounds;
  bounds[0] = -10;
  bounds[1] = 10;
  bounds[2] = 0;
  bounds[3] = 10; */

  /* Starting point */
  b->x[0] = -0.7;
  b->x[1] = 3.2;

  /* PBR_T1_3 */
  pbl->objective = pbr_t1_3;
  pbr_t1_3_set_bounds(bounds);
  pbl->bounds = bounds;
  //  btr_set_proj_data(b, bounds);
  b->dist = nlp_bounds_dist;
  b->stop = nlp_btr_primal_stopping_test;
  b->tol = 1e-4;
  //  b->feasible = nlp_bounds_strict_feasibility;
  b->feasible = nlp_bounds_feasibility;

  b->x[0] = 10.0;
  b->x[1] = 10.0;

  /* PBR_T1_1 */
  /* nlp_bounds_reset(bounds, n);
  pbl->objective = pbr_t1_1;
  pbr_t1_1_set_bounds(bounds);
  pbl->bounds = bounds;

  b->x[0] = -2.0;
  b->x[1] = 1.0; */

  H = nt_matrix_new(n, n);
  nt_matrix_identity(n, n, *H, n);

  no = nlp_outer_new(pbl, b);
  nlp_outer_set_verbose(no);
  no->penalty = nlp_log_barrier_bounds;
  /*  no->penalty_d1 = nlp_log_barrier_bounds_d1;
      no->penalty_d2 = nlp_log_barrier_bounds_d2; */

  nlp_outer_iteration(no, H, work);

  printf("Outer iterations: %d\n", no->k);

  nlp_outer_free(no);
  nlp_problem_free(pbl);

  nt_matrix_free(H);

  return 0;
}
