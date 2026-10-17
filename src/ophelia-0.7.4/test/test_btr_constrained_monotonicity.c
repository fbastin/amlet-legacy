/**
 * test_btr_constrained_monotonicity.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "test_btr_constrained_monotonicity.h"

#include <ophelia/collection.h>

/**
 * @param b is an initialized BTR structure.
 */
int btr_monotonicity_constrained_tests(NTLog *log, BTR *b)
{
  double **H;
  double work[200];
  NLP_Constraint *c = NULL;
  int cidx[4];

  int i, n;

  int error_code = 0;

  /* NLP-UNC2 */
  nt_log_section(log, "NLP-UNC2 problem");

  n = 8;

  op_btr_init(b, n, 0);
  b->max_iter = 30;

  for (i = 0; i < n; i++) b->x[i] = i;
  cidx[0] = 1; cidx[1] = 5;
  cidx[2] = 6; cidx[3] = 8;
  /*  b->x[0] = 0.557242;
  b->x[1] = b->x[2] = 0.863590;
  b->x[3] = 2.28749; */

  H = NULL;

  b->stop = btr_constrained_projection_stopping_test;
  b->step_function = nlp_spcgm_st;
  b->projection_wrapper = nlp_btr_proj_standard;

  /* Set the monotonicity constraints */
  c = nlp_constraint_new(NLP_MONOTONICITY, 2, cidx, 1);
  nlp_btr_constraints_register(b, c);

  //  b->max_iter = 20;

  b->printer = btr_print_iteration;
  btr_set_hessian_approximation(b, NLP_SECANT_SR1, 0);
  /* b->verbose = HIGHLY_VERBOSE; */
  b->verbose = MODERATE_VERBOSE;

  nlp_btr(b, (C_GENERIC)unc2, (C_GENERIC)unc2_grad, NULL,
	  stdout, H, work);

  free(c);

  return error_code;
}
