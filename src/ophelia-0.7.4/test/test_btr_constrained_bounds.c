/**
 * test_btr_constrained_bounds.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2004-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "test_btr_constrained_bounds.h"
#include "bounds.h"
#include "test_bounds.h"
#include "test_minpack2.h"

#include <ophelia/collection.h>

int test_btr_constrained_nc1(NTLog *log, BTR *b)
{
  int i, n;
  double work[200], xl[200], xu[200];

  NlpBounds nb;

  /* NLP-UNC1 */
  nt_log_section(log, "NLP-UNC1 problem");

  n = 4;

  op_btr_init(b, n, 0);

  nb.xl = xl;
  nb.xu = xu;

  /* Set the starting point and bounds. */
  for (i = 0; i < 4; i++) {
    b->x[i] = i/4+0.25;
    nb.xl[i] = 0;
    nb.xu[i] = 2;
  }

  btr_set_hessian_approximation(b, NLP_SECANT_SR1, 0);
  b->step_function = nlp_spcg_st;
  //  b->proj_data = &nb;

  b->stop = btr_constrained_projection_stopping_test;
  /*  b->projection_wrapper = nlp_bounds_projection_st; */
  //  b->projection_data = &nb;

  nlp_btr(b, (C_GENERIC)unc1, NULL, NULL,
	  log->f, NULL, work);

  return 0;
}

int btr_bounds_constrained_tests(NTLog *log, BTR *b)
{
  double **H;
  double work[200], xl[200], xu[200];
  double *bounds;
  NlpBounds nb;

  int i, n;
  int error_code = 0;

  nt_log_subsection(log, "Basic bounds constrained problems");

  /* PBR-T1-3 */

  n = 2;
  op_btr_init(b, n, 0);

  H = nt_matrix_new(n, n);
  nt_matrix_identity(n, n, H, n);

  nlp_bounds_init(&bounds, n);
  nb.xl = xl;
  nb.xu = xu;

  pbr_t1_3_set_bounds(bounds);
  pbr_t1_3_set_bounds_st(nb.xl, nb.xu);

  //  b->stop = btr_constrained_projection_stopping_test;
  b->step_function = nlp_spcg_st;
  /* b->step_function = nlp_proj_cauchy_bounds_st; */
  //  b->projection_wrapper = nlp_bounds_projection_st;
  //b->proj_data = &nb;
  /*  b->projection = nlp_bounds_projection;
  btr_set_step(b, BTR_STEP_PROJ_CAUCHY);
  b->proj_data = bounds;*/

  /* Starting point */
  b->x[0] = 20.0;
  b->x[1] = 20.0;
  /* b->radius_init = NULL;*/  /* otherwise, the algorithm converges... in
			     just one iteration */

  b->printer = btr_print_iteration;
  btr_set_hessian_approximation(b, NLP_SECANT_SR1, 0);

  nlp_btr(b, (C_GENERIC)pbr_t1_3, NULL, NULL, log->f, H, work);

  nt_matrix_free(H);
  free(bounds);

  error_code = test_btr_constrained_nc1(log, b);
  if (error_code) return error_code;

  /*  nt_log_subsection(log, "TRON problems collection");

  error_code = test_tron(log, b); */

  return error_code;
}
