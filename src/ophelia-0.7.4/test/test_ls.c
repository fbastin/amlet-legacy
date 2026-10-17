/**
 * test_ls.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "test_ls.h"
// #include "test_minpack2.h"

#include <ophelia/collection.h>

int test_nlp_least_squares(NTLog *log, BTR *btr)
{
  int error_code = 0;
  double **H;
  double work[200];

  int n;

  nt_log_subsection(log, "Basic least-squares problems");

  //  n = 2;
  n = 1;

  op_btr_init(btr, n, 3);
  btr->verbose = 1;

  H = nt_matrix_new(n, n);
  nt_matrix_identity(n, n, *H, n);

  /* Starting point */
  /*  btr->x[0] = 3.0;
      btr->x[1] = 1.0; */
  btr->x[0] = 1.2;

  btr_set_hessian_approximation(btr, NLP_GAUSS_NEWTON, 1);
  //btr_set_hessian_approximation(btr, NLP_SECANT_BFGS, 1);
  //btr_set_hessian_approximation(btr, NLP_SECANT_SR1, 1);
  btr_set_hessian_approximation(btr, NLP_SECANT_LS, 1);

  btr->printer = btr_print_iteration;

  nlp_btrx(btr, (C_GENERIC)nonlinearls1_fg, NlpFG,
	   NULL, 0, NULL, NlpApprox, NULL, 0.0, stdout, H,
	   CblasRowMajor, CblasUpper, work);
  /* nlp_btrx(btr, (C_GENERIC)nonlinearls1, NlpF,
	   NULL, 0, NULL, NlpApprox, NULL, 0.0, stdout, H,
  	   CblasRowMajor, CblasUpper, work); */
  /* nlp_btrx(btr, (C_GENERIC)linearls1_fg, NlpFG,
	   NULL, 0, NULL, NlpApprox, NULL, 0.0, stdout, H,
	   CblasRowMajor, CblasUpper, work); */

  /*  nlp_btrx(btr, (C_GENERIC)madsen_fg_values, NlpMFG,
	   NULL, 0, NULL, NlpApprox, NULL, 0.0, stdout, H,
	   CblasRowMajor, CblasUpper, work); */

  nt_matrix_free(H);

  if (error_code) return error_code;

  nt_log_subsection(log, "MINPACK-2 problems collection");

  //  error_code = test_minpack2_ls(log, btr);

  return error_code;
}
