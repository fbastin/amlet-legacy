/**
 * test_btr_unconstrained.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2004-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <oratio/ntools.h>

#include <ophelia/nlp.h>
#include <ophelia/collection.h>

int btr_unconstrained_test(NTLog *log, BTR *b)
{
  double **H;
  double work[200];

  int i, n;

  int error_code = 0;

  b->retro = 1;

  /* CGT-1 */

  n = 2;

  H = nt_matrix_new(n, n);
  nt_matrix_identity(n, n, *H, n);

  nt_log_subsection(log, "CGT-1 problem");
  fprintf(log->f,
	  "\nFinite-difference first- and second-order derivatives\n\n");

  op_btr_init(b, n, 0);

  /* Starting point */
  b->x[0] = -0.7;
  b->x[1] = 3.2;

  b->printer = btr_print_iteration;

  nlp_rtr(b, (C_GENERIC)cgt_1, NULL, NULL, log->f, H, work);

  fprintf(log->f, "\nAnalytical first-order derivatives; approximations of second-order derivatives\n\n");

  op_btr_init(b, n, 0);

  /* Starting point */
  b->x[0] = -0.7;
  b->x[1] = 3.2;

  nt_matrix_identity(n, n, *H, n);

  btr_set_hessian_approximation(b, NLP_SECANT_SR1, 0);

  b->printer = btr_print_iteration;
  nlp_rtr(b, (C_GENERIC)cgt_1, (C_GENERIC)cgt_1_deriv1, NULL,
	  log->f, H, work);

  fprintf(log->f, "\nAnalytical first- and second-order derivatives\n\n");

  op_btr_init(b, n, 0);

  /* Starting point */
  b->x[0] = -0.7;
  b->x[1] = 3.2;

  nt_matrix_identity(n, n, *H, n);

  b->printer = btr_print_iteration;

  nlp_btrx(b, (C_GENERIC)cgt_1, NlpF,
	   (C_GENERIC)cgt_1_deriv1, NlpFval,
	   (C_GENERIC)cgt_1_deriv2, NlpExact,
	   NULL, 0.0, log->f, H, CblasRowMajor, CblasUpper, work);

  nt_matrix_free(H);

  /* Rosenbrock */

  n = 2;

  H = nt_matrix_new(n, n);
  nt_matrix_identity(n, n, *H, n);

  nt_log_subsection(log, "Rosenbrock problem");

  btr_set_hessian_approximation(b, NLP_SECANT_BFGS, 0);

  op_btr_init(b, n, 0);

  /* Starting point */
  b->x[0] = 0.0;
  b->x[1] = 0.0;

  b->printer = btr_print_iteration;

  nlp_rtr(b, (C_GENERIC)rosenbrock, NULL, NULL, log->f, H, work);

  nt_matrix_free(H);

  /* NLP-UNC1 */
  nt_log_section(log, "NLP-UNC1 problem");

  n = 4;

  op_btr_init(b, n, 0);

  for (i = 0; i < 4; i++) b->x[i] = i;

  H = NULL;

  btr_set_hessian_approximation(b, NLP_SECANT_SR1, 0);

  nlp_btr(b, (C_GENERIC)unc1, NULL, NULL,
	  log->f, H, work);

  return error_code;
}
