/**
 * test_bfgs.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <ophelia/bfgs.h>
#include <ophelia/collection.h>

int test_bfgs_unconstrained(NTLog *log, BFGS *b)
{
  double **H;
  double *work;

  int n, m;

  /* CGT-1 */

  n = 2;
  m = 0;
  nlp_bfgs_init(b, n, m);

  H = nt_matrix_new(n, n);
  nt_matrix_identity(n, n, *H, n);
  work = _mallocs(4*n, double);

  /* Starting point */
  b->x[0] = -0.7;
  b->x[1] = 3.2;

  b->d1 = cgt_1_deriv1;

  b->printer = bfgs_print_iteration;

  nlp_bfgsx(b, CblasRowMajor, CblasUpper, (C_GENERIC)cgt_1, NlpF,
	    NULL, NULL, stdout, H, b->n, b->x, work);
  bfgs_print_results(stdout, b, NULL);

  nt_matrix_free(H);
  free(work);

  nlp_bfgs_free(b);

  return 0;
}
