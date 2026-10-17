/**
 * test_nlp.c
 *
 * OPHELIA Library
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <ophelia/nlp.h>
#include <ophelia/collection.h>

#include <stdio.h>

#include <oratio/ntools.h>

#include "ophelia_test.h"

#include "test_bfgs.h"
#include "test_btr_unconstrained.h"
#include "test_btr_constrained_bounds.h"
#include "test_bounds.h"
#include "test_isotonic.h"
#include "test_ls.h"
// #include "test_siman.h"
#include "test_outer.h"

int bfgs_bounds_constrained_tests(BFGS *b)
{
  return 0;
}

int main(int argc, char **argv)
{
  BTR *b = op_btr_new();
  //  BFGS bf;
  NTLog *log;

  int info = 0;

  op_btr_init(b, 2, 0);

  log = nt_log_new(NULL);

  nt_log_section(log, "Level 3 Matrix operations");

  test_matrix();

  /*
  nt_log_section(log, "Unconstrained optimization");

  nt_log_section(log, "Unidimensional optimization");

  op_gss_test();

  nt_log_section(log, "BTR optimization");

  btr_unconstrained_test(log, b);
  */
  /* 
  nt_log_section(log, "BFGS optimization");

  test_bfgs_unconstrained(log, &bf); */

  /*
  nt_log_section(log, "Constrained optimization");

  nt_log_section(log, "Bounds constrained optimization");
  
  info = btr_bounds_constrained_tests(log, &b);
  */

  /*
  nt_log_section(log, "Optimization with monotonicity constraints");

  info = btr_monotonicity_constrained_tests(log, &b);
  */

  nt_log_section(log, "Simulated annealing");

  siman_test();

  /*
  nt_log_section(log, "Primal-dual algorithms");

  nlp_test_barrier_parameter(&b);

  nlp_outer_tests(&b);
  */

  //  isotonic_test(log);

  /*
  nt_log_section(log, "Least_squares");

  test_nlp_least_squares(log, &b);

  nlp_btr_free(&b);
  */

  nt_log_free(log);

  op_btr_free(b);

  return info;
}
