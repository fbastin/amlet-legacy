/**
 * test_isotonic.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2005-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "test_isotonic.h"
#include "test_errors.h"

#include <ophelia/isotonic.h>

/**
 * @param log is an opened NTools logfile.
 */
int isotonic_test(NTLog *log)
{
  double y[] = {5, 2, 2, 4, 5, 2, 0, 13, 1, 12, 2, 7, 1, 2, 2, 0, 15, 1};
  double w[] = {1, 3, 4, 2, 6, 3, 5, 4, 3, 5, 1, 1, 2, 3, 4, 5, 6, 7};
  int idx[] = {0, 10, 11, 18};
  int nidx = 2;

  double *x;

  int n, m;

  n = 18;

  x = _mallocs(n, double);

  nt_log_section(log, "Isotonic regression tests");

  nt_log_subsection(log, "Simple test");

  n = 10;

  m = nlp_isotonic_regression(n, y, w, x);

  fprintf(log->f, "\n");
  fprintf(log->f, "Problem data\n");
  nt_vect_print(log->f, "y", y, n);
  nt_vect_print(log->f, "w", w, n);
  fprintf(log->f, "\n");
  fprintf(log->f, "Solution\n");
  nt_vect_print(log->f, "x", x, n);

  if (m != 4) {
    fprintf(log->f, "Error on the number of active sets\n");
    return ERROR_ISOTONIC_ACTIVE_SETS;
  }

  nt_log_subsection(log, "Advanced test");

  n = 18;

  m = nlp_isotonics_wrapper(n, y, w, x, nidx, idx);

  fprintf(log->f, "\n");
  fprintf(log->f, "Problem data\n");
  nt_vect_print(log->f, "y", y, n);
  nt_vect_print(log->f, "w", w, n);
  fprintf(log->f, "\n");
  fprintf(log->f, "Solution\n");
  nt_vect_print(log->f, "x", x, n);

  if (m != 6) {
    fprintf(log->f, "Error on the number of active sets\n");
    return ERROR_ISOTONIC_ACTIVE_SETS;
  }
 
  free(x);

  return ERROR_NONE;
}
