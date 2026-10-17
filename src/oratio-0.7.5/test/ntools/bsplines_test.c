/**
 * bsplines_test.c
 *
 * Random Numbers Library
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "ntools_test.h"

#include <string.h>

/**
 * B-Splines test routine.
 */
int nt_bsplines_check()
{
  double x[] = {0, 1, 2, 3, 4, 4, 5};
  double *U, *N;
  double u, val;
  int n, m;
  int p;
  int span;

  /* Check basis functions with Example 2.2 of Piegl and Tiller (1997) */
  p = 2;
  m = n = 7;

  /* ntools_bsplines_knots_create(p, n, x, &n, &U); */
  U = nt_bsplines_knots_allocate(p, n);

  memcpy(U, x, m*sizeof(double));

  n = m;
  nt_bsplines_knots_initialize(p, &n, &U, 1);

  N = _mallocs(p+1, double);
  u = 2.5;

  val = nt_bsplines_curve_point(n, p, u, U, x, N);

  /* The values of N should be 1/8, 6/8, 1/8. */
  nt_vect_print(stdout, "N", N, p+1);
  printf("Value: %g\n", val);

  span = nt_bsplines_find_span(n, p, u, U);
  printf("span: %d\n", span);

  free(U);
  free(N);

  return 0;
}
