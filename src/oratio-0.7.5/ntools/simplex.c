/**
 * simplex.c
 *
 * ORATIO Library
 * Copyright (C) 2004-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/ctools.h"

/**
 * Check if a point is in the standard simplex.
 *
 * @param n is the dimension of x.
 * @param x is the point to test.
 * @param incx is the increment of x.
 * @param tol is the tolerance to subtract from the sum.
 *
 * @return 1 if x is in the simplex, 0 else.
 */
int nt_simplex_std_inside(int n, double *x, int incx, double tol)
{
  int i, sum = 0.0;

  for (i = 0; i < n; i++) {
    if (x[i*incx] < 0) return 0;
    sum += x[i*incx];
  }

  if (sum >= 1.0+tol) return 0;

  return 1;
}
