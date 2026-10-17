/**
 * nt_simplex.h
 *
 * ORATIO Library
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the BSD License.
 */

#ifndef __ORATIO_NTOOLS_SIMPLEX_H__
#define __ORATIO_NTOOLS_SIMPLEX_H__

#ifdef __cplusplus
extern "C" {
#endif

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
  int nt_simplex_std_inside(int n, double *x, int incx, double tol);

#ifdef __cplusplus
}
#endif

#endif /* __ORATIO_NT_SIMPLEX_H__ */
