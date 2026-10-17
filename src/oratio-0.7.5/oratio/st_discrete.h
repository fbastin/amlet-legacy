/**
 * st_discrete.h
 *
 * ORATIO Library
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the BSD License.
 */

#ifndef __ORATIO_STAT_DISCRETE_H__
#define __ORATIO_STAT_DISCRETE_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Check if a vector can represent discrete probabilities.
   *
   * @param n is the dimension of p.
   * @param p is the vector or probabilities.
   * @param incp is the increment of p.
   * @param tol is the tolerance to use on the sum (taken on absolute value);
   *   if tol <= 0.0, the tolerance n*DBL_EPSILON is used.
   *
   * @return 1 if p properly represent discrete probabilities, 0 else.
   */
  int st_discrete_check(int n, double *p, int incp, double tol);

#ifdef __cplusplus
}
#endif

#endif /* __ORATIO_ST_DISCRETE_H__ */
