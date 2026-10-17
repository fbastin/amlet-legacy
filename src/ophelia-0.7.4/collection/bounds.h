/**
 * bounds.h
 *
 * OPHELIA Library - test problems with bound constraints
 * Copyright (C) 2006-2015 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __OP_COLLECTION_BOUNDS_H__
#define __OP_COLLECTION_BOUNDS_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Willi Hock and Klaus Schittkowski, "Test Examples for Nonlinear
   * Programming Codes", Springer-Verlag, 1981.
   */

  /**
   * PBR_T1_3
   */

  double pbr_t1_3(double *x, int n, void *data);

  void pbr_t1_3_set_bounds(double *b);

  /**
   * Set the lower and upper bounds.
   *
   * @param xl is a pointer to the lower bound.
   * @param xu is a pointer to the upper bound.
   */
  void pbr_t1_3_set_bounds_st(double *xl, double *xu);

#ifdef __cplusplus
}
#endif

#endif  /* __OP_COLLECTION_BOUNDS_H__ */
