/**
 * ran_walk.h
 *
 * ORATIO Library - Random Numbers Generation module
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License
 */

/**
 * Random numbers generator.
 */

#ifndef __ORATIO_RAN_WALK_H__
#define __ORATIO_RAN_WALK_H__

#include "ran_rng.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Produce the next element in a random walk.
   *
   * @param r is an initialized random numbers generator.
   * @param x is the current value in the random walk.
   *
   * @return the next element of the random walk.
   */
  double ran_walk_next(Random *r, double x);

  /**
   * Generate an unidimensional random walk with T+1 values.
   *
   * @param r is an initialized random numbers generator.
   * @param x is a vector where the values along the walk have to be stored
   *   on ouput; x must have a dimension greater or equal to T*incx+1.
   * @param incx is the increment in x.
   * @param x0 is the initial value of the walk.
   * @param T is the horizon of the walk; more specifically,
   *   T values will be generated.
   */
  void ran_walk(Random *r, double *x, int incx,
		double x0, int T);

  /**
   * @param r is an initialized random numbers generator.
   */
  void ran_normal_walk(Random *r, double *x, int incx,
		       double x0, int T, double mu, double sigma);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_RAN_WALK_H__ */
