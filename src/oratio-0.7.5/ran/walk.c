/**
 * walk.c
 *
 * ORATIO Library - Random Numbers Generation module
 * Copyright (C) 2014 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License
 */

/**
 * Random walk 
 */

#include "oratio/ran.h"

/**
 * Produce the next element in a random walk.
 *
 * @param r is an initialized random numbers generator.
 * @param x is the current value in the random walk.
 *
 * @return the next element of the random walk.
 */
double ran_walk_next(Random *r, double x)
{
  return x + ran_random_next(r);
}

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
void ran_walk(Random *r, double *x, int incx, double x0, int T)
{
  int i;

  x[0] = x0;

  for (i = 1; i <= T; i++) x[i*incx] = ran_walk_next(r, x[(i-1)*incx]);
}

/**
 * @param r is an initialized random numbers generator.
 */
void ran_normal_walk(Random *r, double *x, int incx,
		     double x0, int T, double mu, double sigma)
{
  ran_normal_set(r, mu, sigma);

  ran_walk(r, x, incx, x0, T);
}
