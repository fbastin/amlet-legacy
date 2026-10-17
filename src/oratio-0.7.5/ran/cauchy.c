/**
 * cauchy.c
 *
 * ORATIO Library
 * Copyright (C) 1999-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Generation of points following a Cauchy distribution
 */

#include "oratio/ran.h"
#include "oratio/stat.h"

#include <math.h>

/**
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 * @param n is the dimension of the vector seed.
 */
Random *ran_cauchy_new(unsigned long *seed, int n,
		       double a, double b)
{
  Random *random = ran_random_new_with_size(seed, n, 2*sizeof(double));

  ((double *)(random->parameters))[0] = a;
  ((double *)(random->parameters))[1] = b;

  random->type = RAN_CAUCHY;

  return random;
}

void ran_cauchy_get_parameters(const Random *random, double *a, double *b)
{
  *a = ((double *)(random->parameters))[0];
  *b = ((double *)(random->parameters))[1];
}

void ran_cauchy_set_parameters(Random *random, double a, double b)
{
  double temp[2];

  temp[0] = a;
  temp[1] = b;

  ran_random_set_parameters(random, temp, 2*sizeof(double));
}

/**
 * Generate points following a standard Cauchy distribution, i.e. a
 * Cauchy distribution of parameters (0,1).
 */
double ran_cauchy_standard_get_val(const Random *random)
{
  return st_cauchy_icdf(ran_random_get_val(random), 0, 1);
}

/**
 * Draw from a Cauchy distribution.
 *
 * @param random is an initialized random numbers generator.
 * @param a is the location parameter, the median.
 * @param b is scale parameter; b must be stricty greater than 0.
 *
 * @return a (pseudo)-random number from a Cauchy distribution of
 *   location parameter a and scale parameter b.
 */
double ran_cauchy(const Random *random, double a, double b)
{
  return st_cauchy_icdf(ran_random_get_val(random), a, b);
}

double ran_cauchy_get_val(const Random *random)
{
  double a, b;

  ran_cauchy_get_parameters(random, &a, &b);
  return ran_cauchy(random, a, b);
}
