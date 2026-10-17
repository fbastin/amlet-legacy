/**
 * triangular.c
 *
 * ORATIO Library
 * Copyright (C) 2004-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/ran.h"
#include "oratio/stat.h"

/**
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 * @param n is the dimension of the seed vector.
 */
Random *ran_triangular_new(unsigned long *seed, int n,
			   double x_min, double x_max, double mode)
{
  Random *random = ran_random_new_with_size(seed, n, 3*sizeof(double));

  ((double *)(random->parameters))[0] = x_min;
  ((double *)(random->parameters))[1] = x_max;
  ((double *)(random->parameters))[2] = mode;

  random->type = RAN_TRIANGULAR;

  return random;
}

double ran_triangular(const Random *random,
		      double x_min, double x_max, double mode)
{
  double p;

  p = random->rng_get_val(random->rng);

  return st_triangular_icdf(p, x_min, x_max, mode);
}

double ran_triangular_get_val(const Random *random)
{
  double x_min, x_max, mode;

  ran_triangular_get_parameters(random, &x_min, &x_max, &mode);
  return ran_triangular(random, x_min, x_max, mode);
}

void ran_triangular_get_parameters(const Random *random,
				   double *x_min, double *x_max, double *mode)
{
  *x_min = ((double *)(random->parameters))[0];
  *x_max = ((double *)(random->parameters))[1];
  *mode  = ((double *)(random->parameters))[2];
}

void ran_triangular_set_parameters(Random *random,
				   double x_min, double x_max, double mode)
{
  double temp[3];

  temp[0] = x_min;
  temp[1] = x_max;
  temp[2] = mode;

  ran_random_set_parameters(random, temp, 3*sizeof(double));
}
