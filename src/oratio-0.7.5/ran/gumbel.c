/**
 * gumbel.c
 *
 * Random Numbers Library
 * Copyright (C) 2002-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Generation from a Gumbel distribution
 */

#include "oratio/ran.h"
#include "oratio/stat.h"

#include <math.h>

/**
 * Create a Gumbel random variate generator.
 *
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 * @param n is the dimension of the vector seed.
 * @param a is the location parameter.
 * @param b is the scale parameter.
 *
 * @return a pointer to the Gumbel random numbers generator.
 */
Random *ran_gumbel_new(unsigned long *seed, int n,
		       double a, double b)
{
  Random *random = ran_random_new_with_size(seed, n, 2*sizeof(double));

  ((double *)(random->parameters))[0] = a;
  ((double *)(random->parameters))[1] = b;

  random->type = RAN_GUMBEL;

  return random;
}

/**
 * Return a draw from a Gumbel distribution.
 */
double ran_gumbel(const Random *random, double a, double b)
{
  return st_gumbel_icdf(random->rng_get_val(random->rng), a, b);
}

/**
 * Return a draw from a Gumbel distribution.
 */
double ran_gumbel_get_val(const Random *random)
{
  double a, b;

  ran_gumbel_get_parameters(random, &a, &b);
  return ran_gumbel(random, a, b);
}

/**
 * Get the parameters from a Gumbel distribution.
 *
 * @param a is a pointer to a real containing on output the value of the
 *   location parameter.
 * @param b is a pointer to a real containing on output the value of the
 *   scale parameter.
 */
void ran_gumbel_get_parameters(const Random *random,
				      double *a, double *b)
{
  *a = ((double *)(random->parameters))[0];
  *b = ((double *)(random->parameters))[1];
}

/**
 * Set the parameters of a Gumbel distribution.
 *
 * @param a is the value of the location parameter.
 * @param b is the value of the scale parameter.
 */
void ran_gumbel_set_parameters(Random *random, double a, double b)
{
  double temp[2];

  temp[0] = a;
  temp[1] = b;

  ran_random_set_parameters(random, temp, 2*sizeof(double));
}
