/**
 * laplace.c
 *
 * ORATIO Library
 * Copyright (C) 1999-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Random numbers generation for Laplace distribution.
 */

#include "oratio/ran.h"
#include "oratio/stat.h"

/**
 * Return a value drawn from a Laplace distribution of parameter mu and b.
 *
 * @param ran is an initialized random number generator.
 * @param mu is the location parameter.
 * @param b is the scale parameter.
 *
 * @return a number drawn from an Laplace variate of parameter mu and b.
 */
double ran_laplace(const Random *ran, double mu, double b)
{
  return st_laplace_icdf(ran_random_get_val(ran), mu, b);
}

/**
 * Create a Laplace random variate generator.
 *
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 * @param n is the dimension of the vector seed.
 * @param mu is the location parameter.
 * @param b is the scale parameter.
 *
 * @return a pointer to the Laplace random numbers generator.
 */
Random *ran_laplace_new(unsigned long *seed, int n,
			double mu, double b)
{
  Random *random = ran_random_new_with_size(seed, n, 2*sizeof(double));

  ((double *)(random->parameters))[0] = mu;
  ((double *)(random->parameters))[1] = b;

  random->type = RAN_LAPLACE;

  return random;
}

/**
 * Draw a number from a Laplace distribution.
 *
 * @return a draw from a Laplace distribution.
 */
double ran_laplace_get_val(const Random *random)
{
  double mu, b;

  ran_laplace_get_parameters(random, &mu, &b);
  return ran_laplace(random, mu, b);
}

/**
 * Get the parameters of a Laplace distribution.
 *
 * @param mu is a pointer to a real containing on output the value of the
 *   location parameter.
 * @param b is a pointer to a real containing on output the value of the
 *   scale parameter.
 */
void ran_laplace_get_parameters(const Random *random,
				       double *mu, double *b)
{
  *mu = ((double *)(random->parameters))[0];
  *b = ((double *)(random->parameters))[1];
}

/**
 * Set the parameters of a Laplace distribution.
 *
 * @param mu is the value of the location parameter.
 * @param b is the value of the scale parameter.
 */
void ran_laplace_set_parameters(Random *random, double mu, double b)
{
  double temp[2];

  temp[0] = mu;
  temp[1] = b;

  ran_random_set_parameters(random, temp, 2*sizeof(double));
}
