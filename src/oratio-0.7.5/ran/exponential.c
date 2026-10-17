/**
 * exponential.c
 *
 * ORATIO Library
 * Copyright (C) 1999-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/ran.h"
#include "oratio/stat.h"

#include <string.h>

/**
 * Return a value drawn from an exponential distribution of parameter lambda.
 *
 * @param ran is an initialized random number generator.
 * @param lambda is the scale parameter of the exponential distribution;
 *   it corresponds to the inverse of the mean.
 *
 * @return a number drawn from an exponential variate of parameter lambda.
 */
double ran_exponential(const Random *ran, double lambda)
{
  return st_exponential_icdf(ran_random_get_val(ran), lambda);
}

/**
 * @param lambda is the scale parameter of the exponential distribution;
 *   it corresponds to the inverse of the mean.
 */
Random *ran_exponential_new(unsigned long *seed, int n, double lambda)
{
  Random *random;

  random = ran_random_new_with_size(seed, n, sizeof(double));
  memcpy(random->parameters, &lambda, sizeof(double));

  random->type = RAN_EXPONENTIAL;

  return random;
}

/**
 * Generate an exponential variate realization, given an exponential random
 * number generator.
 *
 * @param ran is an initialized random number generator, of type exponential.
 *
 * @return a realization from the exponential random number generator.
 */
double ran_exponential_get_val(const Random *random)
{
  return ran_exponential(random, *((double *)(random->parameters)));
}

/**
 */
double ran_exponential_get_parameters(const Random *random)
{
  return *((double *)(random->parameters));
}

/**
 */
void ran_exponential_set_parameters(Random *random, double lambda)
{
  ran_random_set_parameters(random, &lambda, sizeof(double));
}
