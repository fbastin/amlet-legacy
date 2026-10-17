/**
 * bernouilli.c
 *
 * ORATIO Library - Random Numbers Generation module
 * Copyright (C) 2004-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Generation of binary decisions, using a Bernouilli distribution of
 * parameter p, i.e. P[B = 1] = p.
 */

#include "oratio/ran.h"

#include <string.h>

/**
 * Draw one value from a Bernouilli variate of parameter p.
 *
 * @param p is the Bernouilli parameter, i.e. P[B = 1] = p.
 *
 * @return one value drawn from a Bernouilli variate of parameter p.
 */
int ran_bernouilli(const Random *random, double p)
{
  return (ran_random_get_val(random) <= p) ? 1 : 0;
}

/**
 * Create a Bernouilli generator.
 *
 * @param seed is the seed vector used to initialize the random
 *   numbers generator. A copy of this vector is stored inside the
 *  generator.
 * @param n is the dimension of seed.
 * @param p is the Bernouilli parameter, i.e. P[B = 1] = p.
 */
Random *ran_bernouilli_new(unsigned long *seed, int n,
			   double p)
{
  Random *random;

  random = ran_random_new_with_size(seed, n, sizeof(double));
  memcpy(random->parameters, &p, sizeof(double));

  random->type = RAN_BERNOUILLI;

  return random;
}

/**
 * Compute the valeur of a Bernouilli variable of parameter p for a
 * realisation of a random variable X~U[0,1] equal to u.
 *
 * @param u is a real in [0,1].
 * @param p is the Bernouilli parameter, i.e. P[B = 1] = p.
 *
 * @return 1 if u <= p, 0 otherwise.
 */
int ran_bernouilli_from_U(double u, double p)
{
  return (u <= p) ? 1 : 0;
}

int ran_bernouilli_get_val(const Random *random)
{
  return (random->rng_get_val(random->rng)
	  <= (*(double *)(random->parameters))) ? 1 : 0;
}

double ran_bernouilli_get_param(const Random *random)
{
  return (*(double *)(random->parameters));
}

/**
 * @param p is the Bernouilli parameter, i.e. P[B = 1] = p.
 */
void ran_bernouilli_set_param(Random *random, double p)
{
  ran_random_set_parameters(random, &p, sizeof(double));
}
