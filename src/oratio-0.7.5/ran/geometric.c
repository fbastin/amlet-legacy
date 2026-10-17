/**
 * geometric.c
 *
 * Random Numbers Library
 * Copyright (C) 2005-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Geometric distribution (Bernoulli trial with probability p) 
 *
 * prob(k) =  p (1 - p)^(k) for k = 0, 2, 3, ...
 *
 * Ref.: Luc Devroye, "Non-Uniform Random Variate Generation",
 * Springer-Verlag, 1986, pp. 499-500.
 */

#include "oratio/ran.h"
#include "oratio/st_geometric.h"

#include <string.h>

/**
 * Return a value drawn from a Geometric variate of parameter p.
 */
int ran_geometric(const Random *random, double p)
{
  if (p == 0 || p == 1) return p;
  return st_geometric_icdf(ran_random_get_val(random), p);
}

/**
 * Return a value drawn from a Geometric variate of parameter p.
 *
 * @param seed is the seed vector used to initialize the random
 *   numbers generator. A copy of this vector is stored inside the
 *  generator.
 * @param n is the dimension of seed.
 * @param p is the geometric parameter.
 */
Random *ran_geometric_new(unsigned long *seed, int n, double p)
{
  Random *random;

  random = ran_random_new_with_size(seed, n, sizeof(double));
  memcpy(random->parameters, &p, sizeof(double));

  random->type = RAN_GEOMETRIC;

  return random;
}

double ran_geometric_get_param(const Random *random)
{
  return (*(double *)(random->parameters));
}

int ran_geometric_get_val(const Random *random)
{
  return ran_geometric(random, ran_geometric_get_param(random));
}

void ran_geometric_set_param(Random *random, double p)
{
  ran_random_set_parameters(random, &p, sizeof(double));
}
