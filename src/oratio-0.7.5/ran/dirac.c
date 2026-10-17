/**
 * Random Numbers Library
 * Copyright (C) 2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Generation of points following a Dirac distribution.
 */

#include "oratio/ran.h"

#include <math.h>

/**
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 * @param n is the dimension of the vector seed.
 */
Random *ran_dirac_new(unsigned long *seed, int n, double a)
{
  Random *random = ran_random_new_with_size(seed, n, sizeof(double));

  *((double *)(random->parameters)) = a;

  random->type = RAN_DIRAC;

  return random;
}

void ran_dirac_get_parameters(const Random *random, double *a)
{
  *a = *((double *)(random->parameters));
}

void ran_dirac_set_parameters(Random *random, double a)
{
  ran_random_set_parameters(random, &a, sizeof(double));
}

/**
 * Generate points following a Dirac distribution, i.e. a constant.
 */
double ran_dirac(const Random *random, double a)
{
  return a;
}

double ran_dirac_get_val(const Random *random)
{
  double a;

  ran_dirac_get_parameters(random, &a);

  return ran_dirac(random, a);
}
