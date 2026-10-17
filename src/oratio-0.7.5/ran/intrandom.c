/**
 * intrandom.c
 *
 * Random Numbers Library
 * Copyright (C) 1999-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/ran.h"

#include <stdlib.h>
#include <string.h>

/**
 * Create a new random numbers generator dedicated to the production
 * of integer random numbers.
 *
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 * @param n is the dimension of seed.
 */
Random *ran_intrandom_new(unsigned long *seed, int n,
			  unsigned long max)
{
  Random *random;

  random = ran_random_new_with_size(seed, n, sizeof(unsigned long));
  memcpy(random->parameters, &max, sizeof(unsigned long));

  random->type = RAN_INTEGER;

  return random;
}

/**
 * Generate a pseudo-random random integer in [0, max[
 *
 * @param r is a initialized random numbers generator.
 * @param max is a strictly positive integer.
 *
 * @return an integer from an uniform integer distribution on [0, max[.
 */
unsigned long ran_intrandom(const Random *r, unsigned long max)
{
  return (unsigned long)(ran_random_get_val(r)*max);
}

/**
 * Generate a pseudo-random random integer in [0, max[
 *
 * @param ird is an initialized integer random numbers generator.
 *
 * @return an integer from an uniform integer distribution on [0, max[,
 *   where max is the parameter of the random number generator.
 */
unsigned long ran_intrandom_get_val(const Random *ird)
{
  return ran_intrandom(ird, *((unsigned long *)(ird->parameters)));
}

unsigned long ran_intrandom_get_max(const Random *ird)
{
  return *((unsigned long *)(ird->parameters));
}

void ran_intrandom_set_max(Random *ird, unsigned long max)
{
  ran_random_set_parameters(ird, &max, sizeof(unsigned long)); 
}
