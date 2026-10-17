/**
 * arcsine.c
 *
 * Random Numbers Library
 * Copyright (C) 1999-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/ran.h"

/**
 * Arcsine distribution
 */

/**
 * Create a new arcine generator.
 *
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 * @param n is the dimension of seed.
 */
Random *ran_arcsine_new(unsigned long *seed, int n)
{
  Random *r;

  r = ran_random_new(0, seed, n);
  r->type = RAN_ARCSINE;

  return r;
}

double ran_arcsine(const Random *random)
{
  double q;

  q = sin(M_PI_2*ran_random_get_val(random));

  return q*q;
}

double ran_arcsine_get_val(const Random *random)
{
  return ran_arcsine(random);
}
