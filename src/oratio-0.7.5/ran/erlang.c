/**
 * erlang.c
 *
 * ORATIO Library
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * NOT YET IMPLEMENTED!
 */

#include "oratio/ran.h"

#include <string.h>

/**
 * Return a value drawn from an erlang of parameter b.
 */
double ran_erlang(const Random *random, double alpha, double lambda)
{
  return 0.0;
}

Random *ran_erlang_new(unsigned long *seed, int n,
		       double alpha, double lambda)
{
  Random *random;

  random = ran_random_new_with_size(seed, n, 2*sizeof(double));

  /*  memcpy(random->parameters, &b, sizeof(double));
      memcpy(random->parameters, &b, sizeof(double)); */

  random->type = RAN_ERLANG;

  return random;
}

double ran_erlang_get_val(const Random *random)
{
  return 0.0;
}

double ran_erlang_get_parameters(const Random *random)
{
  return 0.0;
}

void ran_erlang_set_parameters(Random *random, double b)
{
}
