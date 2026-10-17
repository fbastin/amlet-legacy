/**
 * poisson.c
 *
 * ORATIO Library
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/ran.h"
#include "oratio/stat.h"

#include <string.h>

/**
 * Return a value drawn from an Poisson distribution of parameter lambda.
 *
 * @param ran is an initialized random number generator.
 * @param lambda is the scale parameter of the poisson distribution;
 *   it corresponds to the inverse of the mean.
 *
 * @return a number drawn from an exponential variate of parameter lambda.
 */
double ran_poisson(const Random *ran, double lambda)
{
  return st_poisson_icdf(ran_random_get_val(ran), lambda);
}
