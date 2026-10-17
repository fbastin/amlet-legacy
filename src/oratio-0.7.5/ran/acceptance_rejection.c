/**
 * acceptance_rejection.c
 *
 * ORATIO Library
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License.
 */

#include "oratio/ran.h"

/**
 * Generate a realization of the random variable X by mean of
 * acceptance-rejection technique.
 *
 * @param ran1 is an initialized random numbers generator.
 * @param ran2 is an initialized random numbers generator.
 * @param params_f is the parameters vector associated to the density of X.
 * @param f is the density of X.
 * @param params_t is the parameters vector for the majorant function.
 * @param t is the majorant function.
 *
 * @return a randomly generated realization of X.
 */
double ran_acceptance_rejection(Random *ran1, Random *ran2,
				double *params_f, double *params_t,
				double (*f)(double x, double *params),
				double (*t)(double x, double *params))
{
  double y;

  do {
    y = ran_random_next(ran1);
  } while (ran_random_get_val(ran2) > f(y, params_f)/t(y, params_t));

  return y;
}
