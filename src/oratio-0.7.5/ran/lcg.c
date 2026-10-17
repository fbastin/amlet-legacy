/**
 * lcg.c
 *
 * ORATIO Library
 * Copyright (C) 2009-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License
 */

#include "oratio/ran.h"
#include "oratio/memory.h"

/**
 * Return a U(0,1) pseudo-random number using a Linear Congruential Generator.
 *
 * @param lcg is an initialized LCG generator.
 *
 * @return a U(0,1) pseudo-random number using the Standard Minimal algorithm.
 */
double ran_lcg_get_val(void *lcg)
{
  double x;
  unsigned long k;
  double m = ((LCG *)lcg)->m;

  x = ((LCG *)lcg)->a*((LCG *)lcg)->state+((LCG *)lcg)->c;
  k = x/m;
  x -= k*m;
  if (x < 0.0) x += m;

  ((LCG *)lcg)->state = x;

  return x*((LCG *)lcg)->norm;
}

/**
 * Create a new LCG with a given initial state.
 *
 * @param state is the initial state.
 */
LCG *ran_lcg_new_with_state(unsigned long *state)
{
  LCG *lcg;

  lcg = _malloc(LCG);
  ran_lcg_set_state(lcg, state);
  /* We set the generator to be the Standard Minimal by default. */
  lcg->a = 16807;
  lcg->c = 0;
  lcg->m = 2147483647;
  lcg->norm = 1.0/lcg->m;

  return lcg;
}

/**
 * Set the parameters of the LCG, of the form (ax+c) mod m.
 *
 * @param a is the multiplicator of x.
 * @param c is the term added to a*x.
 */
void ran_lcg_set_parameters(LCG *lcg, double a, double c, double m)
{
  lcg->a = a;
  lcg->c = c;
  lcg->m = m;
  lcg->norm = 1.0/m;
}

/**
 * Create a new random number generator
 *
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 *
 * @return a pointer to the newly created random numbers generator.
 */
LCG *ran_lcg_new_with_seed(unsigned long seed)
{
  return ran_lcg_new_with_state(&seed);
}

/**
 * Return the current state of a standard minimal generator.
 *
 * @param lcg is an initialized LCG generator.
 */
unsigned long *ran_lcg_get_state(void *lcg)
{
  return &(((LCG *)lcg)->state);
}

/**
 * Define the generator state.
 *
 * @param lcg is an initialized LCG generator.
 */
void ran_lcg_set_state(void *lcg, void *state)
{
  if (state) {
    unsigned long seed = *((unsigned long *)state);

    if (seed > 0)
      ((LCG *)lcg)->state = seed;
    else
      seed = ((LCG *)lcg)->state = 12345;
  }
  else
    ((LCG *)lcg)->state = 12345;

  state = &(((LCG *)lcg)->state);
}

/**
 * Initialize a LCG random numbers generator.
 */
void ran_lcg_init(Random *random, unsigned long seed)
{
  ran_lcg_set_state(random, &seed);
}

/**
 * Create a new random number generator based on the Standard Minimal LCG.
 *
 * @param seed is the seed used to initialize the random numbers generator.
 */
LCG *ran_lcg_new(unsigned long seed)
{
  return ran_lcg_new_with_state(&seed);
}

/**
 * Free memory associated to the LCG random numbers generator.
 *
 * @param lcg must be an initialized pointer to a LCG structure.
 */
void ran_lcg_free(void *lcg)
{
  free((LCG *)lcg);
}
