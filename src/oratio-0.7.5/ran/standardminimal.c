/**
 * standardminimal.c
 *
 * ORATIO Library
 * Copyright (C) 2002-2007 Jordan Goblet, Gabriel Schwanen,
 *                         Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/ran.h"
#include "standardminimal.h"

#include <oratio/ctools.h>

/**
 * Standard Minimal algorithm constants.
 */
#define a      16807
#define m      2147483647    /* 2^31-1 */
#define AM_MIL (1.0/m)
#define q1     127773
#define r1     2836

/**
 * Return a U(0,1) pseudo-random number using the Standard Minimal algorithm.
 *
 * @param standardminimal is an initialized Standard Minimal generator.
 *
 * @return a U(0,1) pseudo-random number using the Standard Minimal algorithm.
 */
double ran_standardminimal_get_val(void *standardminimal)
{
  long x, y;

  x=((StandardMinimal *)standardminimal)->state;

  y=x/q1;
  x=a*(x-y*q1)-y*r1;
  if (x < 0) x += m;

  /* We update the state of the Standard Minimal generator */
  ((StandardMinimal *)standardminimal)->state = x;

  return AM_MIL*x;
}

StandardMinimal *ran_standardminimal_new_with_state(unsigned long *state)
{
  StandardMinimal *standardminimal;

  standardminimal = _malloc(StandardMinimal);
  ran_standardminimal_set_state(standardminimal, state);

  return standardminimal;
}

/**
 * Create a new random number generator
 *
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 */
StandardMinimal *ran_standardminimal_new_with_seed(unsigned long seed)
{
  return ran_standardminimal_new_with_state(&seed);
}

/**
 * Return the current state of a standard minimal generator.
 */
unsigned long *ran_standardminimal_get_state(void *standardminimal)
{
  return &(((StandardMinimal*) standardminimal)->state);
}

/**
 * Define the generator state.
 */
void ran_standardminimal_set_state(void *standardminimal, void *state)
{ 
  if (state) {
    unsigned long seed = *((unsigned long *)state);

    if (seed > 0)
      ((StandardMinimal *)standardminimal)->state = seed;
    else
      ((StandardMinimal *)standardminimal)->state = 12345;
  }
  else
    ((StandardMinimal *)standardminimal)->state = 12345;

  state = &(((StandardMinimal *)standardminimal)->state);
}

void ran_standardminimal_init(Random *random, unsigned long seed)
{
  ran_standardminimal_set_state(random, &seed);
}

/**
 * Create a new random number generator based on the Standard Minimal LCG.
 *
 * @param seed is the seed used to initialize the random numbers generator.
 */
StandardMinimal *ran_standardminimal_new(unsigned long seed)
{ 
  return ran_standardminimal_new_with_state(&seed);
}

/**
 * Free memory associated to the Park and Miller random numbers generator.
 *
 * @param mrg must be an initialized pointer to a StandardMinimal structure.
 */
void ran_standardminimal_free(void *standardminimal)
{
  free((StandardMinimal *)standardminimal);
}
