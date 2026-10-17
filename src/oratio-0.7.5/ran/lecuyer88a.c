/**
 * lecuyer88a.c
 *
 * ORATIO Library
 * Copyright (C) 2007-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/ran.h"
#include "lecuyer88a.h"

#include <string.h>

#include <oratio/memory.h>

#define m1  2147483563
#define a1  40014
#define q1  53668
#define r1  12211
#define m2  2147483399
#define a2  40692
#define q2  52774
#define r2  3791
#define mm  2147483562
#define imm 4.656613e-10  /* 1/mm */

/**
 * Generation of a pseudo uniform U(0,1) number, using one the 1988
 * L'Ecuyer's algorithm.
 *
 * @param r is a point to a LEcuyer88b structure. It is declared as
 *   null to conform to the common definition of rng_get_val.
 *
 * @return a real number in (0,1), from a pseudo-uniform distribution U(0,1).
 */
double ran_lecuyer88a_get_val(void *r)
{
  long k, z;
  long s1, s2;

  LEcuyer88a *ecu = (LEcuyer88a *)r;

  s1 = ecu->s[0];
  s2 = ecu->s[1];

  /* First LCG */
  k = s1/q1;
  s1 = a1*(s1-k*q1)-k*r1;
  if (s1 < 0) s1 += m1;

  /* Second LCG */
  k = s2/q2;
  s2 = a2*(s2-k*q2)-k*r2;
  if (s2 < 0) s2 += m2;

  /* Combine the LCGs */
  z = s1-s2;
  if (z < 1) z += mm;

  /* Save the new state. */
  ecu->s[0] = s1;
  ecu->s[1] = s2;

  return z*imm;
}

/**
 */
unsigned long *ran_lecuyer88a_get_state(void *rng)
{
  return ((LEcuyer88a *)rng)->s;
}

/**
 * Set the state of a random numbers generator using 1988 L'Ecuyer's
 * algorithm.
 *
 * @param ran is an initialized 1988 L'Ecuyers's algorithm.
 * @param state is an array of at least two long numbers.
 */
void ran_lecuyer88a_set_state(LEcuyer88a *rng, unsigned long *state)
{
  memcpy(rng->s, state, 2*sizeof(double));
}

void ran_lecuyer88a_init(LEcuyer88a *rng, unsigned long *seed, int n)
{
  if (n == 1)
    rng->s[0] = rng->s[1] = *seed;
  else /* n == 2 */
    ran_lecuyer88a_set_state(rng, seed);
}

LEcuyer88a *ran_lecuyer88a_new(unsigned long *seed, int n)
{
  LEcuyer88a *rng;

  if (n < 1 || n > 2) return NULL;

  rng = _malloc(LEcuyer88a);

  if (seed) {
    int i;

    for (i = 0; i < n; i++) if (seed[i] == 0) seed[i] = ran_init_random();
    ran_lecuyer88a_init(rng, seed, n);
  }
  else {
    unsigned long s;

    s = ran_init_random();
    ran_lecuyer88a_init(rng, &s, 1);
  }

  return rng;
}

/**
 * Free memory associated to the LEcuyer88a random number generator.
 */
void ran_lecuyer88a_free(void *rng)
{
  if (rng) free((LEcuyer88a *)rng);
}
