/**
 * lecuyer88a.h
 *
 * ORATIO Library
 * Copyright (C) 2002-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 *
 */

#ifndef __RAN_LECUYER88A_H__
#define __RAN_LECUYER88A_H__

/**
 * Ref.: Pierre L'Ecuyer, Efficient and portable combined random
 * numbers generators, Communications of the ACM, 31(6), pp 742-749,
 * 1988.
 */

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct  {
    /* State of the combined LCG generator.
     * The first number correspond to the state associated to the
     * first LCG, the second one to the second MRG. */
    unsigned long s[2];
  } LEcuyer88a;

  /**
   * Generation of a pseudo uniform U(0,1) number, using one the 1988
   * L'Ecuyer's algorithm.
   *
   * @param r is a point to a LEcuyer88b structure. It is declared as
   *   null to conform to the common definition of rng_get_val.
   *
   * @return a real number in (0,1), from a pseudo-uniform distribution U(0,1).
   */
  double ran_lecuyer88a_get_val(void *r);

  /**
   */
  unsigned long *ran_lecuyer88a_get_state(void *random);

  /**
   * Set the state of a random numbers generator using 1988 L'Ecuyer's
   * algorithm.
   *
   * @param rng is an initialized 1988 L'Ecuyers's algorithm.
   * @param state is an array of at least two long numbers.
   */
  void ran_lecuyer88a_set_state(LEcuyer88a *rng, unsigned long *state);

  void ran_lecuyer88a_init(LEcuyer88a *rng, unsigned long *seed, int n);

  LEcuyer88a *ran_lecuyer88a_new(unsigned long *seed, int n);

  /**
   * Free memory associated to the LEcuyer88a random number generator.
   *
   * @param rng is an initialized 1988 L'Ecuyers's algorithm.
   */
  void ran_lecuyer88a_free(void *rng);

#ifdef __cplusplus
}
#endif

#endif  /* __RAN_LECUYER88A_H__ */
