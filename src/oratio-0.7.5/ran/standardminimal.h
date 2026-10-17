/**
 * standardminimal.h
 *
 * ORATIO Library
 * Copyright (C) 2002-2009 Jordan Goblet, Gabriel Schwanen,
 *                         Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 *
 */

#ifndef __RAN_STANDARDMINIMAL_H__
#define __RAN_STANDARDMINIMAL_H__

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct  {
    unsigned long state;
  } StandardMinimal;

  /**
   * Generation of a pseudo-random number between 0 and 1 using a
   * generic linear congruential algorithm.
   */
  double ran_standardminimal_get_val(void *standardminimal);

  /**
   * Create a new random number generator
   *
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   */
  StandardMinimal *ran_standardminimal_new(unsigned long seed);

  void ran_standardminimal_init(Random *random, unsigned long seed);

  /**
   * Define the generator state.
   */
  void ran_standardminimal_set_state(void *standardminimal, void *state);

  /**
   * Return the current state of a standard minimal generator.
   */
  unsigned long *ran_standardminimal_get_state(void *standardminimal);

  /**
   * Make a new generator of random number by precising the seed of the
   * sequence who will be generated
   * strictement positif; dans le cas contraire, elle ne sera pas prise en
   * compte.
   */
  StandardMinimal *ran_standardminimal_new_with_state(unsigned long *state);

  StandardMinimal *ran_standardminimal_new_with_seed(unsigned long seed);

  /**
   * Free memory associated to the Park and Miller random numbers generator.
   */
  void ran_standardminimal_free(void *standardminimal);

#ifdef __cplusplus
}
#endif

#endif  /* __RAN_STANDARDMINIMAL_H__ */
