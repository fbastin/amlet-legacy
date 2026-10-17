/**
 * mt19937.h
 *
 * Random Numbers Library
 * Copyright (C) 2004-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __RAN_MT19937_H__
#define __RAN_MT19937_H__

#ifdef __cplusplus
extern "C" {
#endif

#define MT19937_STATE_LENGTH 624

  typedef struct {
    /* the array for the state vector  */
    unsigned long mt[MT19937_STATE_LENGTH];
  } MT19937;

  /**
   * Initialize the Mersenne twister generator with a specified seed vector.
   *
   * @param seed is the seed vector.
   * @param n is the dimension of seed.
   */
  void ran_mt19937_init(MT19937 *mt, unsigned long *seed, int n);

  /**
   * Initialize the MT19937 generator with a seed defined as a single integer.
   *
   * @param mt is an allocated MT19937 random numbers generator.
   * @param s is a long unsigned integer, different from 0.
   */
  void ran_mt19937_init_genrand(MT19937 *mt, unsigned long s);

  /**
   * Initialize the MT19937 generator with a seed defined as a vector
   * of integers.
   *
   * @param mt is an allocated MT19937 random numbers generator.
   * @param s is a long unsigned integer, different from 0.
   */
  void ran_mt19937_init_by_array(MT19937 *mt,
				 unsigned long init_key[], int key_length);

  /* generates a random number on [0,0xffffffff]-interval */
  unsigned long ran_mt19937_genrand_int32(MT19937 *mt);

  /* generates a random number on [0,0x7fffffff]-interval */
  long ran_mt19937_genrand_int31(MT19937 *mt);

  /* generates a random number on [0,1]-real-interval */
  double ran_mt19937_genrand_real1(MT19937 *mt);

  /* generates a random number on [0,1)-real-interval */
  double ran_mt19937_genrand_real2(MT19937 *mt);

  /* generates a random number on (0,1)-real-interval */
  double ran_mt19937_genrand_real3(MT19937 *mt);

  /* generates a random number on [0,1) with 53-bit resolution*/
  double ran_mt19937_genrand_res53(MT19937 *mt);

  double ran_mt19937_get_val(void *mt);

  /**
   * Retrieve the generator state.
   */
  unsigned long *ran_mt19937_get_state(void *mt);

  /**
   * UNTESTED FUNCTION!
   */
  void ran_mt19937_set_state(MT19937 *mt, void *state);

  MT19937 *ran_mt19937_new(unsigned long *seed, int n);

  MT19937 *ran_mt19937_new_with_state(unsigned long *state);

  void ran_mt19937_free(void *mt);

#ifdef __cplusplus
}
#endif

#endif  /* __RAN_MT19937_H__ */
