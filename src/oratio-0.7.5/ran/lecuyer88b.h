/**
 * lecuyer88b.h
 *
 * Random Numbers Library
 * Copyright (C) 1999-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __RAN_LECUYER88B_H__
#define __RAN_LECUYER88B_H__

#ifdef __cplusplus
extern "C" {
#endif

#define ntab 32

  typedef struct  {
    /* State of the combined LCG generator.
     * The first number correspond to the state associated to the
     * first LCG, the second one to the second MRG. */
    unsigned long s[2];

    long iy;
    long iv[ntab];
  } LEcuyer88b;

  /**
   * Return a pseudo U(0,1)-random number following the L'Ecuyer's
   * 1988 algorithm and a postponed output technique.
   *
   * @param r is a point to a LEcuyer88b structure. It is declared as
   *   null to conform to the common definition of rng_get_val.
   */
  double ran_lecuyer88b_get_val(void *r);

  /**
   */
  unsigned long *ran_lecuyer88b_get_state(void *random);

  void ran_lecuyer88b_set_state(LEcuyer88b *rng, unsigned long *state);

  void ran_lecuyer88b_init(LEcuyer88b *rng, unsigned long *seed, int n);

  LEcuyer88b *ran_lecuyer88b_new(unsigned long *seed, int n);

  void ran_lecuyer88b_free(LEcuyer88b *rng);

#ifdef __cplusplus
}
#endif

#endif  /* __RAN_ECUYER88B_H__ */
