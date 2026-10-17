/**
 * halton.h
 *
 * ORATIO Library - MOCA module
 * Copyright (C) 2005-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_MOCA_HALTON_H__
#define __ORATIO_MOCA_HALTON_H__

#include <stdio.h>

#include "oratio/ran.h"

#ifdef __cplusplus
extern "C" {
#endif

  enum {
    MOCA_HALTON_NO_SCRAMBLING,
    MOCA_HALTON_BW_SCRAMBLING
  };

  enum {
    MOCA_HALTON_RANDOMIZATION_NONE = 0,
    MOCA_HALTON_RANDOMIZATION_SCRAMBLING,
    MOCA_HALTON_RANDOMIZATION_SHIFTING,
    MOCA_HALTON_RANDOMIZATION_SHIFTING_SCRAMBLING,
    MOCA_HALTON_RANDOMIZATION_SHUFFLING,
    MOCA_HALTON_RANDOMIZATION_SHIFTING_SHUFFLING
  };

  typedef struct {
    int *scrambling;      /* scrambling matrix */
    int lds;
    int maxdim;           /* maximum dimension allowed by the sampling */
    int in;               /* counter */
    unsigned int s;       /* user-defined dimension */
    unsigned int *primes; /* s first prime numbers */
  } MOCA_Halton;

  int moca_halton_std_permutations_max_dim_get();

  int moca_halton_std_permutations_prime_max_dim_get();

  int *moca_halton_std_permutations();

  /**
   * Reset a Halton draws generator.
   *
   * @param h is an initialized Halton generator.
   */
  void moca_halton_reset(MOCA_Halton *h);

  /**
   * @param prime is the prime used to construct the Halton draw.
   * @param n is the integer to transform in a Halton draw.
   */
  double moca_halton(int prime, int s, int n, int *scrambling, int lds);

  /**
   * Generate a Halton draws sample.
   *
   * @param h is a Halton draws generator.
   * @param n is the numbers of vectors forming the sample.
   * @param S contains on output the halton sample.
   * @param ldS is the leading dimension of S; ldS must be grater or equal to
   *   sample dimension, and is typically equal to this dimesion.
   */
  void moca_halton_sample(MOCA_Halton *h, int n, double *S, int ldS);

  /**
   * Generate a sample based on Halton sequences
   *
   * @param ran is random number generator.
   * @param loose is equal to 1 if for some errors can be managed by
   *   only raising a warning, 0 else.
   */
  void moca_halton_generate_sample(Random *ran, void *sample,
                                   unsigned int s, int n,
                                   int randomization, int loose);

  /**
   * Create a new Halton sequence generator.
   *
   * @param s is the dimension of the vectors to generate.
   * @param sid an identifier of the scrambling method to use, 0 meaning
   *   that no scrambling has to be performed.
   *
   * @return a newly allocate MOCA_Halton structure.
   */
  MOCA_Halton *moca_halton_new(unsigned int s, int sid);

  /**
   * Free the memory associated to a MOCA_Halton structure.
   *
   * @param h is an initialized MOCA_Halton structure.
   */
  void moca_halton_free(MOCA_Halton *h);

  /**
   * Generate the next point of a Halton sequence.
   *
   * @param h is a MOCA_Halton initialized generator.
   * @param x is a double vector that contains on output the next point
   *   of the Halton sequence.
   */
  void moca_halton_next(MOCA_Halton *h, double *x);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_MOCA_HALTON_H__ */
