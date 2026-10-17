/**
 * mrg32k3a.h
 *
 * Random Numbers Library
 * Copyright (C) 2002-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 *
 */

#ifndef __RAN_MRG32K3A_H__
#define __RAN_MRG32K3A_H__

/**
 * Ref.: Pierre L'Ecuyer, Good parameters and implementations for
 * combined multiple recursive random number generators, Operations
 * Research, 47(1), pp. 159-164, 1999.
 */

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct  {
    /* State of the combined MRG generator.
       The three first numbers correspond to the state associated to
       the first MRG, the three next ones to the second MRG. */
    double s[6];
    double bs[6];
    double is[6];
  } MRG32K3a;

  /**
   * Return a pseudo U(0,1)-random number following the MRG32KA algorithm.
   */
  double ran_mrg32k3a_get_val(void *mrg);

  /**
   * Get the state of the random number generator.
   *
   * @param mrg is an initialized MRG32K3a random number generator.
   *
   * @return the current of the generator mrg.
   */
  double *ran_mrg32k3a_get_state(MRG32K3a *mrg);

  /**
   */
  void ran_mrg32k3a_set_state(MRG32K3a *mrg, void *state);

  /**
   */
  void ran_mrg32k3a_init(MRG32K3a *mrg, unsigned long *seed, int n);

  /**
   * Create a new mrg based on the MRG32K3a method.
   *
   * @param seed is the seed of the generater. It can be set to NULL, 
   * @param n is the size of the vector pointed by seed. n must be equal
   *   to 1 or 6.
   *
   * @return a pointer to an initialized MRG32K3a structure if n is
   *   equal to 1 or 6, NULL otherwise.
   */
  MRG32K3a *ran_mrg32k3a_new(unsigned long *seed, int n);

  /**
   * Free memory associated to the MRG32K3a random numbers generator.
   *
   * @param mrg must be an initialized pointer to a MRG32K3a structure.
   */
  void ran_mrg32k3a_free(void *mrg);

  /**
   * Reset the current stream to its start.
   *
   * @param mrg is an initialized MRG32K3a random numbers.
   */
  void ran_mrg32k3a_resetstartstream (MRG32K3a *mrg);

  void ran_mrg32k3a_resetstartsubstream (MRG32K3a *mrg);

  void ran_mrg32k3a_resetnextsubstream(MRG32K3a *mrg);

#ifdef __cplusplus
}
#endif

#endif  /* __RAN_MRG32K3A_H__ */
