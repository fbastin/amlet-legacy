/**
 * ran_bernouilli.h
 *
 * ORATIO Library - Random Numbers Generation module
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_RAN_BERNOUILLI_H__
#define __ORATIO_RAN_BERNOUILLI_H__

/**
 * Bernouilli distribution
 */

#include "ran_rng.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Draw one value from a Bernouilli variate of parameter p.
   *
   * @param p is the Bernouilli parameter, i.e. P[B = 1] = p.
   *
   * @return one value drawn from a Bernouilli variate of parameter p.
   */
  int ran_bernouilli(const Random *random, double p);

  /**
   * Create a Bernouilli generator.
   *
   * @param seed is the seed vector used to initialize the random
   *   numbers generator. A copy of this vector is stored inside the
   *  generator.
   * @param n is the dimension of seed.
   * @param p is the Bernouilli parameter, i.e. P[B = 1] = p.
   */
  Random *ran_bernouilli_new(unsigned long *seed, int n, double p);

  /**
   * Compute the valeur of a Bernouilli variable of parameter p for a
   * realisation of a random variable X~U[0,1] equal to u.
   *
   * @param u is a real in [0,1].
   * @param p is the Bernouilli parameter, i.e. P[B = 1] = p.
   *
   * @return 1 if u <= p, 0 otherwise.
   */
  int ran_bernouilli_from_U(double u, double p);

  int ran_bernouilli_get_val(const Random *random);

  double ran_bernouilli_get_param(const Random *random);

  /**
   * @param p is the Bernouilli parameter, i.e. P[B = 1] = p.
   */
  void ran_bernouilli_set_param(Random *random, double p);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_RAN_BERNOUILLI_H__ */
