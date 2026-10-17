/**
 * ran_rng.h
 *
 * ORATIO Library - Random Numbers Generation module
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

/**
 * Random numbers generator.
 */

#ifndef __ORATIO_RAN_RNG_H__
#define __ORATIO_RAN_RNG_H__

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Uniform random numbers generator identifiers
   */
  enum {
    RAN_LECUYER88A = 1,
    RAN_LECUYER88B,
    RAN_MRG32K3A,
    RAN_MT19937,
    RAN_STANDARDMINIMAL,
    END_RANDOM_GENERATORS
  };

#define RAN_DEFAULT RAN_MRG32K3A
#define RAN_DEFAULT_SEED_LENGTH 6

  /**
   * Reproduce, or not, the MRG32k3a generator of ORATIO <= 0.7.5.
   *
   * Up to ORATIO 0.7.5, the first component of MRG32k3a used x(n-1)
   * instead of x(n-2), and the jumps to the next stream or substream were
   * not reduced modulo m (see nt_modular_set_legacy()). The sequences
   * differed from L'Ecuyer's MRG32k3a.
   *
   * @param legacy is 0 (default) for L'Ecuyer's MRG32k3a, or 1 to
   *   reproduce the random draws of ORATIO <= 0.7.5. It must be set
   *   before any random generator is created.
   */
  void ran_mrg32k3a_set_legacy(int legacy);

  /**
   * Uniform random generator structure
   */

  typedef struct _RandomGenerator RandomGenerator;

  /**
   * Random variables generation
   */

  typedef struct _Random Random;

  struct _Random {
    void *rng;

    int type;            /* type of the random numbers generator */
    int id;              /* uniform random numbers generator identificator */

    void *parameters;

    unsigned long *seed; /* initial seed. */

    double l, h;         /* truncature on the cdf; by default, l=0, h=1. */
    double p;            /* value of P[X = a], where a is a real such
			    that l = F(a). It is set to 0.0 by default. */

    double (*rng_get_val)(void *r);  /* U(0,1) random-numbers generators */

    /* TODO */
    /* int (*rng_get_parameters_number)(void *r) */
  };

  /**
   * Create a new generator using default values.
   *
   * @return the newly created random object.
   */
  Random *ran_random();

  /**
   * Create a new generator with an initial seed.
   *
   * @param id is the generator identifier; if set to 0, RAN_DEFAULT is used.
   * @param seed is the seed vector used to initialize the random
   *   numbers generator. A copy of this vector is stored inside the
   *   generator. seed can set to NULL, in which case a default seed is used.
   * @param n is the dimension of seed (irrelevant is seed is set to NULL).
   *
   * @return the newly created random object.
   */
  Random *ran_random_new(int id, unsigned long *seed, int n);

  /**
   * State gestion.
   */

  /**
   * Return the current state of the generator.
   *
   * @param random is an initialized random numbers generator.
   *
   * @return the state associated to the random numbers generators, or
   *   NULL if the generator type is unidentified.
   */
  void *ran_random_get_state(Random *random);

  /**
   * Set the current state of the generator.
   *
   * @param random is an initialized random numbers generator.
   * @param state is the new state to associate to the random numbers
   *   generator. The precise type must correspond to the associated
   *   random object.
   */
  void ran_random_set_state(Random *random, void *state);

  /**
   * Copy the generator current state.
   *
   * @param random is an initialized random numbers generator.
   * @param n contains on output the state length.
   *
   * @return a newly allocated pointer to a copy of the state of the
   *   generator; the state is store as a unsigned long vector.
   */
  unsigned long *ran_random_copy_state(Random *random, int *n);

  /**
   * Set the uniform random numbers generator associated to a random object.
   *
   * @param random is an initialized random numbers generator.
   * @param id is the generator identifier.
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of seed.
   */
  void ran_random_set_generator(Random *random, const int id,
				unsigned long *seed, int n);

  /**
   * Return the generator name.
   *
   * @param random is an initialized random numbers generator.
   *
   * @return a pointer to the name of the generator. There is no new
   *   memory allocation, so the returned pointer has not to be freed.
   */
  char *ran_random_generator_name(Random *random);

  /**
   * Return the number of parameters needed for a generator of type 'type'
   *
   * @param random is an initialized random numbers generator.
   */
  int ran_random_get_parameters_number(Random *random);

  /**
   * Create a new generator with an initial seed and the parameters
   * vector size.
   *
   * @param seed is the seed vector used to initialize the random
   *   numbers generator. A copy of this vector is stored inside the
   *  generator.
   * @param n is the dimension of seed.
   * @param size is the size in bytes of the parameters vector guiding
   *   the generator behavior.
   *
   * @return the newly created random object.
   */
  Random *ran_random_new_with_size(unsigned long *seed, int n,
				   size_t size);

  /**
   * Draw an uniform random number in (F(a), F(b)).
   *
   * @param random is an initialized random numbers generator.
   */
  double ran_random_get_val(const Random *random);

  /**
   * Return a draw from an univariate random number generator, using the type
   * associated to the generator.
   *
   * @param random is an initialized random numbers generator.
   *
   * @return a draw from the univariate random number generator.
   */
  double ran_random_next(const Random *random);

  /**
   * Free the memory associated to a the underlying U(0,1)-random numbers
   * generator.
   *
   * @param random is a previously initialized random generator on
   *   input, and is undefined on output, the associated memory having
   *   been freed.
   */
  void ran_random_free(Random *random);

  void ran_random_set_parameters(Random *random,
				 void *parameters, size_t size);

  /**
   * Define the parameters allowing truncation when inversion
   * technique is used.
   *
   * @param random is an initialized random numbers generator.
   * @param l is a lower bound on the untruncated cdf.
   * @param h is an upper bound on the untruncated cdf.
   * @param p is the value of P[X=a], where a is a real such that l =
   *   F(a).
   */
  void ran_random_set_truncation(Random *random, double l,
				 double h, double p);

  /**
   * Reset the current stream to its start.
   *
   * @param mrg is an initialized random numbers generator whose
   *   implementaion supports streams.
   */
  void ran_random_resetstartstream(Random *random);

  /**
   * Jump to the start of the next substream.
   *
   * @param rantom is an initialized random numbers generator whose
   *   implementation supports streams.
   */
  void ran_random_resetnextsubstream(Random *random);

  /**
   * Jump to the state to the start of the current substream.
   *
   * @param mrg is an initialized random numbers generator whose
   *   implementaion supports streams.
   */
  void ran_random_resetstartsubstream(Random *random);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_RAN_RNG_H__ */
