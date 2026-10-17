/**
 * random_step.h
 *
 * OPHELIA Library
 * Copyright (C) 1999-2015 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __OPHELIA_RANDOM_STEP_H__
#define __OPHELIA_RANDOM_STEP_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Compute a random step.
   *
   * @param ran is an initialized Random numbers generator.
   * @param s is a vector that will contain the step on output.
   * @param n is the dimension of the step to generate.
   * @param incs is the increment of s.
   * @param max_step is the maximum step size.
   */
  void op_random_step(Random *ran, int n, double *s, int incs,
		      double max_step);

#ifdef __cplusplus
}
#endif

#endif  /* __OPHELIA_RANDOM_STEP_H__ */
