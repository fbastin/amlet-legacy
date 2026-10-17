/**
 * random_step.c
 *
 * OPHELIA Library - random step
 * Copyright (C) 1999-2004 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the MIT License
 */

#include <sa_step.h>

#include <bounds.h>

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
			   double max_step)
{
  int i;

  for (i = 0; i < n; i++)
    s[i*incs] += (2.0*ran_random_get_val(ran) - 1.0)*max_step;
}

void nlp_sa_step_bounds(Random *ran, double *x, int n, void *data,
			double step_size)
{
  int i;

  for (i = 0; i < n; i++)
    x[i] += (ran_random_get_val(ran)*2.0-1.0)*step_size;

  nlp_bounds_projection(n, x, NULL, (double *)data, x);
}

void nlp_sa_step_default(Random *ran, double *x, int n, void *data,
			 double step_size)
{
  int i;

  for (i = 0; i < n; i++)
    x[i] += (ran_random_get_val(ran)*2.0-1.0)*step_size;
}
