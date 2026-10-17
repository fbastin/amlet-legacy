/**
 * siman.c
 *
 * OPHELIA Library 
 * Copyright (C) 2004 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "test_siman.h"

#include "bounds.h"

#include <oratio/ran.h>

#include <math.h>
#include <stdio.h>

#include <ophelia/nlp.h>
#include <ophelia/collection.h>

/* Set up parameters for this simulated annealing run */
#define N_TRIES 2000
#define ITERS_FIXED_T 10000
#define STEP_SIZE 1.0
#define K 1.0
#define T_INITIAL 0.008
#define MU_T 1.003
#define T_MIN 2.0e-6

int siman_test()
{
  NlpSa *sa;
  Random *ran;
  unsigned long seed;
  double x;
  double y[2];
  double bounds[4];

  sa = nlp_sa_new();

  seed = 0;
  ran = ran_random_new(0, &seed, 1);

  /* One dimensional test */
  x = 0.4;
  printf("Value: %g\n", nlp_simulated_annealing(sa, ran, gsl_e1, &x, 1, NULL));
  printf("Solution: %g\n", x);
  nlp_sa_free(sa);

  /* Two dimensionals test */
  /*  sa = nlp_sa_new();
  y[0] = y[1] = 5.0;
  sa->step = nlp_sa_step_bounds;
  nlp_bounds_reset(bounds, 2);
  pbr_t1_3_set_bounds(bounds);
  printf("Value: %g\n",
	 nlp_simulated_annealing(sa, ran, pbr_t1_3, y, 2, bounds));
	 printf("Solution: %g %g\n", y[0], y[1]);

  nlp_sa_free(sa);
  */

  ran_random_free(ran);

  return 0;
}
