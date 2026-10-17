/**
 * sa.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 1999-2004 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "sa.h"
#include "sa_step.h"

#include <math.h>
#include <string.h>

#include <oratio/ctools.h>

#define NLP_SA_T_INITIAL           1.0
#define NLP_SA_T_MIN               1.0e-5
#define NLP_SA_MU_T                1.003
#define NLP_SA_BOLTZMANN_CONSTANT  1.0

NlpSa *nlp_sa_new()
{
  NlpSa *sa;

  sa = _malloc(NlpSa);

  /* Some default parameters */
  sa->t = NLP_SA_T_INITIAL;
  sa->t_min = NLP_SA_T_MIN;
  sa->mu_t = NLP_SA_MU_T;
  sa->k = NLP_SA_BOLTZMANN_CONSTANT;
  sa->step = nlp_sa_step_default;
  sa->step_size = 1.0;

  sa->output = NULL;
  sa->printer = NULL;

  sa->data = NULL;

  sa->verbose = 1;

  return sa;
}

void nlp_sa_free(NlpSa *sa)
{
  if (sa) {
    free(sa);
  }
}


int nlp_simulating_annealing_stopping_test()
{
  return 0;
}

/**
 * Simulated annealing (basic implementation).
 */
double nlp_simulated_annealing(NlpSa *sa, Random *ran,
			       double (*f)(double *, int, void *),
			       double *x, int n, void *data)
{
  double e, new_e, best_e;  /* energy, new energy, best energy */
  double t;                 /* temperature */
  double *best_x, *new_x;
  int k;

  int done = 0;

  /* Initialisation */
  t = sa->t;
  best_e = e = f(x, n, data); /* initial energy */

  best_x = _mallocs(n, double);
  new_x = _mallocs(n, double);

  memcpy(new_x, x, n*sizeof(double));
  memcpy(best_x, x, n*sizeof(double));

  k = 0;
  sa->n_evals = 1;

  while (!done) {
    k++;

    sa->step(ran, new_x, n, data, sa->step_size);

    new_e = f(new_x, n, data);
    (sa->n_evals)++;

    if (new_e < e) {
      memcpy(x, new_x, n*sizeof(double));
      e = new_e;
      if (new_e < best_e) {
	memcpy(best_x, new_x, n*sizeof(double));
	best_e = new_e;
      }
    }
    else {
      /* accept the point with the boltzman probability */
      if (ran_random_get_val(ran) < exp((e-new_e)/(sa->k*t))) {
	e = new_e;
	memcpy(x, new_x, n*sizeof(double));
      }
    }

    /* cooling */
    t /= sa->mu_t;
    if (t < sa->t_min) done = 1;
  }

  memcpy(x, best_x, n*sizeof(double));

  free(best_x);
  free(new_x);

  return best_e;
}
