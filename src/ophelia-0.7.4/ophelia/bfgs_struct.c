/**
 * bfgs_struct.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006 Fabian Bastin <bastin@cerfacs.fr>
 *
 * Licensed under the MIT License
 */


#include "bfgs_struct.h"

#include "defaults.h"
#include "stopping.h"
#include "wrapper.h"

/**
 * @param bfgs is a memory allocated BFGS structure.
 */
void nlp_bfgs_wrappers_default(BFGS *bf)
{
  bf->wrapper_fct = nlp_wrapper_fct_dflt;
  bf->wrapper_grad = nlp_wrapper_grad_dflt;
}

/**
 * BFGS method initialization.
 *
 * @param bfgs is a memory allocated BFGS structure.
 * @param n is the main problem dimension.
 * @param m is the secondary problem dimension, usually set to 0.
 */
void nlp_bfgs_init(BFGS *bf, int n, int m)
{
  bf->n = n;
  bf->m = m;

  /* Iterate */
  bf->x = _mallocs(n, double);
  bf->grad = _mallocs(n*(m+1), double);

  /* Values */
  bf->oldgradsqrnrm = bf->gradsqrnrm = DBL_MAX;
  bf->value = DBL_MAX;

  bf->d1 = NULL;
  nlp_bfgs_wrappers_default(bf);

  bf->output = stdout;
  bf->printer = NULL;

  bf->stop = nlp_bfgs_stopping_test;

  bf->info = 0;

  bf->tol = NLP_TOL;
  bf->max_iter = NLP_MAX_ITER;

  bf->verbose = 1;

  bf->data = NULL;
}

/**
 * Free memory associated to a BFGS structure.
 */
void nlp_bfgs_free(BFGS *bf)
{
  free(bf->x);
  free(bf->grad);
}
