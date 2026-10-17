/**
 * solve.c
 *
 * Discrete Choice Library
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <amlet/dcl.h>
#include "likelihood.h"

#include <stdio.h>

#include <ophelia/nlp.h>

/**
 * MEMORY GESTION
 */

void dcl_memory_allocation(int nobs, int nalt, int nvar, int censor,
			   int **choices, double ****obs,
			   int ***availability);

/**
 */

void dcl_ml_set_hessian_approx(BTR *b, int hessian_method)
{
  switch (hessian_method) {
  case HESSIAN_BHHH:
    //    b->d2approx = (C_GENERIC)bhhh_update;
    break;
  case HESSIAN_SR1:
    btr_set_hessian_approximation(b, NLP_SECANT_SR1, 0);
    break;
  case HESSIAN_BFGS:
  default:
    btr_set_hessian_approximation(b, NLP_SECANT_BFGS, 0);
    break;
  }
}

/**
 * @param n is the number of parameters to estimate.
 */
void dcl_ml_init_solver(NTLog *mclog, int n,
			int method, DclOptiMethod *opti,
			double tol, int maxiter, int hessian_method)
{
  double *x;

  /* Optimization algorithm */
  switch(method) {
  case MCBTR:
  case MCBTRC:
  case MCBTRDA:
    /* Trust-region method */
    opti->b = _malloc(BTR);

    op_btr_init(opti->b, n/*+dr.var[NNESTS].i-1*/, 0); 
    opti->b->stop = nlp_btr_stopping_test;

    x = opti->b->x;

    opti->b->tol = nt_max(tol, NLP_MINIMUM_TOL);
    opti->b->max_iter = maxiter;

    btr_set_constants(opti->b, 0.01, 0.75, 0.5, 0.5);
    opti->b->delta = 1.0;

    opti->b->d1 = NULL;
    opti->b->data = NULL;
    opti->b->output = mclog->f;

    dcl_ml_set_hessian_approx(opti->b, hessian_method);
    break;
  case MCBFGS:
  case MCBFGSDA:
    /* BFGS linesearch method */

    opti->bf = _malloc(BFGS);
    nlp_bfgs_init(opti->bf, n/*+dr.var[NNESTS].i-1*/, 0);
    opti->bf->stop = nlp_bfgs_stopping_test;

    x = opti->bf->x;

    opti->b->tol = nt_max(tol, NLP_MINIMUM_TOL);
    opti->b->max_iter = maxiter;

    opti->b->d1 = NULL;
    break;
  default:
#ifdef MPI
    MPI_Finalize();
#endif
    opti->b = NULL;
    c_error("unsupported method invoked for optimization process");
  }
}

/**
 * Very basic starting point initialization. You should consider to
 * replace this function by your own if you have a better strategy for
 * your particular problem.
 *
 * @param n is the number of parameters to estimate.
 * @param x is the parameters vectors, set to an initial value on output.
 * @param incx is the increment to use in x.
 */
void amlet_ml_starting_point(int n, double *x, int incx)
{
  int i;

  for (i = 0; i < n; i++) x[i*incx] = 0.0;
}

/**
 * Solve the maximum likelihood problem.
 *
 * @param log is a ntools log file.
 * @param n is the number of parameters to estimate.
 * @param method is the optimisation methods that has to be used.
 * @param n_random_coeff is the number of random coefficients that are
 *   present in the model.
 * @param verbose is set to 1 if the function is verbose, displaying
 *   various comments on screen, 0 else.
 */
double dcl_ml_solve(NTLog *log, int n,
		    int method, DclOptiMethod *opti,
		    int n_random_coeff,
		    double **H, void *data,
		    int verbose)
{
  double ll = 0.0;
  double *work;

  work = _mallocs(4*n, double);

  if (!n_random_coeff) {
    if (verbose) printf("** Multinomial logit ** \n");
    nt_log_section(log, "Multinomial logit");

    switch (method) {
    case MCBTR:
    case MCBTRC:
    case MCBTRDA:
#ifdef MPI
#else
      ll = nlp_btrx(opti->b,
		    (C_GENERIC)wrapper_grad, NlpFG,
		    (C_GENERIC)wrapper_grad, NlpFval,
		    NULL, 0, data, 0.0, stdout, H,
		    CblasRowMajor, CblasUpper, work);
#endif
      break;
    case MCBFGS:
    case MCBFGSDA:
      break;
    default:
      break;
    }
  }
  else {
    if (verbose) printf("** Mixed logit ** \n");
    nt_log_section(log, "Mixed logit");

    switch (method) {
    default:
      break;
    }
  }

  /* Free allocated memory. */
  free(work);

  return ll;
}
