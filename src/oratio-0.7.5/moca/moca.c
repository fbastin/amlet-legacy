/**
 * moca.c
 *
 * Monte Carlo and Quasi-Monte Carlo Library
 * Copyright (C) 2005-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/moca.h"
#include "oratio/oratio.h"

#define MAX_VALUE 1e25  /* too small values lead to bad behaviour of the BFGS,
			   a well as too large :( */

void *moca_generator_new(int method, int s)
{
  switch (method) {
  case MOCA_SOBOL:
    {
      MOCA_Sobol *sobol = moca_sobol_new(SOBOL_BRATLEY_FOX, s);
      return sobol;
    }
    break;
  default:
    return NULL;
    break;
  }
}

void moca_generator_free(int method, void *gen)
{
  switch (method) {
  case MOCA_SOBOL:
      moca_sobol_free((MOCA_Sobol *)gen);
    break;
  default:
    break;
  }
}

void moca_generator_reset(void *gen, int method)
{
  switch (method) {
  case MOCA_SOBOL:
    moca_sobol_reset(gen);
    break;
  }
}

/**
 * Count the parameters associated to a random distribution law.
 *
 * @param law is an identifier of the random distribution law.
 * @param data is a pointer to specific law information.
 *
 * @return the number of parameters associated to a random distribution law.
 */
int moca_count_parameters(int law, void *data)
{
  switch (law) {
    /* Distributions */
  case CONSTANT:
    return 1;
  case BERNOUILLI:
    return RAN_BERNOUILLI_NPAR;
  case UNIFORM:
    return 2;
  case LOGNORMAL:
    return RAN_LOGNORMAL_NPAR;
  case NORMAL:
    return RAN_NORMAL_NPAR;
  case TRIANGULAR:
    return 2;
  case TRUNCATED_NORMAL:
    return RAN_TRUNCATED_NORMAL_NPAR;
  case TRUNCATED_NORMAL_B:
    return 2;
  case ERRC:
    return 1;
    /* Nonparametric */
  case RAN_SP_SPLINE: {
    NT_BSpline *bs = (NT_BSpline *)data;
    return nt_bsplines_space_dimension(bs->p, bs->n);
  }
  default:
    c_warning("unknown number of parameters; return 1.\n");
    return 1;
  }
}

/**
 * Determine the number of parameters in a vector and the associated
 * number of random coefficients.
 *
 * @param n is the vector dimension.
 * @param pars is the parameters type vector.
 * @param n_random_coeff points to a integer representing the number
 *   of random coefficients.
 * @param idx_rnd is a vector of indexes indicating the positions of
 *   the random coefficients on output, if idx_rnd is not NULL,
 *   otherwise it is not filled.
 */
int moca_count_total_parameters(int n, int *pars, int idx,
				int *npars, int *n_random_coeff,
				int *idx_rnd)
{
  int i, t;

  *npars = 0;
  *n_random_coeff = 0;
  for (i = 0; i < n; i++) {
    t = moca_count_parameters(pars[i*idx], NULL);
    if (t && pars[i*idx] != CONSTANT) {
      if (idx_rnd) idx_rnd[*npars] = *n_random_coeff;
      (*n_random_coeff)++;
    }
    *npars += t;
  }

  return *npars;
}

/**
 * Construct an array that summarizes the basic types for random variables
 */
void moca_set_sampling_types(int *var_types, int *sampling_types,
			     int n_random_coeff)
{
  int i, j;

  i = j = 0;
  while (i < n_random_coeff) {
    if (var_types[j] != CONSTANT) sampling_types[i++] = var_types[j];
    j++;
  }
}

/**
 * Set a particular realisation inside a realizations vector.
 *
 * @param theta is the parameters vector.
 * @param type is the distribution law of nonparametric method to use.
 */
int moca_set_realisation(double *beta, double *theta, double draw,
			 int type, int pos, void *data)
{
  int offset = 1;

  switch(type) {
  case ERRC:
    *beta = theta[pos]*draw;
    break;
  case NORMAL:
  case TRIANGULAR:
  case TRUNCATED_NORMAL_B:
  case UNIFORM:
    *beta = theta[pos]+theta[pos+1]*draw;
    offset = 2;
    break;
  case LOGNORMAL:
    *beta = exp(theta[pos]+theta[pos+1]*draw);
    offset = 2;
    break;
  case RAN_SP_SPLINE: {
    NT_BSpline *bs = (NT_BSpline *)data;
    *beta = nt_bsplines_curve_point(bs->n, bs->p, draw,
				    bs->U, theta+pos, bs->N);
    offset = nt_bsplines_space_dimension(bs->p, bs->n);
    break;
  }
  case TRUNCATED_NORMAL: {
    *beta = ran_truncatednormal_from_U(draw, theta[pos], theta[pos+1],
				       theta[pos+2], theta[pos+3]);
    offset = 4;
    break;
  }
  case CONSTANT:
  default:
    *beta = theta[pos];
    break;
  }

  *beta = nt_min(*beta, MAX_VALUE);
  *beta = nt_max(*beta, -MAX_VALUE);

  return offset;
}
