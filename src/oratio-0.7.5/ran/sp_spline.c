/**
 * sp_splines.c
 *
 * Random Numbers Library
 * Copyright (C) 2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

/**
 * Spline distribution.
 */

#include "oratio/ran.h"
#include "oratio/memory.h"

void *ran_sp_spline_struct_new(int nk, double *knots, double *P)
{
  ran_sp_spline_params *sp;

  sp = _malloc(ran_sp_spline_params);

  sp->P = P;

  sp->bs = nt_bsplines_new(3, nk, knots, 1);

  return sp;
}

void ran_sp_spline_struct_free(Random *ran)
{
  ran_sp_spline_params *sp;

  sp = (ran_sp_spline_params *)(ran->parameters);

  nt_bsplines_free(sp->bs);
}

/**
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 * @param n is the dimension of seed.
 * @param nk is the number of points present in the knots vector
 *   (without repetitions).
 *
 * @return a random numbers generator for a semi-parametric variate X
 *   based on a cubic B-spline defined on the knots vector knots and
 *   parametrized by the vector P.
 */
Random *ran_sp_spline_new(unsigned long *seed, int n,
			  int nk, double *knots, double *P)
{
  Random *ran = ran_random_new_with_size(seed, n, 0);

  ran->type = RAN_SP_SPLINE;

  ran->parameters = (void *)ran_sp_spline_struct_new(nk, knots, P);

  return ran;
}

ran_sp_spline_params *ran_sp_spline_get_parameters(const Random *ran)
{
  return (ran_sp_spline_params *)(ran->parameters);
}

void ran_sp_spline_set_parameters(Random *ran, ran_sp_spline_params *sp)
{
  ran_sp_spline_struct_free(ran);
  free(ran->parameters);
  ran->parameters = sp;
}

double ran_sp_spline_U(const Random *ran, double u)
{
  NT_BSpline *bs = ((ran_sp_spline_params *)(ran->parameters))->bs;
  double *P = ((ran_sp_spline_params *)(ran->parameters))->P;

  return nt_bsplines_curve_point(bs->n, bs->p, u, bs->U, P, bs->N);
}

double ran_sp_spline_get_val(const Random *ran)
{
  return ran_sp_spline_U(ran, ran_random_get_val(ran));
}

