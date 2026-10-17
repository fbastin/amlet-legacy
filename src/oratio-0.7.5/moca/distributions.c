/**
 * distributions.c
 *
 * ORATIO Library
 * Copyright (C) 2006-2007 Fabian Bastin <bastin@iro.umontreal.car>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/ctools.h"
#include "oratio/moca.h"
#include "oratio/stat.h"

#include <stdio.h>

#define TRUNCNORM_LB -10 /* lower bound for the truncated normal of type B */
#define TRUNCNORM_UB  10 /* upper bound for the truncated normal of type B */

/**
 * @param v is the value to convert.
 */
double moca_convert_random_value(int type, double v)
{
  switch (type) {
  case UNIFORM:
    return 2*v-1;
    break;
  case TRIANGULAR:
    return st_triangular_icdf(v, -1.0, 1.0, 0.0);
    break;
  case NORMAL:
  case LOGNORMAL:
  case ERRC:
    return st_gaussian_icdf(v);
    break;
  case CONSTANT:
  default:
    return v;
    break;
  }
}

/**
 * @param law is the probability law identifier.
 * @param r is an initialized random numbers generator.
 * @param data is a pointer to a structure data corresponding to law,
 *   if necessary, otherwise it is set to NULL.
 */
double moca_get_random_value(int law, Random *r, double *param,
			     void *data)
{
  switch (law) {
  case CONSTANT:
    return param[0];
  case BERNOUILLI:
    return ran_bernouilli(r, param[0]);
  case NORMAL:
    return ran_normal_icdf(r, param[0], param[1]);
  case LOGNORMAL:
    return ran_lognormal_icdf(r, param[0], param[1]);
  case TRIANGULAR:
    return ran_triangular(r, param[0]-param[1], param[0]+param[1], param[0]);
  case TRUNCATED_NORMAL:
    return ran_truncatednormal_icdf(r, param[0], param[1], param[2], param[3]);
  case TRUNCATED_NORMAL_B:
    return ran_truncatednormal_b_icdf(r, TRUNCNORM_LB, TRUNCNORM_UB, 0, 1); 
  case UNIFORM:
    return param[0]+param[1]*(ran_random_get_val(r)*2.0-1.0);
  case ERRC:
    return ran_normal_icdf(r, 0.0, param[0]);
  case RAN_SP_SPLINE: {
    NT_BSpline *bs = (NT_BSpline *)data;
    return nt_bsplines_curve_point(bs->n, bs->p, ran_random_get_val(r),
				   bs->U, param, bs->N);
  }
  default:
    c_warning("unsupported distribution %d %d; using U[0,1]", law, BERNOUILLI);
    return ran_random_get_val(r);
  }
}

/**
 * @param law is the probability law identifier.
 * @param data is a pointer to a structure data corresponding to law,
 *   if necessary, otherwise it is set to NULL.
 */
double moca_affect_value(int law, double u, double *param,
			 void *data)
{
  switch (law) {
  case CONSTANT:
    return u;
  case BERNOUILLI:
    return ran_bernouilli_from_U(u, param[0]);
  case NORMAL:
    return st_normal_icdf(u, param[0], param[1]);
  case LOGNORMAL:
    return st_lognormal_icdf(u, param[0], param[1]);
  case TRIANGULAR:
    return st_triangular_icdf(u, param[0]-param[1],
			      param[0]+param[1], param[0]);
  case TRUNCATED_NORMAL:
    return ran_truncatednormal_from_U(u, param[0], param[1],
				      param[2], param[3]);
  case TRUNCATED_NORMAL_B:
    return ran_truncatednormal_b_from_U(u, TRUNCNORM_LB, TRUNCNORM_UB,
					param[0], param[1]);
  case UNIFORM:
    /* We produce a value uniformly distributed on [-1,1] */
    return param[0]+param[1]*(u*2.0-1.0);
  case ERRC:
    return st_normal_icdf(u, 0.0, param[0]);
  case RAN_SP_SPLINE: {
    NT_BSpline *bs = (NT_BSpline *)data;

    return nt_bsplines_curve_point(bs->n, bs->p, u, bs->U, param, bs->N);
  }
  default:
    c_warning("unsupported distribution %d %d; using U[0,1]", law, BERNOUILLI);
    return 1.0;
  }
}
