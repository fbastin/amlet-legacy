/**
 * madsen.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "ophelia/nlp_least_squares.h"

#include <math.h>
#include <stdio.h>

#include <oratio/cblas.h>

double madsen(double *x, int n, void *data)
{
  double r[3];
  double val;

  r[0] = x[0]*x[0]+x[1]*x[1]+x[0]*x[1];
  r[1] = sin(x[0]);
  r[2] = cos(x[1]);

  val = 0.5*cblas_ddot(3, r, 1, r, 1);

  return val;
}

double madsen_fg(double *x, int n, void *data, double *grad)
{
  double r[3];
  double val;
  int m = 3;

  int shift = 0;

  r[0] = x[0]*x[0]+x[1]*x[1]+x[0]*x[1];
  r[1] = sin(x[0]);
  r[2] = cos(x[1]);

  val = 0.5*cblas_ddot(m, r, 1, r, 1);

  grad[0] = 0.0;
  grad[1] = 0.0;

  shift = n;

  /* First element */
  grad[shift++] = 2*x[0]+x[1];
  grad[shift++] = 2*x[1]+x[0];

  /* Second element */
  grad[shift++] = cos(x[0]);
  grad[shift++] = 0.0;

  /* Third element */
  grad[shift++] = 0.0;
  grad[shift++] = -sin(x[1]);

  nlp_ls_grad(m, n, grad, n, r, 1, 0);

  grad[0] = r[0]*(2*x[0]+x[1])+r[1]*cos(x[0]);
  grad[1] = r[0]*(2*x[1]+x[0])+r[2]*(-sin(x[1]));

  return val;
}

double madsen_fg_values(double *x, int n, void *data, double *values, int m,
			double *grad)
{
  double val;

  int shift = 0;

  values[0] = x[0]*x[0]+x[1]*x[1]+x[0]*x[1];
  values[1] = sin(x[0]);
  values[2] = cos(x[1]);

  val = 0.5*cblas_ddot(m, values, 1, values, 1);

  grad[0] = 0.0;
  grad[1] = 0.0;

  shift = n;

  /* First element */
  grad[shift++] = 2*x[0]+x[1];
  grad[shift++] = 2*x[1]+x[0];

  /* Second element */
  grad[shift++] = cos(x[0]);
  grad[shift++] = 0.0;

  /* Third element */
  grad[shift++] = 0.0;
  grad[shift++] = -sin(x[1]);

  nlp_ls_grad(m, n, grad, n, values, 1, 0);

  grad[0] = values[0]*(2*x[0]+x[1])+values[1]*cos(x[0]);
  grad[1] = values[0]*(2*x[1]+x[0])-values[2]*sin(x[1]);

  return val;
}
