/**
 * linearls.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/cblas.h>

#include "ophelia/nlp_least_squares.h"

/**
 * Extremely basic linear least-squares function.
 */

double linearls1(double *x, int n, void *data)
{
  double r[3];
  double val;
  int m=3;

  r[0] = x[0]+x[1]-2.0;
  r[1] = x[0]-5.0;
  r[2] = x[1];

  val = 0.5*cblas_ddot(m, r, 1, r, 1);

  return val;
}

double linearls1_fval(double *x, int n, void *data, double *values, int m)
{
  double val;

  values[0] = x[0]+x[1]-2.0;
  values[1] = x[0]-5.0;
  values[2] = x[1];

  val = 0.5*cblas_ddot(m, values, 1, values, 1);

  return val;
}

double linearls1_fg(double *x, int n, void *data, double *grad)
{
  double r[3];
  double val;
  int m = 3;

  int shift = 0;

  r[0] = x[0]+x[1]-2.0;
  r[1] = x[0]-5.0;
  r[2] = x[1];

  val = 0.5*cblas_ddot(m, r, 1, r, 1);

  shift = n;

  /* First element */
  grad[shift++] = 1.0;
  grad[shift++] = 1.0;

  /* Second element */
  grad[shift++] = 1.0;
  grad[shift++] = 0.0;

  /* Third element */
  grad[shift++] = 0.0;
  grad[shift++] = 1.0;

  nlp_ls_grad(m, n, grad, n, r, 1, 0);

  return val;
}

double linearls1_fgval(double *x, int n, void *data,
		       double *values, int m, double *grad)
{
  double val;

  int shift = 0;

  values[0] = x[0]+x[1];
  values[1] = x[0];
  values[2] = x[1];

  val = 0.5*cblas_ddot(m, values, 1, values, 1);

  shift = n;

  /* First element */
  grad[shift++] = 1.0;
  grad[shift++] = 1.0;

  /* Second element */
  grad[shift++] = 1.0;
  grad[shift++] = 0.0;

  /* Third element */
  grad[shift++] = 0.0;
  grad[shift++] = 1.0;

  nlp_ls_grad(m, n, grad, n, values, 1, 0);

  return val;
}
