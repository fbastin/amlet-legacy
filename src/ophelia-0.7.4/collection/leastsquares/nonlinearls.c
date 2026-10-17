/**
 * nonlinearls.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/cblas.h>

#include <math.h>

#include "ophelia/nlp_least_squares.h"

#include <stdio.h>

/**
 * Basic nonlinear least-squares function.
 */

double nonlinearls1(double *x, int n, void *data)
{
  double r[3];
  double y[]={2, 4, 3};
  double t[]={1, 2, 3};
  double val;
  int m = 3;
  int i;

  for (i = 0; i < m; i++) r[i] = exp(t[i]*x[0])-y[i];

  val = 0.5*cblas_ddot(m, r, 1, r, 1);

  return val;
}

double nonlinearls1_fg(double *x, int n, void *data, double *grad)
{
  double r[3];
  double y[]={2, 4, 3};
  double t[]={1, 2, 3};
  double val;
  int m = 3;
  int i;

  for (i = 0; i < m; i++) r[i] = exp(t[i]*x[0])-y[i];
  for (i = 0; i < m; i++) grad[i+1] = t[i]*exp(t[i]*x[0]);

  val = 0.5*cblas_ddot(m, r, 1, r, 1);

  nlp_ls_grad(m, n, grad, n, r, 1, 0);

  /*  nt_vect_print(stdout, "g", grad, n*(m+1)); */

  return val;
}
