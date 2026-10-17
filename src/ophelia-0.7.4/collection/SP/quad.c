/**
 * quad.c
 *
 * OPHELIA Library
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/matrix.h>
#include <oratio/moca.h>
#include <ophelia/problem.h>

double quad1_eval(OP_Problem *p, int n, double *x, int incx)
{
  double val, t, mean;
  int i, r;

  OR_Matrix *S = (OR_Matrix *)(p->data);

  mean = 0.0;
  for (r = 0; r < S->m; r++) {
    val = 0.0;
    for (i = 1; i < 21; i++) {
      t = i*(x[(i-1)*incx]-(21-i)*S->A[r][i]);
      val += t*t;
    }
    mean += (val-mean)/(r+1);
  }

  return mean;
}

void quad1_gradient(OP_Problem *p, int n, double *x, int incx,
		    double *g, int incg)
{
  int i, r;
  OR_Matrix *S = (OR_Matrix *)(p->data);

  for (i = 0; i < n; i++) g[i] = 0.0;

  for (r = 0; r < S->m; r++) {
    for (i = 1; i < 21; i++) {
      g[i*incg] += (2*i*(x[(i-1)*incx]-(21-i)*S->A[r][i]) -
		    g[i*incg])/(r+1);
    }
  }
}

void op_probl_sp_quad1_free(OP_Problem *p)
{
  if (p) {
    or_matrix_free((OR_Matrix *)(p->data));
    op_problem_free(p);
  }
}

OP_Problem *op_probl_sp_quad1_new(Random *r, int n)
{
  OP_Problem *p;
  OR_Matrix *S;

  p = op_problem_new("QUAD1_SP", 20);
  S = or_matrix_new(n, 20);
  moca_mc_generate_sample(r, n, 20, S->A);

  p->data = S;

  p->eval = &quad1_eval;
  p->gradient = &quad1_gradient;

  return p;
}
