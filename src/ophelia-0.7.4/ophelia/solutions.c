/**
 * solutions.c
 *
 * OPHELIA Library
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "solutions.h"

#include <string.h>

#include <oratio/ctools.h>
#include <oratio/ntools.h>

/**
 * Create a new solutions structure.
 *
 * @param smax is the maximum number of solutions that will be stored.
 * @param n is the solution dimension.
 *
 * @return a newly memory-allocated OP_Solutions structure.
 */
OP_Solutions *op_solutions_new(int smax, int n)
{
  OP_Solutions *ops = NULL;

  ops = _malloc(OP_Solutions);

  ops->n = n;
  ops->smax = smax;
  ops->s = 0;

  ops->obj = _mallocs(smax, double);
  ops->sol = nt_matrix_new(smax, n);

  return ops;
}

/**
 * Free the memory associated to a OP_Solutions structure.
 *
 * @param s is a previously allocated OP_Solutions structure.
 */
void op_solutions_free(OP_Solutions *s)
{
  if (s) {
    free(s->obj);
    nt_matrix_free(s->sol);
  }
}

/**
 * Add a solution to the set of solutions.
 *
 * @param s is an initialized OP_Solutions structure.
 * @param obj is the value of the objective function at x.
 * @param x is the solution to add.
 */
void op_solutions_add(OP_Solutions *s, double obj, double *x)
{
  while (s->s >= s->smax) {
    s->smax *= 2;
    s->obj = _reallocs(s->obj, s->smax, double);
    nt_matrix_realloc(s->sol, s->smax, s->n);
  }

  s->obj[s->s] = obj;
  memcpy(s->sol[s->s], x, s->n*sizeof(double));
  s->s++;
}

/**
 * Dump recorded solutions and objective values.
 *
 * @param s is an initialized OP_Solutions structure.
 * @param f is the file where to dump the values.
 */
void op_solutions_dump(OP_Solutions *s, FILE *f)
{
  int i, j;

  fprintf(f, "%d %d\n", s->s, s->n);
  for (i = 0; i < s->s; i++) {
    fprintf(f, "%10.8g", s->obj[i]);
    for (j = 0; j < s->n; j++)  fprintf(f, " %10.8g", s->sol[i][j]);
    fprintf(f, "\n");
  }
}
