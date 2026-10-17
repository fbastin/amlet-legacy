/**
 * problem.c
 *
 * OPHELIA Library
 * Copyright (C) 2005-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "problem.h"

#include <string.h>

#include <oratio/ctools.h>

/**
 * Create a new, uninitialized, optimization problem.
 *
 * @return a newly created optimization problem.
 */
OP_Problem *op_problem()
{
  OP_Problem *problem = _malloc(OP_Problem);

  problem ->name = NULL;
  problem->variables = NULL;
  problem->data = problem->g = problem->hess = NULL;
  problem->x0 = problem->x = NULL;

  problem->n = 0;

  return problem;
}

/**
 * Set default names to problem variables.
 *
 * @param p is the problem under consideration.
 */
void op_problem_set_varnames_default(OP_Problem *p)
{
  int i;
  int b, len;

  b = 10;
  len = 3;

  if (p->variables)
    for (i = 0; i < p->n; i++) free(p->variables[i]);

  p->variables = _reallocs(p->variables, p->n, char *);
  for (i = 1; i <= p->n; i++) {
    if (i % b == 0) {
      len++;
      b *= 10;
    }
    p->variables[i-1] = _mallocs(len, char);
    sprintf(p->variables[i-1], "x%d", i);
  }
}

/**
 * Create variables in a given problem.
 *
 * @param p is the problem under consideration.
 * @param n is the problem dimension
 */
void op_problem_set_variables(OP_Problem *p, int n)
{
  p->n = n;
  p->x = _reallocs(p->variables, n, double);
  p->x0 = _reallocs(p->variables, n, double);
  op_problem_set_varnames_default(p);
}

/**
 * Create a new optimization problem.
 *
 * @param name is the name of the problem.
 * @param n is the number of variables.
 *
 * @return a newly created optimization problem.
 */
OP_Problem *op_problem_new(char *name, int n)
{
  OP_Problem *problem = _malloc(OP_Problem);

  problem->name = (name) ? strdup(name) : strdup("Untitled");

  problem->n = n;

  problem->x = problem->x0 = NULL;
  problem->variables = NULL;
  op_problem_set_variables(problem, n);

  return problem;
}

/**
 * Free the memory associated to an optimization problem.
 *
 * @param problem is a pointer to an optimization problem.
 */
void op_problem_free(OP_Problem *problem)
{
  if (problem) {
    int i;

    free(problem->x0);
    free(problem->x);

    for (i = 0; i < problem->n; i++) free(problem->variables[i]);
    free(problem->variables);
    if (problem->name) free(problem->name);
    free(problem);
  }
}

/**
 * Get the name of some optimization problem.
 *
 * @param problem is an initialized OP_Problem instance.
 *
 * @return the problem name.
 */
char *op_problem_get_name(OP_Problem *problem)
{
  return problem->name;
}

/**
 * Set the problem name.
 *
 * @param name is a string containing the name of the problem; the string
 *   is duplicated in memory.
 */
void op_problem_set_name(OP_Problem *problem, char *name)
{
  if (name) problem->name = strdup(name);
}

/**
 * Get the initial solution x0.
 *
 * @param problem is an initialized OP_Problem instance.
 *
 * @return a reference to x0. The memory should not be free independently of
 *   problem.
 */
double *op_problem_get_x0(OP_Problem *problem)
{
  return problem->x0;
}

/**
 * Set the initial solution x0.
 *
 * @param problem is an initialized OP_Problem instance.
 * @param x0 is a pointer to the intial solution; the dimension
 *   should be the same then the one defined in problem.
 * @param incx0 is the increment to use with x0 (usually 1).
 *
 * @return x0
 */
void op_problem_set_x0(OP_Problem *problem, double *x0, int incx0)
{
  if (incx0 == 0)
    memcpy(problem->x0, x0, problem->n*sizeof(double));
  else {
    int i;

    for (i = 0; i < problem->n; i++) problem->x0[i] = x0[i*incx0];
  }
}

/**
 * @param n is the problem dimension.
 */
nlp_problem *nlp_problem_new(int n)
{
  nlp_problem *pbl;

  pbl = _malloc(nlp_problem);

  pbl->n = n;

  pbl->objective = NULL;

  pbl->d1 = NULL;
  pbl->d2 = NULL;

  pbl->inequality_cstrs = NULL;
  pbl->n_ineq = 0;

  pbl->equality_cstrs = NULL;
  pbl->n_eq = 0;

  pbl->bounds = NULL;

  return pbl;
}

void nlp_problem_free(nlp_problem *pbl)
{
  if (pbl) free(pbl);
}
