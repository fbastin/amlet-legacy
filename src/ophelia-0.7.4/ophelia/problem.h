/**
 * problem.h
 *
 * OPHELIA Library
 * Copyright (C) 2005-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __OPHELIA_PROBLEM_H__
#define __OPHELIA_PROBLEM_H__

#include "common_struct.h"

#include <oratio/cblas.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct _OP_Problem OP_Problem;

  struct _OP_Problem {
    char *name;
    void *problem;
    int  type;
    void *data;

    int n;        /* number of variables */
    double *x0;   /* initial solution */
    double *x;    /* current solution */
    double *g;
    double *hess;

    int reuse_x;
    int reuse_g;
    int reuse_hess;

    char **variables;

    /**
     * Evaluate the objective function.
     * @param n is the dimension of x.
     * @param incx is the vector x increment.
     */
    double (*eval)(OP_Problem *problem,
		   int n, double *x, int incx);

    /**
     * @param n is the dimension of x and the gradient.
     * @param incx is the vector x increment.
     * @param g is a vector that contains on output the value of the gradient
     *   at x.
     * @param incg is the gradient vector increment.
     */
    void (*gradient)(OP_Problem *problem,
		     int n, double *x, int incx,
		     double *g, int incg);

    /**
     * @param n is the dimension of x.
     * @param incx is the vector x increment.
     * @param Order: CblasRowMajor or CblasColMajor
     * @param H: before entry, the leading m by n part of the array H must
     *   contain the matrix of coefficients. On exit, the leading m by n
     *   part of the array H contains the identity matrix.
     * @param ldH specifies the first dimension of H as declared in the
     *   calling (sub) program. lda must be at least max( 1, m ).
     */
    void (*hessian)(OP_Problem *problem,
		    int n, double *x, int incx,
		    const enum CBLAS_ORDER Order,
		    double *H, const int ldH);

  };

  /**
   * Create a new, uninitialized, optimization problem.
   *
   * @return a newly created optimization problem.
   */
  OP_Problem *op_problem();

  /**
   * Create a new optimization problem.
   *
   * @param name is the name of the problem.
   * @param nv is the number of variables.
   *
   * @return a newly created optimization problem.
   */
  OP_Problem *op_problem_new(char *name, int nv);

  /**
   * Free the memory associated to an optimization problem.
   *
   * @param problem is a pointer to an optimization problem.
   */
  void op_problem_free(OP_Problem *problem);

  /**
   * Get the name of some optimization problem.
   *
   * @param problem is an initialized OP_Problem instance.
   *
   * @return the problem name.
   */
  char *op_problem_get_name(OP_Problem *problem);

  /**
   * Set the problem name.
   *
   * @param name is a string containing the name of the problem; the string
   *   is duplicated in memory.
   */
  void op_problem_set_name(OP_Problem *problem, char *name);

  /**
   * Get the initial solution x0.
   *
   * @param problem is an initialized OP_Problem instance.
   *
   * @return a reference to x0. The memory should not be free independently of
   *   problem.
   */
  double *op_problem_get_x0(OP_Problem *problem);

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
  void op_problem_set_x0(OP_Problem *problem, double *x0, int incx0);

  /**
   * Set default names to problem variables.
   *
   * @param p is the problem under consideration.
   */
  void op_problem_set_varnames_default(OP_Problem *p);

  /**
   * Create variables in a given problem.
   *
   * @param p is the problem under consideration.
   * @param n is the problem dimension
   */
  void op_problem_set_variables(OP_Problem *p, int n);

  /**
   * The remaining of this file is kept for compatibility issues.
   */

  typedef double (*nlp_function)(double *, int, void *);

  typedef struct _nlp_problem nlp_problem;

  struct _nlp_problem {
    /* Number of variables */
    int n;

    /* Objective function */
    nlp_function objective;

    /* derivation functions */
    deriv1 d1;              /* first derivate */
    deriv2 d2;              /* second derivate */

    /* Inequality constraints */
    nlp_function *inequality_cstrs;
    int n_ineq;

    /* Equality constraints */
    nlp_function *equality_cstrs;
    int n_eq;

    /* Bounds */
    double *bounds;
  };

  /**
   * @param n is the problem dimension
   */
  nlp_problem *nlp_problem_new(int n);

  void nlp_problem_free(nlp_problem *pbl);

#ifdef __cplusplus
}
#endif

#endif  /* __OPHELIA_PROBLEM_H__ */
