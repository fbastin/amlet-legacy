/**
 * solutions.h
 *
 * OPHELIA Library
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __OPHELIA_SOLUTIONS_H__
#define __OPHELIA_SOLUTIONS_H__

#include <stdio.h>

#include <oratio/cblas.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct _OP_Solutions OP_Solutions;

  struct _OP_Solutions {
    double *obj;
    double **sol;

    int n;    /* solutions dimensions */
    int s;    /* number of stored solutions */
    int smax; /* maximum number of solutions */
  };

  /**
   * Create a new solutions structure.
   *
   * @param smax is the maximum number of solutions that will be stored.
   * @param n is the solution dimension.
   *
   * @return a newly memory-allocated OP_Solutions structure.
   */
  OP_Solutions *op_solutions_new(int smax, int n);

  /**
   * Free the memory associated to a OP_Solutions structure.
   *
   * @param s is a previously allocated OP_Solutions structure.
   */
  void op_solutions_free(OP_Solutions *s);

  /**
   * Add a solution to the set of solutions.
   *
   * @param s is an initialized OP_Solutions structure.
   * @param obj is the value of the objective function at x.
   * @param x is the solution to add.
   */
  void op_solutions_add(OP_Solutions *s, double obj, double *x);

  /**
   * Dump recorded solutions and objective values.
   *
   * @param s is an initialized OP_Solutions structure.
   * @param f is the file where to dump the values.
   */
  void op_solutions_dump(OP_Solutions *s, FILE *f);

#ifdef __cplusplus
}
#endif

#endif  /* __OPHELIA_SOLUTIONS_H__ */
