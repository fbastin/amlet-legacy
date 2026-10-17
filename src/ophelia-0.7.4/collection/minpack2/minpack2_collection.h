/**
 * minpack2_collection.h
 *
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * This file is licensed under BSD License.
 */

/**
 * MINPACK-2 test problems collection..
 */

#include <oratio/ctools.h>

#ifndef __NLP_MINPACK2_COLLECTION__
#define __NLP_MINPACK2_COLLECTION__

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
    char *name;   /* problem name */
    int id;       /* problem identifier */
    void *data;
    int n;        /* problem dimension */
    int nx, ny;  
    double par;
    double *w;    /* working vector of size 2*n */
  } Minpack2Problem;

  Minpack2Problem *minpack2_set_problem(const char *name, int n, int nx, int ny);

  int minpack2_starting_point(Minpack2Problem *mp, double *x);

  void minpack2_problem_free(Minpack2Problem *mp);

  double minpack2_wrapper(int n, int nx, int ny, double *x,
			  double *f, double *g, int itask,
			  int prob, double par, double *w);

  double minpack2_wrapper_f(double *x, int n, void *problem);

  void minpack2_wrapper_g(double *x, int n, void *problem,
			  double *grad);

  double minpack2_wrapper_f_values(double *x, int n, void *problem,
				   double *values, int m);

  double minpack2_wrapper_fg_values(double *x, int n, void *problem,
				    double *values, int m, double *grad);

  enum {
    DAER = 1,
    DEPT,
    DGL1,
    DGL2,
    DLJ2,
    DLJ3,
    DMSA,
    DMSO,    
    DODC,
    DPJB,
    DSSC
  };

  static const c_symbol minpack2_collection[] = {
    {"DAER", DAER},
    {"DEPT", DEPT},
    {"DGL1", DGL1},
    {"DGL2", DGL2},
    {"DLJ2", DLJ2},
    {"DLJ3", DLJ3},
    {"DMSA", DMSA},
    {"DMSO", DMSO},
    {"DODC", DODC},
    {"DPJB", DPJB},
    {"DSSC", DSSC}
  };

#ifdef __cplusplus
}
#endif

#endif  /* __NLP_MINPACK2_COLLECTION__ */

