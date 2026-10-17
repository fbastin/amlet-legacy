/**
 * bfgs_struct.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2003-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __NLP2_BFGS_STRUCT_H__
#define __NLP2_BFGS_STRUCT_H__

#include "common_struct.h"

#include <oratio/ntools.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * BFGS Structure definition
   */

  typedef struct _BFGS BFGS;

  struct _BFGS {
    /* BFGS parameters */
    double alpha, gamma;

    int n;             /* problem dimension */
    int m;             /* second dimension, set to 0 if unused */

    double *x;
    double value;  /* function value */
    double *grad;  /* function gradient evaluated at x */

    double oldgradsqrnrm, gradsqrnrm;

    int k;        /* iteration counter */
    int max_iter; /* maximal number of iterations */

    double tol;     /* tolerance */

    /* derivation functions */
    deriv1 d1;     /* first derivate */

    double (*wrapper_fct)(C_GENERIC f, int f_type,
			  double *x, int n, int m, void *data,
			  double *values, double *grad, double **hessian,
			  void *prob_info, int type_prob, int *info);
    void (*wrapper_grad)(C_GENERIC f, int f_type, C_GENERIC g, int g_type,
			 double *x, int n, void *data, double val,
			 double *grad, double *scaling,
			 void *prob_info, int prob_type, int *info,
			 double *h, double *work);

    FILE *output;
    /* printing function */
    void (*printer)(FILE *output, BFGS *bf, void *data);

    /* stopping test */
    int (*stop)(BFGS *bfgs);

    /* info on the behaviour of the BFGS linesearch algorithm */
    int info;

    /* generic data pointer */
    void *data; /* induces a problem... why??? */

    int verbose;
  };

  /**
   * BFGS method initialization.
   *
   * @param bfgs is a memory allocated BFGS structure.
   * @param n is the main problem dimension.
   * @param m is the secondary problem dimension, usually set to 0.
   */
  void nlp_bfgs_init(BFGS *bf, int n, int m);

  /**
   * Free memory associated to a BFGS structure.
   */
  void nlp_bfgs_free(BFGS *bf);

#ifdef __cplusplus
}
#endif

#endif /* __NLP2_BFGS_STRUCT_H__ */
