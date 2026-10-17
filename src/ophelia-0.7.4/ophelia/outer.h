/**
 * outer.c
 *
 * ATHOS Library
 * Copyright (C) 2005 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ATHOS_OUTER_H__
#define __ATHOS_OUTER_H__

#include "btr.h"
#include "problem.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct _nlp_outer nlp_outer;

  typedef double (*nlp_penalty)(double, double *, void *, int);
  typedef void   (*nlp_penalty_d1)(double, double *, void *, int, double *);
  typedef void   (*nlp_penalty_d2)(double, double *, void *, int, double **);

  struct _nlp_outer {
    int k;       /* outer iterations index */
    int t;       /* total number of inner iterations */

    double mu;   /* barrier/penalty parameter */

    nlp_problem *pbl;

    BTR *b;

    /* Penalty function */
    nlp_penalty penalty;
    nlp_penalty_d1 penalty_d1;
    nlp_penalty_d2 penalty_d2;

    void *data;

    FILE *output;
    /* printing function */
    void (*printer)(FILE *output, nlp_outer *btr, void *data);

    int verbose;
  };

  nlp_outer *nlp_outer_new(nlp_problem *pbl, BTR *b);

  void nlp_outer_free(nlp_outer *no);

  void nlp_outer_set_verbose(nlp_outer *no);

  void nlp_outer_init_barrier_parameter(nlp_outer *no, double mu);

  void nlp_outer_iteration(nlp_outer *no, double **H, double *work);

#ifdef __cplusplus
}
#endif

#endif  /* __ATHOS_OUTER_H__ */
