/**
 * sa.h
 *
 * OPHELIA Library
 * Copyright (C) 1999-2004 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

/**
 * Simulating annealing method (basic implementation)
 */

#ifndef __OPHELIA_SA_H__
#define __OPHELIA_SA_H__

#include "sa_step.h"

#include <oratio/ran.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
    double t;     /* temperature */
    double mu_t;  /* temperature cooling factor */
    double k;
    double t_min; /* minimum temperature */
    int n_evals;
    void (*step)(Random *, double *, int, void *, double);
    double step_size;

    FILE *output;
    /* printing function */
    void (*printer)(FILE *output, BTR *btr, void *data);

    void *data;

    int verbose;
  } NlpSa;

  NlpSa *nlp_sa_new();
  void nlp_sa_free(NlpSa *sa);

  /**
   * Simulated annealing (basic implementation)
   */
  double nlp_simulated_annealing(NlpSa *sa, Random *ran,
				 double (*f)(double *, int, void *),
				 double *x, int n, void *data);

#ifdef __cplusplus
}
#endif

#endif  /* __OPHELIA_SA_H__ */
