/**
 * model.h
 *
 * AMLET Library
 * Copyright (C) 2010-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __AMLET_MODEL_H__
#define __AMLET_MODEL_H__

#include <ophelia/ophelia.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
    OP_Problem *p;    /* optimization problem */

    int sampling_id;  /* identifier of the sampling method */

    double **H;       /* Hessian */
    double **Hinv;    /* Hessian inverse */
    double **I;       /* Information matrix */
    double **WCov;    /* Weak variance-Covariance matrix */
    double **RCov;    /* Robust variance-Covariance matrix */

    int verbose;      /* verbose mode */
  } AMLET_Model;

  /**
   * Create a new model.
   *
   * @return a newly memory-allocated model.
   */
  AMLET_Model *amlet_model_new();

  /**
   * Free the memory associated to the model.
   *
   * @param m is a pointer to a defined model.
   */
  void amlet_model_free(AMLET_Model *m);

  /**
   * Initialize the model;
   *
   * @param s is the problem dimension.
   */
  void amlet_model_init(AMLET_Model *m, int s);

  /**
   * Return the optimization problem associated to the model.
   *
   * @param m is an initialized model.
   */
  OP_Problem *amlet_get_problem(AMLET_Model *m);

  /**
   * Return the identifer of the sampling method.
   *
   * @param m is an initialized model.
   *
   * @return the identifier of the sampling method associated to the model.
   */
  int amlet_get_sampling_id(AMLET_Model *m);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_MODEL_H__ */
