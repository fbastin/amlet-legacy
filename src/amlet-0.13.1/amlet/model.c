/**
 * amlet.c
 *
 * AMLET (Another Mixed Logit Estimation Tool) Library
 * Copyright (C) 2002-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "amlet.h"

#include <oratio/ctools.h>
#include <ophelia/ophelia.h>

/**
 * Global variables declared in amlet_config.h and log.h.
 */
int loose = 0;
NTLog *mclog = NULL;
int nfctev = 0;
dcl_data data;
unsigned int twopowers[31];

/**
 * Create a new model.
 *
 * @return a newly memory-allocated model.
 */
AMLET_Model *amlet_model_new()
{
  AMLET_Model *m = _malloc(AMLET_Model);

  m->p = op_problem();
  m->sampling_id = MOCA_MONTE_CARLO;

  m->H = NULL;
  m->Hinv = NULL;
  m->I = NULL;
  m->WCov = NULL;
  m->RCov = NULL;

  m->verbose = 0;

  return m;
}

void amlet_model_free(AMLET_Model *m)
{
  if (m) {
    op_problem_free(m->p);
    if (m->H) nt_matrix_free(m->H);
    if (m->Hinv) nt_matrix_free(m->Hinv);
    if (m->I) nt_matrix_free(m->I);
    if (m->WCov) nt_matrix_free(m->WCov);
    if (m->RCov) nt_matrix_free(m->RCov);
    free(m);
  }
}

/**
 * Return the optimization problem associated to the model.
 *
 * @param m is a pointer to a defined model.
 */
OP_Problem *amlet_get_problem(AMLET_Model *m)
{
  return m->p;
}

/**
 * Initialize the model;
 *
 * @param s is the problem dimension.
 */
void amlet_model_init(AMLET_Model *m, int s)
{
  m->H = nt_matrix_new(s, s);
  m->Hinv = nt_matrix_new(s, s);
  m->I = nt_matrix_new(s, s);
  m->WCov = nt_matrix_new(s, s);
  m->RCov = nt_matrix_new(s, s);
  nt_matrix_identity(s, s, *(m->I), s);
  op_problem_set_variables(amlet_get_problem(m), s);
}

/**
 * Return the identifer of the sampling method.
 *
 * @param m is an initialized model.
 *
 * @return the identifier of the sampling method associated to the model.
 */
int amlet_get_sampling_id(AMLET_Model *m)
{
  return m->sampling_id;
}
