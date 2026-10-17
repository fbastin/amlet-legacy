/**
 * individual.h
 *
 * AMLET Library
 * Copyright (C) 2004-2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the BSD License.
 **/

#include <oratio/stat.h>

#include "alternative.h"

#ifndef __AMLET_INDIVIDUAL_H__
#define __AMLET_INDIVIDUAL_H__

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct _amlet_ind AMLET_Ind;
  typedef struct _amlet_obs AMLET_Obs;
  typedef struct _amlet_pop AMLET_Pop;

  /* Observations */
  struct _amlet_obs {
    int id;            /* observal identifier */
    ST_Obs *obs;       /* statistical observation */
    int nalt;          /* number of alternatives available for this
			 observation */
    AMLET_Alt *choice; /* expressed choice */
    AMLET_Alt **alt;   /* pointer to alternatives */
  };

  /* Individual */
  struct _amlet_ind {
    int id;      /* individual identifier */
    int choice;  /* choice expressed by the individual */
    double *prob_choice; /* vector of choice probabilities */
    double *grad_prob_choice; /* gradient of the choice probability, at the
				 expressed choice */
    double *properties;
    AMLET_Obs *obs;
    int nprop;
    //    AMLET_Alternative *alternatives; /* alternatives available for the individual */
    int nalt;          /* number of alternatives available to the individual */
  };

  /* Population */
  struct _amlet_pop {
    int n;
    AMLET_Ind *amlet_ind; /* array of individuals */
  };

  AMLET_Ind *amlet_ind_new();

  void AMLET_Ind_free(AMLET_Ind *i);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_INDIVIDUAL_H__ */
