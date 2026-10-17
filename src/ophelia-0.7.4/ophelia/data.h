/**
 * data.h
 *
 * Discrete Choice Library
 * Copyright (C) 2005-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __DCL_DATA_H__
#define __DCL_DATA_H__

#include <oratio/dstruct.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Define the type of data used for the samplings.
   */
#ifdef USESINGLE
  typedef float MLTYPE;
#else
  typedef double MLTYPE;
#endif

  typedef struct _nest nest;

  struct _nest {
    int n;        /* number of nests */
    int *alt;     /* alt[i] designs the nest number of alternative i*/
    double *mult; /* multipliers associated to nests */
  };

  /* quite crude approach... */
  typedef struct {
    double ***obs;       /* observations table */
    int *rep;            /* repetition table: the observations
			    corresponding to individual i are given in
			    obs[rep[i]]..obs[rep[i+1]-1] */
    int **availability;  /* availability of alernatives; first index:
			    individual, second index: alternative */
    int *choices;        /* choices[i] represent the alternative
			    choice of individual i */
    MLTYPE **sampling;   /* array of draws used in simulation */
    int size_sampling;   /* current sample size */
    int max_size;        /* maximum allowed sample size */

    double *sprob;       /* simulated probabilities for the
			    individual choices */
    double *som;         /* second-order_moments */
    double **sgrad;      /* multinomial logit: gradients of choice
			    probabibilities for each individual, divided
			    by the choice probabilities
			    mixed logit: gradients of simulated
			    probabilities of each individual */
    double **somgrad;

    double *newval;
    double **newgrad;
    unsigned int *twopowers;

    int nobs;            /* number of observations */
    int np  ;            /* number of people */
    int nalt;            /* number of alternatives */
    int nvar;            /* number of coefficients for linear utility model */
    int npar  ;          /* number of parameters to be estimated */
    int n_random_coeff;  /* number of random coefficients */

    int bias;            /* 1 if the bias has to be corrected, 0 else. */

    DS_List *np_list;    /* nonparamatric coefficient objects list */ 
    int n_nonparametric; /* number of nonparametric objects */

    nest *n;             /* nests structure */
 
    int *btype;          /* parameters types for the calibration process */

    double level;        /* level for error estimation */
    double err_quantile;

    int memcheap;        /* memory mode flag */
    int method;          /* approximation method */
    int qmcrep;          /* number of functions evalutation required for
			    estimating the accuracy */

    int safe;

    int neval;
  } dcl_data;

#ifdef __cplusplus
}
#endif

#endif  /* __DCL_DATA_H_ */
