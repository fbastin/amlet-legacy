/**
 * amlet_analysis.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the MIT License
 */

#include "logit.h"
#include "amlet_analysis.h"
#include "amlet_btr.h"
#include "amlet_config.h"

#include "bhhh.h"
#include "log.h"

#include <oratio/oratio.h>
#include <ophelia/nlp.h>
#include <amlet/dcl.h>

#include <math.h>
#include <string.h>

void constant_md(dcl_data *dest, dcl_data *src)
{
  int i;
  dest->obs = src->obs;
  dest->rep = _mallocs(src->nobs, int);
  for (i = 0; i < src->nobs; i++) dest->rep[i] = i;
  dest->availability = src->availability;

  dest->choices = src->choices;
  dest->sampling = NULL;
  dest->size_sampling = 0;
  dest->max_size = 0;

  dest->sprob = _mallocs(src->nobs, double);
  dest->som = _mallocs(src->nobs, double);
  dest->sgrad = nt_matrix_new(src->nobs, src->nvar);

  dest->nobs = src->nobs;
  dest->np = src->nobs;
  dest->nalt = src->nalt;
  dest->nvar = src->nvar;
  dest->npar = src->nvar;

  dest->np_list = src->np_list;

  dest->n = src->n;
  dest->err_quantile = src->err_quantile;
}

double likelihood_constant(AMLET_Driver *dr, double *work, dcl_data *md)
{
  int *choice;
  double llc = 0.0;
  int i, j, k;

  void *obj;

  DS_List *tmp = md->np_list;
  choice = _mallocs(md->nalt, int);

  if (!dr->var[CENSOR].i) {
    memset(choice, 0, md->nalt*sizeof(int));

    for (i = 0; i < md->nobs; i++)
      choice[md->choices[i]]++;

    for (i = 0; i < md->nalt; i++)
      llc += choice[i]*log((double)choice[i]/md->nobs);

    llc /= -md->np;
  }
  else {
    int is_constant, n;
    int **constant_obs;  /* double entry array indicating if all
			    observations for some alternative are set to
			    1.0., i.e. this array identifies constants
			    in utility functions. */
    double **H;
    BTR *b;
    int *btype;
    dcl_data *md2;

    int nfctev_old = nfctev;

    btype = _mallocs(md->nvar, int);

    /* We must calibrate a multinomial logit function */

    md2 = _malloc(dcl_data);
    /*   dcl_data_clone_btrda(md2, md);
    md2->np = md2->nobs;
    md2->rep = (int *)realloc(md2->rep, md2->nobs*sizeof(int));
    for (i = 0; i < md2->nobs; i++) md2->rep[i] = i; */

    constant_md(md2, md);

    /* We use the BTR algorithm */
    b = _malloc(BTR);
    op_btr_init(b, md2->nvar, 0);
    b->d2approx = (C_GENERIC)bfgs_update;

    H = nt_matrix_new(md2->nvar, md2->nvar);
    nt_matrix_identity(md2->nvar, md2->nvar, *H, md2->nvar);

    /* Determine the alternative specific constants */
    constant_obs = _mallocs(md->nalt, int *);
    for (i = 0; i < md->nalt; i++) constant_obs[i] = _mallocs(md->nvar, int);
    /* Search for independant terms */
    for (i = 0; i < md->nalt; i++)
      for (j = 0; j < md->nvar; j++) constant_obs[i][j] = 1;
    for (i = 0; i < md->nobs; i++) {
      for (j = 0; j < md->nalt; j++) {
	k = 0;
	while (k < md->nvar) {
	  /* We assume that alternative specific constants are
	     characterized by the fact that all observations equal to
	     1.0, and that only such constants has this property. */
	  if (md->obs[i][j][k] != 1.0 || md->btype[k] == ERRC)
	    constant_obs[j][k] = 0.0;
	  k++;
	}
      }
    }

    k = 0;
    for (i = 0; i < 4*md->nvar; i++) work[i] = 0.0;
    for (i = 0; i < md->nvar; i++) {
      /* This code could certainly be refined */
      is_constant = 0;
      for (j = 0; j < md->nalt; j++) {
	if (constant_obs[j][i]) {
	  is_constant = 1;
	  break;
	}
      }

      obj = NULL;

      if (amlet_is_nonparametric(md->btype[i])) {
	obj = tmp->data;
	tmp = tmp->next;
      }

      n = moca_count_parameters(md->btype[i], obj);
      btype[i] = is_constant? CONSTANT : FIXED;
      k += n;
    }

    md2->btype = btype;

    llc = nlp_btrx(b,
		   (C_GENERIC)wrapper, NlpFG,
		   (C_GENERIC)wrapper_grad, NlpFval,
		   NULL, 0,
		   md2, 0.0, NULL, H,
		   CblasRowMajor, CblasUpper, work);
    llc *= -md->nobs;
    llc /= -md->np;

    for (i = 0; i < md->nalt; i++) free(constant_obs[i]);
    free(constant_obs);
    nt_matrix_free(H);

    free(md2->sprob);
    free(md2->som);
    nt_matrix_free(md2->sgrad);
    free(md2);

    op_btr_free(b);
    free(btype);

    nfctev = nfctev_old;
  }

  free(choice);
 
  return llc;
}

/**
 * Compute the log-likelihood value when all parameters are set to
 * zero and when only alternative specific constants are considered.
 *
 * @param ismixed is an integer equal to 0 if we face a classical
 *    multinomial model, another value elsewhere.
 */
void special_likelihood_functions(int ismixed, double *ll0, double *llc,
				  dcl_data *md, AMLET_Driver *dr,
				  double *work)
{
  /* Compute the likelihood function at 0 */
  *ll0 = -amlet_mnl_null(md->nalt, md->nobs, md->availability);

  //  *llc = likelihood_constant(dr, work, md);
}
