/**
 * observations.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "info.h"
#include "observations.h"

#include <amlet/amlet.h>
#include <amlet/dcl.h>
#include "amlet_analysis.h"
#include "amlet_data.h"
#include "amlet_driver.h"
#include "io.h"
#include "log.h"

#include <oratio/oratio.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

double scale = 1.0;

/**
 * Exclude invalid observations from the data set.
 *
 * @param np is the number of individiduals.
 * @param nalt is the number of alternatives.
 *
 * @return the number of rejected observations.
 */
int amlet_reject_observations(int *np, int *nobs, int nalt,
			      int **rep, double ****obs, int **choices,
			      int ***availability)
{
  int i, j, k;     /* loops index */
  int start, pos, newpos, newstart;

  int n = 0;       /* number of rejected observations */

  pos = 0;
  newpos = 0;

  for (i = 1; i <= *np; i++) {
    start = pos;
    newstart = newpos;

    for (j = pos; j < start+(*rep)[i]; j++) {
      if ((*choices)[j] == -1) {
	/* free the memory associated to (*obs)[j] */
	for (k = 0; k < nalt; k++) free((*obs)[j][k]);
	free((*obs)[j]);
	/* free the memory associated to (*availability)[j] */
	free((*availability)[j]);
      }
      else{
	/* update (*obs[newpos]) */
	(*obs)[newpos] = (*obs)[j];
	(*choices)[newpos] = (*choices)[j];
	(*availability)[newpos] = (*availability)[j];
	newpos++;
      }

      pos++;
    }
    /* update the number of observations for individual i-1 */
    (*rep)[i] = newpos-newstart;
  }

  /* Memory reallocation */
  *nobs -= pos-newpos;
  *obs = _reallocs(*obs, *nobs, double **);
  *choices = _reallocs(*choices, *nobs, int);
  *availability = _reallocs(*availability, *nobs, int *);

  /* Check that no people has 0 observations */
  pos = 0;
  for (i = 1; i <= *np; i++) {
    if ((*rep)[i] == 0) {
      n++;
      amlet_warning(mclog->f, 0,
		    "individual %d has no more observations; now rejected.",
		    i);
      pos++;
    }
    (*rep)[i] = (*rep)[i+pos];
  }
  *np -= pos;

  return n;
}

double set_var(int partype, double *par, int *pos, Random *r,
		      void *data)
{
  double val;

  val = moca_get_random_value(partype, r, par+(*pos), data);
  *pos += moca_count_parameters(partype, data);

  return val;
}

/**
 * @param nalt is the number of alternatives.
 */
void set_choices(double ***obs, double *par, int *partype, int *choice,
		 Random *r, int np, int nalt, int ncoeff, int *rep,
		 DS_List *np_list)
{
  int i, j, k, nobs;
  int pos;
  double best, utility;

  double *c;

  DS_List *tmp;
  void *obj;

  c = _mallocs(ncoeff, double);

  nobs = 0;

  for (i = 0; i < np; i++) {
    tmp = np_list;

    pos = 0;
    for (j = 0; j < ncoeff; j++) {
      if (amlet_is_nonparametric(partype[j])) {
	obj = tmp->data;
	tmp = tmp->next;
      }
      else
	obj = NULL;
      c[j] = set_var(partype[j], par, &pos, r, obj);
    }

    for (j = nobs; j < nobs+rep[i+1]; j++) {
      choice[j] = 0;
      best = cblas_ddot(ncoeff, c, 1, obs[j][0], 1);
      // we perturb the utility with a random error
      best += scale*ran_gumbel(r, -M_EULER, 1.0);

      for (k = 1; k < nalt; k++) {

	utility = cblas_ddot(ncoeff, c, 1, obs[j][k], 1);
	utility += scale*ran_gumbel(r, -M_EULER, 1.0);

	if (best < utility) {
	  best = utility;
	  choice[j] = k;
	}
      }
    }

    nobs += rep[i+1];
  }

  free(c);
}

/**
 * Generate observations
 *
 * @param r is an initialized random numbers generator.
 * @param nalt is the number of alternatives.
 */
void generate_observations(double ***obs, Random *r,
			   double **alt_values, int **alt_types,
			   int np, int nobs, int nalt, int ncoeff, int *rep)
{
  int i, j, k;
  int pos;

  for (i = 0; i < nobs; i++) {
    for (j = 0; j < nalt; j++) {
      pos = 0;
      for (k = 0; k < ncoeff; k++) {
	obs[i][j][k] = set_var(alt_types[j][k], alt_values[j],
			       &pos, r, NULL);
      }
    }
  }
}

/**
 * Simulate observations, corresponding to the alternatives specifications.
 *
 * @param ran is an initialized random numbers generator.
 * @param obs is three-dimensional array; first dimension correspond
 *   to the individu, second to the alternative and the last to the
 *   observations.
 * @param nalt is the number of alternatives.
 * @param rep is the repetition vector: individual i has
 *   rep[i+1]-rep[i] observations, and rep index the first dimension of
 *   obs.
 */
void simulate_observations(Random *ran, int *choices, double ***obs,
			   int **availability, double *par, int *partype,
			   double **alt_values, int **alt_types,
			   int np, int nobs, int nalt, int ncoeff, int *rep,
			   DS_List *np_list)
{
  int i, j;

  int *prefs;

  //#ifdef DEBUG
  nt_log_section(mclog, "Simulate observations");

  prefs = _mallocs(nalt, int);
  for (i = 0; i < nalt; i++) prefs[i] = 0;
  //#endif

  /* Set observations */
  generate_observations(obs, ran, alt_values, alt_types,
			np, nobs, nalt, ncoeff, rep);

  /* Set choices */
  set_choices(obs, par, partype, choices, ran, np, nalt, ncoeff, rep,
	      np_list);

  /* Set alteratives availability */
  /* Currently all alternatives are assumed to be available */
  for (i = 0; i < nobs; i++)
    for (j = 0; j < nalt; j++)
      availability[i][j] = 1;

  /* We determine the choices proportions */
  for (i = 0; i < nobs; i++) prefs[choices[i]]++;

  nt_log_write(mclog, "Proportions: ");
  for (i = 0; i < nalt; i++)
    nt_log_write(mclog, "%g%c", ((double)prefs[i])/np,
		 (i == nalt-1) ? '\n' : ' ');
  nt_log_write(mclog, "\n");

  free(prefs);
}

/**
 * Generate observations for a synthetic population.
 *
 * @param dr is the driver structure.
 * @param ran is an initialized random numbers generator.
 * @param obs is three-dimensional array; first dimension correspond
 *   to the individu, second to the alternative and the last to the
 *   observations.
 * @param rep is the repetition vector: individual i has
 *   rep[i+1]-rep[i] observations, and rep index the first dimension of
 *   obs.
 */
int amlet_generate_observations(amlet_config *ac,
				Random *ran, AMLET_Driver *dr, int nvar,
				int *choices, double ***obs,
				int **availability, int *rep)
{
  int npargen = 0;      /* number of parameters for the generation process */
  int *bgentype = NULL;
  void *np_object;      /* nonparametric object */

  double *par;          /* model parameters */ /* could we use a working vector? */

  DS_List *np_list, *tmp;

  int i, t;

  //  if (verbose) printf("Data generation.\n");

  np_list = ds_list_new();

  /* Initialize the synthetic population parameters. */
  if (dr->array[BGENTYPE].n != nvar)
    return amlet_error(ac, mclog->f, 1, AMLET_ERROR_PARSING,
		       "incorrect number of parameters in bgentype declaration");
  bgentype = _mallocs(nvar, int);
  driver_affect_integer(bgentype, dr->array+BGENTYPE);

  /* Count how many parameters are present. */
  for (i = 0; i < nvar; i++) {
    np_object = NULL;
    t = dr->array[BGENTYPE].array[i].d.i;
    if (amlet_is_nonparametric(t)) {
      np_object = amlet_nonparametric_set(t);
      np_list = ds_list_prepend(np_list, np_object);
    }
    npargen += moca_count_parameters(t, np_object);
  }
  
  np_list = ds_list_reverse(np_list);

  /* Initialize the parameters. */
  if (dr->array[BGEN].n != npargen)
    return amlet_error(ac, mclog->f, 1, AMLET_ERROR_PARSING,
		       "incorrect number of parameters in bgen declaration");
  par = _mallocs(npargen, double);
  driver_affect_real(par, dr->array+BGEN);

  /* Create synthetic data. */
  simulate_observations(ran, choices, obs, availability,
			par, bgentype, dr->alt_values, dr->alt_types,
			data.np, data.nobs, data.nalt, data.nvar, rep,
			np_list);

  if (dr->obsname)
    /* export the observations */
    export(dr->var[FORMAT].i,
	   dr->obsname, obs, availability,
	   choices, data.np, data.nobs, data.nalt, data.nvar,
	   dr->var[IDCENSOR].i, rep);

  //  if (verbose) printf("End of data generation.\n");

  /* Free memory. */

  tmp = np_list;
  for (i = 0; i < nvar; i++) {  
    t = dr->array[BGENTYPE].array[i].d.i;
    if (amlet_is_nonparametric(t)) {
      amlet_nonparametric_free(tmp->data, t);
      tmp = tmp->next;
    }
  }
  ds_list_free(np_list, NULL);

  driver_free_array(dr, BGENTYPE);
  driver_free_array(dr, BGEN);

  free(bgentype);

  return 0;
}
