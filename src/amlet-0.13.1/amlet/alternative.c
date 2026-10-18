/**
 * alternative.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2004 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the MIT License
 */

/* This code is for preliminary tests only; it should not been used on
   critical applications! In consequence it is not offically
   documented. */

/* QUID DES PANELS??? */
/* KNOWN BUGS: availability not yet supported; temporary file not well
   defined (imposed name) */

#include <alternative.h>
#include <optimization.h>

#include <oratio/ran.h>

#include <amlet_data.h>
#include "amlet_driver.h"

int a_in_alt(int a, int n, int *alt)
{
  int i;
  for (i = 0; i < n; i++) {
    if (a == alt[i]) return 1;
    if (alt[i] == -1) return 0;
  }

  return 0;
}

/***************************************************************************
 *
 ***************************************************************************
 *
 * #PRE:  n: number of alternatives subsamples
 *
 ***************************************************************************/
int amlet_alternatives_sampling(AMLET_Driver *dr, double ***obs,
				int **availability, int *choices,
				int np, int nobs, int nalt, int nvar,
				int *rep, int n)
{
  int i, j, k;
  unsigned long seed = 0;
  int a;
  Random *r;
  double ***new_obs;
  int *new_choices;
  int *new_rep;

  int *chosen_alternatives;

  /* AMLET <= 0.13.1 looped forever */
  if (n < 2 || n > nalt)
    c_error("alternatives_sampling = %d: from 2 to nalt = %d alternatives can be sampled", n, nalt);

  chosen_alternatives = _mallocs(n, int);

  /* integers in [0, nalt), from seed1; AMLET <= 0.13.1 swapped the
     arguments (nalt seeds, read past seed, and integers in [0, 1)) and
     took a seed from the clock */
  seed = dr->var[SEED1].i;
  r = ran_intrandom_new(&seed, 1, nalt);

  /* New observations table */
  new_obs = _mallocs(nobs, double **);
  for (i = 0; i < nobs; i++) {
    new_obs[i] = _mallocs(n, double *);
    for (j = 0; j < n; j++) new_obs[i][j] = _mallocs(nvar, double);
  }

  new_choices = _mallocs(nobs, int);
  new_rep = _mallocs(np+1, int);

  for (i = 1; i <= np; i++) new_rep[i] = rep[i]-rep[i-1];

  for (i = 0; i < nobs; i++) {
    chosen_alternatives[0] = choices[i];
    for (j = 1; j < n; j++) chosen_alternatives[j] = -1;
    new_choices[i] = 0;
    /* Copy observations from the chosen alternatives */
    /* des memcpy seraient plus efficaces */
    for (k = 0; k < nvar; k++)
      new_obs[i][0][k] = obs[i][choices[i]][k];
    for (j = 1; j < n; j++) {
      do {
	a = ran_intrandom_get_val(r);
      } while (a_in_alt(a, n, chosen_alternatives));
      chosen_alternatives[j] = a;
      for (k = 0; k < nvar; k++)
	new_obs[i][j][k] = obs[i][a][k];
    }
  }

  export(dr->var[FORMAT].i,
	 "intermediaire", new_obs, availability,
	 new_choices, np, nobs, n, nvar,
	 dr->var[IDCENSOR].i, new_rep);

  /* Free observations */
  for (i = 0; i < nobs; i++) {
    for (j = 0; j < n; j++) free(new_obs[i][j]);
    free(new_obs[i]);
  }
  free(new_obs);   /* AMLET <= 0.13.1 freed new_obs[0] twice */

  ran_random_free(r);

  free(new_choices);
  free(new_rep);
  free(chosen_alternatives);

  return 1;
}
