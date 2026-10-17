/**
 * memory.c
 *
 * Discrete Choice Library
 * Copyright (C) 2005-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/ctools.h>

/**
 * Some memory allocation for data structures used in discrete choice
 * models estimation.
 *
 * @param nobs is the number of observations.
 * @param nalt is the number of alternatives.
 * @param censor is equal to 1 if all alternatives are available for
 *   all users, 0 else.
 */
void dcl_memory_allocate(int nobs, int nalt, int nvar, int censor,
			 int **choices, double ****obs, int ***availability)
{
  int i, j;

  /* Memory allocation */
  *choices = _mallocs(nobs, int);

  /* observations */
  *obs = _mallocs(nobs, double **);
  for (i = 0; i < nobs; i++) {
    (*obs)[i] = _mallocs(nalt, double *);
    for (j = 0; j < nalt; j++)
      (*obs)[i][j] = _mallocs(nvar, double);
  }

  /* alternatives availability */
  *availability = _mallocs(nobs, int *);
  for (i = 0; i < nobs; i++)
    (*availability)[i] = _mallocs(nalt, int);
  if (!censor) {
    /* all alternatives are always available */
    for (i = 0; i < nobs; i++)
      for (j = 0; j < nalt; j++) (*availability)[i][j] = 1;
  }
}

/**
 * @param nobs is the number of observations.
 * @param nalt is the number of alternatives.
 */
void dcl_memory_deallocate(int nobs, int nalt,
			   int *choices, double ***obs, int **availability)
{
  int i, j;

  if (choices) free(choices);
  
  if (availability) {
    for (i = 0; i < nobs; i++) free(availability[i]);
    free(availability);
  }

  if (obs) {
    /* Why is it problematic??? */
    for (i = 0; i < nobs; i++) {
      for (j = 0; j < nalt; j++) free(obs[i][j]);
      free(obs[i]);
    }
    free(obs);
  }
}
