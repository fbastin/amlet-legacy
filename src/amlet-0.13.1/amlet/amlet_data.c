/**
 * amlet_data.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2004 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "info.h"
#include "amlet_data.h"

#include <amlet/amlet.h>

#include "Train96.h"
#include "Train99.h"

#include <stdio.h>

#include <oratio/ctools.h>

/**
 * Check that expressed choices correspond to available alternatives
 *
 * @return 0 if no error has been found, AMLET_ERROR_CHOICES otherwise.
 */
int amlet_error_choice(FILE *f, int nobs, int *choices, int **availability)
{
  int i;
  int error = 0;

  /* We first tag the invalid observations */
  for (i = 0; i < nobs; i++) {
    if (!availability[i][choices[i]]) {
      choices[i] = -1;
      error = AMLET_ERROR_CHOICES;
      amlet_warning(f, 0, "observation %d is invalid.", i+1);
    }
  }

  return error;
}

/**
 * Check that the number of observations is coherent, i.e. that the sum of the
 * components of the repetition vectors is equal to the number of observations.
 *
 * @param np is number of individuals.
 * @param nobs is the number of observations.
 * @param rep is the vector of repetitions; rep[i] should give the number of
 *   observations delivered by individual i.
 *
 * @return 0 is no error has been found, AMLET_ERROR_REPETITIONS otherwise.
 */
int amlet_check_repetitions(int np, int nobs, int *rep)
{
  int i, j;

  j = rep[0] = 0;
  for (i = 1; i <= np; i++) j += rep[i];

  if (j != nobs) {
    fprintf(stderr, "** FATAL ERROR **: number of observations does not correspond to sum of repetitions.\n");
    return AMLET_ERROR_REPETITIONS;
  }

  return 0;
}

/**
 * Data input
 */

/**
 * Read observations from an external file.
 *
 * @param format is a valid data file format identifier.
 */
int import(int format, char *obsname, double ***obs,
	   int **availability, int *choices,
	   int np, int nind, int nalt, int nobs,
	   int idcensor, int iddep, int *times)
{
  switch (format) {
  case TRAIN96:
    return import_train_1996(obsname, obs, availability, choices,
			     np, nind, nalt, nobs, idcensor, iddep);
    break;
  case TRAIN99:
    return import_train_1999(obsname, obs, availability, choices,
			     np, nind, nalt, nobs, idcensor, iddep, times);
    break;
  default:
    return AMLET_ERROR_UNKNOWN_FORMAT;
    break;
  }
}

/**
 * Export observation in an external file.
 *
 * @param format is a valid data file format identifier.
 */
int export(int format, char *obsname, double ***obs,
	   int **availability, const int *choices,
	   int np, int nind, int nalt, int nobs,
	   int idcensor, int *times)
{
  switch (format) {
  case TRAIN96:
    return export_train_1996(obsname, obs, availability,
			     choices, np, nind, nalt, nobs,
			     idcensor);
    break;
  case TRAIN99:
    return export_train_1999(obsname, obs, availability,
			     choices, np, nind, nalt, nobs,
			     idcensor, times);
    break;
  default:
    return AMLET_ERROR_UNKNOWN_FORMAT;
    break;
  }
}
