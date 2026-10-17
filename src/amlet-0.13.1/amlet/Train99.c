/**
 * Train99.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "info.h"
#include "Train99.h"

#include <stdio.h>
#include <string.h>

#include <amlet/amlet.h>

#include "amlet_data.h"
#include "io.h"

#include "oratio/ctools.h"

/**
 * Read observations written using the Train's format (1999).
 *
 * @param obsname is the name of the file containing the data.
 * @param obs is a three-dimensional array containing the observations
 *   delivered by all the individuals. obs[i][j][k] contains the kth
 *   parameter for the alternative j, in the ith observation.
 * @param times: times[i] (i > 0) represents the number of observations for
 *   individual i-1.
 *
 * @return 0 is the importations has been successfully performed, an
 *   appropriate error code otherwise.
 */
int import_train_1999(const char *obsname, double ***obs,
		      int **availability, int *choices,
		      int np, int nobs, int nalt, int nvar,
		      int idcensor, int iddep, int *times)
{
  FILE *f;

  int i, j, k, l, n;

  char *name;

  name = _mallocs(strlen(obsname)+7, char);
  strcpy(name, obsname);

  strcat(name, "_times");
  if (!(f = c_fopen(name, "r"))) {
    amlet_warning(stderr, 1,
		  "unable to open the repetitions file %s for importation\n",
		  name);
    if (np == nobs) {
      amlet_warning(stderr, 1, "I assume one observation per individual");
      for (i = 1; i <= np; i++) times[i] = 1;
    }
    else {
      free(name);
      return AMLET_ERROR_DATA_READING;
    }
  } else {
    i = amlet_read_repetitions(name, times, np);
    if (i) {
      free(name);
      return i;
    }
  }

  name[strlen(obsname)] = '\0';

  /* Observations */
  strcat(name, "_x");
  if (!(f = c_fopen(name, "r"))) {
    amlet_warning(stderr, 1,
		  "unable to open file %s for observations importation\n",
		  name);
    free(name);
    return AMLET_ERROR_FILE_OPEN;
  }
  name[strlen(obsname)] = '\0';

  n = (idcensor) ? nvar+1 : nvar;

  for (i = 0; i < nobs; i++) {
    k = 0;
    for (l = 0; l < n; l++) {
      if (l == idcensor-1) {
	/* Read the alternatives  availability */
	for (j = 0; j < nalt; j++) fscanf(f, "%d", &availability[i][j]);
      }
      else {
	for (j = 0; j < nalt; j++)
	  fscanf(f, "%lg", &(obs[i][j][k]));
	k++;
      }
    }
  }
 
  /* Choice */
  strcat(name, "_y");
  if (!(f = c_fopen(name, "r"))) {
    amlet_warning(stderr, 1,
		  "unable to open choices file %s for importation\n",
		  name);
    free(name);
    return AMLET_ERROR_FILE_OPEN;
  }
  name[strlen(obsname)] = '\0';

  for (i = 0; i < nobs; i++) {
    fscanf(f, "%d\n", choices+i);
    choices[i]--;
  }
  free(name);

  return 0;
}

/**
 * Save observations written using the Train's format (1999).
 */
int export_train_1999(const char *obsname, double ***obs,
		      int **availability, const int *choices,
		      int np, int nobs, int nalt, int nvar,
		      int idcensor, int *times)
{
  FILE *f;
  char *name;

  int i, j, k;
  int censor, line;

  name = _mallocs(strlen(obsname)+7, char);
  strcpy(name, obsname);

  censor = idcensor ? 1 : 0;
  idcensor--;   /* idcensor can then be directly used as an index */

  /* Observations */
  strcat(name, "_x");
  if (!(f = c_fopen(name, "w"))) {
    amlet_warning(stderr, 1,
		  "unable to open file %s for observations exportation\n",
		  name);
    free(name);
    return AMLET_ERROR_FILE_OPEN;
  }
  name[strlen(obsname)] = '\0';

  for (i = 0; i < nobs; i++) {
    line = 0;
    for (j = 0; j < nvar+censor; j++) {
      if (j == idcensor) {
	/* Alternatives availability */
	for (k = 0; k < nalt-1; k++) fprintf(f, "%d ", availability[i][k]);
	fprintf(f, "%d\n", availability[i][nalt-1]);
      }
      else {
	for (k = 0; k < nalt-1; k++) fprintf(f, "%g ", obs[i][k][line]);
	fprintf(f, "%g\n", obs[i][nalt-1][line++]);
      }
    }
  }
  fclose(f);

  /* Choice */
  strcat(name, "_y");
  if (!(f = c_fopen(name, "w"))) {
    amlet_warning(stderr, 1,
		  "unable to open file %s for observations exportation\n",
		  name);
    free(name);
    return AMLET_ERROR_FILE_OPEN;
  }
  name[strlen(obsname)] = '\0';

  for (i = 0; i < nobs; i++)
    fprintf(f, "%d\n", choices[i]+1);

  fclose(f);

  strcat(name, "_times");
  if (!(f = c_fopen(name, "w"))) {
    amlet_warning(stderr, 1,
		  "unable to open file %s for observations exportation\n",
		  name);
    free(name);
    return AMLET_ERROR_FILE_OPEN;
  }
  name[strlen(obsname)] = '\0';

  for (i = 1; i <= np; i++)
    fprintf(f, "%d\n", times[i]);

  fclose(f);

  free(name);

  return 0;
}
