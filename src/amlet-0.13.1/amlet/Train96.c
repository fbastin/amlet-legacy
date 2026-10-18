/**
 * Train96.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "Train96.h"
#include "info.h"

#include <stdio.h>

#include <amlet/amlet.h>

#include "io.h"

#include <oratio/ctools.h>

/**
 * Train's format file importation (1996 Gauss code).
 *
 * @param idcensor must be set to 0 if there is no line giving the
 *   availibility of the alternatives
 */
int import_train_1996(char *obsname, double ***obs,
		      int **availability, int *choices,
		      int np, int nobs, int nalt, int nvar,
		      int idcensor, int iddep)
{
  FILE *f;
  int i, j, k, l;

  int n;

  int itmp;

  if (!(f = c_fopen(obsname, "r"))) {
    amlet_warning(mclog->f, 1,
		  "unable to open file %s for observations importation\n",
		  obsname);
    return AMLET_ERROR_FILE_OPEN;
  }

  n = (idcensor) ? nvar+2 : nvar + 1;

  /* AMLET <= 0.13.1 read the choice row as one more variable, past the
     arrays, when iddep was not set. */
  if (iddep < 1 || iddep > n || iddep == idcensor) {
    amlet_warning(stderr, 1,
		  "iddep = %d: with the Train96 format, iddep must give the row of the choice in each record (1 to %d)", iddep, n);
    fclose(f);
    return AMLET_ERROR_DATA_READING;
  }

  for (i = 0; i < nobs; i++) {
    k = 0;
    for (l = 0; l < n; l++) {
      if (l == iddep-1) {
	/* Read the choice: one 1, and 0 for the other alternatives. AMLET
	   <= 0.13.1 kept the last nonzero value, and an undefined choice
	   when there was none. */
	int nchosen = 0;

	for (j = 0; j < nalt; j++) {
	  if (fscanf(f, "%d", &itmp) != 1) goto Short;
	  if (itmp != 0 && itmp != 1) nchosen = -nalt;
	  if (itmp) {
	    choices[i] = j;
	    nchosen++;
	  }
	}
	if (nchosen != 1) {
	  amlet_warning(stderr, 1, "%s: record %d: the choice row (row %d) must contain one 1, the other values being 0", obsname, i+1, iddep);
	  fclose(f);
	  return AMLET_ERROR_DATA_READING;
	}
	if (fscanf(f, "\n")) {}
      }
      else if (l == idcensor-1) {
	/* Read the alternatives  availability */
	for (j = 0; j < nalt; j++)
	  if (fscanf(f, "%d", &availability[i][j]) != 1) goto Short;
	if (fscanf(f, "\n")) {}
      }
      else {
	for (j = 0; j < nalt; j++)
	  if (fscanf(f, "%lg", &(obs[i][j][k])) != 1) goto Short;
	k++;
      }

    }
  }

  fclose(f);
  return 0;

 Short:
  /* AMLET <= 0.13.1 went on with undefined values. */
  amlet_warning(stderr, 1, "%s: missing or invalid value in record %d of %d", obsname, i+1, nobs);
  fclose(f);
  return AMLET_ERROR_DATA_READING;
}

/**
 * Save observations at the Gauss format.
 *
 * @param idcensor indexes the line where the alternatives accessibility
 *   is expressed.
 */
int export_train_1996(char *obsname, double ***obs,
		      int **availability, const int *choices,
		      int np, int nobs, int nalt, int nvar,
		      int idcensor)
{
  FILE *f;

  int i, j, k;
  int censor, line;

  if (!(f = c_fopen(obsname, "w"))) {
    amlet_warning(mclog->f, 1,
		  "unable to open file %s for observations exportation\n",
		  obsname);
    return AMLET_ERROR_FILE_OPEN;
  }

  censor = idcensor-- ? 1 : 0;

  for (i = 0; i < nobs; i++) {
    line = 0;
    for (j = 0; j < nvar+censor; j++) {
      if (j == idcensor) {
	/* Availability of alternatives */
	for (k = 0; k < nalt-1; k++) fprintf(f, "%d ", availability[i][k]);
	fprintf(f, "%d\n", availability[i][nalt-1]);
      }
      else {
	/* %.17g writes the values exactly; AMLET <= 0.13.1 wrote 6 digits */
	for (k = 0; k < nalt-1; k++) fprintf(f, "%.17g ", obs[i][k][line]);
	fprintf(f, "%.17g\n", obs[i][nalt-1][line++]);
      }
    }
  //      else if (j == iddep) {
	/* Choice */
	for (k = 0; k < nalt-1; k++) {
	  if (k != choices[i])
	    fprintf(f, "0 ");
	  else
	    fprintf(f, "1 ");
	}
	if (nalt-1 != choices[i])
	  fprintf(f, "0\n");
	else
	  fprintf(f, "1\n");
	//q      }

  }

  fclose(f);

  return 0;
}
