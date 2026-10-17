/**
 * Train96.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
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

  for (i = 0; i < nobs; i++) {
    k = 0;
    for (l = 0; l < n; l++) {
      if (l == iddep-1) {
	/* Read the choice */
	for (j = 0; j < nalt; j++) {
	  fscanf(f, "%d", &itmp);
	  if (itmp) choices[i] = j;
	}
	fscanf(f, "\n");
      }
      else if (l == idcensor-1) {
	/* Read the alternatives  availability */
	for (j = 0; j < nalt; j++) fscanf(f, "%d", &availability[i][j]);
	fscanf(f, "\n");
      }
      else {
	for (j = 0; j < nalt; j++)
	  fscanf(f, "%lg", &(obs[i][j][k]));
	k++;
      }

    }
  }

  return 0;
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
	for (k = 0; k < nalt-1; k++) fprintf(f, "%g ", obs[i][k][line]);
	fprintf(f, "%g\n", obs[i][nalt-1][line++]);
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
