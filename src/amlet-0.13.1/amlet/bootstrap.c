/**
 * bootstrap.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "bootstrap.h"

#include <oratio/ran.h>
#include <oratio/oratio.h>

/**
 * @param idx if 0, in practice, no resampling is performed.
 */
void bootstrap_resample(dcl_data *md, int idx, int *n,
			int nobs, double ***obs, int *rep,
			int *choices, int **availability,
			double ****n_obs, int **nrep,
			int **nchoices, int ***navailability)
{
  Random *r = NULL;
  unsigned long state = 12345;

  int *p;

  int i, j, k;

  p = _mallocs(md->np, int);

  *nrep = _reallocs(*nrep, md->np+1, int);
  (*nrep)[0] = 0;

  if (idx) {
    r = ran_random_new(0, &state, 1);
    ran_random_set_generator(r, RAN_MRG32K3A, &state, 1);

    for (i = 0; i < idx; i++) ran_random_resetnextsubstream(r);

    for (i = 0; i < md->np; i++) p[i] = ran_intrandom(r, md->np);

    /* Create the new and repetitions vector. */
    *n = 0;
    for (i = 1; i <= md->np; i++) {
      (*nrep)[i] = rep[p[i-1]+1];
      *n += (*nrep)[i];
    }
    *n_obs = _reallocs(*n_obs, *n, double **);
  }
  else {
    *n_obs = _reallocs(*n_obs, nobs, double **);
    *nchoices = _reallocs(*nchoices, nobs, int);

    for (i = 1; i <= md->np; i++) (*nrep)[i] = rep[i];    

    *n = nobs;
  }

  *navailability = _reallocs(*navailability, *n, int *);

  if (idx) {
    /*    for (i = 0; i <= md->np; i++) printf("%d ", rep[i]);
	  printf("\n"); */
    /* Reorganize the original repetitions array. */
    for (i = 1; i <= md->np; i++) rep[i] += rep[i-1];
    /*    for (i = 0; i <= md->np; i++) printf("%d ", rep[i]);
	  printf("\n"); */

    *nchoices = _reallocs(*nchoices, *n, int);
    k = 0;
    for (i = 0; i < md->np; i++) {
      for (j = 0; j < (*nrep)[i+1]; j++) {
	(*n_obs)[k+j] = obs[rep[p[i]]+j];
	(*nchoices)[k+j] = choices[rep[p[i]]+j];
	(*navailability)[k+j] = availability[rep[p[i]]+j];
      }
      k += (*nrep)[i+1];
    }
    md->nobs = *n;
    ran_random_free(r);

    /* Recover the original repetitions array. */
    for (i = md->np; i > 0; i--) rep[i] -= rep[i-1];
  }
  else {
    for (i = 0; i < nobs; i++) {
      (*n_obs)[i] = obs[i];
      (*nchoices)[i] = choices[i];
      (*navailability)[i] = availability[i];
    }
  }

  /* Reorganize the repetitions array. */
  for (i = 1; i <= md->np; i++) (*nrep)[i] += (*nrep)[i-1];
  /*    for (i = 0; i <= md->np; i++) printf("%d ", rep[i]);
	printf("\n"); */

  free(p);
}
