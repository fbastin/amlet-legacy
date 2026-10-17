/**
 * bootstrap.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __AMLET_BOOTSTRAP_H__
#define __AMLET_BOOTSTRAP_H__

#include <ophelia/data.h>

#ifdef __cplusplus
extern "C" {
#endif

  void bootstrap_resample(dcl_data *md, int idx, int *n,
			  int nobs, double ***obs, int *rep,
			  int *choices, int **availability,
			  double ****n_obs, int **nrep,
			  int **nchoices, int ***navailability);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_BOOTSTRAP_H__ */
