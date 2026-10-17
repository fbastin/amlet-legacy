/**
 * st_obs.h
 *
 * ORATIO Library
 * Statistical routines
 * Copyright (C) 1999-2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_STAT_OBS_H__
#define __ORATIO_STAT_OBS_H__

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct _st_obs {
    double prob;   /* mass probability or density at the observation */
    double *grad;  /* gradient of the probability with respect to the model
		      parameters */
  } ST_Obs;

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_STAT_OBS_H__ */
