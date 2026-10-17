/**
 * observations.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2006 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __AMLET_OBSERVATIONS_H__
#define __AMLET_OBSERVATIONS_H__

#include "amlet_config.h"
#include "amlet_driver.h"

#include <oratio/ran.h>

#ifdef __cplusplus
extern "C" {
#endif


  /**
   * Exclude invalid observations from the data set.
   *
   * @param np is the number of individiduals.
   *
   * @return the number of rejected observations
   */
  int amlet_reject_observations(int *np, int *nobs, int nalt,
				int **rep, double ****obs, int **choices,
				int ***availability);

  /**
   * Generate observations for a synthetic population.
   *
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
				  int **availability, int *rep);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_OBSERVATIONS_H__ */
