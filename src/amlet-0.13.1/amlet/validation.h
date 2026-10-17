/**
 * validation.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2003-2004 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __AMLET_VALIDATION_H__
#define __AMLET_VALIDATION_H__

#include "amlet.h"
#include "amlet_driver.h"

#include <oratio/ran.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Validate estimation results.
   *
   * @param x is the estimated parameters vectors.
   * @param nrand is the number of random coefficients in the model.
   * @param nsim is the number of likelihood simulations to perform.
   */
  void amlet_validate_results(AMLET_Driver *dr, int n, double *x, dcl_data *md,
			      int nrand, Random *ran, int *sampling_types,
			      int nsim);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_VALIDATION_H__ */
