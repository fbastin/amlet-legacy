/**
 * amlet_analysis.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2004 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __AMLET_ANALYSIS_H__
#define __AMLET_ANALYSIS_H__

#include <amlet/amlet.h>
#include <amlet/amlet_driver.h>

#include <oratio/cblas.h>
#include <oratio/f2c.h>

/**
 * Analysis of AMLET results
 */

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Compute the log-likelihood value when all parameters are set to
   * zero and when only alternative specific constants are considered.
   *
   * @param ismixed is an integer equal to 0 if we face a classical
   *    multinomial model, another value elsewhere.
   */
  void special_likelihood_functions(int ismixed, double *ll0, double *llc,
				    dcl_data *md, AMLET_Driver *dr,
				    double *work);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_ANALYSIS_H__ */
