/**
 * external_strategy.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2004 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * External variable sample size strategy
 */

#ifndef __AMLET_EXTERNAL_STRATEGY_H__
#define __AMLET_EXTERNAL_STRATEGY_H__

#include "optimization.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * External strategy identifiers
   */

  enum {
    AMLET_EXTERNAL_STRATEGY_NONE,
    AMLET_EXTERNAL_STRATEGY_BIERLAIRE
  };

  typedef struct {
    int init;                 /* initial sample size */
    int current;              /* current sample size */
    int final;                /* final sample size */
    int factor;               /* factor multiplier */
    int method;               /* external strategy identifier */
    double tol;               /* tolerance for optimization */
  } AmletExternalStrategy;

  enum {
    AMLET_EXTERNAL_STRATEGY_NO_STOPPING = 0,
    AMLET_EXTERNAL_STRATEGY_FINAL_SIZE
  };

  AmletExternalStrategy *amlet_external_strategy_new(int init, int final,
						     double tol);

  void amlet_external_strategy_free(AmletExternalStrategy* aes);

  int amlet_external_strategy_update(AmletExternalStrategy* aes,
				     DclOptiMethod *opti);

  /**
   * Stopping test for the external variable sample size strategy.
   */
  int amlet_external_strategy_stopping_test(AmletExternalStrategy* aes);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_EXTERNAL_STRATEGY_H__ */
