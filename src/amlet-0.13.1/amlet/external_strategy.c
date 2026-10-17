/**
 * external_strategy.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2004 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * External variable sample size strategy
 */

#include "external_strategy.h"

#include <math.h>

#include <oratio/oratio.h>

//      (dr.var[VARREP].i) ? nt_max(36, rint(dr.var[NREP].i / 40.0)) :

/**
 * External strategies
 */

/**
 *
 * Ref.: Biogeme 1.0
 */
void external_strategy_bierlaire(AmletExternalStrategy *aes,
				 DclOptiMethod *opti)
{
  double optimality_ratio, draw_ratio;
  double gmax;

  gmax = nlp_relative_gradient((opti->b)->value, (opti->b)->x,
			       (opti->b)->grad, NULL, (opti->b)->n, 1.0);
  optimality_ratio = log(gmax)/(aes->tol);
  draw_ratio = (double)aes->current / (double)aes->final;

  if (optimality_ratio >= draw_ratio) {
    aes->current *= aes->factor;
    aes->current = nt_min(aes->current, aes->final);
  }
}

/**
 */

AmletExternalStrategy *amlet_external_strategy_new(int init, int final,
						   double tol)
{
  AmletExternalStrategy *aes = NULL;

  aes = _malloc(AmletExternalStrategy);
  aes->current = aes->init = init;
  aes->final = final;
  aes->factor = 2.0;
  aes->tol = tol;
  aes->method = AMLET_EXTERNAL_STRATEGY_NO_STOPPING;

  return aes;
}

void amlet_external_strategy_free(AmletExternalStrategy* aes)
{
  if (aes) free(aes);
}

int amlet_external_strategy_update(AmletExternalStrategy* aes,
				   DclOptiMethod *opti)
{
  switch (aes->method) {
  case AMLET_EXTERNAL_STRATEGY_BIERLAIRE:
    external_strategy_bierlaire(aes, opti);
    break;
  default:
    aes->current = aes->final;
    break;
  }

  return aes->current;
}

/**
 * Stopping test for the external variable sample size strategy.
 */
int amlet_external_strategy_stopping_test(AmletExternalStrategy* aes)
{
  if (aes->current == aes->final)
    return AMLET_EXTERNAL_STRATEGY_FINAL_SIZE;
  else
    return AMLET_EXTERNAL_STRATEGY_NO_STOPPING;
}
