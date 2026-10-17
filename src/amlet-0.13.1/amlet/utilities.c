/**
 * utilities.c
 *
 * AMLET Library
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <oratio/cblas.h>

#include "utilities.h"

/**
 * Compute the derivative of a (deterministic) linear utility.
 *
 * @param n is the number of parameters.
 * @param ind is the individual index.
 * @param alt is the alternative index.
 */
void derive_utility(int n, int ind, int alt,
			   double *grad, double ***obs)
{
  cblas_dcopy(n, obs[ind][alt], 1, grad, 1);
}
