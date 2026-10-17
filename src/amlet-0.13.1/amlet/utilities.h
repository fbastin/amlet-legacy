/**
 * utilities.h
 *
 * AMLET Library
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __AMLET_UTILITIES_H__
#define __AMLET_UTILITIES_H__

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Compute the derivative of a (deterministic) linear utility.
 *
 * @param n is the number of parameters.
 * @param ind is the individual index.
 * @param alt is the alternative index.
 */
void derive_utility(int n, int ind, int alt,
		    double *grad, double ***obs);

#ifdef __cplusplus
}
#endif

#endif /* __AMLET_UTILITIES_H__ */
