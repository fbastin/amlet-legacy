/**
 * moca.h
 *
 * Monte Carlo and Quasi-Monte Carlo Library
 * Copyright (C) 2005 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_MOCA_MLHS_H__
#define __ORATIO_MOCA_MLHS_H__

#include <stdio.h>

#include <oratio/oratio.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * @param ran is an initialized random numbers generator.
   */
  void mlhs_generate_sample(Random *ran, double **sample,
			    int n, int s);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_MOCA_MLHS_H__ */
