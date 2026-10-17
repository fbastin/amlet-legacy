/**
 * caffeine.h
 *
 * ORATIO Library
 * Copyright (C) 2006-2008 <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_CAFFEINE_H__
#define __ORATIO_CAFFEINE_H__

#include <oratio/oratio.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Lattice tests
   *
   * @param ran is an initialized random numbers generator.
   *
   * @return 0 if the test has been successful, an error code otherwise.
   */
  int moca_lattice_check(Random *ran);

  /**
   * Sobol tests
   *
   * @param ran is an initialized random numbers generator.
   *
   * @return 0 if the test has been successful, an error code otherwise.
   */
  int moca_sobol_check(Random *ran);

  int moca_test_distributions(double **sample);

  int moca_test_rqmc(Random *r);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_CAFFEINE_H__ */
