/**
 * ntools_test.h
 *
 * ORATIO Library
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_NTOOLS_TEST_H__
#define __ORATIO_NTOOLS_TEST_H__

#include <oratio/oratio.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Error codes.
   */
  enum {
    NT_TEST_START = 0,
    NT_TEST_ALGEBRA,
    NT_TEST_BSPLINES,
    NT_TEST_MODULAR,
    NT_TEST_POLYNOMIAL,
    NT_TEST_PRIMES,
    NT_TEST_UTILITIES,
    NT_TEST_END
  };

  /**
   * Test functions.
   */

  /**
   * Algebra test routines.
   */
  int nt_algebra_test();

  /**
   * B-Splines test routine.
   */
  int nt_bsplines_check();

  /**
   * Machine characteristics test.
   */
  int nt_machar_test();

  /**
   * Tests on functions on natural numbers.
   */
  int nt_natural_test();

  /**
   * Modular arithmetic test routine.
   *
   * @return 0 if the test has been successful, some error code otherwise.
   */
  int nt_modular_check();

  /**
   * Polynomial routines test.
   *
   * @return 0 if the test has been successful, some error code otherwise.
   */
  int nt_polynomial_test();

  /**
   * Prime numbers test routine.
   */
  int primes_test();

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_NTOOLS_TEST_H__ */
