/**
 * test_ran.h
 *
 * ORATIO Library.
 * Copyright (C) 2007-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_TEST_RAN_H__
#define __ORATIO_TEST_RAN_H__

#include <oratio/ran.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Perform some tests on uniform U(0,1) distributions.
   *
   * @param f is a write-open file.
   * @param r is an initialized random numbers generator.
   */
  void uniforms_check(FILE *f, Random *r);

  /**
   * Generic methods.
   */

  void test_acceptance_rejection();

  /**
   * Univariate distributions.
   */

  void test_exponential();

  void test_normal();

  void test_truncated_normal();

  /**
   * Multivariate distributions.
   */

  /**
   * Test of multivariate normal distribution generation routines.
   */
  void test_mvn();

#ifdef __cplusplus
}
#endif

#endif
