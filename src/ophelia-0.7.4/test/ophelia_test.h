/**
 * ophelia_test.h
 *
 * OPHELIA Library
 * Copyright (C) 2004-2011 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __OPHELIA_TEST_H__
#define __OPHELIA_TEST_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Test of the Golden Search algorithm.
   *
   * @return 0.
   */
  int op_gss_test();

  void test_matrix();

#ifdef __cplusplus
}
#endif

#endif  /* __OPHELIA_TEST_H__ */
