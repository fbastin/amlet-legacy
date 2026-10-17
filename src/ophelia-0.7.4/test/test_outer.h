/**
 * test_btr_outer.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __OPHELIA_TEST_OUTER_H__
#define __OPHELIA_TEST_OUTER_H__

#include <ophelia/nlp.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Test the barrier parameter evolution.
   */
  int nlp_test_barrier_parameter();

  int nlp_outer_tests(BTR *b);

#ifdef __cplusplus
}
#endif

#endif /* __OPHELIA_TEST_OUTER_H__ */
