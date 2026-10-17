/**
 * test_minpack2.h
 *
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * This file is licensed under BSD License.
 */

#ifndef __OPHELIA_TESTING_MINPACK2_H__
#define __OPHELIA_TESTING_MINPACK2_H__

#include <ophelia/nlp.h>

#include <oratio/ntools.h>

#ifdef __cplusplus
extern "C" {
#endif

  int test_minpack2_ls(NTLog *log, BTR *b);

  /**
   * Problem tests proposed in TRON 1.2, selected from the MINPACK-2
   * problems collection.
   */
  int test_tron(NTLog *log, BTR *b);

#ifdef __cplusplus
}
#endif

#endif  /* __OPHELIA_TESTING_MINPACK2_H__ */
