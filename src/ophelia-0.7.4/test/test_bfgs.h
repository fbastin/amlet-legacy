/**
 * test_btr_unconstrained.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2004-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __OPHELIA_TEST_BFGS_UNCONSTRAINED_H__
#define __OPHELIA_TEST_BFGS_UNCONSTRAINED_H__

#include <ophelia/nlp.h>

#include <oratio/ntools.h>

#ifdef __cplusplus
extern "C" {
#endif

  int test_bfgs_unconstrained(NTLog *log, BFGS *b);

#ifdef __cplusplus
}
#endif

#endif /* __OPHELIA_TEST_BFGS_UNCONSTRAINED_H__ */
