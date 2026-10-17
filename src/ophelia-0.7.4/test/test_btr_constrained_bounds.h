/**
 * test_btr_constrained_bounds.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2004-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __OPHELIA_TESTING_BTR_CONSTRAINED_BOUNDS_H__
#define __OPHELIA_TESTING_BTR_CONSTRAINED_BOUNDS_H__

#include <ophelia/nlp.h>

#include <stdio.h>

#include <oratio/ntools.h>

#ifdef __cplusplus
extern "C" {
#endif

  int btr_bounds_constrained_tests(NTLog *log, BTR *b);

#ifdef __cplusplus
}
#endif

#endif /* __OPHELIA_TESTING_BTR_CONSTRAINED_BOUNDS_H__ */
