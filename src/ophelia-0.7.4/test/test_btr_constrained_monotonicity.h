/**
 * test_btr_constrained_bounds.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __OPHELIA_TESTING_BTR_CONSTRAINED_MONOTONICITY_H__
#define __OPHELIA_TESTING_BTR_CONSTRAINED_MONOTONICITY_H__

#include <ophelia/nlp.h>

#include <stdio.h>

#include <oratio/ntools.h>

#ifdef __cplusplus
extern "C" {
#endif

  int btr_monotonicity_constrained_tests(NTLog *log, BTR *b);

#ifdef __cplusplus
}
#endif

#endif /* __OPHELIA_TESTING_BTR_CONSTRAINED_MONOTONICITY_H__ */
