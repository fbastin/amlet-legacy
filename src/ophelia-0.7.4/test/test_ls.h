/**
 * test_ls.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006 Fabian Bastin <bastin@cerfacs.fr>
 *
 * Licensed under the MIT License
 */

#ifndef __OPHELIA_TEST_LS_H__
#define __OPHELIA_TEST_LS_H__

#include <ophelia/nlp.h>

#include <oratio/ntools.h>

#ifdef __cplusplus
extern "C" {
#endif

  int test_nlp_least_squares(NTLog *log, BTR *btr);

#ifdef __cplusplus
}
#endif

#endif  /* __OPHELIA_TEST_LS_H__ */
