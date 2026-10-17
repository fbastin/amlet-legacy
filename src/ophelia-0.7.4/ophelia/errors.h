/**
 * errors.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 1999-2004 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __NLP2_ERRORS_H__
#define __NLP2_ERRORS_H__

#ifdef __cplusplus
extern "C" {
#endif

  enum {
    NLP_ERROR_TOO_MANY_ITERATIONS = 1
  };

  static const char *btr_error_msg[] ={
    "",
    "Warning: algorithm stopped due to excess of iterations\n"
  };

#ifdef __cplusplus
}
#endif

#endif
