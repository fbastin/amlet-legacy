/**
 * error.h
 *
 * Monte Carlo and Quasi-Monte Carlo Library
 * Copyright (C) 2006 Fabian Bastin <bastin@cerfacs.fr>
 *
 * Licensed under the MIT License
 */

#ifndef __MOCA_ERROR_H__
#define __MOCA_ERROR_H__

#ifdef __cplusplus
extern "C" {
#endif

  enum {
    MOCA_OK = 0,
    MOCA_ERROR_DIMENSION,
    MOCA_ERROR_MAX_CALLS,
    MOCA_ERROR_DIGITS
  };

#ifdef __cplusplus
}
#endif

#endif  /* __MOCA_ERROR_H__ */
