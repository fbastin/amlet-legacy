/**
 * option.h
 *
 * OPHELIA Library
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __OPHELIA_OPTION_H__
#define __OPHELIA_OPTION_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Structure definition
   */
  typedef struct _Option Option;

  /**
   * Financial option structure
   */
  struct _Option {
    double strike;   /* strike price */
    double r;        /* risk-free interest rate */
    double s0;       /* initial capital */
    double sigma;
  };
 
#ifdef __cplusplus
}
#endif

#endif  /* __OPHELIA_BTR_STRUCT_H__ */
