/**
 * moca_io.h
 *
 * Monte Carlo and Quasi-Monte Carlo Library
 * Copyright (C) 2005 Fabian Bastin <bastin@cerfacs.fr>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __MOCA_IO_H__
#define __MOCA_IO_H__

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Write a vector of draws in the file f.
   *
   * @param f non null; it represents the file where the draws have to be saved
   * @param v: vector to write in f
   * @param s: dimension
   */
  void moca_print_vector(FILE *f, double *v, int s);

#ifdef __cplusplus
}
#endif

#endif  /* __MOCA_IO_H__ */
