/**
 * stat_empirical.h
 *
 * STAT Library
 * Copyright (C) 2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_STAT_EMPIRICAL_H__
#define __ORATIO_STAT_EMPIRICAL_H__

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Export in a file empiral cdf coordinates.
   *
   * @param f is the file where the data have to be written.
   * @param n is the number of points to consider.
   * @param v is the observations vector.
   * @param sorted has to bet set to 1 if the data are sorted on input, 0 else.
   */
  void st_export_empirical_cdf(FILE *f, int n, double *v, int sorted);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_STAT_EMPIRICAL_H__ */
