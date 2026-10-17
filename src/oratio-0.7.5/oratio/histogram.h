/**
 * histogram.h
 *
 * STAT Library
 * Copyright (C) 2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __STAT_HISTOGRAM_H__
#define __STAT_HISTOGRAM_H__

#ifdef __cplusplus
extern "C" {
#endif

  double st_histogram_init(int n, int b, double min, double max,
			   double *centers);

  void st_histogram_fill_class(double x, int n, int b, double min,
			       double length, double *classes);

  void st_histogram_output(FILE *f, int n, double *centers, int *classes);

#ifdef __cplusplus
}
#endif

#endif  /* __STAT_HISTOGRAM_H__ */
