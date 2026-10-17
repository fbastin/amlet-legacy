/**
 * empirical.c
 *
 * ORATIO Library
 * Copyright (C) 2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <stdio.h>

#include <oratio/dstruct.h>

/**
 * Export in a file empiral cdf coordinates.
 * Assumptions: all the realizations are different.
 *
 * @param f is the file where the data have to be written.
 * @param n is the number of points to consider.
 * @param v is the observations vector.
 * @param sorted has to be set to 1 if the data are sorted on input, 0 else.
 */
void st_export_empirical_cdf(FILE *f, int n, double *v, int sorted)
{
  int i, j, k;
  double u, factor;

  if (!sorted) ds_sort_dheapsort(n, v);

  i = 0;
  factor = 1.0/n;
  while (i < n) {
    j = 1;
    k = 0;
    while (!k) {
      if (i + j < n)
	if (v[i] == v[i+j]) j++;
	else k = 1;
      else k = 1;
    }
    u = factor*(i+j);
    fprintf(f, "%g %g\n", v[i], u);
    i += j;
  }
}
