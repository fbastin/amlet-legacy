/**
 * histogram.c
 *
 * STAT Library
 * Copyright (C) 2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <math.h>
#include <stdio.h>

double st_histogram_init(int n, int b, double min, double max,
			 double *centers)
{
  double length;
  int s = 0, i;

  if (!b) s = 2;

  /* Compute the center of each classes. */
  length = (max-min)/(n-s);
  min += (b) ? -length/2.0 : length/2.0;

  for (i = 0; i < n; i++) centers[i] = min + i*length;

  return length;
}

void st_histogram_fill_class(double x, int n, int b, double min,
			     double length, double *classes)
{
  int s = 0, i;

  if (!b) s = 1;
  i = floor((x-min)/length)+s;
  if (i >= n) i = n-1;

  classes[i]++;
}

void st_histogram_output(FILE *f, int n, double *centers, int *classes)
{
  int i;

  for (i = 0; i < n; i++) fprintf(f, "%g %d\n", centers[i], classes[i]);
}
