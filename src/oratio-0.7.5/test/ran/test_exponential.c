/**
 * test_exponential.c
 *
 * ORATIO Library
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "test_ran.h"

#include <oratio/stat.h>

#include <stdio.h>

void test_exponential()
{
  int i;
  const int n = 100000;

  double values[n];

  double mean, var;

  Random *r;

  FILE *f;

  printf("*****************************************\n");
  printf("* Exponential test\n");
  printf("*****************************************\n");

  r = ran_random();

  for (i = 0; i < n; i++) values[i] = ran_exponential(r, 2);

  ran_random_free(r);

  f = fopen("exponential.txt", "w");

  for (i = 0; i < n; i++) fprintf(f, "%g\n", values[i]);

  st_mean_variance(values, 1, 0, n, &mean, &var);

  printf("Mean: %g, Variance: %g\n\n", mean, var);

  fclose(f);
}

