/**
 * test_truncatednormal.c
 *
 * Random Numbers Library
 * Copyright (C) 2007 Fabian Bastin <bastin@iro.umontreal.ca>,
 *   Cinzia Cirillo <ccirillo@umd.edu>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "test_ran.h"

#include <stdio.h>

void test_truncated_normal()
{
  int i;
  const int n = 100000;
  double normal[n], truncnormal[n];

  Random *r;

  FILE *f1, *f2;

  r = ran_random();

  for (i = 0; i < n; i++) {
    normal[i] = ran_normal_icdf(r, 2, 2);
    truncnormal[i] = ran_truncatednormal_icdf(r, 0, 4, 2, 2);
  }

  ran_random_free(r);

  f1 = fopen("normal_t1.txt", "w");
  f2 = fopen("normal_t2.txt", "w");

  for (i = 0; i < n; i++) {
    fprintf(f1, "%g\n", normal[i]);
    fprintf(f2, "%g\n", truncnormal[i]);
  }

  fclose(f1);
  fclose(f2);
}

