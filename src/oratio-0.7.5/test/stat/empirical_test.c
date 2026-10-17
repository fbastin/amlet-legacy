/**
 * empirical_test.c
 *
 * Stat Library
 * Copyright (C) 2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <stdio.h>

#include <oratio/oratio.h>

int st_empirical_test()
{
  FILE *f;
  Random *r;
  unsigned long s = 100;
  int n  = 200000;

  double *v;

  int i;

  f = fopen("normal.plot", "w");
  r = ran_normal_new(&s, 1, 0.0, 1.0);
  v = _mallocs(n, double);

  for (i = 0; i < n; i++) v[i] = ran_normal_get_val(r);

  /* Generate the data. */
  st_export_empirical_cdf(f, n, v, 0);

  fclose(f);
  ran_random_free(r);
  free(v);

  return 0;
}
