/**
 * cauchy_test.c
 *
 * Stat Library
 * Copyright (C) 2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <stdio.h>

#include <oratio/oratio.h>

int st_cauchy_test()
{
  FILE *f;
  Random *r1, *r2;
  unsigned long s1 = 100;
  unsigned long s2 = 100;
  int n  = 200000;

  double *v;

  int i;

  f = fopen("cauchy_direct.plot", "w");
  r1 = ran_cauchy_new(&s1, 1, 0.0, 1.0);
  v = _mallocs(n, double);

  for (i = 0; i < n; i++) v[i] = ran_cauchy_get_val(r1);

  /* Generate the data. */
  st_export_empirical_cdf(f, n, v, 0);

  fclose(f);
  ran_random_free(r1);
  free(v);


  f = fopen("cauchy_ratio.plot", "w");
  s1 = 1;
  r1 = ran_normal_new(&s1, 1, 0.0, 1.0);
  s2 = 200;
  r2 = ran_normal_new(&s2, 1, 0.0, 1.0);
  v = _mallocs(n, double);

  for (i = 0; i < n; i++) v[i] = ran_random_next(r1)/ran_random_next(r2);

  /* Generate the data. */
  st_export_empirical_cdf(f, n, v, 0);

  fclose(f);
  ran_random_free(r1);
  ran_random_free(r2);
  free(v);

  return 0;
}
