/**
 * test_truncated_normal.h
 *
 * Random Numbers Library
 * Copyright (C) 2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <oratio/oratio.h>

#include <stdio.h>

void test_ran_sp()
{
  Random *r;
  double knots[] = {0.0, 1.0/3.0, 2.0/3.0,  1.0};
  double P[] = { 1, 1, 1, 2, 2, 2 };
  //double P[] = { -10.317211, -0.79378037, -0.68347875, 0.12695979, 0.12695979, 0.12695979 };

  unsigned long seed = 1;

  double pts[100000];

  FILE *f;
  int i;

  double l, u, a, b, c, m;

  r = ran_sp_spline_new(&seed, 1, 4, knots, P);

  /* Dichotomic search for 0 */
  a = 0.0;
  b = 1.0;
  l = ran_sp_spline_U(r, a);
  u = ran_sp_spline_U(r, b);
  while (b-a >= 1e-10) {
    c = a+0.5*(b-a);
    m = ran_sp_spline_U(r, c);
    if (m < 0.0)
      a = c;
    else if (m > 0.0)
      b = c;
    else break;
  }

  /* Truncate the distribution at 0.0. */
  //  ran_random_set_truncation(r, 0.0, c, 0.0);

  for (i = 0; i < 100000; i++) pts[i] = ran_random_next(r);
  ds_sort_dheapsort(100000, pts);

  f = fopen("ran_sp_spline.txt", "w");
  for (i = 0; i < 100000; i++)
    fprintf(f, "%g %g\n", pts[i], (double)(i+1)/100000.0);
  fclose(f);

  ran_sp_spline_struct_free(r);

  ran_random_free(r);
}
