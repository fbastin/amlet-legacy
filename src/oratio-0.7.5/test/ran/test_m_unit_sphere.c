/**
 * test_.c
 *
 * Random Numbers Library
 * Copyright (C) 2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/oratio.h>

#include <stdio.h>

#include "test_m_unit_sphere.h"

/**
 * @param r is an initialized uniform random numbers generator.
 */
void test_m_unit_sphere(Random *r)
{
  int i, n;
  double *x;

  /* Dimension of the unit_sphere */
  n = 3;

  x = _mallocs(n, double);

  /* repeat the draws */
  for (i = 0; i < 10; i++) {
    ran_m_unit_sphere(r, n, x);
    nt_vect_print(stdout, "x", x, n);
  }

  free(x);
}
