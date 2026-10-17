/**
 * test_gss.c
 *
 * OPHELIA Library 
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the BSD License
 */

#include <ophelia/gss.h>
#include <ophelia/collection.h>

#include <stdio.h>

/**
 * Test of the Golden Search algorithm.
 *
 * @return 0.
 */
int op_gss_test()
{
  double x = 0.0;

  double val = op_gss(op_parabola_1, &x, NULL, -10.0, 10.0, 1e-10, 0);

  printf("Application of Golden Search to x²+4x-10: %g, found at x=%g.\n",
	 val, x);

  return 0;
}
