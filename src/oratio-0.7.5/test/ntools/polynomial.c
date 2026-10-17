/**
 * polynomial.c
 *
 * ORATIO Library
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "ntools_test.h"

#include <stdio.h>

int nt_polynomial_test()
{
  double a[] = {-1.0, 3.0, 2.0};

  printf("Value of 2x^2+3x-1, at x=2: %g\n", 
	 nt_polynomial_eval(2.0, 3, a, 1));

  return 0;
}
