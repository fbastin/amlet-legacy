/**
 * modular_test.c
 *
 * NTools Library
 * Copyright (C) 1999-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "ntools_test.h"

/**
 * Modular arithmetic test routine.
 *
 * @return 0 if the test has been successful, some error code otherwise.
 */
int nt_modular_check()
{
  double x;
  long a;

  a = nt_modular_schrage(9, 3, 10);
  x = nt_modular_daxpy(9, 3, 1, 10);

  printf("9*3 mod 10 = %ld\n9*3 + 1 mod 10 = %g\n", a, x);

  return (a == 7 && x == 8.0) ? 0 : 1;
}
