/**
 * natural_test.c
 *
 * ORATIO Library
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <oratio/natural.h>

#include <stdio.h>

/**
 * Tests on functions on natural numbers.
 */
int nt_natural_test()
{
  int a, b;

  a = 5;
  b = 2;
  printf("GCD of %d %d: %d\n", a, b, nt_gcd(a, b));

  a = 50;
  b = 20;
  printf("GCD of %d %d: %d\n", a, b, nt_gcd(a, b));

  a = 49;
  b = 7;
  printf("GCD of %d %d: %d\n", a, b, nt_gcd(a, b));

  return 0;
}
