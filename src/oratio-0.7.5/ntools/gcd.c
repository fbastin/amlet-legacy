/**
 * gcd.h
 *
 * ORATIO Library
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License
 */

/**
 * Compute the greatest common divisor of two natural numbers.
 *
 * @param a is a natural number.
 * @param b is a natural number.
 *
 * @return the greatest common divisor of a and b.
 */
int nt_gcd(int a, int b)
{
  if (!a) return b;
  while (b) {
    if (a > b)
      a -= b;
    else
      b -= a;
  }
  return a;
}

/**
 * Compute the greatest common divisor of two natural numbers.
 *
 * @param a is a natural number.
 * @param b is a natural number.
 *
 * @return the greatest common divisor of a and b.
 */
long int nt_gcd_long(long int a, long int b)
{
  if (!a) return b;
  while (b) {
    if (a > b)
      a -= b;
    else
      b -= a;
  }
  return a;
}
