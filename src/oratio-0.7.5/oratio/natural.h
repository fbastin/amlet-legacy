/**
 * natural.h
 *
 * ORATIO Library
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_NATURAL_H__
#define __ORATIO_NATURAL_H__

#ifdef __cplusplus
extern "C"{
#endif

/**
 * Compute the greatest common divisor of two natural numbers.
 *
 * @param a is a natural number.
 * @param b is a natural number.
 *
 * @return the greatest common divisor of a and b.
 */
int nt_gcd(int a, int b);

/**
 * Compute the greatest common divisor of two natural numbers.
 *
 * @param a is a natural number.
 * @param b is a natural number.
 *
 * @return the greatest common divisor of a and b.
 */
  long int nt_gcd_long(long int a, long int b);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_NATURAL_H__ */
