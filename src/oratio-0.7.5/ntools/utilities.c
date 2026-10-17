/**
 * utilities.c
 *
 * ORATIO Library
 * Copyright (C) 2002-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <math.h>
#include <stdlib.h>

/**
 * Compute the sum of elements of an array of doubles, using Kahan's method.
 *
 * @param n is the number of terms in the sum.
 * @param x is the vector of observations.
 * @param incx is the increment between observations.
 *
 * @return the sum of elements in x (taking account of increments between
 *   them).
 */
double nt_kahan_sum(int n, double *x, int incx)
{
  /* W. Kahan, "Pracniques: further remarks on reducing truncation errors",
     Communications of the ACM 8(1), p. 40, 1965. */
  double sum = 0.0, y, t;
  double c = 0.0;   /* a running componesation for lost low-order bits */
  int i;

  for (i = 0; i < n; i++) {
    y = x[i*incx]-c;
    t = sum+y;
    c = (t-sum)-y;
    sum = t;
  }

  return sum;
}

/**
 * Compute the index of the maximum element in a vector.
 *
 * @param n is the vector dimension.
 * @param x is the vector to consider.
 * @param incx is the increment in x.
 *
 * @return the index of the greatest element in x.
 */
int nt_max_idx(int n, double *x, int incx)
{
  int max = 0;
  int i;

  for (i = 1; i < n; i++) if (x[i] > x[max]) max = i;

  return max;
}

/**
 * Compute the index of the minimum element in a vector.
 *
 * @param n is the vector dimension.
 * @param x is the vector to consider.
 * @param incx is the increment in x.
 *
 * @return the index of the smallest element in x.
 */
int nt_min_idx(int n, double *x, int incx)
{
  int min = 0;
  int i;

  for (i = 1; i < n; i++) if (x[i] < x[min]) min = i;

  return min;
}

/**
 * Compute consecutive powers of 2, starting from 2^0 = 1.
 *
 * @param n is the number of powers to compute (including 0).
 * @param p is an allocated integer array. On output, this array
 *   contains the successive powers of 2.
 */
void nt_two_powers_array(int n, unsigned int *p)
{
  int i;

  p[0] = 1;
  for (i = 1; i < n; i++) p[i] = p[i-1] << 1;
}

/**
 * Compute the closest power of two to some integer x, but less or equal to x.
 *
 * @param x is the considered integer.
 *
 * @return i, where 2^i <= x, and for all j > i, 2^j > x.
 */
int nt_two_power_floor(int x)
{
  int i = 0;

  while (x /= 2) i++;

  return i;
}

/**
 * Compute the closest power of two to some integer x, but greater or
 * equal to x.
 *
 * @param x is the considered integer.
 *
 * @return i, where 2^i >= x, and for all j < i, 2^j < x.
 */
int nt_two_power_ceil(int x)
{
  int i = 0;

  if (!(x % 2)) i = -1;
  while (x /= 2) i++;

  return i+1;
}

/**
 * Remove the last figure of an integer.
 *
 * @param number is the integer to consider.
 *
 * @return the integer without the last figure.
 */
int nt_remove_last_figure(int number)
{
  div_t result;

  result = div(number, 10);

  return result.quot;
}

/**
 * General purpose routines
 */

int nt_inc_case(int *array, int len, int base)
{
  int i = 0;

  for (;;) {
    (array[i])++;
    if (array[i] != base) break;
    if (i == len-1) {
      array[i] = base-1;
      return 0;
    }
    array[i] = 0;
    i++;
  }

  return 1;
}

/**
 * Truncature functions
 */

/**
 * Tuncate a double precision real into an integer, and return it in
 * double precision.
 *
 * @param a is the real that will be truncated
 */
double nt_trunc_double(double a)
{
  return (double)((int)(a));
}

/**
 * Truncate a double into a long integer
 *
 * @param a is the real that will be truncated
 */
long nt_trunc_long(double a)
{
  if (a < 1.0) return (long)0;
  else return (long)a;
}

/**
 * Compute the machine precision.
 *
 * This subroutine is intended to determine the parameters of the
 * floating-point arithmetic system specified below.  The
 * determination of the first three uses an extension of an algorithm
 * due to M. Malcolm, CACM 15 (1972), pp. 949-951, incorporating some,
 * but not all, of the improvements suggested by M. Gentleman and S.
 * Marovich, CACM 17 (1974), pp. 276-277.  An earlier version of this
 * program was published in the book Software Manual for the
 * Elementary Functions by W. J. Cody and W. Waite, Prentice-Hall,
 * Englewood Cliffs, NJ, 1980.  The present program is a
 * translation of the Fortran 77 program in W. J. Cody, "MACHAR:
 * A subroutine to dynamically determine machine parameters".
 * TOMS (14), 1988.
 *
 * Parameter values reported are as follows:
 * 
 * ibeta   - the radix for the floating-point representation
 * it      - the number of base ibeta digits in the floating-point
 *   significand
 * irnd    - 0 if floating-point addition chops
 *   1 if floating-point addition rounds, but not in the IEEE style
 *   2 if floating-point addition rounds in the IEEE style
 *   3 if floating-point addition chops, and there is partial underflow
 *   4 if floating-point addition rounds, but not in the IEEE style,
 *     and there is partial underflow
 *   5 if floating-point addition rounds in the IEEE style, and there is
 *     partial underflow
 * ngrd    - the number of guard digits for multiplication with
 *   truncating arithmetic.  It is 0 if floating-point arithmetic
 *   rounds, or if it truncates and only  it  base  ibeta digits
 *   participate in the post-normalization shift of the floating-point
 *   significand in multiplication; 1 if floating-point arithmetic
 *   truncates and more than  it  base  ibeta  digits participate in
 *   the post-normalization shift of the floating-point significand in
 *   multiplication
 * machep  - the largest negative integer such that
 *   1.0+FLOAT(ibeta)**machep .NE. 1.0, except that machep is bounded
 *   below by  -(it+3) 
 * negeps  - the largest negative integer such that
 *   1.0-FLOAT(ibeta)**negeps .NE. 1.0, except that
 *   negeps is bounded below by  -(it+3)
 * iexp    - the number of bits (decimal places if ibeta = 10)
 *   reserved for the representation of the exponent (including the
 *   bias or sign) of a floating-point number
 * minexp  - the largest in magnitude negative integer such that
 *   FLOAT(ibeta)**minexp is positive and normalized
 * maxexp  - the smallest positive power of  BETA  that overflows
 * eps     - the smallest positive floating-point number such that
 *   1.0+eps .NE. 1.0. In particular, if either ibeta = 2  or  IRND =
 *   0, eps = FLOAT(ibeta)**machep. Otherwise,  eps =
 *   (FLOAT(ibeta)**machep)/2
 * epsneg  - A small positive floating-point number such that
 *   1.0-epsneg .NE. 1.0. In particular, if ibeta = 2
 *   or  IRND = 0, epsneg = FLOAT(ibeta)**negeps.
 *   Otherwise,  epsneg = (ibeta**negeps)/2.  Because negeps is
 *   bounded below by -(it+3), epsneg may not be the smallest number
 *   that can alter 1.0 by subtraction.
 * xmin    - the smallest non-vanishing normalized floating-point
 *   power of the radix, i.e.,  xmin = FLOAT(ibeta)**minexp
 * xmax    - the largest finite floating-point number.  In particular
 *   xmax = (1.0-epsneg)*FLOAT(ibeta)**maxexp
 * Note - on some machines  xmax  will be only the second, or perhaps
 * third, largest number, being too small by 1 or 2 units in the last
 * digit of the significand.
 *
 * Latest revision - August 4, 1988
 *
 * Initial author - W. J. Cody
 * Argonne National Laboratory
 *
 * @return the machine precision, conditionally to no compiler code
 *   optimisation.
 */
double nt_macheps()
{
  double t = 1.0;
  double eps = 1.0;

  do {
    /* This does not work if any optimization code if performed */
    eps *= 0.5;
    t = 1.0 + eps;
    t = t-1.0;
    t = t+1.0;
  } while (t-1.0 != 0.0);
  eps *= 2;

  return eps;
}

/**
 * Compute the radical inverse of integer i in base b.
 *
 * @param b is the base.
 * @param i is the integer to develop.
 *
 * @return the radical inverse of i in base b.
 */
double nt_radical_inverse(int b, int i)
{
  double digit, radical, inverse;

  digit = radical = 1.0/(double)b;

  inverse = 0.0;
  for (; i > 0; i /= b) {
    inverse += digit*(double)(i%b);
    digit *= radical;
  }

  return inverse;
}
