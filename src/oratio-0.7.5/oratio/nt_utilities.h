/**
 * nt_utilities.h
 *
 * ORATIO Library
 * Copyright (C) 2002-2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_NTOOLS_UTILITIES__
#define __ORATIO_NTOOLS_UTILITIES__

#include <oratio/cblas.h>
#include <oratio/ctools.h>

/**
 * Constants
 */

/**
 * The following constants can be employed to specify if the matrix is
 * stocked by row of by column.
 * They mimics the CBLAS conventions.
 */
#define NT_MATRIX_ROW    CblasRowMajor
#define NT_MATRIX_COLUMN CblasColMajor

#define nt_min(a, b) (((a) < (b)) ? a : b)
#define nt_max(a, b) (((a) < (b)) ? b : a)

#define nt_sign(a) (((a) >= 0) ? 1 : -1)

/**
 * Functions type definition
 */

/**
 * Real functions
 */

typedef double (*NT_f)(double *, int, void *);
typedef void   (*NT_fval)(double *, int, void *, double *);
typedef double (*NT_fg)(double *, int, void *, double *);

/**
 * Real multi-functions
 */

/**
 * The parameters are x, size of x, arbitrary data, image vector
 * and its size.
 */
typedef double (*NT_mf)(double *, int, void *, double *, int);

/**
 * The parameters are x, size of x, arbitrary data, image vector
 * and its size, and the function Jacobian
 */
typedef double (*NT_mfJ)(double *, int, void *, double *, int, double *);

/**
 * Least-squares adapted functions
 */

/**
 * The parameters are x, size of x, arbitrary data, (mean) gradient,
 * the outer product Hessian approximation.
 */
typedef double (*NT_lsfgH)(double *, int, void *, double *, double **);

/**
 * The parameters are x, size of x, arbitrary data,
 * indivuals dependent variables, associated vector size,
 * (mean) gradient, outer product Hessian approximation.
 */
typedef double (*NT_lsfrgH)(double *, int, void *, double *, int,
			    double *, double **);
/**
 * The parameters are x, size of x, arbitrary data,
 * indivuals dependent variables, associated vector size,
 * (mean) gradient, function Jacobian,
 * outer product Hessian approximation.
 */
typedef double (*NT_lsfrgJH)(double *, int, void *, double *, int,
			     double *, double **, double **);

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

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
  double nt_kahan_sum(int n, double *x, int incx);

  /**
   * Incrémente le tableau array suivant la base base. Les nombres sont
   * réprésentés non signés, de la gauche vers la droite.
   * Renvoie 0 en cas d'overflow, 1 sinon.
   */
  int nt_inc_case(int *array, int len, int base);

  /**
   * Compute the index of the maximum element in a vector.
   *
   * @param n is the vector dimension.
   * @param x is the vector to consider.
   * @param incx is the increment in x.
   *
   * @return the index of the greatest element in x.
   */
  int nt_max_idx(int n, double *x, int incx);

  /**
   * Compute the index of the minimum element in a vector.
   *
   * @param n is the vector dimension.
   * @param x is the vector to consider.
   * @param incx is the increment in x.
   *
   * @return the index of the smallest element in x.
   */
  int nt_min_idx(int n, double *x, int incx);

  /**
   * Truncature functions
   */

  /**
   * Tuncate a double precision real into an integer, and return it in
   * double precision.
   *
   * @param a is the real that will be truncated
   */
  double nt_trunc_double(double a);

  /**
   * Truncate a double into a long integer
   *
   * @param a is the real that will be truncated
   */
  long nt_trunc_long(double a);

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
  void nt_machar(int *ibeta, int *it, int *irnd, int *ngrd,
		 int *machep, int *negep, int *iexp, int *minexp,
		 int *maxexp, double *eps, double *epsneg, double *xmin,
		 double *xmax);

  /**
   * Remove the last figure of an integer.
   *
   * @param number is the integer to consider.
   *
   * @return the integer without the last figure.
   */
  int nt_remove_last_figure(int number);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_NTOOLS_UTILITIES__ */
