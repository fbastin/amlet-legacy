/**
 * vector.c
*
 * ORATIO Library
 * Copyright (C) 2002-2005 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <math.h>
#include <string.h>

#include <oratio/error.h>
#include "oratio/nt_utilities.h"

/**
 * Set all vector components to zero.
 *
 * @param n is the vector dimension.
 * @param x is a real vector. On output, all components of x are equal to 0.
 * @param incx is the vector increment.
 */
void nt_vect_zeros(int n, double *x, int incx)
{
  int i;

  /* The use of memset seems ok for IEEE represenation, but I have to
     implement a way to test it. Perhaps in the configure.in, where I
     could compare (int)0 to FP_ZERO constant */
  /*  memset(x, 0, n*sizeof(double)); */
  for (i = 0; i < n; i++) x[i*incx] = 0.0;
}

/**
 * Set all vector components to one.
 *
 * @param n is the vector dimension.
 * @param x is a real vector. On output, all components of x are equal to 1.
 * @param incx is the vector increment.
 */
void nt_vect_ones(int n, double *x, int incx)
{
  int i;

  for (i = 0; i < n; i++) x[i*incx] = 1.0;
}

/**
 * Insert a vector y into a vector x
 *
 * @param x is a vector of dimension nx
 * @param y is a vector of dimension ny
 * @param a is the position where the vector y is to be inserted; the
 *   indexes are supposed to start from 0
 */
int nt_vect_insert_vect(double x[], int nx, double y[], int ny, int a)
{
  if (a + ny > nx) return NT_VECTOR_BAD_DIM;

  memcpy(&x[a], y, ny*sizeof(double));

  return 0;
}

/**
 * Fill selected elements of the vector x with the content of y.
 */
void nt_vect_select(double x[], double y[], int *selector, int n)
{
  int i;

  for (i = 0; i < n; i++) x[selector[i]] = y[i];
}

/**
 * Artithmetic operations
 */

double nt_vect_dist(double x[], double y[], int n)
{
  int i;
  double v = 0.0;

  for (i = 0; i < n; i++) v += (x[i]-y[i])*(x[i]-y[i]);

  return sqrt(v);
}

/**
 * Hadamard product
 *
 * @return y = alpha*x.*y
 */
void nt_vect_dhad(int n, double alpha,
		  double x[], int incx, double y[], int incy)
{
  int i;
  int posx = 0, posy = 0;

  for (i = 0; i < n; i++) {
    *(y+posy) = alpha*(*(x+posx))*(*(y+posy));
    posx += incx;
    posy += incy;
  }
}

/**
 * Inverse Hadamard product
 *
 * @return y = alpha*x^(-1).*y
 */
void nt_vect_dhad_inv(int n, double alpha,
		      double x[], int incx, double y[], int incy)
{
  int i;
  int posx = 0, posy = 0;

  for (i = 0; i < n; i++) {
    *(y+posy) = alpha*(*(y+posy))/(*(x+posx));
    posx += incx;
    posy += incy;
  }
}

/**
 * Comparison functions
 */

int nt_vect_is_zero(double x[], int n)
{
  int i;

  for (i = 0; i < n; i++)
    if (x[i] != 0.0) return 0;

  return 1;
}

int nt_vect_is_one(double x[], int n)
{
  int i;

  for (i = 0; i < n; i++)
    if (x[i] != 1.0) return 0;

  return 1;
}
