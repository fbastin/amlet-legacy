/**
 * modular.c
 *
 * NTools Library
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/ntools.h>

#include "string.h"

/* If not zero, nt_modular_vecvec() reproduces ORATIO <= 0.7.5. */
static int nt_modular_legacy = 0;

void nt_modular_set_legacy(int legacy)
{
  nt_modular_legacy = legacy;
}

/**
 * Compute ax+y mod m.
 *
 * @param a is an integer satisfying 0 <= a*m < 2^53.
 * @param x is an integer.
 * @param y is an integer.
 * @param m is a strictly positive integer satisfying m < 2^35.
 *
 * @return ax+y mod m.
 */
double nt_modular_daxpy(double a, double x, double y, double m)
{
   double v;
   long a1;

   v = a*x+y;

   if ((v >= M_TWO53) || (v <= -M_TWO53)) {
      a1 = (long)(a/M_TWO17);
      a -= a1*M_TWO17;
      v = a1*x;
      a1 = (long)(v/m);
      v -= a1*m;
      v = v*M_TWO17 + a*x + y;
   }

   a1 = (long) (v/m);

   if ((v -= a1*m) < 0.0)
      return v += m;
   else
      return v;
}

/**
 * Recursive method that computes a^e mod m.
 *
 * @param a is a long integer, possibly greater than m.
 * @param e is a long integer.
 * @param m is a strictly positive long integer.
 *
 * @return a^e mod m.
 */
long nt_modular_power(long a, long e, long m)
{
  if (e == 0)
    return 1;
  else if (e == 1)
    return a % m;
  else if ((e & 1) == 1)
    return (a * nt_modular_power(a, e-1, m)) % m;
  else {
    long p = nt_modular_power(a, e/2, m);
    return (p*p) % m;
  }
}

/**
 * Compute ax mod m using Schrage's method.
 *
 * @param a is an integer such that a*a<m or a = floor(m/i), where i*i<m.
 * @param x is a positive integer between 0 and m-1 (included).
 * @param m is a positive integer.
 *
 * @return ax mod m.
 */
long nt_modular_schrage(long a, long x, long m)
{
  long q, r, y;

  q = m/a;
  r = m%a;

  y = x/q;
  x = a*(x-y*q)-y*r;

  if (x < 0) x += m;

  return x;
}

/**
 * Compute (a'x + y) mod m.
 *
 * @param n is the dimension of the vectors a and x.
 * @param a is a vector of dimension n.
 * @param x is a vector of dimension n.
 * @param y is a double with some initialized value on input.
 * @param m is a positive integer such that -m < x[i] < m.
 *
 * @return the value (a'x + y) mod m.
 */
double nt_modular_vecvec(int n, double *a, double *x,
			 double y, double m)
{
  int i;

  for (i = 0; i < n; i++)
    if (nt_modular_legacy)
      /* ORATIO <= 0.7.5 added the running sum again: the result is not
	 reduced modulo m. */
      y += nt_modular_daxpy(a[i], x[i], y, m);
    else
      y = nt_modular_daxpy(a[i], x[i], y, m);

  return y;
}

/**
 * Returns v = A*x mod m.
 *
 * @param A is a matrix of leading dimension incA, stored by rows. In
 *   other terms incA is the number of columns of A.
 * @param x is a vector of dimension n. Its content can be replaced on
 *   output if v is a pointer to the same vector.
 * @param v is a vector whose value on output will be equal to Ax mod m.
 *   On input is content is undefined. v can point to the same pointer
 *   than x.
 * @param m is an integer such that -m < x[i] < m.
 *
 * @return v = A*x mod m.
 */
void nt_modular_matvec(int n, double *A, int incA, double *x,
		       double *v, double m)
{
   int i;

   double t[n];

   for (i = 0; i < n; i++)
     t[i] = nt_modular_vecvec(n, A+incA*i, x, 0.0, m);
 
   cblas_dcopy(n, t, 1, v, 1);
}

/**
 * Compute C = A*B mod m.
 *
 * @param is the effective dimension of A, B, C, i.e. they are viewed
 *   as nxn matrices.
 * @param A is a matrix of dimension at least nxn.
 * @param incA is the leading dimension of A.
 * @param B is a matrix of dimension at least nxn.
 * @param incB is the leading dimension of A.
 * @param C is a matrix of dimension at least nxn.
 *   On output, C = A*B mod m. Work even if A = C or B = C or A = B = C.
 * @param incC is the leading dimension of A.
 * @param m is a strictly positive integer.
 */
void nt_modular_matmat(int n, double *A, int incA, double *B, int incB,
		       double *C, int incC, double m)
{
   int i, j;
   double v[n], W[n*n];

   for (i = 0; i < n; i++) {
     /* Copy the ith column of B. */
     for (j = 0; j < n; j++) v[j] = B[j*incB+i];

     nt_modular_matvec(n, A, n, v, v, m);

     /* Copy v in the ith column of W. */
     for (j = 0; j < 3; j++) W[j*n+i] = v[j];
   }

   /* Copy W in C. */
   for (i = 0; i < n; ++i) memcpy(C+i*incC, W+i*n, n*sizeof(double));
}

/**
 * Compute matrix B = A^p mod m.
 *
 * @param n is the dimensions of matrices A and B, which are assumed
 *   to be square.
 * @param A is a matrix of dimension nxn.
 * @param incA is the leading dimension of A.
 * @param B is a matrix of dimension nxn. B may point to A.
 * @param incB is the leading dimension of B.
 * @param m is a strictly positive integer.
 * @param p is a positive integer.
 *
 * @returns B = A^p mod m.
 */
void nt_modular_matpow(int n, double *A, int ldA,
		       double *B, int ldB, double m, long p)
{
   int i;
   double W[n*n];

   /* Initialize: W = A; B = I */
   for (i = 0; i < n; i++)
     memcpy(W+i*n, A+i*ldA, n*sizeof(double));
   nt_matrix_identity(n, n, B, ldB);

   /* Compute B = A^p % m using the binary decomposition of p */
   while (p > 0) {
     if (p % 2)
       nt_modular_matmat(n, W, n, B, n, B, n, m);
     nt_modular_matmat(n, W, n, W, n, W, n, m);
     p /= 2;
   }
}

/**
 * Compute matrix B = (A^(2^p) % m).
 *
 * @param n is the dimensions of matrices A and B, which are assumed
 *   to be square.
 * @param A is a matrix of dimension nxn.
 * @param incA is the leading dimension of A.
 * @param B is a matrix of dimension nxn. B may point to A.
 * @param incB is the leading dimension of B.
 * @param m is a strictly positive integer.
 * @param p is a positive integer.
 *
 * @returns B = A^(2^p) mod m.
 */
void nt_modular_mattwopow(int n, double *A, int incA,
			  double *B, int incB, double m, long p)
{
   int i;

   /* Initialize: B = A */
   if (A != B) {
     for (i = 0; i < n; i++)
       memcpy(B+incB*i, A+incA*i, n*sizeof(double));
   }

   /* Compute B = A^{2^p} mod m*/
   for (i = 0; i < p; i++) nt_modular_matmat(n, B, n, B, n, B, n, m);
}
