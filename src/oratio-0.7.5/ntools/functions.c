/**
 * functions.c
 *
 * ORATIO Library
 * Copyright (C) 2002 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/ctools.h"

#include <math.h>

#define ST_EPS (DBL_EPSILON*100)
#define FPMIN  (ST_EPS*ST_EPS)
#define ITMAX  10000  /* maximum number of iterations  */

/**
 * Evaluate the function B(z, w)
 */
double nt_beta(double z, double w)
{
  return exp(lgamma(z)+lgamma(w)-lgamma(z+w));
}


/**
 * Evaluates continued fraction for incomplete beta function
 *
 * Ref.: "Numerical Recipes in C", William H. Press, Saul A. Teukolsky,
 * William T. Vettering, Brian P. Flannery, 1992, p 227

 * @param x is the point where to evaluate the continued fraction; x must be
 *   in [0,1].
 */
double nt_beta_cf(double a, double b, double x)
{
  int m, m2;
  double aa, c, d, del, h, qab, qam, qap;
	
  qab = a+b;
  qap = a + 1.0;
  qam = a-1.0;
  c = 1.0;
  d=1.0-qab*x/qap;
  if(fabs(d) < FPMIN) d = FPMIN;
	
  d = 1.0 / d;
  h = d;

  for (m=1; m <= ITMAX; ++m) {
    m2 = 2*m;
    aa = m*(b-m) * x / ( (qam+m2) * (a+m2));
    d = 1.0 + aa *d;
    if(fabs(d) < FPMIN) d = FPMIN;

    c = 1.0 + aa / c;
    if (fabs(c) < FPMIN) c = FPMIN;
		
    d = 1.0 / d;
		
    h *= d * c;
		
    aa= -(a+m) * (qab+m) * x / ( (a+m2) * (qap+m2));

    d = 1.0 + aa * d;
    if (fabs(d) < FPMIN) d = FPMIN;

    c = 1.0 + aa / c;
    if (fabs(c) < FPMIN) c = FPMIN;
		
    d = 1.0 / d;
    del = d * c;
		
    h *= del;
	
    if (fabs(del-1.0) < ST_EPS) break;
		
  }

  if (m > ITMAX)
    c_warning("a or b too big, or MAXIT too small in betacf");

  return h;			
}

/**
 * Returns the incomplete beta function of parameters a, b, at x
 *
 * @param x is the point whete the incomplete beta function has to be
 * evaluated; x must be in [0,1].
 */
double nt_betai(double a, double b, double x)
{
  /* "Numerical Recipes in C", William H. Press, Saul A. Teukolsky,
   * William T. Vettering, Brian P. Flannery, 1992, p 227. */
  double bt;
	
  if(x == 0.0 || x == 1.0)
    bt = 0.0;
  else
    bt = exp(lgamma(a+b)-lgamma(a)-lgamma(b) +
	     a*log(x) + b*log(1.0-x));

  /* Use continued fraction */
  if( x < (a+1.0)/(a+b+2.0))
    return bt*nt_beta_cf(a,b,x)/a;
  else
    return 1.0-bt*nt_beta_cf(b,a,1.0-x)/b;
}

/**
 * Compute n!
 *
 * @param n is a positive integer (n >= 0).
 *
 * @return n! as a double.
 */
double nt_factorial(int n)
{
  static int ntop = 4;
  static double a[33] = {1.0, 1.0, 2.0, 6.0, 24.0};
  int j;

  /* If a larger value than the tabulated ones, we use the property
   * n! = Gamma(n+1). This usually procude a large value, meaning that
   * we can face an overflow, but we choose to try. */
  if (n > 32) return exp(lgamma(n+1.0));

  /* We complete the array a. */
  while (ntop < n) {
    j = ntop++;
    a[ntop] = a[j]*ntop;
  }

  return a[n];
}

/**
 * Compute ln(n!). This can be useful if n is high, since we obtain
 * then a better precision.
 *
 * @param n is a positive integer.
 *
 * @return len(n!) as a double.
 */
double nt_factorial_ln(int n)
{
  static float a[101];  /* a statically allocated array is always
			   initialized to 0. */

  if (n < 0) c_error("** NTOOLS error: negative factorial in routine factln");
  if (n <= 1) return 0.0;
  if (n <= 100)
    return a[n] ? a[n] : (a[n] = lgamma(n+1.0));
  else
    return lgamma(n+1.0);
}

double nt_binomial_coef(int n, int k)
{
  /* The function rint allows us to remove rounding errors for small
     values of n and k. */

  return rint(exp(nt_factorial_ln(n)-nt_factorial_ln(k)-
		  nt_factorial_ln(n-k)));
}

/**
 * Compute x^n, where x is a real and n is a positive integer.
 *
 * @param x is a real number.
 * @param y is an integer.
 *
 * @return x^n.
 */
double nt_ipow(double x, int y)
{
  double t = 1.0;

  if (y < 0) {
    x = 1.0/x;
    y = -y;
  }

  while (y > 0) {
    /* check if y is odd */
    if (y & 1) {
      t *= x;
      if (y == 1) return t;
    }
    x *= x;
    y /= 2;
  }

  return t;
}
