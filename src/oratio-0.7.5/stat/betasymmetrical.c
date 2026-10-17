/**
 * betasymmetrical.c
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Symmetrical Beta distribution
 */

#include "oratio/constants.h"
#include "oratio/ntools.h"
#include "oratio/stat.h"

#include <float.h>
#include <math.h>

/**
 * Code inspired by SSJ (Stochastic Simulation in Java), Pierre L'Ecuyer and
 * Richard Simar.
 */

#define EPSBETA    1.0e-10  /* < 0.75*sqrt(DBL_EPSILON) */
#define EPSSINGLE  1.0e-5
#define MAXI       11       /* maximum number of Newton iterations */
#define MAXIB      50       /* Max number of bisection iterations */
#define MAXJ       2000     /* maximum number of terms in series */

// EPSARRAY[j]: Epsilon required for j decimal degits of precision
static const double EPSARRAY[] = {
  0.5, 0.5E-1, 0.5E-2, 0.5E-3, 0.5E-4, 0.5E-5, 0.5E-6, 0.5E-7, 0.5E-8,
  0.5E-9, 0.5E-10, 0.5E-11, 0.5E-12, 0.5E-13, 0.5E-14, 0.5E-15, 0.5E-16,
  0.5E-17, 0.5E-18, 0.5E-19, 0.5E-20, 0.5E-21, 0.5E-22, 0.5E-23, 0.5E-24,
  0.5E-25, 0.5E-26, 0.5E-27, 0.5E-28, 0.5E-29, 0.5E-30, 0.5E-31, 0.5E-32,
  0.5E-33, 0.5E-34, 0.5E-35 };

/**
 * Compute the mean of a symmetrical beta distribution of parameter alpha.
 *
 * @param alpha is the (shape) parameter of the symmetrical beta
 *   distribution; alpha must be positive.
 *
 * @return the mean of the distribution.
 */
double st_betasymmetrical_mean(double alpha)
{
  if (alpha <= 0) return nt_nan();

  return 0.5;
}

/**
 * Function   (1 - x*x + 2*x*log(x)) / ((1 - x)*(1 - x))
 *
 * @param x is the point where to evaluate the function.
 *
 * @return the function value at x.
 */
double belog(double x)
{
  if (x > 1.0)
    return -belog(1.0/x);
  if (x < 1.0e-200)
    return 1.0;
  if (x < 0.9)
    return (1.0 - x*x + 2.0*x*log(x)) / ((1.0 - x)*(1.0 - x));
  if (x == 1.0)
    return 0.0;

  /* For x near 1, use a series expansion to avoid loss of precision. */
  double term;
  double EPS = 1.0e-12;
  double Y = 1.0 - x;
  double ypow = 1.0;
  double sum = 0.0;
  int j = 2;
  do {
    ypow *= Y;
    term = ypow / (j * (j + 1));
    sum += term;
    j++;
  } while (fabs(term / sum) > EPS);

  return 2.0 * sum;
}

/**
 * Compute the inverse cumulative distribution function of a symmetrical
 * beta distribution using the Peizer and Pratt approximation
 *
 * @param u is the point where the icdf has to be evaluated; u must be in
 *   [0,1].
 * @param alpha is the (shape) parameter of the symmetrical beta
 *   distribution; alpha must be positive.
 */
double peizer_pratt_icdf(double u, double alpha)
{
  /* Inverse of the normal approximation of D. B. Peizer and J. W. Pratt,
   * Journal of the American Statistical Association 63(324), pp. 1416-1456,
   * 1968. */

  double t1, t3, xprev;
  double C2 = alpha - 1.0 / 3.0 + 0.025/alpha;
  double z = st_gaussian_icdf(u);
  double x = 0.5;
  double y = 1.0 - x;
  int i = 0;

  do {
    i++;
    t1 = (2.0*alpha - 5.0/6.0)*x*y;
    t3 = 1.0 - y*belog(2.0*x) - x*belog(2.0*y);
    xprev = x;
    x = 0.5 + 0.5*z*sqrt(t1/t3)/C2;
    y = 1.0 - x;
  } while (i <= MAXI && fabs(x - xprev) > EPSBETA);

  return x;
}

/*
 * This method is used for alpha > 1 and u very close to 0. It will almost
 * never be called, if at all.
 *
 * @param alpha: shape parameter
 * @param logBua: ln(alpha * u * Beta(alpha, alpha))
 * @param a: lower bound on x
 * @param b: upper bound on x
 * @param d: digits of precision
 */
double bisect(double alpha, double logBua, double a, double b, int d)
{
  int i, j;
  double z, sum, term;
  double x, xprev;
  double EPSILON = EPSARRAY[d];

  if (a >= 0.5 || a > b) {
    a = 0.0;
    b = 0.5;
  }

  x = 0.5 * (a + b);
  i = 0;
  do {
    ++i;
    z = -x / (1 - x);

    /* Compute the series for F(x) */
    sum = term = 1.0;
    j = 1;
    do {
      term *= z * (j - alpha) / (j + alpha);
      sum += term;
      ++j;
    } while ((fabs (term/sum) >  EPSILON) && (j < MAXJ));
    sum *= x;

    /* Bisection method */
    term = log(x * (1.0 - x));
    z = logBua - (alpha - 1.0) * term;
    if (z > log(sum))
      a = x;
    else
      b = x;
    xprev = x;
    x = 0.5 * (a + b);
  } while ((fabs(xprev - x) > EPSILON) && (i < MAXIB));

  return x;
}

/**
 * Method used for alpha < 1 and x close to 0.
 *
 * @param alpha is the (shape) parameter.
 * @param bu: u * Beta(alpha, alpha).
 * @param d is the number of digits of precision.
 */
double inverse1(double alpha, double bu, int d)
{
  int i, j;
  double x, xnew, poc, sum, term;
  double EPSILON = EPSARRAY[d];

  // First term of series
  x = pow(bu * alpha, 1.0 / alpha);

  // If T1/T0 is very small, neglect all terms of series except T0
  term = alpha * (1.0 - alpha) * x / (1.0 + alpha);
  if (term < EPSILON)
    return x;

  x = bu * alpha / (1.0 + term);
  xnew = pow (x, 1.0/alpha);  // Starting point of Newton's iterates

  i = 0;
  do {
    ++i;
    x = xnew;

    /* Compute the series for F(x) */
    poc = 1.0;
    sum = 1.0 / alpha;
    j = 1;
    do {
      poc *= x * (j - alpha) / j;
      term = poc / (j + alpha);
      sum += term;
      ++j;
    } while ((term > sum * EPSILON) && (j < MAXJ));
    sum *= pow(x, alpha);

    /* Newton's method */
    term = (sum - bu)*pow(x*(1.0 - x), 1.0 - alpha);
    xnew = x - term;
  } while ((fabs(term) > EPSILON) && (i <= MAXI));

  return xnew;
}

/**
 * Method is used for alpha < 1 and x close to 1/2.
 *
 * @param alpha is the (shape) parameter.
 * @param w: (0.5 - u)B/pow(4, 1 - alpha)
 * @param d is the number of digits of precision.
 */
double inverse2(double alpha, double w, int d)
{
  int i, j;
  double term, y, ynew, z, sum;
  double poc;
  double EPSILON = EPSARRAY[d];

  term = (1.0 - alpha) * w * w * 4.0 / 3.0;
  /* If T1/T0 is very small, neglect all terms of series except T0 */
  if (term < EPSILON)
    return 0.5 - w;

  ynew = w / (1 + term);     /* Starting point of Newton's iterates */
  i = 0;
  do {
    ++i;
    y = ynew;
    z = 4.0 * y * y;

    // Compute the series for G(y)
    poc = sum = 1.0;
    j = 1;
    do {
      poc *= z * (j - alpha) / j;
      term = poc / (2 * j + 1);
      sum += term;
      ++j;
    } while ((term > sum * EPSILON) && (j < MAXJ));

    sum *= y;

    // Newton's method
    ynew = y - (sum - w)*pow(1.0 - z, 1.0 - alpha);
  } while ((fabs(ynew - y) > EPSILON) && (i <= MAXI));

  return 0.5 - ynew;
}

/**
 * Method is used for alpha > 1 and x close to 0.
 *
 * @param alpha is the (shape) parameter.
 * @param logBua: ln(alpha * u * Beta(alpha, alpha))
 * @param d is the number of digits of precision.
 */
double inverse3(double alpha, double logBua, int d)
{
  int i, j;
  double z, x, w, xnew, sum = 0., term;
  double eps = EPSSINGLE;
  double EPSILON = EPSARRAY[d];
  // For alpha <= 100000 and u < 1.0/(2.5 + 2.25*sqrt(alpha)), X0 is always
  // to the right of the solution, so Newton is certain to converge.
  double X0 = 0.497;

  /* Compute starting point of Newton's iterates */
  w = logBua / alpha;
  x = exp(w);
  term = (log1p(-x) + logBua) / alpha;
  z = exp (term);
  if (z >= 0.25)
    xnew = X0;
  else if (z > 1.0e-6)
    xnew = (1.0 - sqrt(1.0 - 4.0*z)) / 2.0;
  else
    xnew = z;

  i = 0;
  do {
    ++i;
    if (xnew >= 0.5)
      xnew = X0;
    x = xnew;

    sum = log(x * (1.0 - x));
    w = logBua - (alpha - 1.0) * sum;
    if (fabs(w) >= DBL_MAX_EXP*M_LN2) {
      xnew = X0;
      continue;
    }
    w = exp(w);
    z = -x/(1 - x);

    /* Compute the series for F(x) */
    sum = term = 1.0;
    j = 1;
    do {
      term *= z * (j - alpha) / (j + alpha);
      sum += term;
      ++j;
    } while ((fabs(term/sum) > eps) && (j < MAXJ));
    sum *= x;

    /* Newton's method */
    term = (sum - w) / alpha;
    xnew = x - term;
    if (fabs(term) < 32.0*EPSSINGLE)
      eps = EPSILON;
  } while ((fabs(xnew - x) > sum * EPSILON) &&
	   (fabs(xnew - x) > EPSILON) && (i <= MAXI));

  /* If Newton has not converged with enough precision, call bisection
     method. It is very slow, but will be called very rarely. */
  if (i >= MAXI && fabs(xnew - x) > 10.0 * EPSILON)
    return bisect (alpha, logBua, 0.0, 0.5, d);
  return xnew;
}

/**
 * Method is used for alpha > 1 and x close to 1/2.
 *
 * @param alpha is the (shape) parameter.
 * @param logBva: ln(B) + ln(1/2 - u) + (alpha - 1)*ln(4)
 * @param d is the number of digits of precision.
 */
double inverse4(double alpha, double logBva, int d)
{
  int i, j;
  double term, sum, y, ynew, z;
  double eps = EPSSINGLE;
  double EPSILON = EPSARRAY[d];

  ynew = exp(logBva);    // Starting point of Newton's iterates
  i = 0;
  do {
    ++i;
    y = ynew;

    /* Compute the series for G(y) */
    z = 4.0 * y * y;
    term = sum = 1.0;
    j = 1;
    do {
      term *= z * (j + alpha - 0.5) / (0.5 + j);
      sum += term;
      ++j;
    } while ((term > sum * eps) && (j < MAXJ));
    sum *= y * (1.0 - z);

    /* Newton's method */
    term = log1p (-z);
    term = sum - exp(logBva - (alpha - 1.0) * term);
    ynew = y - term;
    if (fabs(term) < 32.0*EPSSINGLE)
      eps = EPSILON;
  } while ((fabs (ynew - y) > EPSILON) &&
	   (fabs (ynew - y) > sum*EPSILON) && (i <= MAXI));

  return 0.5 - ynew;
}

/* Compute Beta(alpha, alpha) or Beta(alpha, alpha)*4^(alpha-1).\
 *
 * @param is a vector of at least three elements.
 */
void calc_b4(double alpha, double *bc, double epsilon)
{
  double temp;
  double pB = 0.0;
  double plogB = 0.0;
  double plogC = 0.0;

  if (alpha <= EPSBETA) {
    /* For a -> 0, B(a,a) = (2/a)*(1 - 1.645*a^2 + O(a^3)) */
    pB = 2.0 / alpha;
    plogB = log(pB);
    plogC = plogB + (alpha - 1.0)*M_LN4;
  } else if (alpha <= 1.0) {
    temp = lgamma(alpha);
    plogB = 2.0 * temp - lgamma(2.0*alpha);
    plogC = plogB + (alpha - 1.0)*M_LN4;
    pB = exp(plogB);
  } else if (alpha <= 10.0) {
    plogC = lgamma(alpha) - lgamma(0.5 + alpha) + M_LN_SQPI_2;
    plogB = plogC - (alpha - 1.0)*M_LN4;
    pB = exp(plogB);
  } else if (alpha <= 200.0) {
    /* Convergent series for Gamma(x + 0.5) / Gamma(x) */
    double term = 1.0;
    double sum = 1.0;
    int i = 1;
    while (term > epsilon*sum) {
      term *= (i - 1.5)*(i - 1.5) /(i*(alpha + i - 1.5));
      sum += term;
      i++;
    }
    temp = M_SQPI_2/sqrt((alpha - 0.5)*sum);
    plogC = log(temp);
    plogB = plogC - (alpha - 1.0)*M_LN4;
    pB = exp(plogB);
  } else {
    /* Asymptotic series for Gamma(a + 0.5) / (Gamma(a) * Sqrt(a)) */
    double z = 1.0 / (8.0*alpha);
    temp = 1.0 + z*(-1.0 + z*(0.5 + z*(2.5 - z*(2.625 + 49.875*z))));
    /* This is 4^(alpha - 1)*B(alpha, alpha) */
    temp = M_SQPI_2/(sqrt(alpha)*temp);
    plogC = log(temp);
    plogB = plogC - (alpha - 1.0)*M_LN4;
    pB = exp(plogB);
  }
  bc[0] = pB;
  bc[1] = plogB;
  bc[2] = plogC;
}

/**
 * Compute the inverse cumulative distribution function of a symmetrical
 * beta distribution.
 *
 * @param u is the point where the icdf has to be evaluated; u must be in
 *   [0,1].
 * @param alpha is the (shape) parameter of the symmetrical beta
 *   distribution; alpha must be positive.
 * @param d is the number of digits of precision.
 */
double compute_icdf(double u, double alpha, int d,
		    double logFact, double logBeta, double logCeta, double Ceta) {
  if (alpha <= 0.0 || u > 1.0 || u < 0.0) return nt_nan();
  if (u == 0.0) return 0.0;
  if (u == 1.0) return 1.0;
  if (u == 0.5) return 0.5;

  /* Case alpha = 1 is the uniform law */
  if (alpha == 1.0) return u;
  
  /* Case alpha = 1/2 is the arcsin law */
  if (alpha == 0.5) {
    double temp = sin(u*M_PI_2);
    return temp*temp;
  }

  if (alpha > 100000)
    return peizer_pratt_icdf(u, alpha);

  int up;  /* 1 is u > 0.5, 0 else */
  if (u > 0.5) {
    up = 1;
    u = 1.0 - u;
  } else
    up = 0;

  double x;
  double C = 0.0, B = 0.0, logB = 0.0, logC = 0.0;
  
  if (logFact == DBL_MIN) {
    double bc[] = {0.0, 0.0, 0.0};
    calc_b4(alpha, bc, EPSARRAY[d]);
    B = bc[0]; logB = bc[1]; logC = bc[2];
    C = exp(logC);
  } else {
    B = 1.0/exp(logFact);
    logB = logBeta;
    logC = logCeta;
    C = Ceta;
  }
  
  if (alpha <= 1.0) {
    // First term of integrated series around 1/2
    double y0 = C * (0.5 - u);
    if (y0 > 0.25)
      x = inverse1 (alpha, B * u, d);
    else
      x = inverse2 (alpha, y0, d);
    
  } else {
    if (u < 1.0 / (2.5 + 2.25*sqrt(alpha))) {
      double logBua = logB + log(u * alpha);
      x = inverse3 (alpha, logBua, d);
    } else {
      // logBva = Ln(Beta(a,a) * (0.5 - u)*pow(4, a - 1)
      double logBva = logC - M_LN2 + log1p(-2.0*u);
      
      x = inverse4(alpha, logBva, d);
    }
  }
  
  if (up)
    return 1.0 - x - DBL_EPSILON;
  else
    return x;
}

/**
 * Compute the inverse cumulative distribution function of a symmetrical
 * beta distribution of parameter alpha.
 *
 * @param u is the point where the icdf has to be evaluated; u must be in
 *   [0,1].
 * @param alpha is the (shape) parameter of the symmetrical beta
 *   distribution; alpha must be positive.
 *
 * @return the value of the icdf at u.
 */
double st_betasymmetrical_icdf(double u, double alpha)
{
  double t;
  t = compute_icdf(u, alpha, 14, DBL_MIN, 0.0, 0.0, 0.0);
  return t;
}
