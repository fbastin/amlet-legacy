/**
 * student.c
 *
 * STAT Library
 * Copyright (C) 2001-2007 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/constants.h"
#include "oratio/ntools.h"
#include "oratio/st_betasymmetrical.h"
#include "oratio/st_cauchy.h"
#include "oratio/st_normal.h"

#include <math.h>

/**
 * Compute the Student probability distribution function with nu
 * degrees of freedom, at x.
 *
 * @param x is the point where we want to evaluate the PDF.
 * @param nu is the number of degrees of freedom; nu must be positive.
 *
 * @return the value of the pdf of a Student distribution with nu degrees of
 *   freedom, at x.
 */
double st_student_pdf(double x, double nu)
{
  double factor = nt_gamma_ratio_half(nu/2.0)/sqrt(nu*M_PI);

  return factor*pow(1.0 / (1.0 + x*x/nu), (nu + 1)/2.0);
}

/**
 * Compute the Student cumulative distribution function with nu
 * degrees of freedom, at x.
 *
 * @param x is the point where we want to evaluate the CDF.
 * @param nu is the number of degrees of freedom; nu must be positive.
 *
 * @return the value of the cdf of a Student distribution with nu degrees of
 *   freedom, at x.
 */
double st_student_cdf(double x, double nu)
{
  if (x <= 0.0) return nt_betai(nu/2, 0.5, nu/(x*x+nu))/2;
  return 1-nt_betai(nu/2, 0.5, nu/(x*x+nu))/2;
}

/**
 * Gaver-Kafadar normal approximation for the Student-t inverse cdf.
 *
 * @param u is the point where we want to evaluate the inverse CDF; u must be
 *   in [0,1].
 * @param nu is the number of degrees of freedom; nu must be positive.
 *
 * @return the value of the inverse cumulative distribution function of a
 *   Student with nu degrees of freedom, at u, using the Gaver-Kafadar
 *   approximation.
 */
double inverse_gaver_kafadar(double u, double nu)
{
  /* Ref: Donald P. Gavera and Karen Kafadarb, "A Retrievable Recipe for
   *   Inverse t", The American Statistician, 38(4), pp. 308--311, 1984. */

  double z = st_gaussian_icdf(u);
  double q = z/(nu - 1.0);
  double v = q*q*(nu - 1.5);
  double t = sqrt(nu*expm1(v));
  return (u >= 0.5) ? t : -t;
}

/**
 * Compute the Student inverse cumulative distribution function with nu
 * degrees of freedom, at u.
 *
 * @param u is the point where we want to evaluate the inverse CDF; u must be
 *   in [0,1].
 * @param nu is the number of degrees of freedom; nu must be positive.
 *
 * @return the value of the inverse cumulative distribution function of a
 *   Student with nu degrees of freedom, at u.
 */
double st_student_icdf(double u, double nu)
{
  if (nu < 1 || u < 0.0 || u > 1.0) return nt_nan();
  if (u == 0.0) return -nt_infinity();
  if (u == 1.0) return nt_infinity();

  if (nu == 1) return st_cauchy_icdf(u, 0, 1);

  if (nu == 2) return (2.0*u - 1.0)/sqrt(2.0*u*(1.0 - u));

  if (nu > 100000) return inverse_gaver_kafadar(u, nu);

  double z = st_betasymmetrical_icdf(u, 0.5*nu);
  return (z - 0.5)*sqrt(nu/(z*(1.0 - z)));
}
