/**
 * utils.c
 *
 * STAT Library
 * Copyright (C) 2000-2002 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/ctools.h>

#include <math.h>
#include <float.h>

#define ST_EPS (DBL_EPSILON*100)
#define FPMIN  (ST_EPS*ST_EPS)
#define ITMAX  10000  /* maximum number of iterations  */

/**
 * Evaluates continued fraction for incomplete gamma function
 *
 * Ref.: "Numerical Recipes in C", William H. Press, Saul A. Teukolsky,
 * William T. Vettering, Brian P. Flannery, 1992, p 219.
 */
void st_gamma_cf(double *gammcf, double a, double x, double *gln)
{
  /* L'évaluation de la fraction continuée se fait par la méthode de Lenz
   * modifiée, avec b0 nul, et en utilisant la partie paire. */

  int i;
  double an, b, c, d, frac, delta;
  /* an et b reprendront les termes de la fraction continuée
   * c et d sont les facteurs intervenant dans l'algorithme de Lenz
   * frac contiendra la valeur de la fraction continue
   * delta contiendra le facteur par lequel doit être multiplié la fraction
   * entre deux itérations */

  *gln = lgamma(a);

  /* calcul du premier terme de la fraction continuée */
  b = x + 1.0 - a;
  c = 1.0/FPMIN;   /* c0 = 0; on met c0 petit (et a0 = 0, implicitement) */
  d = 1.0/b;
  frac = d;        /* frac = FPMIN+1/FPMIN*d */

  /* i = 1 correspond au second terme de la fraction continuée */
  for (i = 1; i <= ITMAX; i++) {
    an = -i*(i-a);
    b += 2.0;
    d = an*d+b;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = b+an/c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1.0/d;
    delta = d*c;
    frac *= delta;

    /* test fonctionnant généralement bien: on suppose que l'on peut s'arrêter
     * quand la différence entre la valeur calculée entre deux itérations est
     * suffisamment petite */
    if (fabs(delta-1.0) < ST_EPS) break;
  }

  if (i > ITMAX) c_warning("a too large, ITMAX too small in st_gamma_cf");

  /* multiplie la fraction continuée par les facteurs communs */
  *gammcf = exp(-x+a*log(x)-(*gln))*frac;
}

/**
 * Compute the incomplete function Gamma P(a,x), using the serial
 * representation.
 *
 * @return gamser = P(a,x)
 *         gln = ln Gamma(a)
 */
void st_gamma_ser(double *gamser, double a, double x, double *gln)
{
  double sum, term, ap;
  int n;

  *gln = lgamma(a);
  if (x <= 0.0) {
    if (x < 0.0) c_warning("x less than 0 in routine st_gamma_ser");
    *gamser = 0.0;
    return;
  }
  else {
    ap = a;
    sum = term = 1.0/a;

    for (n = 1; n <= ITMAX; n++) {

      ++ap;
      term *= x/ap;
      sum += term;

      if (fabs(term) < fabs(sum)*ST_EPS) {
	*gamser = sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    c_warning("a too large, ITMAX too small in routine st_gamma_ser");
    return;
  }
}

double st_gamma_p_gln(double a, double x, double *gln)
{
  double gamma;

  if (x < 0.0 || a <= 0.0)
    c_error("invalid arguments in routine st_gamma_p_gln");

  if (x < (a + 1.0)) {
    /* Use the serial representation. */
    st_gamma_ser(&gamma, a, x, gln);
    return gamma;
  }
  else {
    /* Utilise la représentation en fractions continuées */
    st_gamma_cf(&gamma, a, x, gln);
    return 1.0 - gamma;
  }
}

double st_gamma_q_gln(double a, double x, double *gln)
{
  double gamma;

  if (x < 0.0 || a <= 0.0)
    c_error("invalid arguments in routine st_gamma_q_gln");

  if (x < (a + 1.0)) {
    /* Utilise la représentation en série */
    st_gamma_ser(&gamma, a, x, gln);
    return 1.0 - gamma;
  }
  else {
    /* Utilise la représentation en fractions continuées */
    st_gamma_cf(&gamma, a, x, gln);
    return gamma;
  }
}

/**
 * Compute the incomplete gamma function P(a,x).
 */
double st_gamma_p(double a, double x)
{
  double tmp;

  return st_gamma_p_gln(a, x, &tmp);
}

double st_gamma_q(double a, double x)
{
  double tmp;

  return st_gamma_q_gln(a, x, &tmp);
}

/**
 * Error function.
 * This function should have its equivalent in the standard library under
 * the name erf(x).
 *
 * @param x is the point where the function has to be evaluated.
 *
 * @return the value of the error function at x.
 */
double st_erf(double x)
{
  /**
   * Calcule la fonction d'erreur erf(x) = sgn(x)*P(0.5, x²) (gamma incomplète)
   * Explicitement, 
   *   st_erf(x) = 2/sqrt(pi)* integral from 0 to x of exp(-t*t) dt
   */
  return x < 0.0 ? -st_gamma_p(0.5, x*x) : st_gamma_p(0.5, x*x);
}

double st_erfc(double x)
{
  return x < 0.0 ? 1.0 + st_gamma_p(0.5, x*x) : st_gamma_q(0.5, x*x);
}

/**
 * Variable transformation
 */

/**
 * Box-Cox transformation.
 *
 * Reference: George E. P. Box and David R. Cox, "An Analysis of
 *   Transformations", Journal of the Royal Statistical Society, 26(2),
 *   pp. 211--252, 1964.
 *
 * @param x is the variable to transform.
 * @param lambda is the transformation parameter.
 *
 * @return the Box-Cox transformation of x with the parameter lambda.
 */
double st_bct(double x, double lambda)
{
  if (lambda == 0.0) return log(x);

  return (pow(x, lambda) - 1)/lambda;
}
