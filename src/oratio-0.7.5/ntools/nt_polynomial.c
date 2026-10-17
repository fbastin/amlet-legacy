/**
 * nt_polynomial.c
 *
 * NTools Library
 * Copyright (C) 2002-2010 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Polynomial functions
 */

/**
 * Evaluation of a polynomial of one dimension.
 *
 * @param x is the point where the polynomial has to be evaluated.
 * @param int d is the polynomial degree.
 * @param a is the coefficient vector, such tuch a[i*incx] is the factor
 *   of the term of degree i.
 * @param inca is the increment of a.
 *
 * @return the value of the polynomial at x.
 */
double nt_polynomial_eval(double x, int d, double *a, int inca)
{
  double val;

  val = a[d];
  while (d > 0) val = val*x+a[(--d)*inca];

  return val;
}

/**
 * @param x is the point where the polynomial has to be evaluated.
 */
void nt_poly_eval_der(double x, double *coeff, int n,
		      double *dpoly, int nd)
{
  int i, j, nnd;
  double fact = 1.0;

  for (i = 1; i <= nd; i++) dpoly[i] = 0.0;
  dpoly[0] = coeff[n];

  for (i = n-1; i >= 0; i--) {
    nnd = (nd < (n-i) ? nd : n-i);

    /* Evaluation des dérivées */
    for (j = nnd; j >= 1; j--)
      dpoly[j] = dpoly[j]*x+dpoly[j-1];
    /* Evaluation du polynôme */
    dpoly[0] = dpoly[0]*x+coeff[i];
  }

  /* Pour les dérivées d'ordre supérieur à 2, il reste encore à ajuster
   * les coefficients (de type factoriel).
   */
  for (i = 2; i <= nd; i++) {
    fact *= i;
    dpoly[i] *= fact;
  }
}

/*****************************************************************************
 * Détermination d'un polynôme de degré 1.
 *****************************************************************************
 * Note pour les polynômes de dimension supérieure, il est conseillé de
 * procédérer par résolution matricielle.
 *****************************************************************************/
void nt_poly_line(double *a, double *b, double x1, double x2, double y1,
		  double y2)
{
  *a = (y2-y1)/(x2-x1);
  *b = y1-*a*x1;
}
