/**
 * nt_polynomial.h
 *
 * ORATIO Library - NTools module
 *
 * Copyright (C) 2002-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_NTOOLS_POLYNOMIAL_H__
#define __ORATIO_NTOOLS_POLYNOMIAL_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Polynomial functions
   */

  /**
   * Evaluation of a polynomial of one dimension.
   *
   * @param x is the point where the polynomial has to be evaluated.
   * @param int n is the polynomial degree.
   * @param a is the coefficient vector, such tuch a[i*incx] is the factor
   *   of the term of degree i.
   * @param inca is the increment of a.
   *
   * @return the value of the polynomial at x.
   */
  double nt_polynomial_eval(double x, int n, double *a, int inca);

  /***************************************************************************
   * Soit le polynôme de degré n, dont les coefficients sont stockés dans le
   * vecteur coeff, coeff[i] correspondant au coefficient de la puissance i
   * de x. La routine suivante évalue le polynôme en x dans dpoly[0], et les
   * nd dérivées, comme dpoly[1..nd].
   * @param dpoly est un tableau de réels en double précision, dont la mémoire
   *        a déjà été allouée, et contient au moins nd+1 entrées.
   ***************************************************************************/
  /**
   * @param x is the point where the polynomial has to be evaluated.
   */
  void nt_poly_eval_der(double x, double *coeff, int n,
			double *dpoly, int nd);

  /***************************************************************************
   * Détermination d'un polynôme de degré 1.
   ***************************************************************************
   * Note pour les polynômes de dimension supérieure, il est conseillé de
   * procédérer par résolution matricielle.
   ***************************************************************************/
  void nt_poly_line(double *a, double *b, double x1, double x2, double y1,
		    double y2);

#ifdef __cplusplus
}
#endif

#endif  /* __NTOOLS_POLY_H__ */
