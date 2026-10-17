/**
 * st_utils.h
 *
 * ORATIO library (statistical module)
 * Copyright (C) 1999-2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the MIT License
 */

/**
 * Statistical utility functions
 */

#ifndef __ORATIO_ST_UTILS_H__
#define __ORATIO_ST_UTILS_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Incomplete gamma function P(a,x).
   */
  double st_gamma_p(double a, double x);

  /**
   * Fonction semblable à st_gamma_p, mais donne également ln Gamma(a)
   */
  double st_gamma_p_gln(double a, double x, double *gln);

  /**
   * Calcule la fonction gamma incomplète Q(a,x) = 1-P(a,x)
   */
  double st_gamma_q(double a, double x);

  /**
   * Fonction semblable à st_gamma_q, mais donne également ln Gamma(a)
   */
  double st_gamma_q_gln(double a, double x, double *gln);

  /**
   * Error function.
   * Explicitement, 
   *   st_erf(x) = 2/sqrt(pi)* integral from 0 to x of exp(-t*t) dt
   * NOTE: cette fonction a normalement sont équivalent dans la librairie
   * mathématique standard, sous le nom erf(x).
   */
  double st_erf(double x);

  /**
   * Calcule la fonction d'erreur erfc(x)
   * si x >= 0, erfc = Q(0.5, x²) (gamma incomplète), sinon erfc = 2-Q(0.5, x²)
   */
  double st_erfc(double x);

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
  double st_bct(double x, double lambda);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __ORATIO_ST_UTILS__ */
