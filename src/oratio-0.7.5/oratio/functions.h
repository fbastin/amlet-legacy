/**
 * functions.h
 *
 * ORATIO Library
 * Copyright (C) 2002-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_FUNCTIONS_H__
#define __ORATIO_FUNCTIONS_H__

#ifdef __cplusplus
extern "C" {
#endif  /* __cplusplus */

  /**
   * Evaluate the function B(z, w)
   */
  double nt_beta(double z, double w);
  
  /**
   * Returns the incomplete beta function of parameters a, b, at x
   *
   * @param x is the point whete the incomplete beta function has to be
   * evaluated; x must be in [0,1].
   */
  double nt_betai(double a, double b, double x);

  /**
   * Implementation of special functions
   * TODO: check if there is no redefinition of functions implemented
   * in standard librairies (probably coded in a better way)
   */

  /**
   * Compute n!
   *
   * @param n is a positive integer (n >= 0).
   *
   * @return n! as a double.
   */
  double nt_factorial(int n);

  /**
   * Compute the binomial coefficient with parameters n and k (k < n).
   */
  double nt_binomial_coef(int n, int k);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
