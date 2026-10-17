/**
 * percentile.h
 *
 * STAT Library
 * Copyright (C) 1999-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __STAT_PERCENTILE_H__
#define __STAT_PERCENTILE_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Type des fonctions de répartions:
   * La fonction renvoie un double, et prend comme argument, dans l'ordre,
   * le point où elle doit être évaluée, la moyenne de la distribution,
   * son écart-type, un pointeur vers un double qui contiendra la valeur
   * de la fonction de répartition en sortie, et un pointeur vers un double
   * qui contiendra son complément.
   */
  typedef double (*STCDF)(double, double *);

  /**
   * Computation of percentiles of arbitrary continuous cumulative
   * distribution functions.
   *
   * @param alpha is the target qunatile.
   * @param pars is the vector of parameters associated to the random variate.
   * @param tol is the used tolerance; if tol <= 0.0, a default
   *   accuracy will be set.
   */
  double st_percentile(double alpha, STCDF f, double *pars,
		       double median, double sigma, double tol);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
