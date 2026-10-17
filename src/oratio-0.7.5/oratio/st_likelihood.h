/**
 * st_likelihood.h
 *
 * ORATIO Library
 * Copyright (C) 2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_ST_LIKELIHOOD_H__
#define __ORATIO_ST_LIKELIHOOD_H__

#include "st_obs.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Likelihood ratio test: chi-square value.
   * Source: Ben-Akiva and Lerman, "Discrete choice analysis", 1985, p28.
   *
   * @param ll is the log-likelihood value of the model.
   * @param ll0 is the log-likelihood value of the null model.
   *
   * @return the chi-square value associated to the likelihood ratio test.
   */
  double st_likelihood_chi2(double ll, double ll0);

  /**
   * Likelihood ratio test: chi-square level.
   * Source: Ben-Akiva and Lerman, "Discrete choice analysis", 1985, p28.
   *
   * @param npar is the number of parameters in the model.
   * @param ll is the log-likelihood value of the model.
   * @param ll0 is the log-likelihood value of the null model.
   *
   * @return the chi-square level associated to the likelihood ratio test.
   */
  double st_likelihood_chi2_level(int npar, double ll, double ll0);

  /**
   * Compute the log-likelihood of a set of observations.
   *
   * @param n is the number of observations to consider.
   * @param obs is a pointer to an array of observations.
   *
   * @return the log-likelihood of a set of observations.
   */
  double st_loglikelihood(int n, ST_Obs *obs);

  /**
   * Compute the gradient of the log-likelihood function, assuming that the
   * gradients of individual observations are available.
   *
   * @param n is the number of observations to consider.
   * @param m is the parameters vector dimension.
   * @param obs is a pointer to an array of observations.
   * @param ll_grad is a pointer to the gradient of the log-likelihood
   *   function.
   * @param incl is the increment used in the array ll_grad.
   */
  void st_loglikelihood_grad(int n, int m, ST_Obs *obs,
			     double *ll_grad, int incl);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_ST_LIKELIHOOD_H__ */
