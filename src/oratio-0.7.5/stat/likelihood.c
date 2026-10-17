/**
 * likelihood.c
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2005-2015 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <math.h>

#include "oratio/st_chi2.h"
#include "oratio/st_obs.h"

/**
 * Likelihood ratio test: chi-square value.
 * Source: Ben-Akiva and Lerman, "Discrete choice analysis", 1985, p28.
 *
 * @param ll is the log-likelihood value of the model.
 * @param ll0 is the log-likelihood value of the null model.
 *
 * @return the chi-square value associated to the likelihood ratio test.
 */
double st_likelihood_chi2(double ll, double ll0)
{
  return 2*(ll-ll0);
}

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
double st_likelihood_chi2_level(int npar, double ll, double ll0)
{
  return st_chi2_cdf(2*(ll-ll0), npar);
}

/**
 * Compute the log-likelihood of a set of observations.
 *
 * @param n is the number of observations to consider.
 * @param obs is a pointer to an array of observations.
 *
 * @return the log-likelihood of a set of observations.
 */
double st_loglikelihood(int n, ST_Obs *obs)
{
  int i;
  double ll = 0;

  for (i = 0; i < n; i++) ll += (log(obs[i].prob)-ll)/(i+1);

  return ll;
}

/**
 * Update a log-likelihood value of a set of additional observations.
 *
 * @param n0 is the number of observations associated to the original
 *   log-likelihood.
 * @param n1 is the number of additional observations.
 * @param obs is a pointer to the array of additional observations.
 * @param ll is a pointer to the log-likelihood value, containing on
 *   input the initial value associated to the n0 first observations, and
 *   containing on output the log-likelihood value associated to the
 *   complete set of observations.
 *
 * @return the log-likelihood of the complete set of observations.
 */
double st_loglikelihood_update(int n0, int n1, ST_Obs *obs,
				      double *ll)
{
  int i;

  for (i = 0; i < n1; i++) *ll += (log(obs[i].prob)-*ll)/(n0+i+1);

  return *ll;
}

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
			   double *ll_grad, int incl)
{
  int i, j;

  for (j = 0; j < m; j++) ll_grad[j*incl] = 0.0;

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      ll_grad[j*incl] +=
	((obs[i].grad[j])/(obs[i].prob)-ll_grad[j*incl])/(i+1);
}

/**
 * Compute the gradient of the log-likelihood function, assuming that the
 * gradients of individual observations are available.
 *
 * @param n0 is the number of observations associated to the original
 *   log-likelihood.
 * @param n1 is the number of additional observations.
 * @param m is the parameters vector dimension.
 * @param obs is a pointer to an array of observations.
 * @param ll_grad is a pointer to the gradient of the log-likelihood
 *   function.
 * @param incl is the increment used in the array ll_grad.
 */
void st_loglikelihood_grad_update(int n0, int n1, int m, ST_Obs *obs,
				    double *ll_grad, int incl)
{
  int i , j;

  for (i = 0; i < n1; i++)
    for (j = 0; j < m; j++)
      ll_grad[j*incl] +=
	((obs[i].grad[j])/(obs[i].prob)-ll_grad[j*incl])/(n0+i+1);
}

/**
st_BHHH
*/
