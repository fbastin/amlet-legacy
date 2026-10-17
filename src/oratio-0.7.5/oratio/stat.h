/**
 * stat.h
 *
 * ORATIO Library
 * Copyright (C) 1999-2015 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_STAT_H__
#define __ORATIO_STAT_H__

/**
 * THIS LIBRARY CONTAINS A SET OF TOOLS FOR STATISTICAL COMPUTATIONS.
 */

#include <math.h>
#include <sys/types.h>

#include "st_descriptive.h"
#include "st_discrete.h"
#include "st_likelihood.h"
#include "st_obs.h"

#include "st_utils.h"

#include "st_cauchy.h"
#include "st_chi2.h"
#include "st_exponential.h"
#include "st_fisher.h"
#include "st_geometric.h"
#include "st_laplace.h"
#include "st_logistic.h"
#include "st_normal.h"
#include "st_pareto.h"
#include "st_poisson.h"
#include "st_rayleigh.h"
#include "st_student.h"
#include "st_uniform.h"
#include "st_weibull.h"

/* Statistical collectors */
#include "st_tally.h"

#include "percentile.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Compute the bias of an estimator.
   *
   * @param estimator is the estimator value.
   * @param value is the true value.
   *
   * @return the corresponding bias.
   */
  double st_bias(double estimator, double value);

  /**
   * Compute the relative bias of an estimator.
   *
   * @param estimator is the estimator value.
   * @param value is the true value.
   *
   * @return the corresponding relative bias.
   */
  double st_relative_bias(double estimator, double value);

  /**
   * Compute the mean square error of an estimator.
   *
   * @param bias is the estimator bias.
   * @param variance is the estimator variance.
   *
   * @return the mean square error of the estimator.
   */
  double st_mse(double bias, double variance);

  /**
   * Compute the absolute error of an estimator.
   *
   * @param bias is the estimator bias.
   * @param variance is the estimator variance.
   *
   * @return the absolute error of the estimator.
   */
  double st_absolute_error(double bias, double variance);

  /**
   * Compute the relative error of an estimator.
   *
   * @param bias is the estimator bias.
   * @param variance is the estimator variance.
   * @param value is the true value.
   *
   * @return the relative error of the estimator if value is non-null,
   *   otherwise return the absolute error.
   */
  double st_relative_error(double bias, double variance, double value);

  /**
   * Erlang distribution.
   */

  /**
   * Compute the Erlang probability distribution function.
   *
   * @param x is the point where to evaluate the probability
   *   distribution.
   * @param alpha is the shape parameter of the Erlang distribution.
   * @param lambda is the scale parameter of the Erlang distribution.
   *
   * @return the probability distribution of an Erlang random
   *   variable, evaluated at x.
   */
  double st_erlang_pdf(double x, double alpha, double lambda);

  /**
   * Compute the Erlang cumulative distribution function.
   *
   * @param x is the point where to evaluate the distribution.
   * @param alpha is the shape parameter of the Erlang distribution.
   * @param lambda is the scale parameter of the Erlang distribution.
   *
   * @return the cumulative distribution of an Erlang random
   *   variable, evaluated at x.
   */
  //  double st_erlang_cdf(double x, double alpha, double lambda);

  /**
   * Compute the inverse distribution function of an erlang random
   * variable.
   *
   * @param u is the point where to evaluated the inverse distribution
   *   function.
   * @param alpha is the shape parameter of the Erlang distribution.
   * @param lambda is the scale parameter of the Erlang distribution.
   *
   * @return the inverse distribution function of an erlang
   *   distribution, evaluated at u.
   */
  //  double st_erlang_icdf(double u, double alpha, double lambda);

  /**
   * Compute the mean of an Erlang random variable.
   *
   * @param alpha is the shape parameter of the Erlang distribution.
   * @param lambda is the scale parameter of the Erlang distribution.
   *
   * @return the mean of an erlang random variable of parameters alpha
   *   and lambda.
   */
  double st_erlang_mean(double alpha, double lambda);

  /**
   * Compute the standard deviation of an Erlang random variable.
   *
   * @param alpha is the shape parameter of the Erlang distribution.
   * @param lambda is the scale parameter of the Erlang distribution.
   *
   * @return the standard deviation of an erlang random variable of
   *   parameters alpha and lambda.
   */
  double st_erlang_stdv(double alpha, double lambda);

  /**
   * Compute the variance of an Erlang random variable.
   *
   * @param alpha is the shape parameter of the Erlang distribution.
   * @param lambda is the scale parameter of the Erlang distribution.
   *
   * @return the variance of an erlang random variable of parameters
   *   alpha and lambda.
   */
  double st_erlang_variance(double alpha, double lambda);

  /**
   * Gamma distribution.
   */

  /**
   * Gamma probability distribution function.
   *
   * @param x is the point where the probability distribution is evaluated.
   * @param alpha > 0 is the shape parameter of the gamma distribution.
   * @param lambda > 0 is the scale parameter of the gamma distribution.
   *
   * @return the value at x of the density of a gamma variate of shape
   *   parameter k and scale parameter theta.
   */
  double st_gamma_pdf(double x, double alpha, double lambda);

  /**
   * Mean of a Gamma distribution function of parameters alpha and lambda.
   *
   * @param alpha > 0 is the shape parameter of the gamma distribution.
   * @param lambda > 0 is the scale parameter of the gamma distribution.
   *
   * @return the mean of a Gamma distribution function of parameters
   *   alpha and lambda.
   */
  double st_gamma_mean(double alpha, double lambda);

  /**
   * Variance of a Gamma distribution function of parameters alpha and lambda.
   *
   * @param alpha > 0 is the shape parameter of the gamma distribution.
   * @param lambda > 0 is the scale parameter of the gamma distribution.
   *
   * @return the variance of a Gamma distribution function of parameters
   *   alpha and lambda.
   */
  double st_gamma_variance(double alpha, double lambda);

  /**
   * Standard of a Gamma distribution function of parameters alpha and
   * lambda.
   *
   * @param alpha > 0 is the shape parameter of the gamma distribution.
   * @param lambda > 0 is the scale parameter of the gamma distribution.
   *
   * @return the standard deviation of a Gamma distribution function of
   *   parameters alpha and lambda.
   */
  double st_gamma_stdv(double alpha, double lambda);

  /**
   * Mode of a Gamma distribution function of parameters alpha and lambda.
   *
   * @param alpha >= 1.0 is the shape parameter of the gamma distribution.
   * @param lambda > 0 is the scale parameter of the gamma distribution.
   *
   * @return the mode of a Gamma distribution function of parameters
   *   alpha and lambda.
   */
  double st_gamma_mode(double alpha, double lambda);

#include "st_gumbel.h"
#include "st_lognormal.h"
#include "st_triangular.h"

#include "stat_empirical.h"
#include "stat_error.h"

  /**
   * Multivariate distributions.
   */

#include "st_mvn.h"

  /**
   * Goodness of fit.
   */

  /**
   * Pearson's chi-square test.
   *
   * @param n is the number of observations.
   * @param obs is the observations vector.
   * @param ldo is the increment in obs.
   * @param exp is the expected vector.
   * @param lbe is the increment in exp.
   * @param alpha is the test level.
   * @param x2 is the chi-square probability.
   * @param df is the number of degree of freedom; is set to be less or
   *   equal to 0, n-1 degrees of freedom will be considered.
   *
   * @return 1 if the hypothesis that obs follows the same
   *   distribution the exp, 0 else.
   */
  int st_pearson_chi_test(int n, double *obs, int ldo, double *exp, int lde,
			  double alpha, double *x2, int df);

  /**
   * t-test between two datasets representing independent
   * samples. Tests to see if the difference between means of the
   * samples if different from zero. The samples are assumed to be drawn
   * from populations with the same true variance.
   *
   * Ref.: "Numerical Recipes in C", William H. Press, Saul A. Teukolsky,
   * William T. Vettering, Brian P. Flannery, 1992, p 616.
   *
   * @param *t is the t-value
   * @param *prob is the probability that |*t| could be this large or
   * larger just by chance, for distributions with equal means.
   * @return the value of *t.
   */
  double st_ttest(double *data1, int n1,
		  double *data2, int n2,
		  double *t, double *pref);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_STAT_H__ */
