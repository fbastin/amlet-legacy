/**
 * stat_test.c
 *
 * ORATIO Library
 * Copyright (C) 1999-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "stat_test.h"
#include "descriptive_test.h"
#include "empirical_test.h"
#include "lognormal.h"
#include "poisson.h"

#include <oratio/ntools.h>
#include <oratio/stat.h>

#include <stdio.h>

/**
 * Chi-square test.
 *
 * @return an error code if the evaluation has failed.
 */
int chi_test()
{
  double nu;
  double x;

  printf("Test on the chi-square distribution\n");
  printf("-----------------------------------\n");

  nu = 30;
  x = 15.0;

  x = st_chi2_cdf(x, nu);
  printf("P[X2,30 < 15.0] = %g\n", x);
  if (fabs(x-0.01026) > 1e-5) {
    c_warning("incorrect value in chi squre test:\nwe obtained %g instead of an approximate value of 0.01026.", x);
    return 1;
  }

  return 0;
}

double normal_cdf(double x, double *p)
{
  return st_normal_cdf(x, p[0], p[1]);
}

/**
 * Test on generic quantile functions.
 */
void quantile_test()
{
  double quantile;
  double pars[2] = {0.0, 1.0};

  quantile = st_percentile(0.45, normal_cdf, pars,
			   0.0, 1.0, 0.0);
  printf("Quantile 0.450 of a normal variable N(0,1): %g\n", quantile);
  quantile = st_percentile(0.65, normal_cdf, pars,
			   0.0, 1.0, 0.0);
  printf("Quantile 0.650 of a normal variable N(0,1): %g\n", quantile);
  quantile = st_percentile(0.85, normal_cdf, pars,
			   0.0, 1.0, 0.0);
  printf("Quantile 0.850 of a normal variable N(0,1): %g\n", quantile);
  quantile = st_percentile(0.975, normal_cdf, pars,
			   0.0, 1.0, 0.0);
  printf("Quantile 0.975 of a normal variable N(0,1): %g\n", quantile);
}

void inv_cum_test()
{
  double q, x;

  q = st_gaussian_cdf(1.0);
  printf("Quantile at %g: %g\n", 1.0, q);

  x = st_gaussian_icdf(0.45);
  printf("Gaussian at quantile 0.45: %g\n", st_gaussian_icdf(0.45));

  q = st_gaussian_cdf(x);
  printf("Quantile at %g: %g\n", x, q);

  x = st_gaussian_icdf(0.65);
  printf("Gaussian at quantile 0.65: %g\n", x);

  q = st_gaussian_cdf(x);
  printf("Quantile at %g: %g\n", x, q);

  x = st_gaussian_icdf(0.85);
  printf("Gaussian at quantile 0.85: %g\n", st_gaussian_icdf(0.85));

  q = st_gaussian_cdf(x);
  printf("Quantile at %g: %g\n", x, q);

  x = st_gaussian_icdf(0.975);
  printf("Gaussian at quantile 0.975: %g\n", st_gaussian_icdf(0.975));

  q = st_gaussian_cdf(x);
  printf("Quantile at %g: %g\n", x, q);

  x = st_gaussian_icdf(0.99);
  printf("Gaussian at quantile 0.99: %g\n", st_gaussian_icdf(0.99));

  q = st_gaussian_cdf(x);
  printf("Quantile at %g: %g\n", x, q);

  x = st_gaussian_icdf(0.995);
  printf("Gaussian at quantile 0.995: %g\n", x);

  q = st_gaussian_cdf(x);
  printf("Quantile at %g: %g\n", x, q);

  x = st_gaussian_icdf(1e-12);
  printf("Gaussian at quantile 1e-12: %g\n", x);

  q = st_gaussian_cdf(x);
  printf("Quantile at %g: %g\n", x, q);

  x = st_gaussian_icdf(0.0000002);
  printf("Gaussian at quantile 0.0000002: %g\n", x);

  q = st_gaussian_cdf(x);
  printf("Quantile at %g: %g\n", x, q);

  x = st_normal_icdf(0.65, 1, 0.5);
  printf("Normal(1,0.5) at quantile 0.65: %g\n", x);

  q = st_normal_cdf(x, 1, 0.5);
  printf("Quantile at %g: %g\n", x, q);

  q = st_normal_cdf(0.0, -8.419, 5.164);
  printf("Quantile at %g: %g\n", 0.0, q);

  x = st_lognormal_icdf(0.65, 1, 0.5);
  printf("Lognormal(1,0.5) at quantile 0.65: %g\n", x);

  q = st_lognormal_cdf(x, 1, 0.5);
  printf("Quantile at %g: %g\n", x, q);
}

int main()
{
  int err, error_code = 0;

  err = st_descriptive_test();
  if (err) printf("Descriptive statistics test successfully passed.\n\n");
  error_code = nt_max(error_code, err);

  err = chi_test();
  if (!err) printf("Chi square test successfully passed.\n\n");
  error_code = nt_max(error_code, err);

  quantile_test();

  inv_cum_test();

  err = st_cauchy_test();
  error_code = nt_max(error_code, err);

  err = st_empirical_test();
  error_code = nt_max(error_code, err);

  err = st_lognormal_test();
  error_code = nt_max(error_code, err);

  err = st_poisson_test();
  error_code = nt_max(error_code, err);

  err = st_mvn_test();
  error_code = nt_max(error_code, err);

  return (error_code) ? 1 : 0;
}
