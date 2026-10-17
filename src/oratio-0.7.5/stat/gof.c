/**
 * gof.c
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2005-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/ntools.h"
#include "oratio/st_chi2.h"
#include "oratio/stat.h"

#include <math.h>

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
			double alpha, double *x2, int df)
{
  double d;

  int i;

  *x2 = 0.0;
  for (i = 0; i < n; i++) {
    d = obs[i*ldo]-exp[i*lde];
    *x2 += d*d/exp[i*lde];
  }

  if (df <= 0) df = n-1;

  return (st_chi2_cdf(*x2, df) > 1-alpha) ? 0 : 1;
}

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
 *   larger just by chance, for distributions with equal means.
 * @return the value of *t.
 */
double st_ttest(double *data1, int n1,
		double *data2, int n2,
		double *t, double *prob)
{
  double mean1, mean2, var1, var2;
  double pv; /* pooled variance */
  int df;    /* number of degrees of freedom */

  st_mean_variance(data1, 1, 0, n1, &mean1, &var1);
  st_mean_variance(data2, 1, 0, n2, &mean2, &var2);

  df = n1+n2-2;

  pv=((n1-1)*var1+(n2-1)*var2)/df;

  *t = (mean1-mean2)/sqrt(pv*(1.0/n1+1.0/n2));
  *prob = nt_betai(0.5*df, 0.5, df/(df*(*t)*(*t)))/2;

  return *t;
}
