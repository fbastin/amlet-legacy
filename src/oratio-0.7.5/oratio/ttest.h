/**
 * ttest.h
 *
 * STAT Library
 * Copyright (C) 2005 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the MIT License
 */

#ifndef __STAT_TTEST_H__
#define __STAT_TTEST_H__

#ifdef __cplusplus
extern "C" {
#endif

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

#endif  /* __STAT_TTEST_H__ */
