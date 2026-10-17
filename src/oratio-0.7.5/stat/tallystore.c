/**
 * tally_store.c
 *
 * ORATIO Library
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the BSD License
 */

#include <oratio/st_tally.h>
#include <oratio/ctools.h>
#include <oratio/ntools.h>
#include <oratio/st_descriptive.h>
#include <oratio/st_normal.h>

/**
 * Create a new statistical collector of type ST_TallyStore.
 *
 * @return a newly created ST_TallyStore collector.
 */
ST_TallyStore *st_tallystore_new()
{
  ST_TallyStore *ts = _malloc(ST_TallyStore);
  ts->tally = st_tally_new();
  ts->obs = _malloc(double);
  ts->nstored_max = 1;

  return ts;
}

/**
 * Initialize the tallystore collector information.
 *
 * @param ts is an initialized tallystore collector.
 */
void st_tallystore_init(ST_TallyStore *ts)
{
  st_tally_init(ts->tally);
}

/**
 * Reset the tallystore collector information.
 *
 * @param ts is an initialized tally collector.
 */
void st_tallystore_reset(ST_TallyStore *ts)
{
  st_tally_reset(ts->tally);
  free(ts->obs);
  ts->obs = _malloc(double);
  ts->nstored_max = 1;
}

/**
 * Free the memory associated to a TallyStore.
 *
 * @param ts is a statistical collector of type ST_TallyStore.
 */
void st_tallystore_free(ST_TallyStore *ts)
{
  if (ts) {
    if (ts->obs) free(ts->obs);
    if (ts->tally) free(ts->tally);
    free(ts);
  }
}

/**
 * Set the name of the TallyStore.
 *
 * @param ts is a statistical collector of type ST_TallyStore.
 * @param name is an initialized string.
 */
void st_tallystore_set_name(ST_TallyStore *ts, const char *name)
{
  st_tally_set_name(ts->tally, name);
}

/**
 * Get the name of the TallyStore.
 *
 * @param ts is a statistical collector of type ST_TallyStore.
 *
 * @return a pointer to the name of the TallyStore.
 */
char *st_tallystore_get_name(ST_TallyStore *ts)
{
  return st_tally_get_name(ts->tally);
}

/**
 * Add an observation to a tallystore collector.
 *
 * @param ts is a statistical collector of type ST_TallyStore.
 * @param x is the observation to add.
 */
void st_tallystore_add(ST_TallyStore *ts, double x)
{
  int nobs;

  st_tally_add(ts->tally, x);
  nobs = st_tally_nobs(ts->tally);

  if (nobs > ts->nstored_max) {
    ts->nstored_max *= 2;
    ts->obs = _reallocs(ts->obs, ts->nstored_max, double);
  }
  ts->obs[nobs-1] = x;
}

/**
 * Return the number of observation of the collector since its last
 * initialization.
 *
 * @param ts is a statistical collector of type ST_TallyStore.
 *
 * @return the number of observations.
 */
int st_tallystore_nobs(ST_TallyStore *ts)
{
  return st_tally_nobs(ts->tally);
}

/**
 * Return the average of the observations, as collected by a
 * ST_TallyStore collector, since its last initialization.
 *
 * @param ts is a statistical collector of type ST_TallyStore.
 *
 * @return the average of the observations.
 */
double st_tallystore_average(ST_TallyStore *ts)
{
  return st_tally_average(ts->tally);
}

/**
 * Return the variance of the observations collected by a
 * ST_TallyStore collector, since its last initialization.
 *
 * @param ts is a statistical collector of type ST_TallyStore.
 *
 * @return the variance of the observations.
 */
double st_tallystore_variance(ST_TallyStore *ts)
{
  return st_tally_variance(ts->tally);
}

/**
 * Return the standard deviation of the observations collected by a
 * ST_TallyStore collector, since its last initialization.
 *
 * @param ts is a statistical collector of type ST_TallyStore.
 *
 * @return the standard deviation of the observations.
 */
double st_tallystore_stdv(ST_TallyStore *ts)
{
  return st_tally_stdv(ts->tally);
}

/**
 * Computes a Normal confidence interval on the mean.
 *
 * @param ts is a statistical collector of type ST_TallyStore.
 * @param level desired probability that the (random) confidence interval
 *   covers the true mean (a constant).
 * @param center contains on output the center of the confidence interval.
 * @param radius contains on output the radius of the confidence interval.
 */
void st_tallystore_ci_normal(ST_TallyStore *ts, double level,
			     double *center, double *radius)
{
  st_tally_ci_normal(ts->tally, level, center, radius);
}

/**
 * Computes a Student confidence interval on the mean.
 *
 * @param ts is a statistical collector of type ST_TallyStore.
 * @param level desired probability that the (random) confidence interval
 *   covers the true mean (a constant).
 * @param center contains on output the center of the confidence interval.
 * @param radius contains on output the radius of the confidence interval.
 */
void st_tallystore_ci_student(ST_TallyStore *ts, double level,
			      double *center, double *radius)
{
  st_tally_ci_student(ts->tally, level, center, radius);
}

/**
 * Computes the confidence interval of the ratio of expectations using
 * Delta theorem.
 *
 * @param ts1 is a statistical collector of type ST_TallyStore.
 * @param ts2 is a statistical collector of type ST_TallyStore.
 * @param level desired probability that the (random) confidence interval
 *   covers the true mean (a constant).
 * @param center contains on output the center of the confidence interval.
 * @param radius contains on output the radius of the confidence interval.
 */
void st_tallystore_ci_delta(ST_TallyStore *ts1, ST_TallyStore *ts2,
			    double level, double *center,
			    double *radius)
{

  double mean_ts1 = st_tallystore_average(ts1);
  double mean_ts2 = st_tallystore_average(ts2);
  double ratio_ts = mean_ts1/mean_ts2;
  double var_ts1 = st_tallystore_variance(ts1);
  double var_ts2 = st_tallystore_variance(ts2);
  double cov = st_covariance(ts1->tally->nobs, ts1->obs, 1, mean_ts1, ts2->obs,
			     1, mean_ts2, 0, ts1->tally->nobs);
  double ratio_var = (var_ts1 + ratio_ts * ratio_ts * var_ts2
		      - 2 * cov * ratio_ts) / (mean_ts2 * mean_ts2);

  double quantile = st_gaussian_icdf((1 + level) / 2);

  *center = ratio_ts;
  *radius = quantile * sqrt(ratio_var / ts1->tally->nobs);
}

/**
 * Basic report of a Tally collector.
 *
 * @param output is the file where the report has to be print.
 * @param ts is a TallyStore collector.
 */
void st_tallystore_report_basic(FILE *output, ST_TallyStore *ts)
{
  st_tally_report_basic(output, ts->tally);
}

/**
 * Normal confidence interval report of a TallyStore collector.
 *
 * @param output is the file where the report has to be print.
 * @param tally is a Tally collector.
 */
void st_tallystore_report_ci_normal(FILE *output, ST_TallyStore *ts,
				    double level)
{
  st_tally_report_ci_normal(output, ts->tally, level);
}

/**
 * Student confidence interval report of a TallyStore collector.
 *
 * @param output is the file where the report has to be print.
 * @param tally is a Tally collector.
 */
void st_tallystore_report_ci_student(FILE *output, ST_TallyStore *ts,
				     double level)
{
  st_tally_report_ci_student(output, ts->tally, level);
}

/**
 * Confidence interval report of the ratio of expectations using
 * Delta theorem. Two TallyStore collectors are needed.
 *
 * @param output is the file where the report has to be print.
 * @param ts1 is a Tally collector.
 * @param ts2 is a Tally collector.
 * @param level desired probability that the (random) confidence interval
 *   covers the true mean (a constant).
 */
void st_tallystore_report_ci_delta(FILE *output, ST_TallyStore *ts1,
				   ST_TallyStore *ts2, double level)
{
  double center, radius;
  st_tallystore_ci_delta(ts1, ts2, level, &center, &radius);
  fprintf(output, "Mean ratio: %g\n", center);
  fprintf(output, "CI (delta): [ %g, %g ]\n",
	  center-radius, center+radius);
}
