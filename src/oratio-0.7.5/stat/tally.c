/**
 * tally.c
 *
 * ORATIO Library
 * Copyright (C) 2012 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the BSD License
 */

#include "oratio/ntools.h"
#include "oratio/st_normal.h"
#include "oratio/st_student.h"
#include "oratio/st_tally.h"

#include <string.h>

/**
 * Create a new statistical collector of type ST_Tally.
 *
 * @return a newly created ST_Tally collector.
 */
ST_Tally *st_tally_new()
{
  ST_Tally *tally = _malloc(ST_Tally);

  st_tally_init(tally);

  tally->name = NULL;

  return tally;
}

/**
 * Initialize the tally collector information.
 *
 * @param tally is a memory-allocated tally collector.
 */
void st_tally_init(ST_Tally *tally)
{
  tally->nobs = 0;
 
  tally->max_value = -nt_infinity();
  tally->min_value = nt_infinity();
  tally->sum_value = tally->sum_squares = tally->cur_average
    = tally->cur_sum2 = 0.0;
}

/**
 * Reset the tally collector information.
 *
 * @param tally is an initialized tally collector.
 */
void st_tally_reset(ST_Tally *tally)
{
  if (tally->name) free(tally->name);
  st_tally_init(tally);
}

/**
 * Free the memory associated to a ST_Tally collector.
 *
 * @param tally is a statistical collector of type ST_Tally.
 */
void st_tally_free(ST_Tally *tally)
{
  if (tally) {
    if (tally->name) free(tally->name);
    free(tally);
  }
}

/**
 * Set the name of the Tally.
 *
 * @param tally is a statistical collector of type ST_Tally.
 * @param name is an string.
 */
void st_tally_set_name(ST_Tally *tally, const char *name)
{
  tally->name = _reallocs(tally->name, strlen(name)+1, char);
  strcpy(tally->name, name);
}

/**
 * Get the name of the Tally.
 *
 * @param tally is a statistical collector of type ST_Tally.
 *
 * @return a pointer to the name of the Tally.
 */
char *st_tally_get_name(ST_Tally *tally)
{
  return tally->name;
}

/**
 * Add an observation to a tally collector.
 *
 * @param tally is a statistical collector of type ST_Tally.
 * @param x is the observation to add.
 */
void st_tally_add(ST_Tally *tally, double x)
{
  double delta;

  if (x < tally->min_value) tally->min_value = x;
  if (x > tally->max_value) tally->max_value = x;
  tally->nobs++;
  tally->sum_value += x;
  tally->sum_squares += x*x;

  /* The single formula to compute variance is due to
   * B.P. Welford (1962), "Note on a method for calculating corrected sums of
   * squares and products", Technometrics 4(3), pp. 419--420. */
  delta = x - tally->cur_average;
  //  tally->cur_sum2 += (tally->nobs-1)*(x-tally->cur_average)*(x-tally->cur_average)/tally->nobs;
  tally->cur_average += delta/tally->nobs;
  tally->cur_sum2 += delta*(x-tally->cur_average);
}

/**
 * Return the number of observation of the collector since its last
 * initialization.
 *
 * @param tally is a statistical collector of type ST_Tally.
 *
 * @return the number of observations.
 */
int st_tally_nobs(ST_Tally *tally)
{
  return tally->nobs;
}

/**
 * Return the average of the observations, as collected by a
 * ST_Tally collector, since its last initialization.
 *
 * @param tally is a statistical collector of type ST_Tally.
 *
 * @return the average of the observations.
 */
double st_tally_average(ST_Tally *st)
{
  return (st->nobs > 0) ? st->cur_average : nt_nan();
}

/**
 * Return the variance of the observations collected by a
 * ST_Tally collector, since its last initialization.
 *
 * @param tally is a statistical collector of type ST_Tally.
 *
 * @return the variance of the observations.
 */
double st_tally_variance(ST_Tally *tally)
{
  return (tally->nobs > 1) ?
    (tally->cur_sum2)/(tally->nobs-1) : nt_nan();
}

/**
 * Return the standard deviation of the observations collected by a
 * ST_Tally collector, since its last initialization.
 *
 * @param tally is a statistical collector of type ST_Tally.
 *
 * @return the standard deviation of the observations.
 */
double st_tally_stdv(ST_Tally *st)
{
  return sqrt(st_tally_variance(st));
}

/**
 * Computes a Normal confidence interval on the mean.
 *
 * @param tally is a statistical collector of type ST_Tally.
 * @param level desired probability that the (random) confidence interval
 *   covers the true mean (a constant).
 * @param center contains on output the center of the confidence interval.
 * @param radius contains on output the radius of the confidence interval.
 */
void st_tally_ci_normal(ST_Tally *st, double level,
			double *center, double *radius)
{
  double z = st_gaussian_icdf(0.5*(level+1.0));

  *center = st_tally_average(st);
  *radius = z*sqrt(st_tally_variance(st)/st->nobs);
}

/**
 * Computes a Student confidence interval on the mean.
 *
 * @param tally is a statistical collector of type ST_Tally.
 * @param level desired probability that the (random) confidence interval
 *   covers the true mean (a constant).
 * @param center contains on output the center of the confidence interval.
 * @param radius contains on output the radius of the confidence interval.
 */
void st_tally_ci_student(ST_Tally *st, double level,
			 double *center, double *radius)
{
  double z = st_student_icdf(0.5*(level+1.0), st->nobs-1);
  *center = st_tally_average(st);
  *radius = z*sqrt(st_tally_variance(st)/st->nobs);
}

/**
 * Basic report of a Tally collector.
 *
 * @param output is the file where the report has to be print.
 * @param tally is a Tally collector.
 */
void st_tally_report_basic(FILE *output, ST_Tally *tally)
{
  fprintf(output, "Report on Tally");
  if (tally->name) fprintf(output, " %s", tally->name);
  fprintf(output, "\n");
  fprintf(output, "Number of observations: %d\n", tally->nobs);
  fprintf(output, "Minimum: %g\tMaximum: %g\n",
	  tally->min_value, tally->max_value);
  fprintf(output, "Average: %g\tStandard deviation: %g\n",
	  st_tally_average(tally), st_tally_stdv(tally));
}

/**
 * Normal confidence interval report of a Tally collector.
 *
 * @param output is the file where the report has to be print.
 * @param tally is a Tally collector.
 */
void st_tally_report_ci_normal(FILE *output, ST_Tally *tally,
			       double level)
{
  double center, radius;

  st_tally_ci_normal(tally, level, &center, &radius);
  fprintf(output, "CI (normal): [ %g, %g ]\n",
	  center-radius, center+radius);
}

/**
 * Student confidence interval report of a Tally collector.
 *
 * @param output is the file where the report has to be print.
 * @param tally is a Tally collector.
 */
void st_tally_report_ci_student(FILE *output, ST_Tally *tally,
				double level)
{
  double center, radius;

  st_tally_ci_student(tally, level, &center, &radius);
  fprintf(output, "CI (student): [ %g, %g ]\n",
	  center-radius, center+radius);
}
