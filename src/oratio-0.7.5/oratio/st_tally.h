/**
 * st_tally.h
 *
 * ORATIO Library
 * Copyright (C) 2011-2012 Marion Hautbois,
 *                         Augustín Ávila Marín,
 *                         Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the BSD License
 */

#ifndef __ORATIO_ST_TALLY_H__
#define __ORATIO_ST_TALLY_H__

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Tally statistical collector.
   */
  typedef struct _ST_Tally{
    char *name;
    int nobs;
    double max_value;   /* maximum of the observations */
    double min_value;   /* minimum of the observations */
    double sum_value;   /* sum of the observations */
    double sum_squares; /* sum of the square of the observations */
    double cur_average; /* current observations average */
    double cur_sum2;    /* current average of the square of the observations */
  } ST_Tally;

  /**
   * Tally statistical collector, with observations storage.
   */
  typedef struct _ST_TallyStore{
    double *obs;        /* array of observations */
    int nstored_max;    /* maximum number of observations that can be stored;
			 * it can be dynamically adjusted */
    ST_Tally *tally;    /* the other operations are treated by means of
			   a ST_Tally instance */
  } ST_TallyStore;

  /* ST_TALLY */

  /**
   * Create a new statistical collector of type ST_Tally.
   *
   * @return a newly created ST_Tally collector.
   */
  ST_Tally *st_tally_new();

  /**
   * Init the tally collector information.
   *
   * @param tally is a memory-allocated tally collector.
   */
  void st_tally_init(ST_Tally *tally);

  /**
   * Reset the tally collector information.
   *
   * @param tally is an initialized tally collector.
   */
  void st_tally_reset(ST_Tally *tally);

  /**
   * Free the memory associated to a ST_Tally collector.
   *
   * @param tally is a statistical collector of type ST_Tally.
   */
  void st_tally_free(ST_Tally *tally);

  /**
   * Set the name of the Tally.
   *
   * @param tally is a statistical collector of type ST_Tally.
   * @param name is an initialized string.
   */
  void st_tally_set_name(ST_Tally *tally, const char *name);

  /**
   * Get the name of the Tally.
   *
   * @param tally is a statistical collector of type ST_Tally.
   *
   * @return a pointer to the name of the Tally.
   */
  char *st_tally_get_name(ST_Tally *tally);

  /**
   * Add an observation to a tally collector.
   *
   * @param tally is a statistical collector of type ST_Tally.
   * @param x is the observation to add.
   */
  void st_tally_add(ST_Tally *tally, double x);

  /**
   * Return the number of observation of the collector since its last
   * initialization.
   *
   * @param tally is a statistical collector of type ST_Tally.
   *
   * @return the number of observations.
   */
  int st_tally_nobs(ST_Tally *st);

  /**
   * Return the average of the observations, as collected by a
   * ST_Tally collector, since its last initialization.
   *
   * @param tally is a statistical collector of type ST_Tally.
   *
   * @return the average of the observations.
   */
  double st_tally_average(ST_Tally *st);

  /**
   * Return the variance of the observations collected by a
   * ST_Tally collector, since its last initialization.
   *
   * @param tally is a statistical collector of type ST_Tally.
   *
   * @return the variance of the observations.
   */
  double st_tally_variance(ST_Tally *tally);

  /**
   * Return the standard deviation of the observations collected by a
   * ST_Tally collector, since its last initialization.
   *
   * @param tally is a statistical collector of type ST_Tally.
   *
   * @return the standard deviation of the observations.
   */
  double st_tally_stdv(ST_Tally *st);

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
			  double *center, double *radius);

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
			  double *center, double *radius);
  
  /**
   * Normal confidence interval report of a Tally collector.
   *
   * @param output is the file where the report has to be print.
   * @param tally is a Tally collector.
   */
  void st_tally_report_ci_normal(FILE *output, ST_Tally *tally,
				 double level);

  /**
   * Student confidence interval report of a Tally collector.
   *
   * @param output is the file where the report has to be print.
   * @param tally is a Tally collector.
   */
  void st_tally_report_ci_student(FILE *output, ST_Tally *tally,
				  double level);

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
			      double *radius);
   /**
   * Basic report of a ST_Tally collector.
   *
   * @param output if the file where the report has to be print.
   * @param tally is a ST_Tally collector.
   */
  void st_tally_report_basic(FILE *output, ST_Tally *tally);

  /* ST_TALLYSTORE */

  /**
   * Create a new statistical collector of type ST_TallyStore.
   *
   * @return a newly created ST_TallyStore collector.
   */
  ST_TallyStore *st_tallystore_new();

  /**
   * Initialize the tallystore collector information.
   *
   * @param ts is a memory-allocated tallystore collector.
   */
  void st_tallystore_init(ST_TallyStore *ts);

  /**
   * Reset the tallystore collector information.
   *
   * @param ts is an initialized tallystore collector.
   */
  void st_tallystore_reset(ST_TallyStore *ts);

  /**
   * Set the name of the TallyStore.
   *
   * @param ts is a statistical collector of type ST_TallyStore.
   * @param name is an initialized string.
   */
  void st_tallystore_set_name(ST_TallyStore *ts, const char *name);

  /**
   * Get the name of the TallyStore.
   *
   * @param ts is a statistical collector of type ST_TallyStore.
   *
   * @return a pointer to the name of the TallyStore.
   */
  char *st_tallystore_get_name(ST_TallyStore *ts);

  /**
   * Free the memory associated to a TallyStore.
   *
   * @param ts is a statistical collector of type ST_TallyStore.
   */
  void st_tallystore_free(ST_TallyStore *ts);

  /**
   * Add an observation to a tallystore collector.
   *
   * @param ts is a statistical collector of type ST_TallyStore.
   * @param x is the observation to add.
   */
  void st_tallystore_add(ST_TallyStore *ts, double x);

  /**
   * Return the number of observation of the collector since its last
   * initialization.
   *
   * @param ts is a statistical collector of type ST_TallyStore.
   *
   * @return the number of observations.
   */
  int st_tallystore_nobs(ST_TallyStore *ts);

  /**
   * Return the average of the observations, as collected by a
   * ST_TallyStore collector, since its last initialization.
   *
   * @param ts is a statistical collector of type ST_TallyStore.
   *
   * @return the average of the observations.
   */
  double st_tallystore_average(ST_TallyStore *ts);

  /**
   * Return the variance of the observations collected by a
   * ST_TallyStore collector, since its last initialization.
   *
   * @param ts is a statistical collector of type ST_TallyStore.
   *
   * @return the variance of the observations.
   */
  double st_tallystore_variance(ST_TallyStore *ts);

  /**
   * Return the standard deviation of the observations collected by a
   * ST_TallyStore collector, since its last initialization.
   *
   * @param ts is a statistical collector of type ST_TallyStore.
   *
   * @return the standard deviation of the observations.
   */
  double st_tallystore_stdv(ST_TallyStore *ts);

  /**
   * Computes a confidence interval on the mean.
   *
   * @param ts is a statistical collector of type ST_TallyStore.
   * @param level desired probability that the (random) confidence interval
   *   covers the true mean (a constant).
   * @param center contains on output the center of the confidence interval.
   * @param radius contains on output the radius of the confidence interval.
   */
  void st_tallystore_ci_normal(ST_TallyStore *ts, double level,
			       double *center, double *radius);

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
				double *center, double *radius);

    /**
   * Basic report of a TallyStore collector.
   *
   * @param output if the file where the report has to be print.
   * @param ts is a ST_TallyStore collector.
   */
  void st_tallystore_report_basic(FILE *output, ST_TallyStore *ts);

  /**
   * Normal confidence interval report of a TallyStore collector.
   *
   * @param output is the file where the report has to be print.
   * @param tally is a Tally collector.
   */
  void st_tallystore_report_ci_normal(FILE *output, ST_TallyStore *ts,
				      double level);

  /**
   * Student confidence interval report of a TallyStore collector.
   *
   * @param output is the file where the report has to be print.
   * @param tally is a Tally collector.
   */
  void st_tallystore_report_ci_student(FILE *output, ST_TallyStore *ts,
				       double level);

  /**
   * Confidence interval of the ratio of expectations using Delta theorem.
   * Two TallyStore collectors are needed.
   *
   * @param output is the file where the report has to be print.
   * @param ts1 is a TallyStore collector.
   * @param ts2 is a TallyStore collector.
   * @param level desired probability that the (random) confidence interval
   *   covers the true mean (a constant).
   */
  void st_tallystore_report_ci_delta(FILE *output, ST_TallyStore *ts1,
				     ST_TallyStore *ts2, double level);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_ST_TALLY_H__ */
