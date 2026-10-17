/**
 * estimator.c
 *
 * ORATIO Library
 * Copyright (C) 1999-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "math.h"

/**
 * Compute the bias of an estimator.
 *
 * @param estimator is the estimator value.
 * @param value is the true value.
 *
 * @return the corresponding bias.
 */
double st_bias(double estimator, double value)
{
  return estimator-value;
}

/**
 * Compute the relative bias of an estimator.
 *
 * @param estimator is the estimator value.
 * @param value is the true value.
 *
 * @return the corresponding relative bias if value is non-null,
 *   otherwise return the standard bias.
 */
double st_relative_bias(double estimator, double value)
{
  return (value == 0.0) ? estimator-value : (estimator-value)/fabs(value);
}

/**
 * Compute the mean square error of an estimator.
 *
 * @param bias is the estimator bias.
 * @param variance is the estimator variance.
 *
 * @return the mean square error of the estimator.
 */
double st_mse(double bias, double variance)
{
  return bias*bias+variance;
}


/**
 * Compute the absolute error of an estimator.
 *
 * @param bias is the estimator bias.
 * @param variance is the estimator variance.
 *
 * @return the absolute error of the estimator.
 */
double st_absolute_error(double bias, double variance)
{
  return sqrt(bias*bias+variance);
}

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
double st_relative_error(double bias, double variance, double value)
{
  return (value == 0.0) ? sqrt(bias*bias+variance) :
    sqrt(bias*bias+variance)/fabs(value);
}
