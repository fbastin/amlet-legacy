/**
 * truncation.c
 *
 * ORATIO Library
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Evaluation of the probability distribution of a truncated
 * distribution.
 *
 * @param Fmin is the value of the CDF at the lower bound of the
 *   interval.
 * @param Fmax is the value of the CDF at the upper bound of the
 *   interval.
 *
 * @return the probability function value of the truncated distribution.
 */
double st_truncated_pdf(double fx, double Fmin, double Fmax)
{
  return fx/(Fmax-Fmin);
}

/**
 * Evaluation of the cumulative distribution function of a truncated
 * distribution.
 *
 * @param Fmin is the value of the CDF at the lower bound of the
 *   interval.
 * @param Fmax is the value of the CDF at the upper bound of the
 *   interval.
 *
 * @return the CDF value of the truncated distribution.
 */
double st_truncated_cdf(double Fx, double Fmin, double Fmax)
{
  return (Fx-Fmin)/(Fmax-Fmin);
}
