/**
 * logit.h
 *
 * AMLET Library
 * Copyright (C) 2011 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __AMLET_LOGIT_H__
#define __AMLET_LOGIT_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Apply a safeguard to the utilities, such that the logit formula is
   * mathematically equivalent to the original one, but is numerically
   * more stable.
   *
   * @param nalt is the number of alternatives.
   */
  void amlet_basic_utilities_safeguard(double *utilities,
				       int *availability, int nalt);

  /**
   * Apply a safeguard to the utilities, such that the logit formula is
   * mathematically equivalent to the original one, but is numerically
   * more stable.
   *
   * @param utilities[i] is the utility value of alternative i; and on
   *   return utilities[i] is the exponential of the original value.
   * @param den is the sum of the exponential of utilities values on
   *   return; its value on input in undefined.
   * @param availability[i] is equal to 1 if the alternative i is
   *   available, 0 else.
   * @param nalt is the number of alternatives.
   * @param choice is the alternative choice if it superior or equal to
   *   0; a negative value can be use to indicate that the choice value is
   *   unknown.
   *
   * @return the probability choice if choice is positive (or null),
   *   -1.0 else.
   */
  double amlet_utilities_safeguard(double *utilities, double *den,
				   int *availability, int nalt, int choice);

  double dcl_logit(double *utils, int *avail, int nalt, int j);

  void dcl_logit_probs(double *prob,
		       double *utils, int *avail, int nalt);

  /**
   * Compute the mean null multinomial logit log-likelihood.
   *
   * @param nalt is the total number of alternatives.
   * @param nobs is the number of observations.
   * @param availability is the alternatives availability array; it can
   *   be set to NULL if all alternatives are available to all
   *   individuals.
   *
   * @return the null multinomial logit log-likelihood, averaged over the number
   *   of individuals.
   */
  double amlet_mnl_null(int nalt, int nobs, int **availability);

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_LOGIT_H__ */
