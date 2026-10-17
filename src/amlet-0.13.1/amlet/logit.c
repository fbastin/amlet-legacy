/**
 * logit.c
 *
 * AMLET Library
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <math.h>

#include <oratio/ntools.h>

/**
 * Apply a safeguard to the utilities, such that the logit formula is
 * mathematically equivalent to the original one, but is numerically
 * more stable.
 *
 * @param nalt is the number of alternatives.
 */
void amlet_basic_utilities_safeguard(double *utilities,
					  int *availability, int nalt)
{
  double trans = -1.0e300;
  double util;

  int k;

  /* Determine the maximum utility. */
  for (k = 0; k < nalt; k++)
      if (availability[k]) trans = nt_max(trans, utilities[k]);

  util = exp(trans);
  if (!finite(util) || util == 0.0) {
    for (k = 0; k < nalt; k++) {
      if (availability[k]) utilities[k] -= trans;
    }
  }
}

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
					int *availability, int nalt, int choice)
{
  const double ln2 = log(2);
  double trans = -1.0e300;
  double util;

  int k;

  /* Determine the maximum utility. */
  for (k = 0; k < nalt; k++)
      if (availability[k]) trans = nt_max(trans, utilities[k]);

  util = exp(trans);
  if (!finite(util) || util == 0.0) {
    for (k = 0; k < nalt; k++)
      if (availability[k]) utilities[k] -= trans-ln2;
  }

  *den = 0.0;
  for (k = 0; k < nalt; k++) {
    if (availability[k]) {
      utilities[k] = exp(utilities[k]);
      *den += utilities[k];
    }
  }

  if (choice >= 0)
    return utilities[choice]/(*den);
  else
    return -1.0;
}

/**
 * Logit formula.
 */
double dcl_logit(double *utils, int *avail, int nalt, int j)
{
  int i;
  double den;

  den = 0.0;
  for (i = 0; i < nalt; i++) if (avail[i]) den += exp(utils[i]);

  return exp(utils[j])/den;
}

void dcl_logit_probs(double *prob,
			    double *utils, int *avail, int nalt)
{
  int i;
  double den;

  den = 0.0;
  for (i = 0; i < nalt; i++) if (avail[i]) den += exp(utils[i]);

  for (i = 0; i < nalt; i++) prob[i] = exp(utils[i])/den;
}

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
double amlet_mnl_null(int nalt, int nobs, int **availability)
{
  int i, k;
  double ll0 = 0.0;
  int na;

  for (i = 0; i < nobs; i++) {
    /* Count the number of available alternatives */
    if (availability) {
      na = 0;
      for (k = 0; k < nalt; k++) if (availability[i][k]) na++;
    }
    else na = nalt;
    ll0 += (log(1.0/na)-ll0)/(i+1);
  }

  return ll0;
}
