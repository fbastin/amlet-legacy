/**
 * bootstrap.c
 *
 * ORATIO Library
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * This piece of code is covered by the BSD License.
 */

#include "oratio/ran.h"

/**
 * Generate a bootstrap sample.
 *
 * @param ran is an initialized Random numbers generators.
 * @param n is the sample size.
 * @param x is the original samble, which is not modified on output.
 * @param xb is a vector of size a least n, and contains a bootstrap
 *   sample on output.
 */
void ran_bootstrap_sampling(const Random *ran, int n,
			    double *x, double *xb)
{
  int i;

  for (i = 0; i < n; i++) xb[i] = x[ran_intrandom(ran, n)];
}

/**
 * Generate a bootstrap sample.
 *
 * @param ran is an initialized Random numbers generators.
 * @param n is the sample size.
 * @param x is the original samble, which is not modified on output.
 * @param xb is a vector of size a least n, and contains a bootstrap
 *   sample on output.
 * @param f is vector with, on output, f[i] denoting the frequency of
 *   element x[i] in the bootstrap sample. 
 */
void ran_bootstrap_sampling_freq(const Random *ran, int n,
				 double *x, double *xb, double *f)
{
  int i, pos;

  for (i = 0; i < n; i++) f[i] = 0;

  for (i = 0; i < n; i++) {
    pos = ran_intrandom(ran, n);
    xb[i] = x[pos];
    f[pos]++;
  }

  for (i = 0; i < n; i++) f[i] /= n;
}
