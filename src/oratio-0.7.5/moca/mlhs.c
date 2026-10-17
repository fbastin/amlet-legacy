/**
 * mlhs.c
 *
 * ORATIO Library
 * Copyright (C) 2004-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "mlhs.h"

#include "oratio/moca.h"

void mlhs_vectors_ind(Random *r, double *pts, int n)
{
  int j;
  double v;

  for (j = 0; j < n; j++) pts[j] = ((double)j)/n;

  /* Shift the sequence */
  v = ran_random_get_val(r)/n;
  for (j = 0; j < n; j++) pts[j] += v;

  /* Shuffle the sequence in order to break correlations
     between dimensions */
  ran_shuffle_double(r, n, pts, 1);
}

/**
 * @param ran is an initialized random numbers generator.
 */
void mlhs_generate_sample(Random *ran, double **sample,
			  int n, int s)
{
  /* Array of draws, for one individual */
  double **draws;

  int j, k;

  draws = nt_matrix_new(s, n);

  for (j = 0; j < s; j++) {
    mlhs_vectors_ind(ran, draws[j], n);
    if (sample) {
      for (k = 0; k < n; k++)
	sample[k][j] = draws[j][k];
    }
  }

  nt_matrix_free(draws);
}
