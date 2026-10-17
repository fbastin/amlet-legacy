/**
 * ratio.c
 *
 * ORATIO Library
 * Copyright (C) 2008-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/ran.h"

/**
 * Draw from a ratio of two random distributions.
 *
 * @param ran1 is an initialized random numbers generator.
 * @param ran2 is an initialized random numbers generator, independent
 *   of ran1.
 *
 * @return a draw from a ratio of two random distributions.
 */
double ran_generate_ratio(Random *ran1, Random *ran2,
			  double *x1, double *x2)
{
    *x1 = ran_random_next(ran1);
    *x2 = ran_random_next(ran2);

    return (*x1)/(*x2);
}

/**
 * Generate realizations from a ratio of two random distributions.
 *
 * @param ran1 is an initialized random numbers generator.
 * @param ran2 is an initialized random numbers generator, independent
 *   of ran1.
 * @param n is the number of realizations to draw.
 * @param x1 is a double array, of size at least equal to n, used to
 *   store draws produced by ran1.
 * @param incx1 is the increment to use in x1.
 * @param x2 is a double array, of size at least equal to n, used to
 *   store draws produced by ran2.
 * @param incx2 is the increment to use in x2.
 * @param r is a double array, of size at least equal to n, that
 *   will contain the ratios between elements of x1 and x2.
 * @param incr is the increment to use in r.
 */
void ran_generate_ratio_distribution(Random *ran1, Random *ran2,
				     int n, double *x1, int incx1,
				     double *x2, int incx2,
				     double *r, int incr)
{
  int i;

  for (i = 0; i < n; i++)
    r[incr*i] = ran_generate_ratio(ran1, ran2, x1+incx1*i, x2+incx2*i);
}
