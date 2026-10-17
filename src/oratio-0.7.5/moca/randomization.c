/**
 * randomization.c
 *
 * Monte Carlo and Quasi-Monte Carlo Library
 * Copyright (C) 2005-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/moca.h"
#include "oratio/cblas.h"

/**
 * Randomly shift a vector of points, modulo one.
 *
 * @param r is a pseudo-random numbers generator (must be initialized).
 * @param pts is the vector of points to shift. On output, each
 *   component of pts has been shifted by some random value (identical
 *   between components), modulo one.
 * @param n is the size of pts.
 */
void moca_shift_modulo_one(Random *r, double *pts, int n)
{
  double v;
  int i;

  v = ran_random_get_val(r);

  for (i = 0; i < n; i++) {
    pts[i] += v;
    if (pts[i] >= 1.0) pts[i] -= 1.0;
  }
}

/**
 * @param s is the dimension of x and the shift vector.
 */
void moca_shift_apply(int s, double *x, double *shift)
{
  int i;

  cblas_daxpy(s, 1.0, shift, 1, x, 1);

  /* Apply the modulo one operation. */
  for (i = 0; i < s; i++) if (x[i] >= 1.0) x[i] -= 1.0;
}

/**
 * Initialize a random shift vector with a padding for 32-bit random
 * numbers generators.
 *
 * @param r is a pseudo-random numbers generator (must be
 *   initialized).
 * @param s is the dimension of the shift vector.
 * @param shift contains on output the random shift vector.
 */
void moca_shift_init32(Random *r, int s, double *shift)
{
  int i;

  for (i = 0; i < s; i++) shift[i] = ran_random_get_val(r)+1.0/(4.0*(1 << 31));
}

/**
 * Initialize a random shift vector with a padding for 32-bit random
 * numbers generators; each component is set to the same number.
 *
 * @param r is a pseudo-random numbers generator (must be
 *   initialized).
 * @param s is the dimension of the shift vector.
 * @param shift contains on output the random shift vector.
 */
void moca_shift_uniform_init32(Random *r, int s, double *shift)
{
  int i;
  double sh = ran_random_get_val(r)+1.0/(4.0*(1 << 31));

  for (i = 0; i < s; i++) shift[i] = sh;
}

/**
 * Initialize a random shift vector.
 *
 * @param r is a pseudo-random numbers generator (must be
 *   initialized).
 * @param s is the dimension of the shift vector.
 * @param shift contains on output the random shift vector.
 */
void moca_shift_init(Random *r, int s, double *shift)
{
  int i;

  for (i = 0; i < s; i++) shift[i] = ran_random_get_val(r);
}
