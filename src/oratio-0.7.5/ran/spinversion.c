/**
 * spinversion.c
 *
 * Random Number Library
 * Copyright (C) 2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

/**
 * Random numbers generation by semi-parametric inversion.
 */

#include "oratio/ran.h"
#include "oratio/ctools.h"

/**
 * Draw pseudo-random numbers by using semi-parametric inversion.
 *
 * @param ran is an initialized random numbers generator.
 * @param n is the number of values to draw.
 * @param icdf is the inverse cumulative distribution function.
 * @param type is the type identifier of idf.
 */
void ran_sp_inversion(Random *ran, int n, double *draws,
		      C_GENERIC *icdf, int type,
		      int npars, int *pars)
{
  int i;

  /* Generate uniform draws on (0,1). */
  for (i = 0; i < n; i++) draws[i] = ran_random_get_val(ran);

  /* Apply the inversion. */
  switch (type) {
  case RAN_SP_SPLINE:
    break;
  default:
    break;
  }
}
