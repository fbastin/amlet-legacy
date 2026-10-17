/**
 * values.c
 *
 * ORATIO Library
 * Copyright (C) 2005-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/ntools.h>

#include <math.h>

/**
 * Return the infinity representation.
 */
double nt_infinity()
{
#ifdef INFINITY
  return INFINITY;
#else
  double a=0.0;
  return 1.0/a;
#endif
}

/**
 * Return the Not a Number representation.
 *
 * This representation is architecture dependent.
 * Be careful: this will produce a floating exception on some
 * architectures, in particular when no direct double representation
 * is available.
 */
double nt_nan()
{
#ifdef NAN
  return NAN;
#else
  double a = 0.0;
  return 0.0/a;
#endif
}
