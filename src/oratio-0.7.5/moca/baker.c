/**
 * baker.c
 *
 * ORATIO Library
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/moca.h"

/**
 * Apply Baker's transformation to the scalar u.
 *
 * @param u is the scalar to transform; u must be in [0,1).
 *
 * @param return the transformed scalar.
 */
double moca_baker(double u)
{
  if (u < 0.5) return 2.0*u;
  else return 2.0*(1-u);
}

/**
 * Apply Baker's transformation to vector u.
 *
 * @param s is the dimension of vector u.
 * @param u is a vector in [0,1)^s, and is the vector after Baker's
 *   transformation on output.
 * @param inc is the increment of u.
 */
void moca_baker_point(int s, double *u, int incx)
{
  int i;

  for (i = 0; i < s; i++) u[i*incx] = moca_baker(u[i*incx]);
}
