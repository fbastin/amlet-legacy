/**
 * norm.c
 *
 * NTools Library
 * Copyright (C) 2002 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <ntools/norm.h>

#include <math.h>

double nt_vect_nrm1(int n, double *x, int inc) 
{
  int i;
  double nrm = 0.0;

  for (i = 0; i < n; i++) {
    nrm += fabs(*x);
    x += inc;
  }

  return nrm;
}
