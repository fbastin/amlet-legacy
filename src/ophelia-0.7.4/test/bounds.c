/**
 * bounds.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2004 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "bounds.h"

/**
 *
 * Willi Hock and Klaus Schittkowski, "Test Examples for Nonlinear
 * Programming Codes", Springer-Verlag, 1981.
 *
 */

/**
 * PBR_T1_1
 */

double pbr_t1_1(double *x, int n, void *data)
{
  double t, t2;

  t = x[1]-x[0]*x[0];
  t2 = 1-x[0];

  return 100.0*t*t + t2*t2;
}

void pbr_t1_1_set_bounds(double *b)
{
  b[2] = -1.5;
}
