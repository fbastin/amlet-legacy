/**
 * rho.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2004 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the MIT License
 */

#include <float.h>

#define NLP_RHO_MULTIPLIER 10

/**
 * Robust rho computation.
 * Ref. Conn, Gould and Toint, "Trust-Region Methods", SIAM, 2000,
 * Section 17.4.2, page 793.
 */
double robust_rho()
{
  //  double shift = NLP_RHO_MULTIPLIER*DBL_EPSILON;

  return DBL_EPSILON;
}
