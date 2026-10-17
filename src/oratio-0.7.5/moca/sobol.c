/**
 * sobol.c
 *
 * ORATIO Library
 * Copyright (C) 2004-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/moca.h"

#include <oratio/oratio.h>

/**
 * Sobol primitive polynomials. Source: Bratley and Fox.
 */

#define SOBOL_BF_MAXDIM   40
#define SOBOL_BF_MAXDIGIT 30
#define SOBOL_BF_MAXDEG   8

/**
 * Primitive polynomials in binary encoding.
 */
static const int primitive_polynomials[SOBOL_BF_MAXDIM] = 
{
  1,     3,   7,  11,  13,  19,  25,  37,  59,  47,
  61,   55,  41,  67,  97,  91, 109, 103, 115, 131,
  193, 137, 145, 143, 241, 157, 185, 167, 229, 171,
  213, 191, 253, 203, 211, 239, 247, 285, 369, 299
};

/**
 * Degrees of the primitive polynomials.
 */
static const int degree_table[SOBOL_BF_MAXDIM] = 
{
  0, 1, 2, 3, 3, 4, 4, 5, 5, 5,
  5, 5, 5, 6, 6, 6, 6, 6, 6, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
  7, 7, 7, 7, 7, 7, 7, 8, 8, 8
};


/**
 * Initial values for direction tables, following
 * Bratley+Fox, taken from [Sobol+Levitan, preprint 1976]
 */
static const int v_init[SOBOL_BF_MAXDEG][SOBOL_BF_MAXDIM] =
{
  {
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  },
  {
    0, 0, 1, 3, 1, 3, 1, 3, 3, 1,
    3, 1, 3, 1, 3, 1, 1, 3, 1, 3,
    1, 3, 1, 3, 3, 1, 3, 1, 3, 1,
    3, 1, 1, 3, 1, 3, 1, 3, 1, 3
  }, 
  {
    0, 0, 0, 7, 5, 1, 3, 3, 7, 5,
    5, 7, 7, 1, 3, 3, 7, 5, 1, 1,
    5, 3, 3, 1, 7, 5, 1, 3, 3, 7,
    5, 1, 1, 5, 7, 7, 5, 1, 3, 3
  }, 
  {
    0,  0,  0,  0,  0,  1,  7,  9, 13, 11,
    1,  3,  7,  9,  5, 13, 13, 11,  3, 15,
    5,  3, 15,  7,  9, 13,  9,  1, 11,  7,
    5, 15,  1, 15, 11,  5,  3,  1,  7,  9
  }, 
  {
     0,  0,  0,  0,  0,  0,  0,  9,  3, 27,
    15, 29, 21, 23, 19, 11, 25,  7, 13, 17,
     1, 25, 29,  3, 31, 11,  5, 23, 27, 19,
    21,  5,  1, 17, 13,  7, 15,  9, 31,  9
  }, 
  {
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0, 37, 33,  7,  5, 11, 39, 63,
    27, 17, 15, 23, 29,  3, 21, 13, 31, 25,
     9, 49, 33, 19, 29, 11, 19, 27, 15, 25
  }, 
  {
     0,   0,  0,  0,  0,  0,    0,  0,  0,   0,
     0,   0,  0,  0,  0,  0,    0,  0,  0,  13,
    33, 115, 41, 79, 17,  29, 119, 75, 73, 105,
     7,  59, 65, 21,  3, 113,  61, 89, 45, 107
  }, 
  {
    0, 0, 0, 0, 0, 0, 0, 0,  0,  0,
    0, 0, 0, 0, 0, 0, 0, 0,  0,  0,
    0, 0, 0, 0, 0, 0, 0, 0,  0,  0,
    0, 0, 0, 0, 0, 0, 0, 7, 23, 39
  }
};

/**
 * Compute the next vector in the sequence.
 *
 * x is the next generated vector.
 */
void moca_sobol_next(MOCA_Sobol *sbl, double *x)
{
  int k;
  int im;
  int ell = 0;

  im = sbl->in++;

  /* Find the rightmost 0 bit. */
  while (1) {
    ++ell;
    if (!(im & 1)) break;
    im >>=1;
  }

  if (ell > sbl->max_digits) {printf("Error in Sobol procedure (not enought digits).\n"); exit(1);};

  for (k = 0; k < sbl->s; k++) {
    sbl->ix[k] ^= sbl->vi[(ell-1)*SOBOL_BF_MAXDIM+k];
    x[k] = sbl->ix[k]*sbl->fac;
  }
}

void moca_sobol_reset(MOCA_Sobol *sbl)
{
  int k;

  sbl->in = 0;
  for(k = 0; k < sbl->max_dim; k++) sbl->ix[k]=0;
}

void moca_sobol_init(MOCA_Sobol *sbl)
{
  int i, j, k, deg;
  int ell;

  /* Initialize direction table in dimension 0. */
  for (k = 0; k < SOBOL_BF_MAXDIGIT; k++) sbl->vi[k*SOBOL_BF_MAXDIM] = 1;

  /* Initialize in remaining dimensions. */
  for (i = 1; i < sbl->s; i++) {

    int includ[sbl->maxdeg];

    /* Expand the polynomial bit pattern to separate
     * components of the logical array includ[]. */
    int p_i = primitive_polynomials[i];

    deg = degree_table[i];

    for (k = deg-1; k >= 0; k--) {
      includ[k] = ((p_i % 2) == 1);
      p_i /= 2;
    }

    /* Leading elements for dimension i come from v_init[][]. */
    for (j = 0; j < deg; j++)
      sbl->vi[j*SOBOL_BF_MAXDIM+i] = v_init[j][i];

    /* Calculate remaining elements for this dimension,
     * as explained in Bratley+Fox, section 2. */
    for (j = deg; j < SOBOL_BF_MAXDIGIT; j++) {
      int newv = sbl->vi[(j-deg)*SOBOL_BF_MAXDIM+i];
      ell = 1;
      for (k = 0; k <deg; k++) {
        ell *= 2;
        if(includ[k])
	  newv ^= (ell * sbl->vi[(j-k-1)*SOBOL_BF_MAXDIM+i]);
      }
      sbl->vi[j*SOBOL_BF_MAXDIM+i] = newv;
    }
  }

  /* Multiply columns of v by appropriate power of 2. */
  ell = 1;
  for (j = SOBOL_BF_MAXDIGIT-2; j >= 0; j--) {
    ell *= 2;
    for (i = 0; i < sbl->s; i++) {
      sbl->vi[j*SOBOL_BF_MAXDIM+i] *= ell;
    }
  }

  /* 1/(common denominator of the elements in vi) */
  sbl->fac = 1.0 / (2.0 * ell);

  //  sbl->fac = 1.0/(1L << sbl->max_digits);

  sbl->in = 0;
  for(k = 0; k < sbl->max_dim; k++) sbl->ix[k]=0;
}

/**
 * Create a new Sobol quasi-random numbers generator.
 *
 * @param s is the quasi-random vector dimension.
 */
MOCA_Sobol *moca_sobol_new(int pid, int s)
{
  MOCA_Sobol *sbl;

  sbl = _malloc(MOCA_Sobol);

  switch (pid) {
  case SOBOL_BRATLEY_FOX:
  default:
    sbl->max_dim = SOBOL_BF_MAXDIM;
    sbl->max_digits = SOBOL_BF_MAXDIGIT;
    sbl->ix = _mallocs(SOBOL_BF_MAXDIM, int);
    sbl->vi = _mallocs(SOBOL_BF_MAXDIGIT*SOBOL_BF_MAXDIM, int);
    sbl->maxdeg = SOBOL_BF_MAXDEG;
    sbl->mdeg = (int *)degree_table;
    break;
  }
  sbl->s = s;
  moca_sobol_init(sbl);

  return sbl;
}

/**
 * Free the memory associated to a Sobol quasi-random numbers generator.
 *
 * @param sbl is an initialized Sobol quasi-random numbers generator.
 */
void moca_sobol_free(MOCA_Sobol *sbl)
{
  if (sbl) {
    free(sbl->ix);
    free(sbl->vi);
    free(sbl);
  }
}

/**
 * Generate a Sobol sample.
 *
 * @param sobol is an initialized Sobol quasi-random numbers generator.
 * @param sample contains on output the Sobol sample.
 * @param n is the numbers of vectors forming the sample.
 * @param s is the quasi-random vector dimension.
 */
void moca_sobol_sample(MOCA_Sobol *sobol, int n, int s,
		       double *S, int ldS)
{
  int i;
  
  for (i = 0; i < n; i++) moca_sobol_next(sobol, S+i*s);
}
