/**
 * permutations.c
 *
 * ORATIO Library
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/ran.h"
#include "oratio/dstruct.h"

#include <string.h>
#define SWAP(a,b) { swap=(a); (a)=(b); (b)=swap; }

/**
 * Shuffle a generic pointers vector using the Fisher-Yates shuffle.
 *
 * @param r is pseudo-random numbers generator (must be initialized).
 * @param x is the vector to shuffle.
 * @param n is the dimension of x.
 */
void ran_shuffle_generic(Random *r, int n, void *x, int incx,
				int size)
{
  /**
   * Ref.:
   * D. E. Knuth: The Art of Computer Programming, Volume 2,
   * Third edition. Section 3.4.2, Algorithm P, pp 145. Reading:
   * Addison-Wesley, 1997. ISBN: 0-201-89684-2.
   *
   * R. A. Fisher and F. Yates: Statistical Tables. London, 1938.
   * Example 12.
   */

  int i, j;
  void *swap = c_malloc(size);

  for (i = n-1; i > 0; i--) {
    j = ran_intrandom(r, i+1);
    memcpy(swap, x+i*size, size);
    memcpy(x+i*size, x+j*size, size);
    memcpy(x+j*size, swap, size);
  }
}

/**
 * Shuffle a double vector using the Fisher-Yates shuffle.
 *
 * @param r is pseudo-random numbers generator (must be initialized).
 * @param x is the vector to shuffle.
 * @param n is the dimension of x.
 */
void ran_shuffle_double(Random *r, int n, double *x, int incx)
{
  int i, j;
  double swap;

  for (i = n-1; i > 0; i--) {
    j = ran_intrandom(r, i+1);
    SWAP(x[i], x[j])
  }
}
