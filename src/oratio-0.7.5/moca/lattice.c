/**
 * lattice.c
 *
 * ORATIO Library
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/moca.h"
#include "oratio/oratio.h"

/**
 * Rank-1 lattice.
 */

/**
 * Create a new rank-1 lattice.
 *
 * @param n is the number of points.
 * @param s dimension of the lattice vector a.
 * @param a the lattice vector; each component must be in [0,n[.
 *
 * @return a newly memory-allocated rank-1 lattice.
 */
MOCA_Rank1Lattice *moca_rank1_lattice_new(int n, int s, int *a)
{
  MOCA_Rank1Lattice *lattice = _malloc(MOCA_Rank1Lattice);

  int i;

  lattice->s = s;
  lattice->n = n;
  lattice->normfactor = 1.0/(double)n;
  lattice->v = _mallocs(s, double);

  for (i = 0; i < s; i++) lattice->v[i] = lattice->normfactor*a[i];

  return lattice;
}

/**
 * Free the memory associated to a rank-1 lattice->
 *
 * @param lattice is a memory-allocated rank-1 lattice->
 */
void moca_rank1_lattice_free(MOCA_Rank1Lattice *lattice)
{
  if (lattice) {
    free(lattice->v);
    free(lattice);
  }
}

/**
 * Compute the jth coordinate of the ith rank-1 lattice point.
 *
 * @param lattice is an initialixed rank-1 lattice.
 * @param i is the point index.
 * @param is the coordinate index.
 *
 * @return the jth coordinate of the ith rank-1 lattice point.
 */
double moca_rank1_lattice_get_ij(MOCA_Rank1Lattice *lattice,
				 int i, int j)
{
  int k;
  double x;

  x = i*lattice->v[j];
  k = x;

  return x-k;
}

/**
 * Compute the current coordinate of the current lattice point.
 *
 * @param lattice is an initialixed rank-1 lattice->
 *
 * @return the next coordinate of the current lattice point.
 */
double moca_rank1_lattice_get_coordinate(MOCA_Rank1Lattice *lattice)
{
  int k;
  double x;

  x = lattice->ptindex*lattice->v[lattice->coordindex];
  k = x;

  return x-k;

}

/**
 * Compute the next coordinate of the current lattice point.
 *
 * @param lattice is an initialixed rank-1 lattice->
 *
 * @return the next coordinate of the current lattice point.
 */
double moca_rank1_lattice_next_coordinate(MOCA_Rank1Lattice *lattice)
{
  int k;
  double x;

  x = lattice->ptindex*lattice->v[(lattice->coordindex)++];
  k = x;

  return x-k;
}

/**
 * Compute the point of a rank-1 lattice at the current point index.
 *
 * @param lattice is an initialixed rank-1 lattice->
 * @param x is the vector that will contain the point.
 * @param incx is the increment of x.
 */
void moca_rank1_lattice_get_point(MOCA_Rank1Lattice *lattice,
				  double *x, int incx)
{
  int i;

  lattice->coordindex = 0;
  for (i = 0; i < lattice->s; i++)
    x[i*incx] = moca_rank1_lattice_next_coordinate(lattice);
}


/**
 * Compute the point of a rank-1 lattice at the next point index,
 * which is incremented (after point computation).
 *
 * @param lattice is an initialixed rank-1 lattice->
 * @param x is the vector that will contain the point.
 * @param incx is the increment of x.
 */
void moca_rank1_lattice_next(MOCA_Rank1Lattice *lattice,
			     double *x, int incx)
{
  moca_rank1_lattice_get_point(lattice, x, incx);
  lattice->ptindex++;
}

/**
 * Increment the point index of a rank-1 lattice->
 *
 * @param lattice is an initialixed rank-1 lattice->
 */
void moca_rank1_lattice_inc_point(MOCA_Rank1Lattice *lattice)
{
  lattice->ptindex++;
}

/**
 * Reset the current rank-1.
 *
 * @param lattice is an initialixed rank-1 lattice.
 */
void moca_rank1_lattice_reset(MOCA_Rank1Lattice *lattice)
{
  lattice->ptindex = lattice->coordindex = 0;
}

/**
 * Generate a rank-1 lattice sample.
 *
 * @param lattice is a rank-1 lattice generator.
 * @param n is the numbers of vectors forming the sample.
 * @param S contains on output the lattice sample.
 * @param ldS is the leading dimension of S; ldS must be grater or equal to
 *   lattice points dimension, and is typically equal to this dimesion.
 */
void moca_rank1_lattice_sample(MOCA_Rank1Lattice *lattice,
			       int n, double *S, int ldS)
{
  int i;
  
  moca_rank1_lattice_reset(lattice);
  for (i = 0; i < n; i++) moca_rank1_lattice_next(lattice, S+i*ldS, 1);
}
