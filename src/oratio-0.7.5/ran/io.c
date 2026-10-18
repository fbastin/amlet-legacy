/**
 * io.c
 *
 * ORATIO Library
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/ran.h"

#include <stdio.h>

/**
 * Write on output file a vector of random draws.
 *
 * @param f is the file where the draws have to be written; if NULL,
 *   the output occurs on stdout.
 * @param v is the realizations vector to write.
 * @param d is the dimension of the random vector.
 */
void ran_write_vector(FILE *f, double *v, int d)
{
  int i;

  if (!f) f = stdout;
  /* 17 significant digits, so that the draws are read back exactly; ORATIO
     <= 0.7.5 wrote 16 */
  for (i = 0; i < d-1; i++) fprintf(f, "%24.17g ", v[i]);
  fprintf(f, "%24.17g\n", v[d-1]);
}

/**
 * Read vectors of draws stored in a file.
 *
 * @param f is the file where the draws have to be read.
 * @param n is the number of vectors to read.
 * @param s it the dimension of each vector to read.
 * @param S is the array where the draws will be saved.
 *
 * @return RAN_ERROR_READ_DRAW is some draw ccould not be read, 0 otherwise.
 */
int ran_import_sample(FILE *f, int n, int s, double **S)
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < s; j++)
      if (fscanf(f, "%lg", S[i]+j) != 1) return RAN_ERROR_READ_DRAW;

  return 0;
}

/**
 * Export a draws sample array.
 *
 * @param f represents the file where the draws have to be saved.
 * @param m is the number of generated random vectors.
 * @param d is the dimension of the random vectors.
 * @param S is the sample array.
 * @param ldS is the leading dimension of the draws sample array S.
 */
void ran_export_sample(FILE *f, int m, int d, double *S, int ldS)
{
  int i;

  for (i = 0; i < m; i++) ran_write_vector(f, S+i*ldS, d);
}
