/**
 * nt_io.h
 *
 * ORATIO Library
 * Copyright (C) 2002-2011 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_NTOOLS_IO_H__
#define __ORATIO_NTOOLS_IO_H__

#include <stdio.h>
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Matrix input and output
   */

  /**
   * Read the matrix A in the file in.
   *
   * @param in is the input file, containing the matrix to read.
   * @param name contains on output the matrix name.
   * @param A is on input a pointer to a NULL matrix or a matrix that will be
   *   destroyed, and on output a pointer to a newly created containing the
   *   read elements.
   * @param m is the number of rows of the matrix.
   * @param n is the number of columns of the matrix.
   */
  int nt_matrix_read(FILE *in, char **name,
		     double ***A, int *nrow, int *ncol);

  /**
   * Print the matrix A in the file out.
   *
   * @param out is the file where the matrix has to be written.
   * @param m is the number of rows of the matrix.
   * @param n is the number of columns of the matrix.
   */
  void nt_matrix_print(FILE *out, char *name,
		       double **A, int n, int m);

  /**
   * Print a vector on a file.
   *
   * @param out is the file where the vector has to be written.
   * @param x is the vector to print.
   * @param n is dimension of x.
   */
  void nt_vect_print(FILE *out, char *name,
		     double *v, int n);

  /**
   *
   * @param out is the file where the matrix has to be written.
   * @param m is the number of rows of the matrix.
   * @param n is the number of columns of the matrix.
   */
  void nt_matrix_print2(FILE *out, char *name,
			double *A, int n, int m,
			int lda);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_NTOOLS_IO_H__ */
