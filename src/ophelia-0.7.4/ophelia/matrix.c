/**
 * matrix.c
 *
 * OPHELIA Library
 * Copyright (C) 2005-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/ctools.h>

#include <oratio/f2c.h>
#include <oratio/cblas.h>

#include "clapack.h"

#include "matrix.h"

#include <stdlib.h>

/**
 * Copy the off-diagonal lower part of a square matrix in the off-diagonal
 * upper-part.
 *
 * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
 * @param n is the matrix order.
 * @param A is a real square matrix nxn.
 * @param ldA is the leading dimension of A.
 */
void op_matrix_copy_lower(const enum CBLAS_ORDER Order,
				 int n, double *A, int ldA)
{
  int i, j;

  if (Order == CblasRowMajor) {
    for (i = 0; i < n; i++)
      for (j = i+1; j < n; j++)
	A[i*ldA+j] = A[j*ldA+i];
  }
  else
    op_matrix_copy_upper(CblasColMajor, n, A, ldA);;
}

/**
 * Copy the off-diagonal lower part of a square matrix in the off-diagonal
 * upper-part.
 *
 * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
 * @param n is the matrix order.
 * @param A is a real square matrix nxn.
 * @param ldA is the leading dimension of A.
 */
void op_matrix_copy_upper(const enum CBLAS_ORDER Order,
				 int n, double *A, int ldA)
{
  int i, j;

  if (Order == CblasRowMajor) {
    for (i = 1; i < n; i++)
      for (j = 0; j < i; j++)
	A[i*ldA+j] = A[j*ldA+i];
  }
  else
    op_matrix_copy_lower(CblasColMajor, n, A, ldA);;
}

/**
 * Compute the inverse of a symmetric indefinite square matrix.
 *
 * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
 * @param Uplo designs which half-matrix is filled (CblasUpper or
 *   CblasLower).
 * @param I is the matrix to invert. The storage on input and output
 *   is given by Order and Uplo.
 * @param npar is the matrix main dimension.
 *
 * @return the information value delivered but Lapack routines used to
 *   compute the inverse.
 */
int op_matrix_syin_inverse(const enum CBLAS_ORDER Order,
			   const enum CBLAS_UPLO Uplo,
			   double *I, int npar)
{
  char uplo;
  integer info, n;

  integer *ipiv;
  double *work;
  integer lwork;

  n = npar;

  switch (Order) {
  case CblasRowMajor:
    uplo =  (Uplo == CblasLower) ? 'U' : 'L';
    break;
  case CblasColMajor:
    uplo =  (Uplo == CblasLower) ? 'L' : 'U';
    break;
  }

  ipiv = _mallocs(npar, integer);
  info = 1;
  ipiv[0] = -1;
  lwork = ilaenv_(&info, "dsytrf", &uplo, &n, ipiv, ipiv, ipiv, 6, 1);
  lwork *= n;
  work = _mallocs(lwork, double);
  dsytrf_(&uplo, &n, I, &n, ipiv, work, &lwork, &info);
  dsytri_(&uplo, &n, I, &n, ipiv, work, &info);
  free(work);
  free(ipiv);

  return info;
}

/**
 * Compute the inverse of a symmetric positive definite square matrix.
 *
 * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
 * @param Uplo designs which half-matrix is filled (CblasUpper or
 *   CblasLower).
 * @param I is the matrix to invert. The storage on input and output
 *   is given by Order and Uplo.
 * @param npar is the matrix main dimension.
 *
 * @return the information value delivered but Lapack routines used to
 *   compute the inverse.
 */
int op_matrix_sydp_inverse(const enum CBLAS_ORDER Order,
			   const enum CBLAS_UPLO Uplo,
			   double *I, int npar)
{
  char uplo;
  integer info, n;

  n = npar;

  switch (Order) {
  case CblasRowMajor:
    uplo =  (Uplo == CblasLower) ? 'U' : 'L';
    break;
  case CblasColMajor:
    uplo =  (Uplo == CblasLower) ? 'L' : 'U';
    break;
  }

  dpotrf_(&uplo, &n, I, &n, &info);
  dpotri_(&uplo, &n, I, &n, &info);

  return info;
}


/**
 * Inverse a real symmetric matrix.
 *
 * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
 * @param Uplo designs which half-matrix is filled (CblasUpper or
 *   CblasLower).
 * @param I is the matrix to invert. The storage on input and output
 *   is given by Order and Uplo.
 * @param npar is the matrix main dimension.
 *
 * @return the information value delivered but Lapack routines used to
 *   compute the inverse.
 */
int op_matrix_sy_inverse(const enum CBLAS_ORDER Order,
			 const enum CBLAS_UPLO Uplo,
			 double *I, int npar)
{
  char uplo;
  integer info, n;
  double diag[npar];

  int i;

  n = npar;  /* n is used only to ensure that we work with the same
		data types */

  /* We follow Fortran column orientation */
  switch (Order) {
  case CblasRowMajor:
    uplo =  (Uplo == CblasLower) ? 'U' : 'L';
    break;
  case CblasColMajor:
    uplo =  (Uplo == CblasLower) ? 'L' : 'U';
    break;
  }

  /* Copy the diagonal. */
  for (i = 0; i < npar; i++) diag[i] = I[i*(npar+1)];

  /* We try first to apply the Cholesky factorization */ 
  dpotrf_(&uplo, &n, I, &n, &info);
  if (!info) {
    dpotri_(&uplo, &n, I, &n, &info);
    return info;
  }
  else {
    /* The matrix is not positive definite. */
    integer *ipiv;
    double *work;
    integer lwork;

    switch (Order) {
    case CblasRowMajor:
      uplo =  (Uplo == CblasLower) ? 'U' : 'L';
      op_matrix_copy_lower(Order, npar, I, npar);
      break;
    case CblasColMajor:
      op_matrix_copy_upper(Order, npar, I, npar);
      break;
    }
    for (i = 0; i < npar; i++)  I[i*(npar+1)] = diag[i];

    ipiv = _mallocs(npar, integer);
    info = 1;
    ipiv[0] = -1;
    lwork = ilaenv_(&info, "dsytrf", &uplo, &n, ipiv, ipiv, ipiv, 6, 1);
    lwork *= n;
    work = _mallocs(lwork, double);
    dsytrf_(&uplo, &n, I, &n, ipiv, work, &lwork, &info);
    dsytri_(&uplo, &n, I, &n, ipiv, work, &info);
    free(work);
    free(ipiv);
  }

  return info;
}
