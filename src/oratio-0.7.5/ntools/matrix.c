/**
 * matrix.c
 *
 * ORATIO Library
 * Copyright (C) 2002-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/ctools.h"
#include "oratio/ntools.h"

#include <stdlib.h>
#include <string.h>

/**
 * Matrix definitions
 */


/**
 * Create a new matrix of dimensions m*n.
 *
 * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
 * @param m is the number of rows of the matrix to create.
 * @param n is the number of rows of the matrix to create.
 *
 * @return a pointer to the newly allocated matrix. The matrix is
 *   stored is contiguous memory, and its orientation is defined by Order.
 */
double **nt_ma_new(const enum CBLAS_ORDER Order, int m, int n)
{
  double **A;
  int ldA, sdA; /* leading and secundary dimensions */
  int i;

  /* Allocate the lines pointer. */
  ldA = (Order == CblasRowMajor) ? m : n;
  sdA = (Order == CblasRowMajor) ? n : m;

  A = _mallocs(ldA, double *);

  A[0] = _mallocs(m*n, double);

  /* Adjust the pointers. */
  for (i = 1; i < ldA; i++) A[i] = A[i-1]+sdA;

  return A;
}

/**
 * Create a new matrix of dimensions m*n.
 *
 * @param m is the number of rows of the matrix to create.
 * @param n is the number of columns of the matrix to create.
 *
 * @return a pointer to the newly allocated matrix. The matrix is
 *   stored is contiguous memory.
 */
double **nt_matrix_new(int m, int n)
{
  return nt_ma_new(CblasRowMajor, m, n);
}

/**
 * Set A to the zeros matrix.
 *
 * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
 * @param m is the number of rows of matrix A.
 * @param n is the number of columns of matrix A.
 * @param  A: before entry, the leading m by n part of the array A must
 *   contain the matrix of coefficients. On exit, the leading m by n
 *   part of the array A contains the zeros matrix.
 * @param ldA specifies the first dimension of A as declared in the
 *   calling (sub) program. ldA must be at least max( 1, m ) if Order
 *   is equal to CblasColMajor, or max (1, n) if Order is equal to
 *   CblasRowMajor.
 */
void nt_ma_zeros(const enum CBLAS_ORDER Order,
		 const int m, const int n, double *A, const int ldA)
{
  int i, j;

  switch(Order) {
  case CblasColMajor:
    for (j = 0; j < n; j++)
      for (i = 0; i < m; i++) A[i+j*ldA] = 0.0;
    break;
  case CblasRowMajor:
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++) A[i*ldA+j] = 0.0;
    break;
  default:
    c_warning("invalid Order argument in nt_ma_zeros");
  }
}

/**
 * Set A to a diagonal matrix.
 *
 * @param Order: CblasRowMajor or CblasColMajor
 * @param m is the number of rows of matrix A.
 * @param n is the number of columns of matrix A.
 * @param  A: before entry, the leading m by n part of the array A must
 *   contain the matrix of coefficients. On exit, the leading m by n
 *   part of the array A contains the diagonal matrix with entries
 *   equal to alpha.
 * @param ldA specifies the first dimension of A as declared in the
 *   calling (sub) program. ldA must be at least max( 1, m ).
 * @param alha is some real number.
 */
void nt_ma_diag(const enum CBLAS_ORDER Order,
		const int m, const int n, double *A, const int ldA,
		double alpha)
{
  int i;

  nt_ma_zeros(Order, m, n, A, ldA);

  for (i = 0; i < nt_min(m,n); i++) A[i*(ldA+1)] = alpha;
}

/**
 * Set A to the identity matrix.
 *
 * @param Order: CblasRowMajor or CblasColMajor
 * @param m is the number of rows of matrix A.
 * @param n is the number of columns of matrix A.
 * @param A: before entry, the leading m by n part of the array A must
 *   contain the matrix of coefficients. On exit, the leading m by n
 *   part of the array A contains the identity matrix.
 * @param ldA specifies the first dimension of A as declared in the
 *   calling (sub) program. lda must be at least max( 1, m ).
 */
void nt_ma_identity(const enum CBLAS_ORDER Order,
		    const int m, const int n, double *A, const int ldA)
{
  nt_ma_diag(Order, m, n, A, ldA, 1.0);
}

/**
 * Scale the matrix A.
 *
 * @param Order: CblasRowMajor or CblasColMajor
 * @param m is the number of rows of matrix A.
 * @param n is the number of columns of matrix A.
 * @param alpha is the scaling factor.
 * @param  A: before entry, the leading m by n part of the array A must
 *   contain the matrix of coefficients. On exit, the leading m by n
 *   part of the array A contains the identity matrix.
 * @param ldA specifies the first dimension of A as declared in the
 *   calling (sub) program. ldA must be at least max( 1, m ).
 */
void nt_ma_scale(const enum CBLAS_ORDER Order, const int m, const int n,
		 double alpha, double *A, int ldA)
{
  int i, j;

  switch(Order) {
  case CblasRowMajor:
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++) A[i*ldA+j] *= alpha;
    break;
  case CblasColMajor:
    for (j = 0; j < n; j++)
      for (i = 0; i < m; i++) A[i+j*ldA] *= alpha;
    break;
  default:
    c_warning("invalid Order argument in nt_ma_scale");
  }
}

/**
 * Set the triangular lower part of a matrix based of a defining vector.
 *
 * @param Order: CblasRowMajor or CblasColMajor
 * @param m is the number of rows of matrix A.
 * @param n is the number of columns of matrix A.
 * @param A is the matrix under consideration.
 * @param ldA is the leading dimension of A.
 * @param x is a vector containing the elements defining the lower
 *   triangular part of a square matrix; it contains the elements row by
 *   row if Order is set to CblasRowMajor, and column by column if Order
 *   is set to CblasColMajor. If some elements of x falls outside the
 *   definition range of A, they are skipped.
 * @param incx is the increment of x.
 */
void nt_ma_set_L(const enum CBLAS_ORDER Order, const int m, const int n,
		 double *A, int ldA, double *x, int incx)
{
  int i, j, k;

  k = 0;

  switch(Order) {
  case CblasRowMajor:
    for (i = 0; i < m; i++)
      for (j = 0; j <= i; j++) {
	if (j < n) A[i*ldA+j] = x[k*incx];
	k++;
      }
    break;
  case CblasColMajor:
    for (i = 0; i < n; i++)
      for (j = i; j < m; j++)
	A[i*ldA+j] = x[(k++)*incx];
    break;
  default:
    c_warning("invalid Order argument in nt_ma_scale");
  }
}

/**
 * Set the triangular lower part of a matrix with unit diagonal,
 * on the basis of a defining vector for the off-diagonal elements.
 *
 * @param Order: CblasRowMajor or CblasColMajor
 * @param m is the number of rows of matrix A.
 * @param n is the number of columns of matrix A.
 * @param A is the matrix under consideration.
 * @param ldA is the leading dimension of A.
 * @param x is a vector containing the elements defining the lower
 *   triangular part of a square matrix; it contains the elements row by
 *   row if Order is set to CblasRowMajor, and column by column if Order
 *   is set to CblasColMajor. If some elements of x falls outside the
 *   definition range of A, they are skipped.
 * @param incx is the increment of x.
 */
void nt_ma_set_L_unitary(const enum CBLAS_ORDER Order,
			 const int m, const int n,
			 double *A, int ldA, double *x, int incx)
{
  int i, j, k;

  k = 0;

  nt_matrix_identity(m, n, A, ldA);
  switch(Order) {
  case CblasRowMajor:
    for (i = 0; i < m; i++)
      for (j = 0; j < i; j++) {
	if (j < n) A[i*ldA+j] = x[k*incx];
	k++;
      }
    break;
  case CblasColMajor:
    for (i = 0; i < n; i++)
      for (j = i+1; j < m; j++)
	A[i*ldA+j] = x[(k++)*incx];
    break;
  default:
    c_warning("invalid Order argument in nt_ma_scale");
  }
}

/**
 * Reallocate the memory associated to the matrix A.
 *
 * @param m is the number of rows of the matrix.
 * @param n is the number of rows of the matrix.
 */
double **nt_matrix_realloc(double **A, int m, int n)
{
  int i;

  if (!A) return nt_matrix_new(m, n);

  A[0] = _reallocs(A[0], m*n, double);
  A = _reallocs(A, m, double *);
  for (i = 1; i < m; i++) A[i] = A[i-1]+n;

  return A;
}

/**
 * Free the memory allocated to the matrix A.
 *
 * @param A is a previously allocated matrix.
 */
void nt_matrix_free(double **A)
{
  if (A) {
    free(A[0]);
    free(A);
  }
}

/**
 * Free the memory allocated to the matrix A.
 *
 * @param A is a previously allocated matrix.
 */
void nt_ma_free(double **A)
{
  if (A) {
    free(A[0]);
    free(A);
  }
}

/**
 * Duplicate a matrix into another place memory
 *
 * @param nrow is the number of rows of the matrix.
 * @param ncol is the number of cols of the matrix.
 */
double **nt_matrix_duplicate(const double **m, int nrow, int ncol)
{
  double **copy;

  copy = nt_matrix_new(nrow, ncol);

  memcpy(copy[0], m[0], nrow*ncol*sizeof(double));

  return copy;
}

/**
 * Set all components of a mxn matrix to the value 0.0.
 *
 * @param m is the number of rows of the matrix A.
 * @param n is the number of columns of the matrix A.
 * @param A is a pointer to a matrix of dimensions (at least) mxn. On
 *   output all elements of the principal submatrix of dimension mxn of
 *   A are set to 0.0.
 * @param ldA is the leading dimensions of A.
 */
void nt_matrix_zeros(int m, int n, double *A, int ldA)
{
  int i, j;

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      A[i*ldA+j] = 0.0;
}

/**
 * Set a mxn matrix to the identity matrix.
 *
 * @param m is the number of rows of the matrix A.
 * @param n is the number of columns of the matrix A.
 * @param A is a pointer to a matrix of dimensions (at least) mxn. On
 * the principal submatrix of dimension mxn of A is the identity matrix.
 * @param ldA is the leading dimensions of A.
 */
void nt_matrix_identity(int m, int n, double *A, int ldA)
{
  int i;

  nt_matrix_zeros(m, n, A, ldA);
  for (i = 0; i < nt_min(n, m); i++) A[i*(ldA+1)] = 1.0;
}

/**
 *
 * @param n is the number of rows of the matrix.
 * @param m is the number of cols of the matrix.
 */
void nt_matrix_insert_vector(double **A, double x[], int n, int m,
			     int orientation, int pos)
{
  if (orientation != NT_MATRIX_ROW && orientation != NT_MATRIX_COLUMN)
    {
      c_warning("** NTOOLS WARNING: incorrect orientation in routine "\
		"nt_matrix_2_vector, using NT_ROW instead");
      orientation = NT_MATRIX_ROW;
    }

  if (orientation == NT_MATRIX_ROW)
    memcpy(A[pos], x, m*sizeof(double));
  else
    memcpy(A[pos], x, n*sizeof(double));
}

/**
 * Insert a matrix B into a matrix A
 *
 * #PRE:  A: matrix of dimension n1xm1
 *        B: matrix of dimension n2xm2
 *        (a,b): position where the matrix B is to inserted; the
 *               indexes are supposed to start from 0
 *        Both matrices are supposed to have the same orientation;
 *        and the indexes reflect this orientation. In other terms
 *        if the matrices are oriented by columns, the first index
 *        indicates the column, while is the orientation is by rows,
 *        it reflects the row
 *
 * #POST: if the dimensions agree, B is inserted in A at the
 *        position (a,b)
 */
int nt_matrix_insert_matrix(double **A, int n1, int m1,
			    double **B, int n2, int m2,
			    int a, int b)
{
  int i;

  if (a + n2 > n1) return NT_MATRIX_BAD_DIM1;
  if (b + m2 > m1) return NT_MATRIX_BAD_DIM2;

  for (i = 0; i < n2; i++)
    memcpy(A[a+i]+b, B[i], m2*sizeof(double));

  return 0;
}

/**
 * CONVERSION ROUTINES
 */

void nt_matrix_2_vector(double **A, double x[],
			int n, int m, int orientation)
{
  int i, j, l;
  if (orientation != NT_MATRIX_ROW && orientation != NT_MATRIX_COLUMN) {
    c_warning("** NTOOLS WARNING: incorrect orientation in routine "\
	      "nt_matrix_2_vector, using NT_ROW instead");
    orientation = NT_MATRIX_ROW;
  }

  l = 0;
  if (orientation == NT_MATRIX_ROW) {
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
	x[l++] = A[i][j];
  }
  else {
    for (j = 0; j < m; j++)
      for (i = 0; i < n; i++)
	x[l++] = A[i][j];
  }
}

/**
 *
 * @param n is the number of rows of the matrix.
 * @param m is the number of cols of the matrix.
 */
void nt_vect_2_matrix(double x[], double **A,
		      int n, int m, int orientation)
{
  int i, j, l;

  if (orientation != NT_MATRIX_ROW && orientation != NT_MATRIX_COLUMN) {
    c_warning("** NTOOLS WARNING: incorrect orientation in routine nt_vect_2_matrix, using NT_MATRIX_ROW instead");
    orientation = NT_MATRIX_ROW;
  }

  l = 0;
  if (orientation == NT_MATRIX_ROW) {
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
	A[i][j] = x[l++];
  }
  else {
    for (j = 0; j < m; j++)
      for (i = 0; i < n; i++)
	A[i][j] = x[l++];
  }
}

/**
 * Compute the transpose of a matrix
 *
 * @param A is the matrix to transpose, of dimensions n*m.
 * @param B is a matrix of dimension m*n, different of A. On output, B
 *   is the transpose of A.
 * @param n is the number of rows of the matrix.
 * @param m is the number of cols of the matrix.
 */
void nt_matrix_transpose(double **A, double **B, int n, int m)
{
  int i, j;

  if (A == B)
    c_warning("** NTOOLS WARNING: invalide use of nt_matrix_transpose: A == B");
  else {
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
	B[j][i] = A[i][j];
  }
}

/**
 * Compute A'A.
 *
 * @param n is the number of rows of the matrix.
 * @param m is the number of cols of the matrix.
 */
void nt_matrix_ata(double **A, double **C, int n, int m)
{
  int i, j, k;

  for (i = 0; i < n; i++)
    for (j = 0; i < n; j++) {
      C[i][j] = 0.0;
      for (k = 0; k < m; k++) {
	C[i][j] += A[k][i]*A[k][j];
      }
      C[j][i] = C[i][j];
    }
}


/**
 * Create an OR_Matrix structure using a preexisting matrix.
 *
 * @param A is a pointer to a double array representing a matrix, stored
 *   by rows.
 * @param m is the number of rows.
 * @param n is the number of columns.
 *
 * @return a pointer to a newly created OR_Matrix structure. The matrix
 *   content is shared with the original double array.
 */
OR_Matrix *or_matrix_create(double **A, int m, int n)
{
  OR_Matrix *om = _malloc(OR_Matrix);

  om->A = A;
  om->m = m;
  om->n = n;

  return om;
}

/**
 * Allocate memory for a new OR_Matrix structure.
 *
 * @param m is the number of rows.
 * @param n is the number of columns.
 *
 * @return a pointer to a newly created OR_Matrix structure, allocating memory
 *   for underlying matrix.
 */
OR_Matrix *or_matrix_new(int m, int n)
{
  OR_Matrix *om = _malloc(OR_Matrix);

  om->A = nt_matrix_new(m, n);
  om->m = m;
  om->n = n;

  return om;
}

/**
 * Free an OR_Matrix structure.
 *
 * @param A is an initialized OR_Matrix. The inner matrix component has
 *   to be allocated the same way as does the function nt_matrix_new.
 */
void or_matrix_free(OR_Matrix *A)
{
  if (A) {
    nt_matrix_free(A->A);
    free(A);
  }
}
