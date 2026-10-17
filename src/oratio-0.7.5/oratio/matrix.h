/**
 * nt_matrix.h
 *
 * ORATIO Library
 * Copyright (C) 2010-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_MATRIX_H__
#define __ORATIO_MATRIX_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Real matrix.
   */
  typedef struct {
    double **A;
    int m, n;   /* m and n represents the number of rows and columns,
		   respectively. */
  } OR_Matrix;

  /**
   * Integer matrix.
   */
  typedef struct {
    int **A;
    int m, n;   /* m and n represents the number of rows and columns,
		   respectively. */
  } OR_iMatrix;

  /**
   * Create an OR_Matrix structure using a preexisting matrix.
   *
   * @param A is a pointer to a matrix, stored by rows.
   * @param m is the number of rows.
   * @param n is the number of columns.
   *
   * @return a pointer to a newly created OR_Matrix structure. The matrix
   *   content is shared with the original double array.
   */
  OR_Matrix *or_matrix_create(double **A, int m, int n);

  /**
   * Allocate memory for a new OR_Matrix structure.   *
   * @param m is the number of rows.
   * @param n is the number of columns.
   *
   * @return a pointer to a newly created OR_Matrix structure, allocating memory
   *   for underlying matrix.
   */
  OR_Matrix *or_matrix_new(int m, int n);

  /**
   * Free an OR_Matrix structure.
   *
   * @param A is an initialized OR_Matrix. The inner matrix component has
   *   to be allocated the same way as does the function nt_matrix_new.
   */
  void or_matrix_free(OR_Matrix *A);

  /**
   * Matrix definitions
   */

  /**
   * Create a new matrix of dimensions m*n.
   *
   * @param m is the number of rows of the matrix to create.
   * @param n is the number of columns of the matrix to create.
   *
   * @return a pointer to the newly allocated matrix. The matrix is
   *   stored is contiguous memory.
   */
  double **nt_matrix_new(int m, int n);

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
  double **nt_ma_new(const enum CBLAS_ORDER Order, int m, int n);

  /* Need to be validated */
  /**
   * Reallocate the memory associated to the matrix A.
   *
   * @param m is the number of rows of the matrix.
   * @param n is the number of rows of the matrix.
   */
  double **nt_matrix_realloc(double **A, int nrow, int ncol);

  /**
   * Free the memory allocated to the matrix A.
   *
   * @param A is a previously allocated matrix.
   */
  void nt_matrix_free(double **A);

  /**
   * Free the memory allocated to the matrix A.
   *
   * @param A is a previously allocated matrix.
   */
  void nt_ma_free(double **A);


  /**
   * Duplicate a matrix into another place memory
   */
  double **nt_matrix_duplicate(const double **m, int nrow, int ncol);

  void nt_matrix_insert_vector(double **A, double x[], int n, int m,
			       int orientation, int pos);

  /**
   * Insert a matrix B into a matrix A
   *
   * @param A is a matrix of dimension (n1,m1).
   * @param B is a matrix of dimension (n2,m2).
   * @param a
   * @param b       (a,b): position where the matrix B is to inserted; the
   *               indexes are supposed to start at 0
   *        Both matrices are supposed to have the same orientation;
   *        and the indexes reflect this orientation. In other terms
   *        if the matrices are oriented by columns, the first index
   *        indicates the column, while is the orientation is by rows,
   *        it reflects the row. If the dimensions agree, B is
   * inserted in A at the position (a,b) on output. 
   */
  int nt_matrix_insert_matrix(double **A, int n1, int m1,
			      double **B, int n2, int m2,
			      int a, int b);

  /**
   * Compute the transpose of a matrix
   *
   * @param A is the matrix to transpose, of dimensions n*m.
   * @param B is a matrix of dimension m*n, different of A. On output, B
   *   is the transpose of A.
   */
  void nt_matrix_transpose(double **A, double **B, int n, int m);
  
  /**
   * Compute  A'A
   */
  void nt_matrix_ata(double **A, double **C, int n, int m);

  /**
   * Set A to the zeros matrix.
   *
   * @param Order: CblasRowMajor or CblasColMajor
   * @param m is the number of rows of matrix A.
   * @param n is the number of columns of matrix A.
   * @param  A: before entry, the leading m by n part of the array A must
   *   contain the matrix of coefficients. On exit, the leading m by n
   *   part of the array A contains the zeros matrix.
   * @param ldA specifies the first dimension of A as declared in the
   *   calling (sub) program. lda must be at least max( 1, m ) if Order
   *   is equal to CblasColMajor, or max (1, n) if Order is equal to
   *   CblasRowMajor.
   */
  void nt_ma_zeros(const enum CBLAS_ORDER Order,
		   const int m, const int n, double *A, const int ldA);

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
		      const int m, const int n, double *A, const int ldA);

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
   *   calling (sub) program. lda must be at least max( 1, m ).
   * @param alpha is some real number.
   */
  void nt_ma_diag(const enum CBLAS_ORDER Order,
		  const int m, const int n, double *A, const int ldA,
		  double alpha);

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
   *   calling (sub) program. lda must be at least max( 1, m ).
   */
  void nt_ma_scale(const enum CBLAS_ORDER Order, const int m, const int n,
		   double alpha, double *A, int ldA);

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
		   double *A, int ldA, double *x, int incx);

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
			   double *A, int ldA, double *x, int incx);

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
  void nt_matrix_zeros(int m, int n, double *A, int ldA);

  /**
   * Set a mxn matrix to the identity matrix.
   *
   * @param m is the number of rows of the matrix A.
   * @param n is the number of columns of the matrix A.
   * @param A is a pointer to a matrix of dimensions (at least) mxn. On
   * the principal submatrix of dimension mxn of A is the identity matrix.
   * @param ldA is the leading dimensions of A.
   */
  void nt_matrix_identity(int m, int n, double *A, int ldA);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_MATRIX_H__ */
