/**
 * cblaswr_cxml.h
 *
 * NTools Library
 * Copyright (C) 2005-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ATHOR_NTOOLS_CBLASWR_CXML_H__
#define __ATHOR_NTOOLS_CBLASWR_CXML_H__

#include "cblas.h"

#ifdef CBLASWR_CXML

#include "ntools_f2c.h"

extern doublecomplex zsum_(int *, doublecomplex *, int *);

extern float sdot_(int* N, float* X, int* incX, float* Y, int* incY);

extern double ddot_(int* N, double* X, int* incX, double* Y, int* incY);

extern void cdotu_(complex* retval, int* N, complex* X, int* incX,
		   complex* Y, int* incY);

extern void cdotc_(complex* retval, int* N, complex* X, int* incX,
		   complex* Y, int* incY);

extern void zdotu_(doublecomplex* retval, int* N, doublecomplex* X, int* incX,
		   doublecomplex* Y, int* incY);

extern void zdotc_(doublecomplex* retval, int* N, doublecomplex* X, int* incX, 
		   doublecomplex* Y, int* incY);

extern float snrm2_(int* N, float* X, int* incX);

extern float sasum_(int* N, float* X, int* incX);

extern double dnrm2_(int* N, double* X, int* incX);

extern double dasum_(int* N, double* X, int* incX);

extern float scnrm2_(int* N, complex* X, int* incX);

extern float scasum_(int* N, complex* X, int* incX);

extern double dznrm2_(int* N, doublecomplex* X, int* incX);

extern double dzasum_(int* N, doublecomplex* X, int* incX);

extern int isamax_(int* N, float* X, int* incX);

extern int idamax_(int* N, double* X, int* incX);

extern int icamax_(int* N, complex* X, int* incX);

extern int izamax_(int* N, doublecomplex* X, int* incX);

extern int sswap_(int* N, float* X, int* incX, float* Y, int* incY);

extern int scopy_(int* N, float* X, int* incX, float* Y, int* incY);

extern int saxpy_(int* N, float* alpha, float* X, int* incX,
		  float* Y, int* incY);

extern int dswap_(int* N, double* X, int* incX, double* Y, int* incY);

extern int dcopy_(int* N, double* X, int* incX, double* Y, int* incY);

extern int daxpy_(int* N, double* alpha, double* X, int* incX,
		  double* Y, int* incY);

extern int cswap_(int* N, complex* X, int* incX, complex* Y, int* incY);

extern int ccopy_(int* N, complex* X, int* incX, complex* Y, int* incY);

extern int caxpy_(int* N, complex* alpha, complex* X, int* incX,
		  complex* Y, int* incY);

extern int zswap_(int* N, doublecomplex* X, int* incX, doublecomplex* Y, int* incY);

extern int zcopy_(int* N, doublecomplex* X, int* incX, doublecomplex* Y, int* incY);

extern int zaxpy_(int* N, doublecomplex* alpha, doublecomplex* X, int* incX,
		  doublecomplex* Y, int* incY);

extern int srotg_(float* a, float* b, float* c, float* s);

extern int srot_(int* N, float* X, int* incX, float* Y, int* incY,
		 float* c, float* s);

extern int drotg_(double* a, double* b, double* c, double* s);

extern int drot_(int* N, double* X, int* incX, double* Y, int* incY,
		 double* c, double* s);

extern int sscal_(int* N, float* alpha, float* X, int* incX);

extern int dscal_(int* N, double* alpha, double* X, int* incX);

extern int cscal_(int* N, complex* alpha, complex* X, int* incX);

extern int zscal_(int* N, doublecomplex* alpha, doublecomplex* X, int* incX);

extern int csscal_(int* N, float* alpha, complex* X, int* incX);

extern int zdscal_(int* N, double* alpha, doublecomplex* X, int* incX);

extern int sgemv_(char* trans, int* M, int* N, float* alpha,
		  float* A, int* lda, float* X, int* incX,
		  float* beta, float* Y, int* incY);

extern int sgbmv_(char *trans, int *M, int *N, int *KL, int *KU, 
		  float *alpha, float *A, int *lda, float *X, int *incX, 
		  float *beta, float *Y, int *incY);

extern int strmv_(char* uplo, char *trans, char* diag, int *N,  
		  float *A, int *lda, float *X, int *incX);

extern int stbmv_(char* uplo, char* trans, char* diag, int* N, int* K,
		  float* A, int* lda, float* X, int* incX);

extern int stpmv_(char* uplo, char* trans, char* diag, int* N, 
		  float* Ap, float* X, int* incX);

extern int strsv_(char* uplo, char* trans, char* diag, int* N,
		  float* A, int* lda, float* X, int* incX);

extern int stbsv_(char* uplo, char* trans, char* diag, int* N, int* K,
		  float* A, int* lda, float* X, int* incX);

extern int stpsv_(char* uplo, char* trans, char* diag, int* N,
		  float* Ap, float* X, int* incX);

extern int dgemv_(char* trans, int* M, int* N, double* alpha,
		  double* A, int* lda, double* X, int* incX, double* beta,
		  double* Y, int* incY);

extern int dgbmv_(char *trans, int *M, int *N, int *KL, int *KU,
		  double *alpha, double *A, int *lda, double *X, int *incX, 
		  double *beta, double *Y, int *incY);

extern int dtrmv_(char* uplo, char *trans, char* diag, int *N,
		  double *A, int *lda, double *X, int *incX);

extern int dtbmv_(char* uplo, char* trans, char* diag, int* N, int* K,
		  double* A, int* lda, double* X, int* incX);

extern int dtpmv_(char* uplo, char* trans, char* diag, int* N, 
		  double* Ap, double* X, int* incX);

extern int dtrsv_(char* uplo, char* trans, char* diag, int* N,
		  double* A, int* lda, double* X, int* incX);

extern int dtbsv_(char* uplo, char* trans, char* diag, int* N, int* K,
		  double* A, int* lda, double* X, int* incX);

extern int dtpsv_(char* uplo, char* trans, char* diag, int* N, 
		  double* Ap, double* X, int* incX);

extern int cgemv_(char* trans, int* M, int* N, complex* alpha,
		  complex* A, int* lda, complex* X, int* incX,
		  complex* beta, complex* Y, int* incY);

extern int cgbmv_(char *trans, int *M, int *N, int *KL, int *KU, 
		  complex *alpha, complex *A, int *lda, 
		  complex *X, int *incX, complex *beta, 
		  complex *Y, int *incY);

extern int ctrmv_(char* uplo, char *trans, char* diag, int *N,
		  complex *A, int *lda, complex *X, int *incX);

extern int ctbmv_(char* uplo, char* trans, char* diag, int* N, int* K,
		  complex* A, int* lda, complex* X, int* incX);

extern int ctpmv_(char* uplo, char* trans, char* diag, int* N, 
		  complex* Ap, complex* X, int* incX);

extern int ctrsv_(char* uplo, char* trans, char* diag, int* N,
		  complex* A, int* lda, complex* X, int* incX);

extern int ctbsv_(char* uplo, char* trans, char* diag, int* N, int* K,
		  complex* A, int* lda, complex* X, int* incX);

extern int ctpsv_(char* uplo, char* trans, char* diag, int* N, 
		  complex* Ap, complex* X, int* incX);

extern int zgemv_(char* trans, int* M, int* N, doublecomplex* alpha,
		  doublecomplex* A, int* lda, doublecomplex* X, int* incX,
		  doublecomplex* beta, doublecomplex* Y, int* incY);

extern int zgbmv_(char *trans, int *M, int *N, int *KL, int *KU, 
		  doublecomplex *alpha, doublecomplex *A, int *lda, 
		  doublecomplex *X, int *incX, doublecomplex *beta, 
		  doublecomplex *Y, int *incY);

extern int  ztrmv_(char* uplo, char *trans, char* diag, int *N,  
		   doublecomplex *A, int *lda, doublecomplex *X, int *incX);

extern int ztbmv_(char* uplo, char* trans, char* diag, int* N, int* K,
		  doublecomplex* A, int* lda, doublecomplex* X, int* incX);

extern void ztpmv_(char* uplo, char* trans, char* diag, int* N, 
		   doublecomplex* Ap, doublecomplex* X, int* incX);

extern int ztrsv_(char* uplo, char* trans, char* diag, int* N,
		  doublecomplex* A, int* lda, doublecomplex* X, int* incX);

extern int ztbsv_(char* uplo, char* trans, char* diag, int* N, int* K,
		  doublecomplex* A, int* lda, doublecomplex* X, int* incX);

extern int ztpsv_(char* uplo, char* trans, char* diag, int* N, 
		  doublecomplex* Ap, doublecomplex* X, int* incX);

extern int ssymv_(char* uplo, int* N, float* alpha, float* A, int* lda,
		  float* X, int* incX, float* beta,
		  float* Y, int* incY);

extern int ssbmv_(char* uplo, int* N, int* K, float* alpha,
		  float* A, int* lda, float* X, int* incX,
		  float* beta, float* Y, int* incY);

extern int sspmv_(char* uplo, int* N, float* alpha, float* Ap,
		  float* X, int* incX, float* beta,
		  float* Y, int* incY);

extern int sger_(int* M, int* N, float* alpha, float* X, int* incX,
		 float* Y, int* incY, float* A, int* lda);

extern int ssyr_(char* uplo, int* N, float* alpha, float* X, int* incX,
		 float* A, int* lda);

extern int sspr_(char* uplo, int* N, float* alpha, float* X, int* incX,
		 float* Ap);

extern int ssyr2_(char* uplo, int* N, float* alpha, float* X, int* incX,
		  float* Y, int* incY, float* A, int* lda);

extern int sspr2_(char* uplo, int* N, float* alpha, float* X, int* incX,
		  float* Y, int* incY, float* A);

extern int dsymv_(char* uplo, int* N, double* alpha, double* A, int* lda,
		  double* X, int* incX, double* beta, double* Y, int* incY);

extern int dsbmv_(char* uplo, int* N, int* K, double* alpha,
		  double* A, int* lda, double* X, int* incX,
		  double* beta, double* Y, int* incY);

extern int dspmv_(char* uplo, int* N, double* alpha, double* Ap,
		  double* X, int* incX, double* beta,
		  double* Y, int* incY);

extern int dger_(int* M, int* N, double* alpha, double* X, int* incX,
		 double* Y, int* incY, double* A, int* lda);

extern int dsyr_(char* uplo, int* N, double* alpha, double* X, int* incX,
		 double* A, int* lda);

extern int dspr_(char* uplo, int* N, double* alpha,
		 double* X, int* incX, double* Ap);

extern int dsyr2_(char* uplo, int* N, double* alpha,
		  double* X, int* incX, double* Y, int* incY,
		  double* A, int* lda);

extern int dspr2_(char* uplo, int* N, double* alpha, 
		  double* X, int* incX, double* Y, int* incY,
		  double* A);

extern int chemv_(char* uplo, int* N, complex* alpha,
		  complex* A, int* lda, complex* X, int* incX,
		  complex* beta, complex* Y, int* incY);

extern int chbmv_(char* uplo, int* N, int* K, complex* alpha,
		  complex* A, int* lda, complex* X, int* incX,
		  complex* beta, complex* Y, int* incY);

extern int chpmv_(char* uplo, int* N, complex* alpha,
		  complex* Ap, complex* X, int* incX,
		  complex* beta, complex* Y, int* incY);

extern int cgeru_(int* M, int* N, complex* alpha,
		  complex* X, int* incX, complex* Y, int* incY,
		  complex* A, int* lda);

extern int cgerc_(int* M, int* N, complex* alpha,
		  complex* X, int* incX, complex* Y, int* incY,
		  complex* A, int* lda);

extern int cher_(char* uplo, int* N, float* alpha,
		 complex* X, int* incX, complex* A, int* lda);

extern int chpr_(char* uplo, int* N, float* alpha,
		 complex* X, int* incX, complex* Ap);

extern int cher2_(char* uplo, int* N, complex* alpha,
		  complex* X, int* incX, complex* Y, int* incY,
		  complex* A, int* lda);

extern int chpr2_(char* uplo, int* N, complex* alpha,
		  complex* X, int* incX, complex* Y, int* incY,
		  complex* Ap);

extern int zhemv_(char* uplo, int* N, doublecomplex* alpha,
		  doublecomplex* A, int* lda, doublecomplex* X, int* incX,
		  doublecomplex* beta, doublecomplex* Y, int* incY);

extern int zhbmv_(char* uplo, int* N, int* K, doublecomplex* alpha,
		  doublecomplex* A, int* lda, doublecomplex* X, int* incX,
		  doublecomplex* beta, doublecomplex* Y, int* incY);

extern int zhpmv_(char* uplo, int* N, doublecomplex* alpha,
		  doublecomplex* Ap, doublecomplex* X, int* incX,
		  doublecomplex* beta, doublecomplex* Y, int* incY);

extern int zgeru_(int* M, int* N, doublecomplex* alpha,
		  doublecomplex* X, int* incX, doublecomplex* Y, int* incY,
		  doublecomplex* A, int* lda);

extern int zgerc_(int* M, int* N, doublecomplex* alpha,
		  doublecomplex* X, int* incX, doublecomplex* Y, int* incY,
		  doublecomplex* A, int* lda);

extern int zher_(char* uplo, int* N, double* alpha,
		 doublecomplex* X, int* incX, doublecomplex* A, int* lda);

extern int zhpr_(char* uplo, int* N, double* alpha,
		 doublecomplex* X, int* incX, doublecomplex* Ap);

extern int zher2_(char* uplo, int* N, doublecomplex* alpha,
		  doublecomplex* X, int* incX, doublecomplex* Y, int* incY,
		  doublecomplex* A, int* lda);

extern int zhpr2_(char* uplo, int* N, doublecomplex* alpha,
		  doublecomplex* X, int* incX, doublecomplex* Y, int* incY,
		  doublecomplex* Ap);

extern int sgemm_(char* transA, char* transB, int* M, int* N, int* K,
		  float* alpha, float* A, int* lda, float* B, int* ldb,
		  float* beta, float* C, int* ldc);

extern int ssymm_(char* side, char* uplo, int* M, int* N,
		  float* alpha, float* A, int* lda, float* B, int* ldb,
		  float* beta, float* C, int* ldc);

extern int ssyrk_(char* uplo, char* trans, int* N, int* K,
		  float* alpha, float* A, int* lda, float* beta,
		  float* C, int* ldc);

extern int ssyr2k_(char* uplo, char* trans, int* N, int* K,
		   float* alpha, float* A, int* lda, float* B, int* ldb,
		   float* beta, float* C, int* ldc);

extern int strmm_(char* side, char* uplo, char* trans, char* diag, 
		  int* M, int* N, float* alpha, float* A, int* lda,
		  float* B, int* ldb);

extern int strsm_(char* side, char* uplo, char* trans, char* diag,
		  int* M, int* N, float* alpha, float* A, int* lda,
		  float* B, int* ldb);

extern int dgemm_(char* transA, char* transB, int* M, int* N, int* K,
		  double* alpha, double* A, int* lda, double* B, int* ldb,
		  double* beta, double* C, int* ldc);

extern int dsymm_(char* side, char* uplo, int* M, int* N, double* alpha,
		  double* A, int* lda, double* B, int* ldb,
		  double* beta, double* C, int* ldc);

extern int dsyrk_(char* uplo, char* trans, int* N, int* K,
		  double* alpha, double* A, int* lda,
		  double* beta, double* C, int* ldc);

extern int dsyr2k_(char* uplo, char* trans, int* N, int* K,
		   double* alpha, double* A, int* lda,
		   double* B, int* ldb, double* beta,
		   double* C, int* ldc);

extern int dtrmm_(char* side, char* uplo, char* trans, char* diag, 
		  int* M, int* N, double* alpha, double* A, int* lda,
		  double* B, int* ldb);

extern int dtrsm_(char* side, char* uplo, char* trans, char* diag,
		  int* M, int* N, double* alpha, double* A, int* lda,
		  double* B, int* ldb);

extern int cgemm_(char* transA, char* transB, int* M, int* N, int* K,
		  complex* alpha, complex* A, int* lda,
		  complex* B, int* ldb, complex* beta,
		  complex* C, int* ldc);

extern int csymm_(char* side, char* uplo, int* M, int* N,
		  complex* alpha, complex* A, int* lda,
		  complex* B, int* ldb, complex* beta,
		  complex* C, int* ldc);

extern int csyrk_(char* uplo, char* trans, int* N, int* K,
		  complex* alpha, complex* A, int* lda,
		  complex* beta, complex* C, int* ldc);

extern int csyr2k_(char* uplo, char* trans, int* N, int* K,
		   complex* alpha, complex* A, int* lda,
		   complex* B, int* ldb, complex* beta,
		   complex* C, int* ldc);

extern int ctrmm_(char* side, char* uplo, char* trans, char* diag, 
		  int* M, int* N, complex* alpha,
		  complex* A, int* lda, complex* B, int* ldb);

extern int ctrsm_(char* side, char* uplo, char* trans, char* diag,
		  int* M, int* N, complex* alpha, complex* A, int* lda,
		  complex* B, int* ldb);

extern int zgemm_(char* transA, char* transB, int* M, int* N, int* K,
		  doublecomplex* alpha, doublecomplex* A, int* lda,
		  doublecomplex* B, int* ldb, doublecomplex* beta,
		  doublecomplex* C, int* ldc);

extern int zsymm_(char* side, char* uplo, int* M, int* N,
		  doublecomplex* alpha, doublecomplex* A, int* lda,
		  doublecomplex* B, int* ldb, doublecomplex* beta,
		  doublecomplex* C, int* ldc);

extern int zsyrk_(char* uplo, char* trans, int* N, int* K,
		  doublecomplex* alpha, doublecomplex* A, int* lda,
		  doublecomplex* beta, doublecomplex* C, int* ldc);

extern int zsyr2k_(char* uplo, char* trans, int* N, int* K,
		   doublecomplex* alpha, doublecomplex* A, int* lda,
		   doublecomplex* B, int* ldb, doublecomplex* beta,
		   doublecomplex* C, int* ldc);

extern int ztrmm_(char* side, char* uplo, char* trans, char* diag, 
		  int* M, int* N, doublecomplex* alpha,
		  doublecomplex* A, int* lda, doublecomplex* B, int* ldb);

extern int ztrsm_(char* side, char* uplo, char* trans, char* diag,
		  int* M, int* N, doublecomplex* alpha,
		  doublecomplex* A, int* lda, doublecomplex* B, int* ldb);

extern int chemm_(char* side, char* uplo, int* M, int* N,
		  complex* alpha, complex* A, int* lda,
		  complex* B, int* ldb, complex* beta,
		  complex* C, int* ldc);

extern int cherk_(char* uplo, char* trans, int* N, int* K,
		  float* alpha, complex* A, int* lda,
		  float* beta, complex* C, int* ldc);

extern int cher2k_(char* uplo, char* trans, int* N, int* K,
		   complex* alpha, complex* A, int* lda,
		   complex* B, int* ldb, float* beta,
		   complex* C, int* ldc);

extern int zhemm_(char* side, char* uplo, int* M, int* N,
		  doublecomplex* alpha, doublecomplex* A, int* lda,
		  doublecomplex* B, int* ldb, doublecomplex* beta,
		  doublecomplex* C, int* ldc);

extern int zherk_(char* uplo, char* trans, int* N, int* K,
		  double* alpha, doublecomplex* A, int* lda,
		  double* beta, doublecomplex* C, int* ldc);

extern int zher2k_(char* uplo, char* trans, int* N, int* K,
		   doublecomplex* alpha, doublecomplex* A, int* lda,
		   doublecomplex* B, int* ldb, double* beta,
		   doublecomplex* C, int* ldc);

#endif

#endif /* __ATHOR_NTOOLS_CBLASWR_CXML_H__ */
