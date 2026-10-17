/**
 * bfgs.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2002-2007 Laurent Martin
 *                         Emilie Wanufelle <ewa@math.fundp.ac.be>,
 *                         Fabian Bastin    <fabian.bastin@cerfacs.fr>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "bfgs.h"

#include "defaults.h"
#include "linesearch.h"
#include "stopping.h"

#include <oratio/ctools.h>
#ifndef CBLASWR_CXML
#include <oratio/f2c.h>
#endif

#include <float.h>
#include <stdio.h>

extern /* Subroutine */ int dposvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	char *equed, doublereal *s, doublereal *b, integer *ldb, doublereal *
	x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *
	berr, doublereal *work, integer *iwork, integer *info);
 

#define LINESEARCH_STEP_MAX 1.e3
#define LINESEARCH_STEP_MIN 1.e-9

/**
 * Positive definite secant update (BFGS), unfactored form.
 *
 * Ref.: J.E. Dennis, Jr. and Robert B. Schnabel, "Numerical Methods
 * for Unconstrained Optimization and Nonlinear Equations", SIAM,
 * 1983, Algorithm 9.4.1, page 355.
 *
 * The implementation in Algorithm 9.4.1 is more conservative than the
 * current, due to the use of an additional parameter eta, described
 * in the guidelines 2, in particular section "fdigits".
 *
 * @param Order designs the orientation of the matrix Bk. The
 *   possible values are CblasRowMajor and CblasColMajor.
 * @param Uplo designs the triangular part that has to be updated. The
 *   possible value CblasUpper and CblasLower.
 * @param ldbk is the leading dimension of matrix Bk.
 * @param n is the problem dimension.
 * @param yk is the gradient difference (g+-gk).
 * @param incyk is the increment in vector yk.
 * @param sk is the current step: sk = x+-xk.
 * @param incsk is the increment in vector sk.
 * @param work is a working vector of size n.
 *
 * @return 1 is the updated has been performed, 0 else.
 */
int bfgs_update(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
		const int n, double *Bk, const int ldbk,
		const double *yk, const int incyk,
		const double *sk, const int incsk,
		double *work)
{
  double snorm, ynorm, temp;

  snorm = cblas_dnrm2(n, sk, incsk);
  ynorm = cblas_dnrm2(n, yk, incyk);
  temp = cblas_ddot(n, sk, 1, yk, 1);

  if (temp > sqrt(DBL_EPSILON)*snorm*ynorm) {
    /* Compute p = Bk*sk */
    cblas_dsymv(Order, Uplo, n, 1, Bk, ldbk, sk, incsk, 0.0, work, 1);

    cblas_dsyr(Order, Uplo, n,  1/cblas_ddot(n, yk, incyk, sk, incsk),
	       yk, incyk, Bk, ldbk);

    cblas_dsyr(Order, Uplo, n, -1/cblas_ddot(n, sk, incsk, work, 1), work, 1,
	       Bk, ldbk);

    return 1;
  }

  return 0;
}

/**
 * Function computing an iterate of the BFGS's method
 *
 * @param Order designs the orientation of the matrix Bk. The
 *   possible values are CblasRowMajor and CblasColMajor.
 * @param Uplo designs the triangular part that has to be updated. The
 *   possible value CblasUpper and CblasLower.
 * @param f is the function to minimize
 * @param fg is function to minimize with its gradient (optional).
 * @param n is the problem dimension.
 * @param x is the initial point. On return,s is the point obtained
 *   after one algo BFGS's iteration.
 * @param data represents the data needed by function f (NULL if no
 *   data is necessary).
 * @param scaling is the  scaling factor; if no scaling is used in the
 *   function evaluation, the scaling factor should be set to 1.0.
 * @param Bk is the initial secant matrix. On return Bk is the matrix
 *   obtained after one algo BFGS's iteration.
 * @param work is a double working vector of dimension 3*n.
 * @param iwork is an integer working vector of dimension n.
 *
 * @return the function value at x after one algo BFGS's iteration.
 */
double bfgs_iterate(BFGS *b,
		    const enum CBLAS_ORDER Order,
		    const enum CBLAS_UPLO Uplo,
		    C_GENERIC f, int f_type,
		    double (*fg)(double *, int, void *, double *),
		    double *scaling,
		    void *data, int n, double *x, double **Bk,
		    double *p, double *work, integer *iwork)
{
  double *gradplus, *invBk, *xplus,*s;
  /* Bk will contain the Cholesky factorization of Bk. */
  double a,rcond,ferr,berr;
  const double alpha = 0.0001, gamma = 0.9 ;
  integer taille, nrhs;
  int i; 
  char uplo, fact, equed;

  double val;

  integer info;

  info = b->info;

  invBk = _mallocs(n*n, double);
  gradplus = _mallocs(n, double);
  xplus = _mallocs(n, double);

  s = xplus;

  taille=n;
  nrhs=1;
  a=1;
  uplo='L';
  fact='N';
  equed='E';

  for (i=0; i<n*n; i++) invBk[i]= (*Bk)[i];  
  
  /* Computation of the descent direction, by solving the system p =
     -B^{-1} grad */
  cblas_dcopy(n, b->grad, 1, gradplus,1); 
  cblas_dscal(n,-1,gradplus,1);
  dposvx_(&fact, &uplo, &taille, &nrhs, *Bk, &taille, invBk, &taille,
	  &equed, s, gradplus, &taille, p, &taille, &rcond, &ferr, &berr,
	  work ,iwork, &info);
  cblas_dscal(n,-1,gradplus,1);
  cblas_dcopy(n,x,1,xplus,1);
  
  /* Linesearch (Wolfe-Thuente algorithm) */

  /* If the search direction isn't a descent direction, we take by default
     -grad */
  val = b->value;
  info = dcsrch(n, gradplus, &val, xplus, b->wrapper_fct,
		f, f_type, fg, data,
		b->d1, scaling, p, &a, alpha, 2.2e-16, gamma,
		LINESEARCH_STEP_MAX, LINESEARCH_STEP_MIN, work);
  if (info == NLP_LINESEARCH_DESCENT_DIRECTION) {
    c_warning("initial step is not a descent direction in BFGS algorithm");
    cblas_dscal(n, -1, b->grad, 1);
    cblas_dcopy(n, b->grad, 1, p, 1);
    info = dcsrch(n, gradplus, &val, xplus, b->wrapper_fct,
		  f, f_type, fg, data, b->d1,
		  scaling, p, &a, alpha, 2.2e-16, gamma,
		  LINESEARCH_STEP_MAX, LINESEARCH_STEP_MIN, work);
  }
  else
    cblas_dscal(n, -1, b->grad, 1);

  if (info != NLP_LINESEARCH_WOLFE_CONDITIONS) {
    /* We try x+s in order to see if it produces some decrease in the
       objective */
    /* Try to avoid to recompute these informations */
    cblas_dcopy(n, x, 1, xplus, 1);
    cblas_daxpy(n, a, p, 1, xplus, 1);

    if (fg)
      val = fg(xplus, n, data, gradplus);
    else {
      val = b->wrapper_fct(f, f_type, xplus, n, b->m, data,
			   NULL, gradplus,
			   Bk, NULL, 0, &(b->info));
      if (b->d1)
	b->d1(xplus, n, data, gradplus);
      else
	nt_derive_grad_fd((NT_f)f, xplus, gradplus, work, 0, scaling, n, 0, 0, data);
    }
  }

  if (info == NLP_LINESEARCH_WOLFE_CONDITIONS || val < b->value) {
    /* The linesearch has succeeded. */

    /* BFGS Update */
    info = 1;
    b->value = val;

    /* Computation of sk (in x) and yk (in grad) */ 
    cblas_dscal(n,-1,x,1);
    cblas_daxpy(n,1,xplus,1,x,1);
    cblas_daxpy(n, 1, gradplus, 1, b->grad, 1);

    bfgs_update(CblasRowMajor, CblasUpper, n, *Bk, n, b->grad, 1, x, 1, p);
 
    /* x=xplus   grad=gradplus */
    cblas_dcopy(n, xplus, 1, x, 1);
    cblas_dcopy(n*(b->m+1), gradplus, 1, b->grad, 1);

    b->oldgradsqrnrm = b->gradsqrnrm;
    b->gradsqrnrm = cblas_ddot(n, b->grad, 1, b->grad, 1);
  }

  free(invBk);
  free(gradplus);
  free(xplus);

  b->info = info;

  return b->value;
}

/**
 * BFGS method. Expert routine.
 * Function computing the minimum of a function of class C1
 *
 * @param bf is an initialized BFGS structure.
 * @param Order designs the orientation of the matrix Bk. The
 *   possible values are CblasRowMajor and CblasColMajor.
 * @param Uplo designs the triangular part that has to be updated. The
 *   possible value CblasUpper and CblasLower.
 * @param f is the function to minimize.
 * @param n is the primary problem dimension.
 * @param x is on input the starting point. It is the solution on output.
 * @param output is the output file for results and iteration informations
 * @param Bk is the initial secant matrix.
 * @param data represents the data needed by function f (NULL if no
 *   data is necessary).
 * @param work is a double precision working vector of size at least 3*n.
 *
 * @return the optimal value.
 */
double nlp_bfgsx(BFGS *bf,
		 const enum CBLAS_ORDER Order,
		 const enum CBLAS_UPLO Uplo,
		 C_GENERIC f, int f_type,
		 double (*fg)(double *, int, void *, double *),
		 void *data, FILE *output,
		 double **Bk, int n, double *x,
		 double *work)
{
  int i;
  double *h, *scaling;
  integer *iwork;
  integer info = 1;

  int has_work = 1;

  bf->k = 0;  /* set the iteration counter */
  bf->info = info;

  if (!work) {
    has_work = 0;
    work = _mallocs(3*n, double);
  }

  h = _mallocs(n, double);
  scaling = _mallocs(n, double);
  iwork = _mallocs(n, integer);

  for(i=0;i<n;i++) h[i]=0;
  for(i=0;i<n;i++) scaling[i]=1.0;
  
  /* gradient initialisation */
  if (fg)
    bf->value = fg(x, n, data, bf->grad);
  else {
    bf->value = bf->wrapper_fct(f, f_type, x, n, bf->m, data,
				NULL, bf->grad,
				Bk, NULL, 0, &(bf->info));
    if (bf->d1)
      bf->d1(x, n, data, bf->grad);
    else
      nt_derive_grad_fd((NT_f)f, x, bf->grad, h, 0, scaling, n, 0, 0, data);
    bf->gradsqrnrm = cblas_ddot(n, bf->grad, 1, bf->grad, 1);
  }

  /* Update the secant initialization. Not OPTIMAL!!! */
  nlp_secant_init_default(Order, Uplo, bf->n, bf->m, *Bk,
			  NULL, NULL, bf->grad, NULL,
			  NULL, &(bf->info), work, NULL);
  while (!bf->stop(bf)) {
    if (bf->printer) bf->printer(output, bf, data);
    bf->k++;
    bfgs_iterate(bf, Order, Uplo,
		 f, f_type, fg, scaling, data, n, x, Bk, h, work, iwork);
  }

  if (bf->printer) bf->printer(output, bf, data);

  free(h);
  free(scaling);
  free(iwork);
  if (!has_work) free(work);
  
  return bf->value;
}

/**
 * BFGS method. Simple routine.
 * Function computing the minimum of a function of class C1
 *
 * @param work is a double precision working vector of size at least 3*n.
 */
double nlp_bfgs(BFGS *bf, C_GENERIC f, C_GENERIC g, void *data,
		FILE *output, double **hessian, double *work)
{
  if (g)
    return nlp_bfgsx(bf, CblasRowMajor, CblasUpper,
		     f, NlpF, (NT_fg)g, data, output, hessian,
		     bf->n, bf->x, work);
  else
    return nlp_bfgsx(bf, CblasRowMajor, CblasUpper,
		     f, NlpFG, (NT_fg)f, data, output, hessian,
		     bf->n, bf->x, work);
}

/**
 * BFGS secant update.
 *
 * @param Order designs the orientation of the matrix Bk. The
 *   possible values are CblasRowMajor and CblasColMajor.
 * @param Uplo designs the triangular part that has to be updated. The
 *   possible value CblasUpper and CblasLower.
 * @param n is the primary problem dimension.
 * @param m is the secondary problem dimension.
 */
int nlp_bfgs_update(const enum CBLAS_ORDER Order,
		    const enum CBLAS_UPLO Uplo,
		    int n, int m, double *H,
		    double *values, double *valuesplus,
		    double *grad, double *gradplus, double *step,
		    int *info, double *work, void *data)
{
  cblas_dcopy(n, gradplus, 1, work, 1);
  cblas_daxpy(n, -1.0, grad, 1, work, 1);

  return bfgs_update(Order, Uplo, n, H, n, work, 1, step, 1, work+n);
}
