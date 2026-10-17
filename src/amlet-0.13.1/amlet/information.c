/**
 * information.c
 *
 * AMLET Library
 * Copyright (C) 2005-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <string.h>

#include <oratio/ctools.h>
#include <oratio/ntools.h>
#include <amlet/dcl.h>
#include <ophelia/matrix.h>

#include "likelihood.h"
#include "model.h"

void amlet_dump_information_matrices(AMLET_Model *m)
{
  OP_Problem *p = amlet_get_problem(m);
  char *filename = strdup(op_problem_get_name(p));
  FILE *f = NULL;

  filename = c_string_append(filename, ".cov");
  f = c_fopen(filename, "w");

  if (!f) f = stdout;
  free(filename);

  nt_matrix_print(f, "H", m->H, p->n, p->n);
  nt_matrix_print(f, "Hinv", m->Hinv, p->n, p->n);
  nt_matrix_print(f, "I", m->I, p->n, p->n);
  nt_matrix_print(f, "WCov", m->WCov, p->n, p->n);
  nt_matrix_print(f, "RCov", m->RCov, p->n, p->n);

  c_fclose(f);
}

/**
 * t-test
 *
 * @param n is the number of parameters.
 * @param theta if the vector of parameters.
 * @param I is a matrix having parameters variances on its diagonal
 *   and zero as other entries.
 * @paran alpha us the level used for the test.
 * @param t is a vector of size n on input, and contains the
 *   t-statistics on output.
 *
 * @return the percentile used for the test.
 */
double amlet_t_statistics(int n, double *theta, double *hypothetical,
			  double **I, double alpha, double *t)
{

  int i;
  double p;
  //  double pars = rint(n-1);
  double pars[2];

  pars[0] = 0.0;
  pars[1] = 1.0;

  /* Compute the percentile needed for the test */
  p = st_gaussian_icdf(1-alpha/2);

  /* We should use a student... */
  for (i = 0; i < n; i++)
    t[i] = fabs((theta[i]-hypothetical[i])/sqrt(fabs(I[i][i])));

  return p;
}

/**
 * Compute the correlation matrix.
 *
 * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
 * @param Uplo designs which half-matrix is filled (CblasUpper or
 *   CblasLower).
 * @param V is the variance/covariance matrix.
 *
 * #POST: The correlation values are stocked in triangular off-diagonal of
 *        of V that is not used for covariance values.
 */
void amlet_correlation(const enum CBLAS_ORDER Order,
		       const enum CBLAS_UPLO Uplo,
		       double **V, int n)
{
  int c = 1;
  int i, j;

  /* We follow C row orientation */
  switch (Order) {
  case CblasRowMajor:
    c =  (Uplo == CblasLower) ? 'L' : 'U';
    break;
  case CblasColMajor:
    c =  (Uplo == CblasLower) ? 'U' : 'L';
    break;
  }

  switch (c) {
  case 'L':
    for (i = 1; i < n; i++) {
      for (j = 0; j < i; j++)
	V[i][j] = V[j][i]/sqrt(fabs(V[i][i]*V[j][j]));
    }
    break;
  case 'U':
    for (i = 0; i < n; i++) {
      for (j = i+1; j < n; j++)
	V[i][j] = V[j][i]/sqrt(fabs(V[i][i]*V[j][j]));
    }
    break;
  }
}

/**
 * Compute the inverse of the likehood Hessian.
 *
 * @param m is an estimated model, including the maximum likelihood Hessian.
 */
void amlet_likelihood_Hinv(AMLET_Model *m)
{
  int i, j;
  OP_Problem *p = amlet_get_problem(m);

  memcpy(m->Hinv[0], m->H[0], p->n*p->n*sizeof(double));

  op_matrix_sy_inverse(CblasRowMajor, CblasUpper, m->Hinv[0], p->n);

  /* Copy the upper part of the hessian into the lower part, as
     op_matrix_inverse only fills half of the matrix. */
  for (i = 1; i < p->n; i++)
    for (j = 0; j < i; j++)
      m->Hinv[i][j] = m->Hinv[j][i];
}

/**
 * Compute the Hessian associated to the likelihood of the model.
 *
 * @param m is an initialized AMLET model.
 */
double *amlet_model_likelihood_hessian(AMLET_Model *m,
				       double *x, double *grad,
				       int n_random_coeff, dcl_data *md,
				       double *work)
{
  int s = m->p->n;
  double *scaling = NULL;

  int i;

  if (m->verbose)
    printf("Compute the Hessian matrix; this may take a while.\n");
  scaling = _mallocs(md->npar, double);
  nt_vect_zeros(2*s, work, 1);

  for (i = 0; i < s; i++) scaling[i] = 1.0;

  if (!n_random_coeff)
    /* Multinomial logit */
    nt_derive_hess_fd(likelihood_gradient, x, work,
		      m->H, work+s, 0.0, scaling, s,
		      NULL, NULL, md);
  else
    /* Mixed logit */
    nt_derive_hess_fd_gradval(wrapper_simulated_grad, x,
			      grad, m->H, work, 0.0, scaling,
			      s, NULL, work+s, md);

  return m->H[0];
}

void amlet_weak_covariance_matrix(AMLET_Model *m, int n,
				  double scale)
{
  memcpy(m->WCov[0], m->Hinv[0], n*n*sizeof(double));

  /* We rescale the matrix to take into account the number of observations */
  cblas_dsymm(CblasRowMajor, CblasLeft, CblasUpper, n, n,
	      0.0, m->WCov[0], n, m->WCov[0], n, scale, m->WCov[0], n);

  amlet_correlation(CblasRowMajor, CblasLower, m->WCov, n);
}

/**
 * Compute the information and correlation matrices.
 */
void amlet_robust_covariance_matrix(AMLET_Model *m, int n,
				    double scale)
{
  /* We need A for cblas_dgemm... but this is a bit strange. */
  double **A = nt_matrix_new(n, n);
  nt_matrix_zeros(n, n, *A, n);

  /* The information is then computed as H^-1 I H^-1, so we need the
     Hessian */
  nt_matrix_zeros(n, n, m->RCov[0], n);

  /* A = HinvI */
  cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, n, n, 1.0,
	      m->I[0], n, m->Hinv[0], n, 0.0, *A, n);

  /* Cov = CovHinv */
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
  	      n, n, n, 1.0, *A, n, m->Hinv[0], n, 0.0, m->RCov[0],  n);

  /* We rescale the matrix to take into account the number of observations */
  cblas_dsymm(CblasRowMajor, CblasLeft, CblasUpper, n, n,
  	      0.0, m->RCov[0], n, m->RCov[0], n, scale, m->RCov[0], n);

  nt_matrix_free(A);

  amlet_correlation(CblasRowMajor, CblasLower, m->RCov, n);
}
