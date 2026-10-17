/**
 * test_mvn.c
 *
 * Random Numbers Library
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License.
 */

#include "test_ran.h"

#include <oratio/f2c.h>
#include <oratio/stat.h>

#include <stdio.h>

/**
 * Test of multivariate normal distribution generation routines.
 */
void test_mvn()
{
  int i, j;

  integer n, m;

  double **x;
  double *mu;
  double **Sigma;

  /* We use upper storage since we call a Fortran routine, so
     implicitely perform a transposition operation. */
  char uplo = 'U';
  integer info;

  Random *ran;

  printf("*****************************************\n");
  printf("* Multivariate normal test\n");
  printf("*****************************************\n");

  n = 3;
  m = 20000;

  ran = ran_normal_new(NULL, 0, 0, 1);

  x = nt_matrix_new(m, n);
  mu = _mallocs(n, double);
  Sigma = nt_matrix_new(n, n);

  mu[0] = 1.0;
  mu[1] = 0.0;
  mu[2] = -1.0;

  /* Store the covariances in the triangular upper part. The variances,
   * that should be present on the diagonal, are ignored. */
  Sigma[0][1] = Sigma[1][0] = 1.0;
  Sigma[0][2] = Sigma[2][0] = 0.5;
  Sigma[1][2] = Sigma[2][1] = 1.0;
  Sigma[0][0] = 2.0;
  Sigma[1][1] = 3.0;
  Sigma[2][2] = 2.0;

  /* Store the Cholesky factor of the variance-covariance matrix in the
   * triangular lower part. */
  /*  Sigma[0][0] = sqrt(2);
  Sigma[1][1] = 1.5811388;
  Sigma[2][2] = 1.2845233;
  Sigma[1][0] = 0.7071068;
  Sigma[2][0] = 0.3535534;
  Sigma[2][1] = 0.4743416; */

  nt_matrix_print(stdout, "Covariance matrix: ", Sigma, n, n);

  /* Cholesky factorization, stored in the lower triangular part. */
  dpotrf_(&uplo, &n, *Sigma, &n, &info);

  nt_matrix_print(stdout, "Triangular lower factorization: ", Sigma, n, n);

  for (i = 0; i < m; i++) {
    /* The storage of all the vectors is not efficient. I should add
       some statitical collectors. */    
    ran_mvn(ran, CblasRowMajor, n, *Sigma, n, mu, 1, x[i], 1);
  }

  for (i = 0; i < n; i++) {
    st_mean_variance(&(x[0][i]), n, 0, m, mu+i, &(Sigma[i][i]));
    for (j = 0; j < i; j++)
      Sigma[i][j] = Sigma[j][i] =
	st_covariance(m, &(x[0][i]), n, mu[i], &(x[0][j]), n, mu[j], 0, m);
  }

  nt_vect_print(stdout, "Empirical mean", mu, n);
  nt_matrix_print(stdout, "Empirical covariance matrix", Sigma, n, n);

  free(mu);
  nt_matrix_free(Sigma);
  nt_matrix_free(x);

  ran_random_free(ran);
}
