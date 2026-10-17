#include "oratio/cblas.h"
#include "oratio/ctools.h"
#include "oratio/f2c.h"
#include "oratio/clapack.h"
#include "oratio/ntools.h"
#include "oratio/ran.h"
#include <string.h>

/**
 * Create a new season (a one time period structure) for a mpar process.
 *
 * @param n is the state vector dimension.
 * @param p is the number of autoregressive steps.
 *
 * @return a newly allocated and initialized RAN_MPAR_season structure,
 *   with memory allocated to store the state vector and the Cholesky factor
 *   of the variance-covariance matrix.
 */
RAN_MPAR_Season *ran_mpar_season_new(int n, int p)
{
  int i;
  RAN_MPAR_Season *s = _malloc(RAN_MPAR_Season);

  s->n = n;
  s->mean = _mallocs(n, double);
  s->var = _mallocs(n, double);

  s->p = p;
  if (p) {
    s->par_matrices = _mallocs(n, double *);
    for (i = 0; i < p; i++) s->par_matrices[i] = _mallocs(n*n, double);
  }

  s->chol = _mallocs(n*n, double);

  s->x = _mallocs(n, double);

  return s;
}

/**
 * Free the memory of a mpar process season.
 *
 * @param s is a pointer to a filled RAN_MPAR_Season on input, and has an
 *   undetermined value on output, the consumed memory being freed.
 */
void ran_mpar_season_free(RAN_MPAR_Season *s)
{
  int i;

  if (s) {
    free(s->mean);
    free(s->var);

    for (i = 0; i < s->p; i++) free(s->par_matrices[i]);
    free(s->par_matrices);

    free(s->chol);

    free(s->x);

    free(s);
  }
}

/**
 * Set the state mean of a season of a mpar process.
 *
 * @param s is a pointer to a RAN_MPAR_Season process.
 * @param n is the state dimension of the mpar process.
 * @param m is the pointer to a double array containing the state mean.
 * @param incm is the array increment used in m.
 *
 * @return 0 is the m vector dimension does not coincide the season state
 *   structure, in which case the season content is not modified, 1 otherwise,
 *   indicating a successful mean affectation.
 */
int ran_mpar_season_set_mean(RAN_MPAR_Season *s, int n,
			     double *m, int incm)
{
  if (s->n != n) return 0;

  cblas_dcopy(n, m, incm, s->mean, 1);

  return 1;
}

/**
 * Set the state vector of a season of a mpar process.
 *
 * @param s is a pointer to a RAN_MPAR_Season process.
 * @param n is the state dimension of the mpar process.
 * @param x is the pointer to a double array containing the state.
 * @param incx is the array increment used in x.
 *
 * @return 0 is the x vector dimension does not coincide the season state
 *   structure, in which case the season content is not modified, 1 otherwise,
 *   indicating a successful state affectation.
 */
int ran_mpar_season_set_x(RAN_MPAR_Season *s, int n,
			  double *x, int incx)
{
  if (s->n != n) return 0;

  cblas_dcopy(n, x, incx, s->x, 1);

  return 1;
}

/**
 * Set the state variance of a season of a mpar process.
 *
 * @param s is a pointer to a RAN_MPAR_Season process.
 * @param n is the state dimension of the mpar process.
 * @param v is the pointer to a double array containing the state variance.
 * @param incv is the array increment used in v.
 *
 * @return 0 is the v vector dimension does not coincide the season state
 *   structure, in which case the season content is not modified, 1 otherwise,
 *   indicating a successful mean affectation.
 */
int ran_mpar_season_set_variance(RAN_MPAR_Season *s, int n,
				 double *v, int incv)
{
  if (s->n != n) return 0;

  cblas_dcopy(n, v, incv, s->var, 1);

  return 1;
}

/**
 * Set the parameters vectors of the autoregression part of a season of a
 * mpar process.
 *
 * @param s is a pointer to a RAN_MPAR_Season process.
 * @param n is the state dimension of the mpar process.
 * @param par_matrices is a vector containing the parameters used for each
 *   period of the autoregressive process, for increasing p.
 * @param p is the number of periods used in the autoregression.
 */
void ran_mpar_season_set_par_matrices(RAN_MPAR_Season *s, int n,
				      double *par_matrices, int p)
{
  int i, ns;

  ns = n*n;

  for (i = 0; i < p; i++)
    memcpy(s->par_matrices[i], par_matrices+i*ns, ns*sizeof(double));
}

/**
 * Set the Cholesky vector of the variance-covariance matrix used in a season
 * of an mpar process.
 *
 * @param s is a pointer to a RAN_MPAR_Season process.
 * @param n is the state dimension of the mpar process.
 * @param chol is a double array of size n by n, containing in the lower
 *   triangular part the Cholesky factor of a variance-covariance matrix,
 *   stored by lines.
 */
void ran_mpar_season_set_cholesky(RAN_MPAR_Season *s, int n, double *chol)
{
  memcpy(s->chol, chol, n*n*sizeof(double));
}

/**
 * Set the Cholesky vector of the variance-covariance matrix used in a season
 * of an mpar process.
 *
 * @param s is a pointer to a RAN_MPAR_Season process.
 * @param n is the state dimension of the mpar process.
 * @param cov is a double array of size n by n, containing a
 *   variance-covariance matrix, stored by lines.
 *
 * @return the info integer produced by the BLAS dpotrf routine, during the
 *   Cholesky factorization of the variance-covariance matrix.
 */
int ran_mpar_season_set_covariance(RAN_MPAR_Season *s,
				   int n, double *cov)
{
  /* We use upper storage since we call a Fortran routine, so
     implicitely perform a transposition operation. */
  char uplo = 'U';
  integer info, dim;

  dim = n;

  memcpy(s->chol, cov, n*n*sizeof(double));

  dpotrf_(&uplo, &dim, s->chol, &dim, &info);

  return info;
}

/**
 * Create a new mpar process.
 *
 * @param n is the state dimension of the mpar process.
 * @param nseasons is the number of seasons contained in the process.
 * @param p is the number of autoregressive steps.
 *
 * @return a pointer the a newly memory allocated mpar process.
 */
RAN_MPAR *ran_mpar_new(int n, int nseasons, int p)
{
  int i;
  RAN_MPAR *m = _malloc(RAN_MPAR);

  m->n = n;
  m->p = p;
  m->t = 0;

  m->x = _mallocs(n, double);

  m->nseasons = nseasons;
  m->seasons = _mallocs(nseasons, RAN_MPAR_Season *);
  for (i = 0; i < m->nseasons; i++)
    m->seasons[i] = ran_mpar_season_new(n, p);

  m->truncate = 1;
  
  return m;
}

/**
 * Free the memory of a mpar process.
 *
 * @param mpar is a pointer to an initialized RAN_MPAR process.
 */
void ran_mpar_free(RAN_MPAR *mpar)
{
  int i;

  if (mpar) {
    free(mpar->x);
    if (mpar->seasons) {
      for (i = 0; i < mpar->nseasons; i++)
	ran_mpar_season_free(mpar->seasons[i]);
      free(mpar->seasons);
    }
    free(mpar);
  }
}

/**
 * Generate a new draw from a mpar process.
 *
 * @param ran is an initialized U(0,1) random numbers generator.
 *
 * @return a pointer to the mpar current state.
 */
double *ran_mpar_next(Random *ran, RAN_MPAR *mpar)
{
  int i, j;

  RAN_MPAR_Season *season;
  double epsilon[mpar->n];

  mpar->t = (mpar->t+1)%(mpar->nseasons);

  season = mpar->seasons[mpar->t];

  for (i = 0; i < mpar->n; i++) season->x[i] = 0.0;

  // Compute the autoregressive part
  for (i = 0; i < season->p; i++) {
    j = mpar->t-i-1;
    if (j < 0) j += mpar->nseasons;
    /*void cblas_dgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY); */
    cblas_dgemv(CblasRowMajor, CblasNoTrans, season->n, season->n, 1.0,
		season->par_matrices[i], season->n, mpar->seasons[j]->x, 1,
		1.0, season->x, 1);
  }

  // White noise computation
  ran_mvn(ran, CblasRowMajor, mpar->n, season->chol, mpar->n, NULL, 1,
	  epsilon, 1);

  cblas_daxpy(mpar->n, 1.0, epsilon, 1, season->x, 1);

  // Truncate if required.
  if (mpar->truncate)
    for (i = 0; i < mpar->n; i++)
    season->x[i] = nt_max(season->x[i], -season->mean[i]);

  cblas_dcopy(mpar->n, season->mean, 1, mpar->x, 1);
  cblas_daxpy(mpar->n, 1.0, season->x, 1, mpar->x, 1);

  return mpar->x;
}

/**
 * Set the mean vector of the t-th season of a mpar process.
 *
 * @param mpar is a pointer to an initialized RAN_MPAR process.
 * @param t is the season index.
 * @param n is the mean vector dimension.
 * @param m is a pointer to the mean vector.
 * @param incm is the increment used in m.
 *
 * @param 1 is the mean affection has been successful, 0 otherwise.
 */
int ran_mpar_set_season_mean(RAN_MPAR *mpar, int t,
			     int n, double *m, int incm)
{
  return ran_mpar_season_set_mean(mpar->seasons[t], n, m, incm);
}

/**
 * Set the state vector of the t-th season of a mpar process.
 *
 * @param mpar is a pointer to an initialized RAN_MPAR process.
 * @param t is the season index.
 * @param n is the mean vector dimension.
 * @param x is a pointer to the state vector.
 * @param incx is the increment used in x.
 *
 * @param 1 is the state affection has been successful, 0 otherwise.
 */
int ran_mpar_set_season_x(RAN_MPAR *mpar, int t,
			  int n, double *x, int incx)
{
  return ran_mpar_season_set_mean(mpar->seasons[t], n, x, incx);
}

/**
 * Proceed to a warmup of a mpar process generation.
 *
 * @param ran is an initialized U(0,1) random numbers generator.
 * @param mpar is a pointer to an initialized RAN_MPAR process.
 * @param start is the season index corresponding to the first season
 *   to generate.
 * @param warmup is the warmup length. 
 */
void ran_mpar_warmup(Random *ran, RAN_MPAR *mpar, int start, int warmup)
{
  int i, j;

  int t = start-(warmup % mpar->nseasons)-1;

  if (t < 0) t += mpar->nseasons;
  mpar->t = t;

  for (i = 0; i < mpar->nseasons; i++)
    for (j = 0; j < mpar->n; j++) mpar->seasons[i]->x[j] = 0.0;

  for (i = 0; i < warmup; i++) ran_mpar_next(ran, mpar);
}
