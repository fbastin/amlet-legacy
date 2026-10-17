/**
 * brownian.c
 *
 * ORATIO Library - Random Numbers Generation module
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License
 */

/**
 * Brownian motion 
 */

#include "oratio/ntools.h"
#include "oratio/ran_brownian.h"

double ran_brownian_inc(const Random *r,
			       double mu, double sigma,
			       double delta)
{
  return mu*delta + sigma*sqrt(delta)*ran_gaussian_icdf(r);
}

double ran_brownian_incx(const Random *r,
				double mu_delta, double sigma_sqrt_delta)
{
  return mu_delta + sigma_sqrt_delta*ran_gaussian_icdf(r);
}

/**
 * @param r is an initialized random numbers generator.
 * @param x is a double array of size greater or equal to n+1.
 * @param incx is the increment of x.
 */
void ran_brownian_generate_path(const Random *r,
				       int n, double *x, int incx,
				       double x0, double *mu,
				       double *sigma, double delta)
{
  int i;

  x[0] = x0;
  for (i = 0; i < n; i++)
    x[(i+1)*incx] = x[i*incx] +
      ran_brownian_inc(r, mu[i], sigma[i], delta);
}

/**
 * @param r is an initialized random numbers generator.
 * @param x is a double array of size greater or equal to n+1.
 * @param incx is the increment of x.
 */
void ran_brownian_generate_pathx(const Random *r,
				 int n, double *x, int incx,
				 double x0, double* mu_delta,
				 double* sigma_sqrt_delta)
{
  int i;

  x[0] = x0;
  for (i = 0; i < n; i++)
    x[(i+1)*incx] = x[i*incx] +
      ran_brownian_incx(r, mu_delta[i], sigma_sqrt_delta[i]);
}

/**
 * @param r is an initialized random numbers generator.
 * @param x is a double array of size greater or equal to n+1.
 * @param incx is the increment in x (usually 1).
 */
void ran_geometricbrownian_generate_pathx(const Random *r,
					  int n, double *x, int incx,
					  double x0, double* mu_delta,
					  double* sigma_sqrt_delta)
{
  int i;

  x[0] = log(x0);

  for (i = 0; i < n; i++) {
    x[(i+1)*incx] = x[i*incx] +
      ran_brownian_incx(r, mu_delta[i], sigma_sqrt_delta[i]);
  }

  x[0] = x0;
  for (i = 1; i <= n; i++) x[i] = exp(x[i]);
}

/**
 * Generate a Brownian path.
 *
 * @param b is an initialized Brownian motion generator.
 * @param x is a double array of size greater or equal to n+1.
 */
void ran_brownian_path(RAN_Brownian *b, double *x)
{
  ran_brownian_generate_pathx(b->r, b->T*b->d, x, 1, b->x0,
			      b->mu_delta, b->sigma_sqrt_delta);
}

/**
 * @param b is an initialized Brownian motion generator.
 */
void ran_geometricbrownian_path(RAN_Brownian *b, double *x)
{
  ran_geometricbrownian_generate_pathx(b->r, b->T*b->d, x, 1, b->x0,
				       b->mu_delta, b->sigma_sqrt_delta);
}

/**
 * @param n is the number of observation times (excluding time 0);
 */
void ran_brownian_setdelta(int n, double mu, double sigma,
			   double *obs_times, double *mu_delta,
			   double *sigma_sqrtdelta) {
  int i;
  double delta;

  for (i = 0; i < n; i++) {
    delta = obs_times[i+1] - obs_times[i];
    mu_delta[i] = mu*delta;
    sigma_sqrtdelta[i] = sigma*sqrt(delta);
  }
}

void ran_brownian_setobstimes(int T, int d, double *obs_times) {
  int i;

  obs_times[0] = 0;
  for (i = 1; i <= T*d; i++) obs_times[i] = (double)i/(double)d;
}

/**
 * Create a new Brownian motion generator.
 *
 * @param r is an initialized random number generator.
 * @param T is the time horizon.
 * @param d is the numbers of overvations per period.
 */
RAN_Brownian *ran_brownian_new(Random *r,
				  int T, int d, double x0,
				  double mu, double sigma)
{
  RAN_Brownian *b = _malloc(RAN_Brownian);

  b->r = r;

  b->T = T;
  b->d = d;

  b->x0 = x0;

  b->mu = mu;
  b->sigma = sigma;
  b->mu_delta = _mallocs(T*d, double);
  b->sigma_sqrt_delta = _mallocs(T*d, double);

  b->obs_times = _mallocs(T*d+1, double);

  ran_brownian_setobstimes(T, d, b->obs_times);
  ran_brownian_setdelta(T*d, mu, sigma, b->obs_times,
			b->mu_delta, b->sigma_sqrt_delta);

  return b;
}

/**
 * Free the memory associated to a Brownian motion generator.
 *
 * @param b is an initialized Brownian motion generator.
 */
void ran_brownian_free(RAN_Brownian *b)
{
  if (b) {
    free(b->mu_delta);
    free(b->sigma_sqrt_delta);
    free(b->obs_times);
    free(b);
  }
}

/**
 * Create a new geometric Brownian motion generator.
 *
 * @param T is the time horizon.
 * @param d is the numbers of observations per period.
 * @param drift is the percentage drift.
 */
RAN_Brownian *ran_geometricbrownian_new(Random *r, 
						  int T, int d, double x0,
						  double drift, double sigma)
{
  return ran_brownian_new(r, T, d, x0, drift - 0.5*sigma*sigma, sigma);
}

/**
 * Free the memory associated to a geometric Brownian motion generator.
 *
 * @param b is an initialized Brownian motion generator.
 */
void ran_geometricbrownian_free(RAN_Brownian *b)
{
  ran_brownian_free(b);
}

/**
 * @param r is an initialized random numbers generator.
 * @param n is the number of observation times.
 * @param d is the dimension of the random vector at each time step.
 * @param x0 is the initial level in each component.
 */
void ran_mbm_generate_pathx(const Random *r,
			    int n, int d, double **x,
			    double *x0, double *mu,
			    const enum CBLAS_UPLO uplo,
			    double *B, double ldB,
			    double *obs_times)
{
  int i, j;
  double delta, sqrtdelta;

  for (j = 0; j < d; j++) x[0][j] = x0[j];

  for (i = 1; i <= n; i++) {
    for (j = 0; j < d; j++)
      x[i][j] = ran_gaussian_icdf(r);
    cblas_dtrmv(CblasRowMajor, uplo, CblasNoTrans, CblasNonUnit,
		d, B, ldB, x[i], 1);
    delta = obs_times[i+1] - obs_times[i];
    sqrtdelta = sqrt(delta);
    for (j = 0; j < d; j++)
      x[i][j] = x[i-1][j] + mu[j]*delta+sqrtdelta*x[i][j]; 
  }
}

/**
 * @param r is an initialized random numbers generator.
 * @param n is the number of observation times.
 * @param d is the dimension of the random vector at each time step.
 * @param x0 is the initial level in each component.
 */
void ran_mgbm_generate_pathx(const Random *r,
			     int n, int d, double **x,
			     double *x0, double *mu,
			     const enum CBLAS_UPLO uplo,
			     double *B, double ldB,
			     double *obs_times)
{
  int i, j;

  for (i = 0; i < d; i++) x[0][i] = log(x0[i]);
  ran_mgbm_generate_pathx(r, n, d, x, x0, mu, uplo, B, ldB, obs_times);

  for (i = 0; i < d; i++) x[0][i] = x0[i];
  for (i = 1; i <= n; i++)
    for (j = 0; j < d; j++) x[i][j] = exp(x[i][j]);
}
