/**
 * brownian.h
 *
 * ORATIO Library - Random Numbers Generation module
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License
 */

#ifndef __ORATIO_RAN_BROWNIAN_H__
#define __ORATIO_RAN_BROWNIAN_H__

#include "oratio/cblas.h"
#include "oratio/ran.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
    Random *r;

    double x0;

    double mu, sigma;
    double *mu_delta, *sigma_sqrt_delta;

    int T;              /* Time horizon */
    int d;              /* Number of observations per time period */

    double *obs_times;  /* Array obs_times[0..s] must contain
			   obs_times[0]=0.0, plus the n observations
			   times. */
  } RAN_Brownian;

  double ran_brownian_inc(const Random *r,
			  double mu, double sigma,
			  double delta);

  double ran_brownian_incx(const Random *r,
			   double mu_delta, double sigma_sqrt_delta);

  /**
   * @param r is an initialized random numbers generator.
   */
  void ran_brownian_generate_path(const Random *r,
				  int n, double *x, int incx,
				  double x0, double *mu,
				  double *sigma, double delta);

  /**
   * @param r is an initialized random numbers generator.
   * @param x is a double array of size greater or equal to n+1.
   */
  void ran_brownian_generate_pathx(const Random *r,
				   int n, double *x, int incx,
				   double x0, double* mu_delta,
				   double* sigma_sqrt_delta);

  /**
   * @param r is an initialized random numbers generator.
   * @param x is a double array of size greater or equal to n+1.
   * @param incx is the increment in x (usually 1).
   */
  void ran_geometricbrownian_generate_pathx(const Random *r,
					    int n, double *x, int incx,
					    double x0, double* mu_delta,
					    double* sigma_sqrt_delta);

  /**
   * Generate a Brownian path.
   *
   * @param b is an initialized Brownian motion generator.
   * @param x is a double array of size greater or equal to n+1.
   */
  void ran_brownian_path(RAN_Brownian *b, double *x);

  /**
   * @param b is an initialized Brownian motion generator.
   */
  void ran_geometricbrownian_path(RAN_Brownian *b, double *x);

  RAN_Brownian *ran_brownian_new(Random *r, int T, int d, double x0,
				    double mu, double sigma);

  /**
   * Free the memory associated to a Brownian motion generator.
   *
   * @param b is an initialized Brownian motion generator.
   */
  void ran_brownian_free(RAN_Brownian *b);

  /**
   * Create a new geometric Brownian motion generator.
   *
   * @param T is the time horizon.
   * @param d is the numbers of observations per period.
   * @param drift is the percentage drift.
   */
  RAN_Brownian *ran_geometricbrownian_new(Random *r, int T, int d,
					     double x0,
					     double drift, double sigma);
  
  /**
   * Free the memory associated to a geometric Brownian motion generator.
   *
   * @param b is an initialized Brownian motion generator.
   */
  void ran_geometricbrownian_free(RAN_Brownian *b);

  /**
   * @param r is an initialized random numbers generator.
   * @param n is the number of observation times.
   * @param d is the dimension of the random vector at each time step.
   * @param x0 is the initial level in each component.
   */
  inline void ran_mbm_generate_pathx(const Random *r,
				     int n, int d, double **x,
				     double *x0, double *mu,
				     const enum CBLAS_UPLO uplo,
				     double *B, double ldB,
				     double *obs_times);

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
			       double *obs_times);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_RAN_BROWNIAN_H__ */
