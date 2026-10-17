/**
 * ran_mpar.h
 *
 * ORATIO Library - Random Numbers Generation module
 * Copyright (C) 2015 Fabian Bastin <fabian.bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License
 */

/**
 * MPAR(p) stochastic process
 */

#ifndef __ORATIO_RAN_MPAR_H__
#define __ORATIO_RAN_MPAR_H__

#include "ran_rng.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct _ran_mpar_season {
    int n;                     /* state vector dimension */
    double *mean;
    double *var;

    int p;
    double **par_matrices;
    double *chol;

    double *x;                 /* state at this season */
  } RAN_MPAR_Season;

  typedef struct _ran_mpar {
    int n;                     /* state vector dimension */
    double *x;                 /* state vector */

    int p;

    int nseasons;              /* number of seasons */
    RAN_MPAR_Season **seasons; /* seasons */

    int t;                     /* current season index */

    int truncate;              /* 1 if the state vector has only nonnegative
				  values, 0 otherwise */
  } RAN_MPAR;

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
  RAN_MPAR_Season *ran_mpar_season_new(int n, int p);

  /**
   * Free the memory of a mpar process season.
   *
   * @param s is a pointer to a filled RAN_MPAR_Season on input, and has an
   *   undetermined value on output, the consumed memory being freed.
   */
  void ran_mpar_season_free(RAN_MPAR_Season *s);

  int ran_mpar_season_set_mean(RAN_MPAR_Season *s, int n,
			       double *m, int incm);

  int ran_mpar_season_set_variance(RAN_MPAR_Season *s, int n,
				   double *v, int incv);

  void ran_mpar_season_set_par_matrices(RAN_MPAR_Season *s, int n,
					double *par_matrices, int p);

  /**
   * @param chol is the cholesky factor matrix. The matrix is assumed to be
   *   lower triangular and stored by rows.
   */
  void ran_mpar_season_set_cholesky(RAN_MPAR_Season *s,
				    int n, double *chol);

  int ran_mpar_season_set_covariance(RAN_MPAR_Season *s,
				      int n, double *cov);

  RAN_MPAR *ran_mpar_new(int n, int nseasons, int p);

  void ran_mpar_free(RAN_MPAR *mpar);

  int ran_mpar_set_season_mean(RAN_MPAR *mpar, int t,
			       int n, double *m, int incm);

  void ran_mpar_warmup(Random *ran, RAN_MPAR *mpar, int start, int warmup);

  double *ran_mpar_next(Random *ran, RAN_MPAR *mpar);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_RAN_WALK_H__ */
