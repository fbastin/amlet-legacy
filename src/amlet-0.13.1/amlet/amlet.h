/**
 * amlet.h
 *
 * AMLET (Another Mixed Logit Estimation Tool) Library
 * Copyright (C) 2002-2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the MIT License
 */

#ifndef __AMLET_LIBRARY_H__
#define __AMLET_LIBRARY_H__

#include <oratio/ctools.h>
#include <oratio/ran.h>

#include <amlet/amlet_driver.h>
#include <amlet/credits.h>
#include <amlet/data.h>
#include <amlet/amlet_driver.h>
#include <amlet/info.h>
#include <amlet/information.h>
#include <amlet/model.h>
#include <amlet/lattice.h>
#include <amlet/likelihood.h>
#include <amlet/nonparametric.h>
#include <amlet/validation.h>
#include <amlet/variable.h>

#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Default file names */
#define AMLET_DEFAULT_LOG          "amlet.log"
#define AMLET_DEFAULT_PROBLEM_NAME "unnamed"

  /**
   * Default values
   */

#define AMLET_DEFAULT_NITER 1000
#define AMLET_DEFAULT_TOL   1e-4
#define AMLET_DEFAULT_LEVEL 0.05

  /**
   * Error and warnings functions
   */

  /**
   * @param f is the file where the warning will be issued.
   * @param err is the error code.
   */
  int amlet_warning_va(FILE *f, int err, char *type,
		       const char *msg, va_list ap);

  /**
   * @param f is the file where the warning will be issued.
   */
  void amlet_warning(FILE *f, int verbose, char *msg, ...);

  /**
   * Allocate memory for (quasi-)random draws.
   *
   * @param pop is the population size.
   * @param n is the number of random vectors to draw per individual.
   * @param s is the size of each vector.
   *
   * @return an allocated tabular of the required size to store the
   *   draws, with each line corresponding to one individual.
   */
  double **amlet_sample_alloc(int pop, int n, int s);
  
  /**
   * Deallocate memory associated to (quasi-)random draws.
   *
   * @param S is the (allocated) array of draws.
   * @param pop is the population size.
   */
  void amlet_sample_free(double **S, int pop);

  /**
   * Sample generation for use in (quasi-)Monte Carlo simulation
   *
   * @param f is the file where the U(0,1) draws will be saved, if not NULL.
   * @param ran is an initialized random numbers generator.
   * @param method is the sampling method.
   * @param n is the number of vectors to draw.
   * @param s is the size of draws vectors.
   * @param sample contains, on output, the generated sample.
   * @param randomization adds scrambling or shuffling to the randomly
   *   shifted Halton sequences (MOCA_HALTON_RANDOMIZATION_*); it is
   *   ignored by the other methods.
   */
  void amlet_sample_generate_01(FILE *f, Random *ran, int np, int method,
				int n, int s, double **sample,
				int randomization);

  /**
   * @param f is the file where the U(0,1) draws will be saved, if not NULL.
   * @param ran is an initialized random numbers generator.
   */
  void amlet_sample_generate(FILE *f, Random *ran, AMLET_Driver *dr, int newly,
			     double **sample, int *types, int n_random_coeff);

  /**
   * Compute the Hessian associated to the likelihood of the model.
   *
   * @param m is an initialized AMLET model.
   */
  double *amlet_model_likelihood_hessian(AMLET_Model *m,
					 double *x, double *grad,
					 int n_random_coeff, dcl_data *md,
					 double *work);

  /**
   * Very basic starting point initialization. You should consider to
   * replace this function by your own if you have a better strategy for
   * your particular problem.
   *
   * @param n is the number of parameters to estimate.
   * @param x is the parameters vectors, set to an initial value on output.
   * @param incx is the increment to use in x.
   */
  void amlet_ml_starting_point(int n, double *x, int incx);

  /**
   * INPUT-OUTPUT functions.
   */

  /**
   * Read the repetitions file.
   *
   * @param repname is the name of the file containing the repetitions. The
   *   ith entry in the file represents the number of observations associated
   *   to individual i.
   * @param rep is a pointer to an integer array of size at least
   *   equal to np+1. On output rep[i+1] will the contain the number of
   *   observations associated to indivual i.
   * @param np is number of individuals.
   *
   * @return 0 is the lecture has encountered no difficulties, else an
   *   error identifier.
   */
  int amlet_read_repetitions(char *repname, int *rep, int np);

  /**
   * MEMORY GESTION
   */

  /**
   * Some memory allocation for data structures used in discrete choice
   * models estimation.
   *
   * @param nobs is the number of observations.
   * @param nalt is the number of alternatives.
   * @param censor is equal to 1 if all alternatives are available for
   *   all users, 0 else.
   */
  void dcl_memory_allocate(int nobs, int nalt, int nvar, int censor,
			   int **choices, double ****obs,
			   int ***availability);
  /**
   * @param nobs is the number of observations.
   * @param nalt is the number of alternatives.
   */
  void dcl_memory_deallocate(int nobs, int nalt,
			     int *choices, double ***obs, int **availability);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_LIBRARY_H__ */
