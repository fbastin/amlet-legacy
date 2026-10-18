/**
 * moca.h
 *
 * ORATIO Library
 * Copyright (C) 2004-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_MOCA_H__
#define __ORATIO_MOCA_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "ctools.h"
#include "ran.h"

#include "faure.h"
#include "halton.h"
#include "mlhs.h"

#include "error.h"

  enum {
    MOCA_MONTE_CARLO,
    MOCA_FAURE,
    MOCA_HALTON,
    MOCA_RANK1_LATTICE,
    MOCA_MLHS,
    MOCA_SOBOL,
    MOCA_METHODS_END
  };

  static const c_symbol moca_methods[] = {
    {"monte_carlo", MOCA_MONTE_CARLO},
    {"faure", MOCA_FAURE},
    {"halton", MOCA_HALTON},
    {"lattice", MOCA_RANK1_LATTICE},
    {"mlhs", MOCA_MLHS},
    {"sobol", MOCA_SOBOL}
  };

  /**
   * Return a pointer to the set of moca methods symbols.
   *
   * @return a pointer to the methods symbols array.
   */
  c_symbol *moca_get_methods();

  /**
   * Give the number of implemented methods.
   *
   * @return the number of implemented methods.
   */
  int moca_methods_number();

  /**
   * Return the name of the sampling method.
   *
   * @param method is a sampling method identifier.
   *
   * @param method the name of the sampling method associated to the
   *   given identifier, or NULL if this identifier is invalid.
   */
  char *moca_get_method_name(int method);

  /**
   * Generate a Monte Carlo sample.
   *
   * @param ran is an initialized random numbers generator.
   * @param n is the number of vectors to generate.
   * @param s is the dimension of each vector.
   * @param sample contains on output the Monte Carlo sample; sample has
   *   to be of size n by s.
   */
   void moca_mc_generate_sample(Random *ran, int n, int s, double **sample);

  /**
   * Count the parameters associated to a random distribution law.
   *
   * @param law is an identifier of the random distribution law.
   * @param data is a pointer to specific law information.
   *
   * @return the number of parameters associated to a random distribution law.
   */
  int moca_count_parameters(int law, void *data);

  /**
   * Determine the number of parameters in a vector and the associated
   * number of random coefficients.
   *
   * @param n is the vector dimension.
   * @param pars is the parameters type vector.
   * @param n_random_coeff points to a integer representing the number
   *   of random coefficients.
   * @param idx_rnd is a vector of indexes indicating the positions of
   *   the random coefficients on output, if idx_rnd is not NULL,
   *   otherwise it is not filled.
   */
  int moca_count_total_parameters(int n, int *pars, int idx,
				  int *npars, int *n_random_coeff,
				  int *idx_rnd);

  /**
   * Construct an array that summarizes the basic types for random variables
   */
  void moca_set_sampling_types(int *var_types, int *sampling_types,
			       int n_random_coeff);

  /**
   * Set a particular realisation inside a realizations vector.
   *
   * @param theta is the parameters vector.
   * @param type is the distribution law of nonparametric method to use.
   */
  int moca_set_realisation(double *beta, double *theta, double draw,
			   int type, int pos, void *data);

  /**
   * Baker's tranformation.
   */

  /**
   * Apply Baker's transformation to the scalar u.
   *
   * @param u is the scalar to transform; u must be in [0,1).
   *
   * @param return the transformed scalar.
   */
  double moca_baker(double u);

  /**
   * Apply Baker's transformation to vector u.
   *
   * @param s is the dimension of vector u.
   * @param u is a vector in [0,1)^s, and is the vector after Baker's
   *   transformation on output.
   * @param inc is the increment of u.
   */
  void moca_baker_point(int s, double *u, int incx);

  /**
   * Input/output routines.
   */

  /**
   * Write a draws vector in output file.
   *
   * @param f is the file where the draws have to be written; if NULL,
   *   the output occurs on stdout.
   * @param v is the vector to write.
   * @param s is the dimension of v.
   */
  void moca_print_vector(FILE *f, double *v, int s);

  /**
   * Randomization techniques.
   */

  /**
   * Randomly shift a vector of points, modulo one.
   *
   * @param r is a pseudo-random numbers generator (must be initialized).
   * @param pts is the vector of points to shift. On output, each
   *   component of pts has been shifted by some random value (identical
   *   between components), modulo one.
   * @param n is the size of pts.
   */
  void moca_shift_modulo_one(Random *r, double *pts, int n);

  /**
   * @param s is the dimension of x and the shift vector.
   */
  void moca_shift_apply(int s, double *x, double *shift);

  /**
   * Initialize a random shift vector with a padding for 32-bit random
   * numbers generators.
   *
   * @param r is a pseudo-random numbers generator (must be
   *   initialized).
   * @param s is the dimension of the shift vector.
   * @param shift contains on output the random shift vector.
   */
  void moca_shift_init32(Random *r, int s, double *shift);

  /**
   * Reproduce, or not, the random shifts of ORATIO <= 0.7.5, whose
   * padding of 2^-33 was subtracted instead of added (an integer
   * overflow).
   *
   * @param legacy is 0 (default) or 1 to reproduce ORATIO <= 0.7.5.
   */
  void moca_set_legacy(int legacy);

  /**
   * Initialize a random shift vector with a padding for 32-bit random
   * numbers generators; each component is set to the same number.
   *
   * @param r is a pseudo-random numbers generator (must be
   *   initialized).
   * @param s is the dimension of the shift vector.
   * @param shift contains on output the random shift vector.
   */
  void moca_shift_uniform_init32(Random *r, int s, double *shift);

  /**
   * Initialize a random shift vector.
   *
   * @param r is a pseudo-random numbers generator (must be
   *   initialized).
   * @param s is the dimension of the shift vector.
   * @param shift contains on output the random shift vector.
   */
  void moca_shift_init(Random *r, int s, double *shift);

  /**
   * RQMC
   */

  /**
   * Rank-1 lattice nets.
   */

  typedef struct {
    int s;                /* lattice point dimension */
    int n;                /* number of points in the lattice */
    double normfactor;    /* scaling factor (1.0/n) */
    double *v;            /* generator vector */
    int ptindex;          /* point index */
    int coordindex;       /* coordinate index */
  } MOCA_Rank1Lattice;

  /**
   * Create a new rank-1 lattice.
   *
   * @param n is the number of points.
   * @param s dimension of the lattice vector a.
   * @param a the lattice vector; each component must be in [0,n[.
   *
   * @return a newly memory-allocated rank-1 lattice.
   */
  MOCA_Rank1Lattice *moca_rank1_lattice_new(int n, int s, int *a);

  /**
   * Free the memory associated to a rank-1 lattice->
   *
   * @param lattice is a memory-allocated rank-1 lattice->
   */
  void moca_rank1_lattice_free(MOCA_Rank1Lattice *lattice);

  /**
   * Compute the jth coordinate of the ith rank-1 lattice point.
   *
   * @param lattice is an initialixed rank-1 lattice.
   * @param i is the point index.
   * @param is the coordinate index.
   *
   * @return the jth coordinate of the ith rank-1 lattice point.
   */
  double moca_rank1_lattice_get_ij(MOCA_Rank1Lattice *lattice,
				   int i, int j);

  /**
   * Compute the current coordinate of the current lattice point.
   *
   * @param lattice is an initialixed rank-1 lattice->
   *
   * @return the next coordinate of the current lattice point.
   */
  double moca_rank1_lattice_get_coordinate(MOCA_Rank1Lattice *lattice);

  /**
   * Compute the next coordinate of the current lattice point.
   *
   * @param lattice is an initialixed rank-1 lattice->
   *
   * @return the next coordinate of the current lattice point.
   */
  double moca_rank1_lattice_next_coordinate(MOCA_Rank1Lattice *lattice);

  /**
   * Compute the point of a rank-1 lattice at the current point index.
   *
   * @param lattice is an initialixed rank-1 lattice->
   * @param x is the vector that will contain the point.
   * @param incx is the increment of x.
   */
  void moca_rank1_lattice_get_point(MOCA_Rank1Lattice *lattice,
				    double *x, int incx);

  /**
   * Compute the point of a rank-1 lattice at the next point index,
   * which is incremented (after point computation).
   *
   * @param lattice is an initialixed rank-1 lattice->
   * @param x is the vector that will contain the point.
   * @param incx is the increment of x.
   */
  void moca_rank1_lattice_next(MOCA_Rank1Lattice *lattice,
			       double *x, int incx);

  /**
   * Increment the point index of a rank-1 lattice->
   *
   * @param lattice is an initialixed rank-1 lattice->
   */
  void moca_rank1_lattice_inc_point(MOCA_Rank1Lattice *lattice);

  /**
   * Reset the current rank-1.
   *
   * @param lattice is an initialixed rank-1 lattice.
   */
  void moca_rank1_lattice_reset(MOCA_Rank1Lattice *lattice);

  /**
   * Generate a rank-1 lattice sample.
   *
   * @param lattice is a rank-1 lattice generator.
   * @param sample contains on output the latticex sample.
   * @param n is the numbers of vectors forming the sample.
   * @param ldS is the leading dimension of S; ldS must be grater or equal to
   *   lattice points dimension, and is typically equal to this dimesion.
   */
  void moca_rank1_lattice_sample(MOCA_Rank1Lattice *lattice,
				 int n, double *S, int ldS);

  /**
   * Sobol sequences.
   */

  typedef struct {
    int max_dim;        /* maximum dimension of generated vectors. */
    int max_digits;
    int *mdeg;          /* degree of primitive polynomials. */
    int maxdeg;         /* maximum degree amongst used primitive polynomials. */
    int *mi;
    int *vi;
    int in;             /* index of the generated vector. */
    int *ix;            /* current point. */
    int *iu;
    double fac;
    int s;              /* user-defined dimension. */
  } MOCA_Sobol;

  enum {
    SOBOL_BRATLEY_FOX
  };

  /**
   * Generate a Sobol sample.
   *
   * @param sobol is an initialized Sobol quasi-random numbers generator.
   * @param n is the numbers of vectors forming the sample.
   * @param s is the quasi-random vector dimension.
   * @param S contains on output the Sobol sample.
   * $param ldS is the leading dimension of S.
   */
  void moca_sobol_sample(MOCA_Sobol *sobol, int n, int s,
			 double *S, int lds);

  void moca_sobol_init(MOCA_Sobol *sbl);

  /**
   * @param s is the dimension of the vectors to generate.
   */
  MOCA_Sobol *moca_sobol_new(int pid, int s);

  void moca_sobol_free(MOCA_Sobol *sbl);

  void moca_sobol_reset(MOCA_Sobol *sbl);

  void moca_sobol_next(MOCA_Sobol *sbl, double *x);

  /**
   * Distributions.
   */

  enum {
    CONSTANT = 0,
    FIXED,
    BERNOUILLI,
    UNIFORM,
    NORMAL,
    LOGNORMAL,
    TRIANGULAR,
    TRUNCATED_NORMAL,
    TRUNCATED_NORMAL_B,
    ERRC,
    MOCA_END_DISTRIBUTIONS
  };

  static const c_symbol distributions[] = {
    {"constant", CONSTANT},
    {"fixed", FIXED},
    {"bernouilli", BERNOUILLI},
    {"uniform", UNIFORM},
    {"normal", NORMAL},
    {"lognormal", LOGNORMAL},
    {"triangular", TRIANGULAR},
    {"truncated_normal", TRUNCATED_NORMAL},
    {"truncated_normal_b", TRUNCATED_NORMAL_B},
    {"errc", ERRC}
  };
  
  /**
   * Nonparametric distributions.
   */
  /*  enum {
    SPLINE = MOCA_END_DISTRIBUTIONS
    }; */

  static const c_symbol nonparametric_methods[] = {
    {"spline", RAN_SP_SPLINE}
  };

  double moca_convert_random_value(int type, double v);

  /**
   * @param law is the probability law identifier.
   * @param r is an initialized random numbers generator.
   * @param data is a pointer to a structure data corresponding to law,
   *   if necessary, otherwise it is set to NULL.
   */
  double moca_get_random_value(int law, Random *r, double *param,
			       void *data);

  /**
   * @param law is the probability law identifier.
   * @param data is a pointer to a structure data corresponding to law,
   *   if necessary, otherwise it is set to NULL.
   */
  double moca_affect_value(int law, double u, double *param, void *data);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_MOCA_H__ */
