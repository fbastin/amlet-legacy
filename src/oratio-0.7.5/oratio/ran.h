/**
 * ran.h
 *
 * ORATIO Library
 * Copyright (C) 1999-2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the MIT License
 */

/**
 * Random numbers generator.
 */

#ifndef __ORATIO_RANDOM_H__
#define __ORATIO_RANDOM_H__

#include <math.h>
#include <stdlib.h>
#include <sys/types.h>

#include "ran_mpar.h"
#include "ran_rng.h"
#include "ran_walk.h"

#include "ntools.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * CONSTANTS DEFINITIONS
   */

  /**
   * Supported distributions
   */
  enum {
    RAN_CONSTANT = 0,
    RAN_ARCSINE,
    RAN_BERNOUILLI,
    RAN_CAUCHY,
    RAN_DIRAC,
    RAN_ERLANG,
    RAN_EXPONENTIAL,
    RAN_GAMMA,
    RAN_GEOMETRIC,
    RAN_GUMBEL,
    RAN_INTEGER,
    RAN_JOHNSONSB,
    RAN_JOHNSONSU,
    RAN_LAPLACE,
    RAN_LOGNORMAL,
    RAN_NORMAL,
    RAN_STUDENT,
    RAN_TRIANGULAR,
    RAN_TRUNCATED_NORMAL,
    RAN_TRUNCATED_NORMAL_B,
    RAN_UNIFORM,
    RAN_M_UNIT_SPHERE,
    RAN_SP_SPLINE,
    RAN_CONSTANT_END
  };

  /* Macros giving the number of parameters for each distribution */
#define RAN_ARCSINE_NPAR             0
#define RAN_BERNOUILLI_NPAR          1
#define RAN_CAUCHY_NPAR              2
#define RAN_DIRAC_NPAR               1
#define RAN_ERLANG_NPAR              2
#define RAN_EXPONENTIAL_NPAR         1
#define RAN_GAMMA_NPAR               2
#define RAN_GEOMETRIC_NPAR           1
#define RAN_GUMBEL_NPAR              2
#define RAN_INTEGER_NPAR             1
#define RAN_JOHNSONSB_NPAR           4
#define RAN_JOHNSONSU_NPAR           4
#define RAN_LAPLACE_NPAR             2
#define RAN_LOGNORMAL_NPAR           2
#define RAN_NORMAL_NPAR              2
#define RAN_STUDENT_NPAR             1
#define RAN_TRIANGULAR_NPAR          0
#define RAN_TRUNCATED_NORMAL_NPAR    4
#define RAN_TRUNCATED_NORMAL_B_NPAR  4
#define RAN_UNIFORM_NPAR             0

  /**
   * Error constants.
   */
  enum {
    RAN_ERROR_READ_DRAW = 1
  };

  /**
   * Produce an integer, using the intern clock
   *
   * To review!!
   */
  unsigned long ran_init_random();

  /**
   * Generic LCG
   */

  typedef struct  {
    unsigned long state;
    double a, c, m, norm;
    long q1, r1;
  } LCG;

  /**
   * Return a draw from a U(0,1) distribution, using a Linear
   * Congruential Generator.
   *
   * @param lcg is an initialized LCG generator.
   *
   * @return a U(0,1) pseudo-random number using the Standard Minimal algorithm.
   */
  double ran_lcg_get_val(void *lcg);

  /**
   * @param state is the initial state.
   */
  LCG *ran_lcg_new_with_state(unsigned long *state);

  /**
   * Create a new random number generator
   *
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   */
  LCG *ran_lcg_new_with_seed(unsigned long seed);

  /**
   * Set the parameters of the LCG, of the form (ax+c) mod m.
   */
  void ran_lcg_set_parameters(LCG *lcg, double a, double c, double m);

  /**
   * Return the current state of a standard minimal generator.
   *
   * @param lcg is an initialized LCG generator.
   */
  unsigned long *ran_lcg_get_state(void *lcg);

  /**
   * Define the generator state.
   */
  void ran_lcg_set_state(void *lcg, void *state);

  void ran_lcg_init(Random *random, unsigned long seed);

  /**
   * Create a new random number generator based on the Standard Minimal LCG.
   *
   * @param seed is the seed used to initialize the random numbers generator.
   */
  LCG *ran_lcg_new(unsigned long seed);

  /**
   * Free memory associated to the Park and Miller random numbers generator.
   *
   * @param mrg must be an initialized pointer to a LCG structure.
   */
  void ran_lcg_free(void *lcg);

  /**
   * Univariate random variables.
   */

  /**
   * Discrete random variables
   */

#include "ran_bernouilli.h"

  /**
   * Generation of points following a Dirac distribution
   */

  /**
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of the vector seed.
   */
  Random *ran_dirac_new(unsigned long *seed, int n, double a);

  void ran_dirac_get_parameters(const Random *random, double *a);

  void ran_dirac_set_parameters(Random *random, double a);
  /**
   * Generate points following a Dirac distribution, i.e. a constant.
   */
  double ran_dirac(const Random *random, double a);

  double ran_dirac_get_val(const Random *random);

  /**
   * Geometric distribution
   */

  int ran_geometric(const Random *random, double p);

  /**
   * Return a value drawn from a Geometric variate of parameter p.
   *
   * @param seed is the seed vector used to initialize the random
   *   numbers generator. A copy of this vector is stored inside the
   *  generator.
   * @param n is the dimension of seed.
   * @param p is the geometric parameter.
   */
  Random *ran_geometric_new(unsigned long *seed, int n, double p);

  int ran_geometric_get_val(const Random *random);

  double ran_geometric_get_param(const Random *random);

  void ran_geometric_set_param(Random *random, double p);

  /**
   * Generate a pseudo-random random integer in [0, max[
   *
   * @param r is a initialized random numbers generator.
   * @param max is a strictly positive integer.
   *
   * @return an integer from an uniform integer distribution on [0, max[.
   */
  unsigned long ran_intrandom(const Random *r, unsigned long max);

  /**
   * Generate a pseudo-random random integer in [0, max[
   *
   * @param ird is an initialized integer random numbers generator.
   *
   * @return an integer from an uniform integer distribution on [0, max[,
   *   where max is the parameter of the random number generator.
   */
  unsigned long ran_intrandom_get_val(const Random *random);

  /**
   * Create a new random numbers generator dedicated to the production
   * of integer random numbers.
   *
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of seed.
   */
  Random *ran_intrandom_new(unsigned long *seed, int n,
			    unsigned long max);

  unsigned long ran_intrandom_get_max(const Random *ird);

  void ran_intrandom_set_max(Random *ird, unsigned long max);

  /**
   * Continuous random variables
   */

  /**
   * Arcsine distribution
   */

  /**
   * Create a new arcine generator.
   *
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of seed.
   */
  Random *ran_arcsine_new(unsigned long *seed, int n);

  double ran_arcsine(const Random *random);

  double ran_arcsine_get_val(const Random *random);


  /**
   * Cauchy Distribution
   */

  /**
   * Draw from a Cauchy distribution.
   *
   * @param random is an initialized random numbers generator.
   * @param a is the location parameter, the median.
   * @param b is scale parameter; b must be stricty greater than 0.
   *
   * @return a (pseudo)-random number from a Cauchy distribution of
   *   location parameter a and scale parameter b.
   */
  double ran_cauchy(const Random *random, double a, double b);

  /**
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of the vector seed.
   */
  Random *ran_cauchy_new(unsigned long *seed, int n,
			 double a, double b);

  /**
   * Generate points following a standard Cauchy distribution, i.e. a
   * Cauchy distribution of parameters (0,1).
   */
  double ran_cauchy_standard_get_val(const Random *random);

  double ran_cauchy_get_val(const Random *random);

  void ran_cauchy_get_parameters(const Random *random, double *a, double *b);

  void ran_cauchy_set_parameters(Random *random, double a, double b);

  /**
   * Erlang distribution
   */

  /**
   * Exponential distribution
   */

  /**
   * Return a value drawn from an exponential of parameter lambda.
   *
   * @param ran is an initialized random number generator.
   * @param lambda is the scale parameter of the exponential distribution;
   *   it corresponds to the inverse of the mean.
   *
   * @return a number drawn from an exponential variate of parameter lambda.
   */
  double ran_exponential(const Random *ran, double lambda);

  /**
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of the vector seed.
   * @param lambda is the parameter of the exponential distribution.
   */
  Random *ran_exponential_new(unsigned long *seed, int n, double lambda);

  /**
   * Generate an exponential variate realization, given an exponential random
   * number generator.
   *
   * @param ran is an initialized random number generator, of type exponential.
   *
   * @return a realization from the exponential random number generator.
   */
  double ran_exponential_get_val(const Random *ran);

  /**
   */
  double ran_exponential_get_parameters(const Random *random);

  /**
   * @param lambda is the parameter of the exponential distribution.
   */
  void ran_exponential_set_parameters(Random *random, double b);

  /**
   * Normal distribution
   */

  /**
   * Generate a pseudo-random from a Normal distribution with mean mu
   * and standard deviation sigma.
   *
   * @param mu is the mean of the distribution.
   * @param sigma is the standard deviation of the distribution.
   *
   * @return a number drawn from a normal of mean mu and of standard
   *   deviation sigma.
   */
  double ran_normal(const Random *random, double mu, double sigma);

  /**
   * Create a new normal random numbers generator.
   *
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of seed.
   * @param mu is the mean of the distribution.
   * @param sigma is the standard deviation of the distribution.
   *
   * @return a random numbers generator for a variate X~N(mu, sigma^2).
   */
  Random *ran_normal_new(unsigned long *seed, int n, double mu, double sigma);

  /**
   * Set a random numbers generator to produce normally distributed numbers,
   * following N(mu, sigma^2).
   *
   * @param random is an initialized normal random numbers generator;
   *   it corresponds on output to a random numbers generator for a variate
   *   X~N(mu, sigma^2).
   * @param mu is the mean of the distribution.
   * @param sigma is the standard deviation of the distribution.
   */
  void ran_normal_set(Random *random, double mu, double sigma);

  double ran_gaussian_icdf(const Random *random);

  double ran_normal_icdf(const Random *random, double mu, double sigma);

  double ran_normal_get_val(const Random *random);

  void ran_normal_get_parameters(const Random *random,
				 double *mu, double *sigma);

  void ran_normal_set_parameters(Random *random, double mu, double sigma);

  /**
   * Draw a number from a standard normal distribution, using the
   * Box-Muller transformation.
   *
   * @param random is an initialized normal random numbers generator.
   *
   * @return one real number drawn from a Gaussian distribution.
   */
  double ran_gaussian_box_muller(const Random *random);

  /**
   * Draw from a N(0,1) using the acceptance-rejection method and a Laplace
   * distribution.
   *
   * @param random is an initialized normal random numbers generator.
   *
   * @return one real number drawn from a Gaussian distribution.
   */
  double ran_gaussian_ar_laplace(const Random *random);

  /**
   * Lognormal distribution
   */

  /**
   * Generate a pseudo-random from a Lognormally-distributed variate X
   * with parameters defined by the associated normal distribution.
   *
   * @param mu is the mean of ln X.
   * @param sigma is the standard deviation of ln X.
   *
   * @return a number drawn from a lognormal variate X such that ln(X)
   *   has mean mu and standard deviation sigma.
   */
  double ran_lognormal(const Random *random, double mu, double sigma);

  /**
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of seed.
   * @param mu is the mean of ln X.
   * @param sigma is the standard deviation of ln X.
   */
  Random *ran_lognormal_new(unsigned long *seed, int n,
			    double mu, double sigma);

  double ran_lognormal_icdf(const Random *random, double mu, double sigma);

  double ran_lognormal_get_val(const Random *random);

  void ran_lognormal_get_parameters(const Random *random,
				    double *mu, double *sigma);

  void ran_lognormal_set_parameters(Random *random, double mu, double sigma);

  /**
   * Gamma distribution
   */

  double ran_gamma(const Random *random, double a, double b);

  /**
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of the vector seed.
   */
  Random *ran_gamma_new(unsigned long *seed, int n,
                        double a, double b);

  double ran_gamma_get_val(const Random *random);

  void ran_gamma_get_parameters(const Random *random, double *a, double *b);

  void ran_gamma_set_parameters(Random *random, double a, double b);

  /**
   * Gumbel distribution
   */

  double ran_gumbel(const Random *random, double a, double b);

  /**
   * Create a Gumbel random variate generator.
   *
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of the vector seed.
   * @param a is the location parameter.
   * @param b is the scale parameter.
   *
   * @return a pointer to the Gumbel random numbers generator.
   */
  Random *ran_gumbel_new(unsigned long *seed, int n,
			 double a, double b);

  double ran_gumbel_get_val(const Random *random);

  void ran_gumbel_get_parameters(const Random *random, double *a, double *b);

  /**
   * Set the parameters of a Gumbel distribution.
   *
   * @param a is the value of the location parameter.
   * @param b is the value of the scale parameter.
   */
  void ran_gumbel_set_parameters(Random *random, double a, double b);

  /**
   * Johnson distribution.
   */

  /**
   * Unbounded Johnson.
   */

  /**
   * Draw from an unbounded Johnson variate.
   *
   * @param random is a generic random number generator.
   * @param xi is the localisation parameter.
   * @param lamdda is the scale parameter.
   * @param gamma and delta are the form parameters.
   *
   * @return the generated draw.
   */
  double ran_johnsonSU(const Random *random,
		       double xi, double lambda, double gamma, double delta);

  /**
   * Create a new unbounded Johnson random variable generator with some
   * initial state.
   *
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of the seed vector.
   * @param xi is the localisation parameter.
   * @param lamdda is the scale parameter.
   * @param gamma and delta are the form parameters.
   *
   * @return the newly allocated generator object.
   */
  Random *ran_johnsonSU_new(unsigned long *seed, int n,
			    double xi, double lambda,
			    double gamma, double delta);


  /**
   * Draw from an unbounded Johnson variate unsing a specialized generator.
   *
   * @param random is a random number generator created specifically for
   *   some unbounded Johnson variate.
   * @param xi is the localisation parameter.
   * @param lamdda is the scale parameter.
   * @param gamma and delta are the form parameters.
   *
   * @return the generated draw.
   */
  double ran_johnsonSU_get_val(const Random *random);

  /**
   * Retrieve the parameters of an unbounded Johnson variate.
   *
   * @param random is an initialized random numbers generator.
   * @param xi contains on output the localisation parameter.
   * @param lamdda contains on ouput the scale parameter.
   * @param gamma and delta contain on output the form parameters.
   */
  void ran_johnsonSU_get_parameters(const Random *random,
				    double *xi, double *lambda,
				    double *gamma, double *delta);

  /**
   * Set the parameters of an unbounded bounded Johnson random variate,
   *
   * @param random is an initialized random numbers generator.
   * @param xi is the localisation parameter.
   * @param lamdda is the scale parameter.
   * @param gamma and delta are the form parameters.
   */
  void ran_johnsonSU_set_parameters(Random *random,
				    double xi, double lambda, 
				    double gamma, double delta);

  /**
   * Johnson SB: bounded form
   */

  /**
   * Draw from an bounded Johnson variate.
   *
   * @param random is a generic random number generator.
   * @param xi is the localisation parameter.
   * @param lamdda is the scale parameter.
   * @param gamma and delta are the form parameters.
   *
   * @return the generated draw.
   */
  double ran_johnsonSB(const Random *random,
		       double xi, double lambda, double gamma, double delta);

  /**
   * Create a new bounded Johnson random variable generator with some
   * initial state.
   *
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of the seed vector.
   * @param xi is the localisation parameter.
   * @param lamdda is the scale parameter.
   * @param gamma and delta are the form parameters.
   *
   * @return the newly allocated generator object.
   */
  Random *ran_johnsonSB_new(unsigned long *seed, int n,
			    double xi, double lambda,
			    double gamma, double delta);

  /**
   * Draw from an bounded Johnson variate unsing a specialized generator.
   *
   * @param random is a random number generator created specifically for
   *   some bounded Johnson variate.
   * @param xi is the localisation parameter.
   * @param lamdda is the scale parameter.
   * @param gamma and delta are the form parameters.
   *
   * @return the generated draw.
   */
  double ran_johnsonSB_get_val(const Random *random);

  /**
   * Retrieve the parameters of a bounded Johnson variate.
   *
   * @param random is an initialized random numbers generator.
   * @param xi contains on output the localisation parameter.
   * @param lamdda contains on ouput the scale parameter.
   * @param gamma and delta contain on output the form parameters.
   */
  void ran_johnsonSB_get_parameters(const Random *random,
				    double *xi, double *lambda,
				    double *gamma, double *delta);
  
  /**
   * Set the parameters of a bounded Johnson random variate,
   *
   * @param random is an initialized random numbers generator.
   * @param xi is the localisation parameter.
   * @param lamdda is the scale parameter.
   * @param gamma and delta are the form parameters.
   */
  void ran_johnsonSB_set_parameters(Random *random,
				    double xi, double lambda, 
				    double gamma, double delta);

  /**
   * Laplace Distribution
   */

  /**
   * Return a value drawn from a Laplace distribution of parameter mu and b.
   *
   * @param ran is an initialized random number generator.
   * @param mu is the location parameter.
   * @param b is the scale parameter.
   *
   * @return a number drawn from an Laplace variate of parameter mu and b.
   */
  double ran_laplace(const Random *ran, double mu, double b);

  /**
   * Create a Laplace random variate generator.
   *
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of the vector seed.
   * @param mu is the location parameter.
   * @param b is the scale parameter.
   *
   * @return a pointer to the Laplace random numbers generator.
   */
  Random *ran_laplace_new(unsigned long *seed, int n,
			  double mu, double b);

  /**
   * Draw a number from a Laplace distribution.
   *
   * @return a draw from a Laplace distribution.
   */
  double ran_laplace_get_val(const Random *random);

  /**
   * Get the parameters of a Laplace distribution.
   *
   * @param mu is a pointer to a real containing on output the value of the
   *   location parameter.
   * @param b is a pointer to a real containing on output the value of the
   *   scale parameter.
   */
  void ran_laplace_get_parameters(const Random *random,
				  double *mu, double *b);

  /**
   * Set the parameters of a Laplace distribution.
   *
   * @param mu is the value of the location parameter.
   * @param b is the value of the scale parameter.
   */
  void ran_laplace_set_parameters(Random *random, double mu, double b);

  /**
   * Student Distribution
   */

  /**
   * Generate a pseudo-random from a Student distribution with the parameter nu.
   *
   * @param nu is the number of degrees of freedom of the Student distribution.
   *
   * @return a number draw from a Student with the number of degrees of
   *   freedom equal to nu.
   */
  double ran_student(const Random *random, double nu);

  /**
   * Create a new generator object dedicated to Student distribution.
   *
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of the seed vector.
   * @param nu is the number of degrees of freedom of the Student distribution.
   *
   * @return the newly allocated generator object.
   */
  Random *ran_student_new(unsigned long *seed, int n, double nu);

  /**
   * Generate a pseudo-random from a Student distribution using a
   * specialized generator.
   *
   * @return a number draw from a Student with the number of degrees of
   *   freedom equal to nu.
 */
  double ran_student_get_val(const Random *random);

  /**
   * Get the number of degrees of freedoms of a Student random variate
   * generator.
   *
   * @param random is an initialized random numbers generator, dedicated
   *   to the generation from a Student random variate with degrees of
   *   freedoms equal to nu.
   * @param nu is on output the number of degrees of freedom of the
   *   Student distribution.
   */
  void ran_student_get_parameters(const Random *random, double *nu);

  /**
   * Set the number of degrees of freedom of a Student random variate
   * generator.
   *
   * @param random is an initialized random numbers generator.
   * @param nu is the number of degrees of freedom of the Student distribution
   */
  void ran_student_set_parameters(Random *random, double nu);

  /**
   * Triangular distribution
   */

  /**
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of the seed vector.
   */
  Random *ran_triangular_new(unsigned long *seed, int n,
			     double x_min, double x_max, double mode);

  double ran_triangular(const Random *random,
			double x_min, double x_max, double mode);

  /**
   * @param u is a real in [0,1].
   */
  double ran_triangular_from_U(double u,
			       double x_min, double x_max, double mode);

  double ran_triangular_get_val(const Random *random);

  void ran_triangular_get_parameters(const Random *random,
				     double *x_min, double *x_max,
				     double *mode);

  void ran_triangular_set_parameters(Random *random,
				     double x_min, double x_max,
				     double mode);

  /**
   * Truncated normal distribution
   */

  double ran_truncatednormal(const Random *random, double a, double b,
			     double mu, double sigma);

  /**
   * @param u is a real in [0,1].
   * @param a is the lower bound on the distribution.
   * @param b is the upper bound on the distribution.
   * @param mu is the mean is the untruncated distribution.
   * @param sigma is the mean is the untruncated distribution.
   */
  double ran_truncatednormal_from_U(double u, double a, double b, 
				    double mu, double sigma);

  /**
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of the seed vector.
   * @param a is the lower bound on the distribution.
   * @param b is the upper bound on the distribution.
   * @param mu is the mean is the untruncated distribution.
   * @param sigma is the mean is the untruncated distribution.
   */
  Random *ran_truncatednormal_new(unsigned long *seed, int n,
				  double a, double b,
				  double mu, double sigma);

  double ran_truncatednormal_icdf(const Random *random, double a, double b,
				  double mu, double sigma);

  double ran_truncatednormal_get_val(const Random *random);

  void ran_truncatednormal_get_parameters(const Random *random,
					  double *a, double *b,
					  double *mu, double *sigma);

  void ran_truncatednormal_set_parameters(Random *random,
					  double a, double b,
					  double mu, double sigma);

  double ran_truncatednormal_b(const Random *random, double a, double b,
			       double mu, double sigma);

  double ran_truncatednormal_b_icdf(const Random *random,
				    double a, double b,
				    double mu, double sigma);

  double ran_truncatednormal_b_get_val(const Random *random);

  /**
   * @param u is a real in [0,1].
   */
  double ran_truncatednormal_b_from_U(double u, double a, double b, 
				      double mu, double sigma);

  void ran_truncatednormal_b_get_parameters(const Random *random,
					    double *a, double *b,
					    double *mu, double *sigma);
  
  void ran_truncatednormal_b_set_parameters(Random *random,
					    double a, double b,
					    double mu, double sigma);

  /**
   * Multivariate random variables.
   */

#include "ran_mvn.h"

  /**
   * Create a new uniform generator on the unit hypersphere with some
   * initial state.
   *
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of the seed vector.
   * @param m is the size of the random vector.
   *
   * @return the newly allocated random generator object.
   */
  Random *ran_m_unit_sphere_new(unsigned long *seed, int n, int m);

  /**
   * Draw a vector uniformly distributed on the unit hypersphere.
   *
   * @param ran is an initialized random numbers generator.
   * @param m is the dimension of the sphere.
   * @param x is on output a realization of the random vector defined on
   *   the unit sphere.
   *
   * @return x, a realization of the random vector defined uniformly on
   *   the unit sphere and the state of ran is updated.
   */
  void ran_m_unit_sphere(Random *ran, int m, double *x);

  /**
   * Draw a vector uniformly distributed on the unit hypersphere.
   *
   * @param ran is an initialized random numbers generator.
   * @param m is the dimension of the sphere.
   * @param x is on output a realization of the random vector defined on
   *   the unit sphere.
   *
   * @return x, a realization of the random vector defined uniformly on
   *   the unit sphere and the state of ran is updated.
   */
  void ran_m_unit_sphere_get_val(Random *ran, int m, double *x);

  /**
   * Multivariate normal distribution.
   */

  /**
   * Draw from a multivariate normal distribution.
   *
   * @param Order is the matrix orientation (CblasRowMajor or CblasColMajor).
   * @param ran is an initialized random numbers generator.
   * @param n is the dimension of the random vector.
   * @param L is the Cholesky factor of the variance-covariance matrix.
   * @param ldL is the leading dimension of L.
   * @param mu is the mean vector.
   * @param incmu is the increment of mu.
   * @param x is on output a realization of the multivariate normal
   *   random vector.
   * @param incx is the increment of x.
   */
  void ran_mvn(Random *ran,
	       const enum CBLAS_ORDER Order, int n, double *L, int ldL,
	       double *mu, int incmu, double *x, int incx);

  /**
   * SEMI-PARAMETRIC GENERATION.
   */

  typedef struct {
    int nk;
    double *knots;
    double *P;
    NT_BSpline *bs;
  } ran_sp_spline_params;

  void *ran_sp_spline_struct_new(int nk, double *knots, double *P);

  void ran_sp_spline_struct_free(Random *ran);

  /**
   * @param seed is the seed vector used to initialize the random
   *   numbers generator.
   * @param n is the dimension of seed.
   * @param nk is the number of points present in the knots vector
   *   (without repetitions).
   *
   * @return a random numbers generator for a semi-parametric variate X
   *   based on a cubic B-spline defined on the knots vector knots and
   *   parametrized by the vector P.
   */
  Random *ran_sp_spline_new(unsigned long *seed, int n,
			    int nk, double *knots, double *P);

  ran_sp_spline_params *ran_sp_spline_get_parameters(const Random *ran);

  void ran_sp_spline_set_parameters(Random *ran, ran_sp_spline_params *sp);

  double ran_sp_spline_U(const Random *ran, double u);

  double ran_sp_spline_get_val(const Random *ran);

  /**
   * Ratio of distributions.
   */

  /**
   * Draw from a ratio of two random distributions.
   *
   * @param ran1 is an initialized random numbers generator.
   * @param ran2 is an initialized random numbers generator, independent
   *   of ran1.
   *
   * @return a draw from a ratio of two random distributions.
   */
  double ran_generate_ratio(Random *ran1, Random *ran2,
			    double *x1, double *x2);

  /**
   * Generate realizations from a ratio of two random distributions.
   *
   * @param ran1 is an initialized random numbers generator.
   * @param ran2 is an initialized random numbers generator, independent
   *   of ran1.
   * @param n is the number of realizations to draw.
   * @param x1 is a double array, of size at least equal to n, used to
   *   store draws produced by ran1.
   * @param incx1 is the increment to use in x1.
   * @param x2 is a double array, of size at least equal to n, used to
   *   store draws produced by ran2.
   * @param incx2 is the increment to use in x2.
   * @param r is a double array, of size at least equal to n, that
   *   will contain the ratios between elements of x1 and x2.
   * @param incr is the increment to use in r.
   */
  void ran_generate_ratio_distribution(Random *ran1, Random *ran2,
				       int n, double *x1, int incx1,
				       double *x2, int incx2,
				       double *r, int incr);

  /**
   * Bootstrap.
   */

  /**
   * Generate a bootstrap sample.
   *
   * @param ran is an initialized Random numbers generators.
   * @param n is the sample size.
   * @param x is the original samble, which is not modified on output.
   * @param xb is a vector of size a least n, and contains a bootstrap
   *   sample on output.
   */
  void ran_bootstrap_sampling(const Random *ran, int n,
			     double *x, double *xb);

  /**
   * Generate a bootstrap sample.
   *
   * @param ran is an initialized Random numbers generators.
   * @param n is the sample size.
   * @param x is the original samble, which is not modified on output.
   * @param xb is a vector of size a least n, and contains a bootstrap
   *   sample on output.
   * @param f is vector with, on output, f[i] denoting the frequency of
   *   element x[i] in the bootstrap sample. 
   */
  void ran_bootstrap_sampling_freq(const Random *ran, int n,
				   double *x, double *xb, double *f);

  /**
   * Input-output functions.
   */

  /**
   * Write on output file a vector of random draws.
   *
   * @param f is the file where the draws have to be written; if NULL,
   *   the output occurs on stdout.
   * @param v is the realizations vector to write.
   * @param d is the dimension of the random vector.
   */
  void ran_write_vector(FILE *f, double *v, int d);

  /**
   * Read vectors of draws stored in a file.
   *
   * @param f is the file where the draws have to be read.
   * @param n is the number of vectors to read.
   * @param s it the dimension of each vector to read.
   * @param S is the array where the draws will be saved.
   *
   * @return RAN_ERROR_READ_DRAW is some draw ccould not be read, 0 otherwise.
   */
  int ran_import_sample(FILE *f, int n, int s, double **S);

  /**
   * Export a draws sample array.
   *
   * @param f represents the file where the draws have to be saved.
   * @param m is the number of generated random vectors.
   * @param d is the dimension of the random vectors.
   * @param S is the sample array.
   * @param ldS is the leading dimension of the draws sample array S.
   */
  void ran_export_sample(FILE *f, int m, int d, double *S, int ldS);

  /**
   * Random permutations.
   */

  /**
   * Shuffle a generic pointers vector using the Fisher-Yates shuffle.
   *
   * @param r is pseudo-random numbers generator (must be initialized).
   * @param n is the dimension of x.
   * @param x is the vector to shuffle.
   * @param incx is the increment in vector x.
   */
  void ran_shuffle_generic(Random *r, int n, void *x, int incx, int size);

  /**
   * Shuffle a double vector using the Fisher-Yates shuffle.
   *
   * @param r is pseudo-random numbers generator (must be initialized).
   * @param n is the dimension of x.
   * @param x is the vector to shuffle.
   * @param incx is the increment in vector x.
 */
  void ran_shuffle_double(Random *r, int n, double *x, int incx);

  /**
   * Generate a realization of the random variable X by mean of
   * acceptance-rejection technique.
   *
   * @param ran1 is an initialized random numbers generator.
   * @param ran2 is an initialized random numbers generator.
   * @param params_f is the parameters vector associated to the density of X.
   * @param f is the density of X.
   * @param params_t is the parameters vector for the majorant function.
   * @param t is the majorant function.
   *
   * @return a randomly generated realization of X.
   */
  double ran_acceptance_rejection(Random *ran1, Random *ran2,
				  double *params_f, double *params_t,
				  double (*f)(double x, double *params),
				  double (*t)(double x, double *params));

#include "ran_copula.h"

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_RANDOM_H__ */
