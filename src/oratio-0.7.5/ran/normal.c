/**
 * normal.c
 *
 * Random Numbers Library
 * Copyright (C) 1999-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

/**
 * Random numbers generation for normal and lognormal distributions
 */

#include "oratio/ran.h"
#include "oratio/ntools.h"
#include "oratio/stat.h"

/**
 * Normal Distribution
 */

typedef struct {
  double m, s;  /* mean and standard deviation */
  /* buff and flag are used with the Box-Muller method. */
  double buff;
  int flag;
} normal_params;

double ran_gaussian_icdf(const Random *random)
{
  return st_gaussian_icdf(ran_random_get_val(random));
}

/**
 * Generate a pseudo-random from a Normal distribution with mean mu
 * and standard deviation sigma, using inverse distribution function.
 *
 * @param mu is the mean of the distribution.
 */
double ran_normal_icdf(const Random *random, double mu, double sigma)
{
  double u;

  u = ran_random_get_val(random);

  return st_normal_icdf(u, mu, sigma);
}

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
double ran_normal(const Random *random, double mu, double sigma)
{
  return st_normal_icdf(ran_random_get_val(random), mu, sigma);
}

double ran_normal_get_val(const Random *random)
{
  return st_normal_icdf(ran_random_get_val(random),
			((normal_params *)(random->parameters))->m,
			((normal_params *)(random->parameters))->s);
}

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
void ran_normal_set(Random *random, double mu, double sigma)
{
  random->type = RAN_NORMAL;

  if (random->parameters) free(random->parameters);
  random->parameters = malloc(sizeof(normal_params));

  ((normal_params *)(random->parameters))->m = mu;
  ((normal_params *)(random->parameters))->s = sigma;
  ((normal_params *)(random->parameters))->flag = 0;
}

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
Random *ran_normal_new(unsigned long *seed, int n, double mu, double sigma)
{
  Random *random = ran_random_new(0, seed, n);

  ran_normal_set(random, mu, sigma);

  return random;
}

/**
 * @param mu is on output the mean of the distribution.
 * @param sigma is ontput the standard deviation of the distribution.
 */
void ran_normal_get_parameters(const Random *random, double *mu, double *sigma)
{
  *mu = ((normal_params *)(random->parameters))->m;
  *sigma = ((normal_params *)(random->parameters))->s;
}

/**
 * @param mu is the mean of the distribution.
 */
void ran_normal_set_parameters(Random *random, double mu, double sigma)
{
  normal_params np;

  np.m = mu;
  np.s = sigma;
  np.flag = 0;

  ran_random_set_parameters(random, &np, sizeof(np));
}

/**
 * Draw a number from a standard normal distribution, using the
 * Box-Muller transformation.
 *
 * @param random is an initialized normal random numbers generator.
 *
 * @return one real number drawn from a Gaussian distribution.
 */
double ran_gaussian_box_muller(const Random *random)
{
  if (((normal_params *)(random->parameters))->flag) {
    ((normal_params *)(random->parameters))->flag = 0;

    return ((normal_params *)(random->parameters))->buff;
  }
  else {
    double r, theta;

    ((normal_params *)(random->parameters))->flag = 1;

    theta = 2*M_PI*ran_random_get_val(random);
    r = sqrt(-2*log(ran_random_get_val(random)));

    ((normal_params *)(random->parameters))->buff = r*sin(theta);

    return r*cos(theta);
  }
}

/**
 * Draw from a N(0,1) using the acceptance-rejection method and a Laplace
 * distribution.
 *
 * @param random is an initialized normal random numbers generator.
 *
 * @return one real number drawn from a Gaussian distribution.
 */
double ran_gaussian_ar_laplace(const Random *random)
{
  double y;
  double r;

  do {
    y = ran_laplace(random, 0.0, 1.0);
    r = (y > 0.0) ? exp(-0.5*(y-1)*(y-1)) : exp(-0.5*(y+1)*(y+1));
  } while (ran_random_get_val(random) > r);

  return y;
}

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
double ran_lognormal(const Random *random, double mu, double sigma)
{
  double x;

  x = ran_gaussian_icdf(random);

  return exp(x*sigma + mu);
}

double ran_lognormal_icdf(const Random *random, double mu, double sigma)
{
  double u;

  u = ran_random_get_val(random);

  return st_lognormal_icdf(u, mu, sigma);
}

/**
 * Create a new random numbers generator dedicated to the lognormal
 * distribution.
 *
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 * @param n is the dimension of seed.
 * @param mu is the mean of ln X.
 * @param sigma is the standard deviation of ln X.
 */
Random *ran_lognormal_new(unsigned long *seed, int n,
			  double mu, double sigma)
{
  Random *r;

  r = ran_normal_new(seed, n, mu, sigma);
  r->type = RAN_LOGNORMAL;

  return r;
}

double ran_lognormal_get_val(const Random *random)
{
  return exp(ran_normal_get_val(random));
}


void ran_lognormal_get_parameters(const Random *random,
				  double *mu, double *sigma)
{
  *mu = ((normal_params *)(random->parameters))->m;
  *sigma = ((normal_params *)(random->parameters))->s;
}

void ran_lognormal_set_parameters(Random *random,
				  double mu, double sigma)
{
  ran_normal_set_parameters(random, mu, sigma);
}
