/**
 * truncatednormal.c
 *
 * Random Numbers Library
 * Copyright (C) 2007-2008 Cinzia Cirillo <ccirillo@umd.edu>,
 *                         Fabian Bastin <bastin@iro.umontreal.ca>
 *                         Pierre Derian <pierre.derian@gmail.com>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Truncated Normal.
 *
 * Reference: "Simulation: aspects stochastiques".
 * Curse given at University of Montreal
 */

#include "oratio/ran.h"
#include "oratio/stat.h"

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
				double mu, double sigma)
{
  Random *random = ran_random_new_with_size(seed, n, 4*sizeof(double));

  ((double *)(random->parameters))[0] = a;
  ((double *)(random->parameters))[1] = b;
  ((double *)(random->parameters))[2] = mu;
  ((double *)(random->parameters))[3] = sigma;

  random->type = RAN_TRUNCATED_NORMAL;

  return random;
}

/**
 * Generate a normal truncated on the interval [a,b].
 */
double ran_truncatednormal(const Random *random, double a, double b,
			   double mu, double sigma)
{
  double u, p;
  double x_min, x_max;

  p = random->rng_get_val(random->rng);
  x_min = st_normal_cdf(a, mu, sigma);
  x_max = st_normal_cdf(b, mu, sigma);

  u = x_min+p*(x_max - x_min);

  return st_normal_icdf(u, mu, sigma);
}

double ran_truncatednormal_icdf(const Random *random,
				       double a, double b,
				       double mu, double sigma)
{
  return ran_truncatednormal(random, a, b, mu, sigma);
}

double ran_truncatednormal_get_val(const Random *random)
{
  double a, b, mu, sigma;

  ran_truncatednormal_get_parameters(random, &a, &b, &mu, &sigma);
  return ran_truncatednormal(random, a, b, mu, sigma);
}

void ran_truncatednormal_get_parameters(const Random *random,
				        double *a, double *b,
                                        double *mu, double *sigma)
{
  *a = ((double *)(random->parameters))[0];
  *b = ((double *)(random->parameters))[1];
  *mu = ((double *)(random->parameters))[2];
  *sigma = ((double *)(random->parameters))[3];
}

void ran_truncatednormal_set_parameters(Random *random,
			                double a, double b,
                                        double mu, double sigma)
{
  double temp[4];

  temp[0] = a;
  temp[1] = b;
  temp[2] = mu;
  temp[3] = sigma;

  ran_random_set_parameters(random, temp, 4*sizeof(double));
}

/**
 * @param a is the lower bound on the distribution.
 * @param b is the upper bound on the distribution.
 * @param mu is the mean is the untruncated distribution.
 * @param sigma is the mean is the untruncated distribution.
 */
double ran_truncatednormal_from_U(double u, double a, double b, 
				  double mu, double sigma)
{  
  double  p;
  double x_min, x_max;

  x_min = st_normal_cdf(a, mu, sigma);
  x_max = st_normal_cdf(b, mu, sigma);

  p = x_min+u*(x_max - x_min);

  return st_normal_icdf(p, mu, sigma);

}

/**
 * Draw a random number from at truncated normal distribution with
 * probability mass for truncature points.
 *
 * @param random is an initialized random number generator.
 * @param a is the left truncature.
 * @param b is the right truncature.
 * @param mu is the distribution mean.
 * @param sigma is the distribution standard deviation.
 *
 * @returns a random draw from at truncated normal distribution with
 * probability mass for truncature points.
 */
double ran_truncatednormal_b(const Random *random, double a, double b,
			     double mu, double sigma)
{
  double p;

  p = random->rng_get_val(random->rng);
  p = st_normal_icdf(p, mu, sigma);
  
  if (p<a) p = a;
  else if (p>b) p = b;

  return p;
}

double ran_truncatednormal_b_icdf(const Random *random,
					 double a, double b,
					 double mu, double sigma)
{
  return ran_truncatednormal_b(random, a, b, mu, sigma);
}

double ran_truncatednormal_b_get_val(const Random *random)
{
  double a, b, mu, sigma;

  ran_truncatednormal_b_get_parameters(random, &a, &b, &mu, &sigma);
  return ran_truncatednormal_b(random, a, b, mu, sigma);
}
 
double ran_truncatednormal_b_from_U(double u, double a, double b, 
				    double mu, double sigma)
{
  double p;

  p = st_normal_icdf(u, mu, sigma);
  
  if (p < a) p = a;
  else if (p > b) p = b;

  return p;
}

void ran_truncatednormal_b_get_parameters(const Random *random,
					  double *a, double *b,
					  double *mu, double *sigma)
{
  *a = ((double *)(random->parameters))[0];
  *b = ((double *)(random->parameters))[1];
  *mu = ((double *)(random->parameters))[2];
  *sigma = ((double *)(random->parameters))[3];
}

void ran_truncatednormal_b_set_parameters(Random *random,
					  double a, double b,
					  double mu, double sigma)
{
  double temp[4];

  temp[0] = a;
  temp[1] = b;
  temp[2] = mu;
  temp[3] = sigma;

  ran_random_set_parameters(random, temp, 4*sizeof(double));
}
