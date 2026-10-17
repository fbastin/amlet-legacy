/**
 * gamma.c
 *
 * ORATIO Library
 * Copyright (C) 2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

/**
 * Gamma distribution
 *
 * The Gamma distribution of order a>0 is defined by:
 *
 *   p(x) = 1/(\Gamma(a)b^a)x^{a-1}e^{-x/b}
 *
 * for x stricly positive.
 * Si X et Y sont des variables aléatoires indépendants de
 * distribution gamma, d'ordre a1 et a2, avec le même paramètre d'échelle b,
 * alors X+Y suit une distribution gamma d'ordre a1+a2.
 */

#include "oratio/ran.h"

#include <math.h>

double ran_gamma_int(const Random *random, const unsigned int a)
{
  double x;

  if (a < 12) {
    /* Direct method */

    unsigned int i;
    double prod = 1;

    for (i = 0; i < a; i++)
      prod *= ran_random_get_val(random);

    x = -log(prod);
  }
  else {
    /* Rejection method */

    double sqa, am, y, e;
    am = a-1;
    sqa = sqrt(2*a-1);

    do {
      do {
	y = tan(M_PI*ran_random_get_val(random));
	x = sqa*y+am;
      } while (x <= 0.0);
      e = (1.0+y*y)*exp(am*log(x/am)-sqa*y);
    } while (ran_random_get_val(random) > e);
  }

  return x;
}

double ran_gamma_frac(const Random *random, const double a)
{
  double p, q, x, u, v;

  p = M_E/(a+M_E);
  do {
    u = ran_random_get_val(random);
    v = ran_random_get_val(random); /* ran_rng_get_val_pos */

    if (u < p) {
      x = exp((1/a)*log(v));
      q = exp(-x);
    }
    else {
      x = 1-log(v);
      q = exp((a-1)*log(x));
    }
  } while (ran_random_get_val(random) >= q);

  return x;
}

/**
 * @param random is an initialized random numbers generator.
 */
double ran_gamma(const Random *random, double a, double b)
{
  /* We assume a > 0 */
  unsigned int na = floor(a);

  if (a == na)
    return b*ran_gamma_int(random, na);
  else if (na == 0)
    return b*ran_gamma_frac(random, a);
  else
    return b*(ran_gamma_int(random, na)+ran_gamma_frac(random, a-na));
}

/**
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 * @param n is the dimension of the vector seed.
 */
Random *ran_gamma_new(unsigned long *seed, int n,
		      double a, double b)
{
  Random *random = ran_random_new_with_size(seed, n, 2*sizeof(double));

  ((double *)(random->parameters))[0] = a;
  ((double *)(random->parameters))[1] = b;

  random->type = RAN_GAMMA;

  return random;
}

double ran_gamma_get_val(const Random *random)
{
  double a, b;

  ran_gamma_get_parameters(random, &a, &b);
  return ran_gamma(random, a, b);
}

void ran_gamma_get_parameters(const Random *random, double *a, double *b)
{
  *a = ((double *)(random->parameters))[0];
  *b = ((double *)(random->parameters))[1];
}

void ran_gamma_set_parameters(Random *random, double a, double b)
{
  double temp[2];

  temp[0] = a;
  temp[1] = b;

  ran_random_set_parameters(random, temp, 2*sizeof(double));
}
