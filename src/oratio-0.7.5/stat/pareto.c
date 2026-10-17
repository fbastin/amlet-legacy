/**
 * pareto.c
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <math.h>
/*
La distribution de Pareto de paramètres $(\alpha,\beta)$,
dont la fonction de distribution est $F(x) = 1 - (\beta/x)^\alpha$ 
pour $x\ge \beta$.
*/
/**
 * Compute the cumulative distribution function of a Pareto variate.
 *
 * @param x is the point where cdf has to be evaluated.
 * @param alpha is the Pareto index.
 * @param beta is the mode of the distribution.
 *
 * @return the cdf of a Pareto variate of mode beta and index alpha.
 */
double st_pareto_cdf(double x, double alpha, double beta)
{
  return 1-pow(beta/x, alpha);
}

/**
 * Compute the inverse cumulative distribution function of a Pareto variate.
 *
 * @param u is the point where inverse cdf has to be evaluated.
 * @param alpha is the Pareto index.
 * @param beta is the mode of the distribution.
 *
 * @return the inverse cdf of a Pareto variate of mode beta and index alpha.
 */
double st_pareto_icdf(double u, double alpha, double beta)
{
  return beta*pow(1-u, -1.0/alpha);
}

