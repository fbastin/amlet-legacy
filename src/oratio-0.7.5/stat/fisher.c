/**
 * fisher.c
 *
 * STAT Library
 * Copyright (C) 2001-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <math.h>

/**
 * Fisher-Snedecor distribution.
 */

/**
 * Evaluation of the Fisher-Snedecor function with parameters a and b,
 * at point x.
 */
double st_fisher_pdf(double x, double a, double b)
{
  double f;

  f = lgamma(0.5*(a+b))+0.5*a*log(a)+0.5*b*log(b)+(0.5*a-1)*log(x)-
    lgamma(0.5*a)-lgamma(0.5*b)-0.5*(a+b)*log(b+a*x);
  
  return exp(f);
}

/**
 * Moyenne de la fonction de Fisher-Snedecor de paramètres a, b
 *
 * #PRE: b > 2
 */
double st_fisher_mean(double x, double a, double b)
{
  return b/(b-2);
}

/**
 * Variance de la fonction de Fisher-Snedecor de paramètres a, b
 *
 * #PRE: b > 4
 */
double st_fisher_variance(double x, double a, double b)
{
  return 2*a*a*(a+b-2)/(a*(b-4)*(b-2)*(b-2));
}
