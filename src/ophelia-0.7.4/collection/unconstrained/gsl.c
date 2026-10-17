#include <oratio/ran.h>

/**
 * Function E1 in simulated annealing GSL test.
 *
 * @param x is the point where the function is evaluated.
 *
 * @return the value of E1 function at x.
 */
double gsl_e1(double x)
{
  return exp(-(x-1)*(x-1))*sin(8*x)-exp(-(x-1000)*(x-1000))*0.89;
}

double gsl_e1_stoch(Random *ran, double x)
{
  double r1 = 1+ran_gaussian_icdf(ran);
  double r2 = 0.5*ran_gaussian_icdf(ran);

  return r1*exp(-(x-1)*(x-1))*sin(8*x)-exp(-(x-1000)*(x-1000))*0.89+r2;
}
