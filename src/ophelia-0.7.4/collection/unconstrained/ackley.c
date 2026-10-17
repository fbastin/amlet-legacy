#include <oratio/constants.h>

/**
 * Ackley test function.
 *
 * References:
 * - 2D form
 *   David H. Ackley, "A connectionist machine for genetic hillclimbing", Kluwer Academics
 *   Publishers, Boston, MA, USA, 1987.
 * - Generalization to n dimensions
 *   Thomas Bäck, "Evolutionary algorithms in theory and practice: evolution strategies,
 *   evolutionary programming, genetic algorithms", Oxford University Press, Oxford, United
 *   Kingdom, 1996.
 *
 * @param n is the dimension of the vector x.
 * @param x is the point where the function is evaluated.
 * @param a is the first parameter of Ackley function.
 * @param b is the second parameter of Ackley function.
 * @param c is the third parameter of Ackley function.
 * Typical variable values are a=20, b=0.2, c=2pi
 *
 * @return the Ackley function value at x, given the parameters a, b, and c.
 */
double ackley(int n, double *x, double a, double b, double c)
{
  int i;
  double s1, s2;

  s1 = s2 = 0.0;

  for (i = 0; i < n; i++) {
    s1 += x[i]*x[i];
    s2 += cos(c*x[i]);
  }

  return -a*exp(-b*sqrt(s1/n))-exp(s2/n)+a+M_E;
}
