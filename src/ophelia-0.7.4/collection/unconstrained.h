/**
 * unconstrained.h
 *
 * OPHELIA Library - unconstrained test problems
 * Copyright (C) 2006-2015 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __OP_COLLECTION_UNCONSTRAINED_H__
#define __OP_COLLECTION_UNCONSTRAINED_H__

#include <oratio/ran.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Unidimensional problems
   */

  /**
   * Evaluate the function x²+4x-10.
   *
   * @param x is the point where the function has to be evaluated.
   * @param data is unused.
   *
   * @return the value of x²+4x-10.
   */
  double op_parabola_1(double x, void *data);

  /**
   * Nonlinear optimization test problems.
   */

  /**
   * Problem CGT-1
   *
   * Source: Conn, Gould, and Toint, "Trust-Region Methods", SIAM, 2000.
   */

  /**
   * Problem CGT-1 objective function.
   *
   * @param data is unused.
   *
   * Source: Conn, Gould, and Toint, "Trust-Region Methods", SIAM, 2000.
   */
  double cgt_1(double *x, int n, void *data);

  /**
   * Problem CGT-1 gradient.
   *
   * @param data is unused.
   */
  void cgt_1_deriv1(double *x, int n, void *data, double *grad);
  
  /**
   * CGT-1 Hessian.
   */
  void cgt_1_deriv2(double *x, int n, void *data, double *H);

  /**
   * Rosenbrock function.
   */

  double rosenbrock(double *x, int n, void *data);

  void rosenbrock_deriv1(double *x, int n, void *data, double *grad);

  /**
   * NLP-UNC1 problem
   */

  double unc1(double *x, int n, void *data);

  void unc1_grad(double *x, int n, void *data, double *grad);

  double unc2(double *x, int n, void *data);

  void unc2_grad(double *x, int n, void *data, double *grad);

  /**
   * Global optimization test problems.
   */

  /**
   * Ackley problem
   */

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
  double ackley(int n, double *x, double a, double b, double c);

  /**
   * Chichinadze problem
   */

  /**
   * Chichinadze test function.
   *
   * Reference: V. K. Chichinadze, "The $\psy$-transform for Solving Linear and Nonlinear Programming Problems", Automatica 5(3), pp. 347--355, 1969.
   *
   * @param x is the vector (of dimension 2) where the function should be evaluated.
   */
  double chichinadze(double *x);

  /**
   * Problems taken from the GNU scientific library.
   */

  /**
   * Function E1 in simulated annealing GSL test.
   *
   * @param x is the point where the function is evaluated.
   *
   * @return the value of E1 function at x.
   */
  double gsl_e1(double x);

  double gsl_e1_stoch(Random *ran, double x);

#ifdef __cplusplus
}
#endif

#endif  /* __OP_COLLECTION_UNCONSTRAINED_H__ */
