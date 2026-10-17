/**
 * leastsquares.h
 *
 * OPHELIA Library - least squares problems
 * Copyright (C) 2006-2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the MIT License
 */

#ifndef __OP_COLLECTION_LEASTSQUARES_H__
#define __OP_COLLECTION_LEASTSQUARES_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Basic linear least-squares problem.
   */
  double linearls1(double *x, int n, void *data);

  /**
   * Basic linear least-squares problem with gradient computation.
   */
  double linearls1_fg(double *x, int n, void *data, double *grad);

  /**
   * This example is taken from NL2SOL.
   */

  /**
   * @param x is the point where the function has to be evaluated.
   */
  double madsen(double *x, int n, void *data);

  /**
   * @param x is the point where the function has to be evaluated.
   */
  double madsen_fg(double *x, int n, void *data, double *grad);

  /**
   * @param x is the point where the function has to be evaluated.
   */
  double madsen_fg_values(double *x, int n, void *data, double *values, int m,
			  double *grad);

  /**
   * Basic nonlinear least-squares problem.
   */
  double nonlinearls1(double *x, int n, void *data);

  /**
   * Basic nonlinear least-squares problem with gradient computation.
   */
  double nonlinearls1_fg(double *x, int n, void *data, double *grad);

#ifdef __cplusplus
}
#endif

#endif  /* __COLLECTION_LEASTSQUARES_H__ */
