/**
 * activesets.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2006-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ATHOS_ACTIVE_SETS_H__
#define __ATHOS_ACTIVE_SETS_H__

/**
 * Active sets routines.
 */

#include "constraints.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Active sets construction prototype.
   *
   * @param n is the problem dimension.
   * @param g is the objective gradient.
   * @param incg the the increment in vector g.
   */
  typedef int (* NlpAsConstruct)(void *as,
				 int n, double *x, int incx,
				 double *g, int incg,
				 int *a, int inca);

  /**
   * Free variables initialization. By default, all variables are free.
   *
   * @param n is the problem dimension.
   */
  void nlp_as_free_init_dflt(int n, int *free, int incfree);

  /**
   * Determine the active set for bounds constraints.
   *
   * @param n is the problem dimension.
   * @param x is the point where the active set is determined.
   * @param incx is the increment of x.
   */
  int nlp_as_bounds(void *as, int n, double *x, int incx,
		    double *g, int incg, int *free, int incfree);

  /**
   * Determine the active set for bounds constraints and step validation.
   *
   * @param n is the problem dimension.
   * @param x is the point where the active set is determined.
   * @param incx is the increment of x.
   */
  int nlp_as_bounds_step(void *as, int n, double *x, int incx,
			 double *step, int incs, double *p, int incp,
			 int *free, int incfree, int *nfree);

  /**
   * Determine the active set under monotonicity constraints.
   *
   * @param n is the problem dimension.
   * @param x is the point where the active set is determined.
   * @param incx is the increment of x.
   * @param p is the search direction.
   * @param sets[i] designs the start of the i-th active set.
   */
  int nlp_as_monotonicity(NLP_Constraint *c, int n, double *x, int incx,
			  double *p, int incp, int *sets, int incsets);

  /**
   * Determine the active set for monotonicity constraints and step validation.
   *
   * @param n is the problem dimension.
   * @param x is the point where the active set is determined.
   * @param incx is the increment of x.
   */
  int nlp_as_monotonicity_step(NLP_Constraint *c, int n, double *x, int incx,
			       double *step, int incs, double *p, int incp,
			       int *nsets, int *sets, int incsets);

  /**
   * Compute the restricted gradient under bounds constraints.
   *
   * @param n is the problem dimension.
   */
  void nlp_as_rg_bounds(int n, double *p, int incp,
			double *g, int incg, int *xfree, int incf, int nfree);

  /**
   * Compute the restricted gradient under monotonicity constraints.
   *
   * We assume that the gradient is additive when different variables
   * are treated as one variable.
   *
   * @param n is the problem dimension.
   * @param incg is the increment in g.
   */
  void nlp_as_rg_monotonicity(int n, double *p, int incp,
			      double *g, int incg, int *sets, int incs,
			      int nsets);

  double nlp_as_rg_nrm2(int n, double *g, int incg, int nsets, int *sets, int incs);

  double nlp_as_rdot(int n, double *x, int incx, double *y, int incy,
		     int nsets, int *sets, int incs);

  void nlp_as_monotonicity_print(int nset, int *sets);

#ifdef __cplusplus
}
#endif

#endif  /* __ATHOS_ACTIVE_SETS_H__ */
