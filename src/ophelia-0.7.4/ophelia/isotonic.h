/**
 * isotonic.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2005 Fabian Bastin, CERFACS <bastin@cerfacs.fr>
 *
 * Licensed under the MIT License
 */

#ifndef __NLP2_ISOTONIC_H__
#define __NLP2_ISOTONIC_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Isotonic regression.
   *
   * @param n is the problem dimension.
   * @param y is the observations vector.
   * @param w is the weighting vector; is w is set to NULL, the weights
   *   are assumed to be equal to 1.0.
   * @param x is unknown on input, and the solution on output; x may be
   *   the same pointer as y.
   * @param workspace is a double working vector of size at least equal
   *   to 3*n+2.
   *
   * @return the number of active sets.
   */
  int nlp_isotonic_regression_workspace(int n, double *y, double *w, double *x,
					double *workspace);

  /**
   * Isotonic regression.
   *
   * @param n is the dimension of the vectors y, w and x.
   * @param y is the observations vector.
   * @param w is the weighting vector; is w is set to NULL, the weights
   *   are assumed to be equal to 1.0.
   * @param x is unknown on input, and the solution on output; x may be
   *   the same pointer as y.
   *
   * @return the number of active sets
   */
  int nlp_isotonic_regression(int n, double *y, double *w, double *x);

  /**
   * Wrapper for the isotonic regression codes. See header file projection.h
   *
   * @param n is the dimension of the vectors y, w and x.
   * @param y is the observations vector.
   * @param w is the weighting vector; is w is set to NULL, the weights
   *   are assumed to be equal to 1.0.
   * @param x is unknown on input, and the solution on output; x may be
   *   the same pointer as y.
   * @param data is a data vector, not used here, but present to match
   *   some functional definitions.
   *
   * @return the number of active sets produced by the regression.
   */
  int nlp_isotonic_wrapper(int n, double *y, double *w, double *x,
			   void *data);

  /**
   * Wrapper for the isotonic regression codes. See header file projection.h
   *
   * @param n is the dimension of the vectors y, w and x.
   * @param y is the observations vector.
   * @param w is the weighting vector; is w is set to NULL, the weights
   *   are assumed to be equal to 1.0.
   * @param x is unknown on input, and the solution on output; x may be
   *   the same pointer as y.
   * @param nidx is the size of the index vector idx, divided by 2,
   *   indicating the number of regressions to perform.
   * @param idx is an index vector defining the subsets that will be
   *   used in the regression. More precisely, each subset defined by
   *   the start index idx[2*i] (inclusive) and end index idx[2*i+1]
   *   (exclusive) will be considered for one (independent) regression.
   *
   * @return the total number of active sets. x is the vector containing
   *   the subsets for which a regression has been applied, at the
   *   corresponding positions compared to y, and a copy of original
   *   elements when no regression was applied.
   */
  int nlp_isotonics_wrapper(int n, double *y, double *w, double *x,
			    int nidx, int *idx);

#ifdef __cplusplus
}
#endif

#endif  /* __NLP2_ISOTONIC_H__ */
