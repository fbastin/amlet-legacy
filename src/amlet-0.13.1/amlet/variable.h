/**
 * variable.h
 *
 * OPHELIA Library
 * Copyright (C) 2006-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __DCL_VARIABLE_H__
#define __DCL_VARIABLE_H__

#include <oratio/oratio.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Recombine fixed and variable parameters.
   *
   * @param nvar is the number of variables.
   * @param type is a vector such that type[i] corresponds to the type
   *   of variable i.
   * @param x is the variable parameters vector.
   * @param np is the list of non-parametric variables.
   * @param fixed is the fixed variables vector.
   * @param f is on output a vector such that f[i] = 0 if i is a fixed
   *   variable.
   * @param p is on output the recombined parameters vectors. p must
   *   have a sufficient size on input.
   */
  void dcl_combine_parameters(int nvar, int *type, double *x,
			      DS_List *np, double *fixed, int *f, double *p);

  /**
   * Permutate the variable and type vectors in order to make the
   * optimization task easier.
   *
   * @param n is the dimension of x.
   * @param x is the variables vector.
   * @param px is a vector containing on output the permutation applied to x.
   * @param ipx is on output the vector containing the inverse permutation
   *   of px.
   * @param ptype is a vector containing on output the permutation
   *   applied to type.
   * @param offset is a vector containing on output the offsets between
   *   types. [offset[i], offset[i+1]-1] designs a group of similar
   *   variables (FIXED, CONSTANT, RANDOM, NONPARAMETRIC) in x. Each
   *   nonparametric variable has its own group, the other variables
   *   sharing their group with similar variables.
   * @param w is a working vector of dimension at least equal to nx.
   */
  void dcl_variable_permutate(int nx, double *x, int nt, int *type,
			      int *px, int *ipx, int *ptype,
			      int *offsets, double *w);

  /**
   * Restore the variable in x in the original order.
   *
   * @param n is the dimension of x.
   * @param x is the variables vector.
   * @param pinv it the permutation vector that has to be applied to x.
   */
  void dcl_variable_restore(int n, double *x, int *pinv);

#ifdef __cplusplus
}
#endif

#endif  /* __DCL_VARIABLE_H__ */
