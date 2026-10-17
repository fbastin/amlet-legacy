/**
 * variable.c
 *
 * OPHELIA Library
 * Copyright (C) 2006-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "variable.h"

#include <ophelia/nonparametric.h>

#include <oratio/oratio.h>

#include <string.h>
#include <unistd.h>

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
			    DS_List *np, double *fixed, int *f, double *p)
{
  int i, j = 0, n;
  int pos = 0, k = 0, l;
  void *obj;

  for (i = 0; i < nvar; i++) {
    if (amlet_is_nonparametric(type[i])) {
      obj = np->data;
      np = np->next;
    }
    else
      obj = NULL;
    if (type[i] == FIXED) {
      f[pos] = 1;
      p[pos++] = fixed[j++];
    }
    else {
      n = moca_count_parameters(type[i], obj);
      for (l = 0; l < n; l++) f[pos+l] = 0;
      memcpy(p+pos, x+k, n*sizeof(double));
      pos += n;
      k += n;
    }
  }
}

/**
 * Permutate the variable and type vectors in order to make the
 * optimization task easier.
 *
 * @param n is the dimension of x.
 * @param x is the variables vector.
 * @param px is a vector containing on output the permutation applied to x.
 * @param ipx is on output the vector containing the inverse
 *   permutation of px.
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
			    int *offsets, double *w)
{
  int i;
  int *s = _mallocs(nt+1, int);

  s[0] = 0;
  for (i = 0; i < nt; i++)
    s[i+1] = s[i]+moca_count_parameters(type[i], NULL);
  /* Generate the permutation for type. The sorting is performed on
     basis of the distributions definition in the library MOCA. */
  ds_sort_gheapsort_with_perm(type, nt, sizeof(int),
			      ds_search_int, ptype);
  /* Extend it to x. */
  for (i = 0; i < nx; i++) px[i] = -1;  
  for (i = 0; i < nt; i++) {
    px[s[i]] = s[ptype[i]];
    memcpy(w+s[i]*sizeof(double), &(x[px[s[i]]]), s[ptype[i]]+1-s[ptype[i]]); 
  }
  /* Regularise the vector px. */
  for (i = 1; i < nx; i++) {
    if (px[i] == -1) px[i] = px[i-1]+1;
  }
  /* Inverse the permutation px. */
  ds_perm_inv(nx, px, ipx);

  /* Set offsets */
  offsets[0] = 0;
  i = 0;
  while (type[i] == FIXED) i++;
  offsets[1] = s[i];
  while (type[i] == CONSTANT) i++;
  offsets[2] = s[i];
  while (type[i] < MOCA_END_DISTRIBUTIONS) i++;
  offsets[3] = s[i];
  offsets[4] = nx;

  free(s);
}

/**
 * Restore the variable in x in the original order.
 *
 * @param n is the dimension of x.
 * @param x is the variables vector.
 * @param pinv it the permutation vector that has to be applied to x.
 */
void dcl_variable_restore(int n, double *x, int *pinv)
{
  ds_perm_inp(n, x, sizeof(double), pinv);
}
