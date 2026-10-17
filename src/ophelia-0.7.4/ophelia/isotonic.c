/**
 * isotonic.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2005-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <oratio/ctools.h>

#include <string.h>

/**
 * Isotonic regression.
 *
 * Ref: Michael J. Best and Nilotpal Chakravarti, "Active set
 *   algorithms for isotonic regression; a unifying framework",
 *   Mathematical Programming 47, pp 425-439, 1990.
 *
 * @param n is the problem dimension
 * @param y is the observations vector
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
				      double *workspace)
{
  double *av, *weights, *avJ;
  double wbm, avbm, wpbm;
  int *J;

  int m = 0;

  int i, j, stop;

  J = _mallocs(n+2, int);
  avJ = workspace;
  av = avJ+n+1;
  weights = av+n+1;

  /* Initialization */
  J[0] = 0;
  av[n] = 0.0;
  av[n-1] = y[n-1];
  weights[n-1] = w ? w[n-1] : 1.0;

  for (i = n-2; i >=0; i--) {
    weights[i] = weights[i+1]+ (w ? w[i] : 1.0);
    av[i] = ((w ? w[i] : 1.0)*y[i]+weights[i+1]*av[i+1])/weights[i];
  }
  avJ[0] = av[0];
  j = 0;

  while (j != n-1) {
    j = J[m];
    while(av[j+1] <= avJ[m] && j < n-1) j++;

    if (j != n-1) {
      wbm = weights[J[m]]-weights[j+1];
      avbm = (weights[J[m]]*avJ[m]-weights[j+1]*av[j+1])/wbm;
      avJ[m] = avbm;

      stop = 0;
      while (m && !stop) {
	if (avJ[m] > avJ[m-1])
	  stop = 1;
	else {
	  /* Fusion of two blocks */
	  wbm = weights[J[m]] - weights[j+1];
	  wpbm = weights[J[m-1]] - weights[J[m]];
	  m--;
	  avJ[m] = (avJ[m+1]*wbm+avJ[m]*wpbm)/(wpbm+wbm);
	}
      }
    }

    /* Add the block U_j */
    J[++m] = j+1;
    avJ[m] = av[j+1];
  }

  J[m+1] = n+1;
  /* Set up the solution. */
  for (i = 0; i < m; i++) {
    for (j = J[i]; j < J[i+1]; j++) x[j] = avJ[i];
  }

  free(J);

  return m;
}


/**
 * Isotonic regression.
 *
 * @param n is the problem dimension.
 * @param y is the observations vector.
 * @param w is the weighting vector; is w is set to NULL, the weights
 *   are assumed to be equal to 1.0.
 * @param x is unknown on input, and the solution on output; x may be
 *   the same pointer as y.
 *
 * @return the number of active sets.
 */
int nlp_isotonic_regression(int n, double *y, double *w, double *x)
{
  double *workspace;
  int as;

  workspace = _mallocs(3*n+2, double);

  as = nlp_isotonic_regression_workspace(n, y, w, x, workspace);

  free(workspace);

  return as;
}

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
				void *data)
{
  return nlp_isotonic_regression(n, y, w, x);
}

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
				 int nidx, int *idx)
{
  double *yp, *wp, *xp;

  int as = 0; /* total number of active sets. */
  int i;

  /* This copy is necessary in order to take account of invariant subsets. */
  if (x != y) ;
  memcpy(x, y, n*sizeof(double));

  for (i = 0; i < nidx; i++) {
    n = idx[2*i+1]-idx[2*i];

    yp = y+idx[2*i];
    wp = (w) ? w+idx[2*i] : 0;
    xp = x+idx[2*i];

    as += nlp_isotonic_regression(n, yp, wp, xp);
  }

  return as;
}
