/**
 * minpack2_wrapper.c
 *
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * This file is licensed under BSD License.
 */

/**
 * Wrapper for MINPACK-2 test problem collection.
 */

#include "minpack2_collection.h"

#include "ophelia/nlp_least_squares.h"

#include <limits.h>
#include <oratio/ctools.h>
#include <oratio/stat.h>

#include <string.h>

/* Declaration of MINPACK-2 functions in C-Style. */

extern void wdaerfj_(int *m, int *n, double *x, double *fvec, double *fjac,
		     int *ldfjac, int *task);
extern void dmsabc_(int *nx, int *ny,
		    double *bottom, double *top,
		    double *left, double *right);
extern void dmsafg_(int *nx, int *ny, double *x, double *f, double *fgrad,
		    char *task,
		    double *bottom, double *top,
		    double *left, double *right);

/**
 * Set a problem from the Minpack 2 collection.
 *
 * @param n is the number of variables.
 * @param nx is the number of grid points in the first coordinate direction.
 * @param ny is the number of grid points in the second coordinate
 *   direction. If the problem is formulated in one spatial direction,
 *   ny=1.
 *
 * @return a newly Minpack2 allocated Minpack2Problem structure with
 *   relevant information.
 */
Minpack2Problem *minpack2_set_problem(const char *name, int n, int nx, int ny)
{
  Minpack2Problem *mp = NULL;
  char prob[5];
  int nb, idx;

  strncpy(prob, name, 4);
  prob[4] = '\0';

  nb = sizeof(minpack2_collection)/sizeof(minpack2_collection[0]);

  idx = c_symbol_get_token(prob, (const c_symbol *)minpack2_collection, nb);

  if (idx != INT_MIN) {
    mp = _malloc(Minpack2Problem);

    mp->name = strdup(name);
    mp->id = idx;
    mp->n = n;
    mp->nx = nx;
    mp->ny = ny;
    mp->w = _mallocs(2*n, double);
  }

  return mp;
}

void minpack2_problem_free(Minpack2Problem *mp)
{
  if (mp) {
    free(mp->name);
    free(mp->w);
    free(mp);
  }
}

/**
 * @param itask
 *   itask = 1 corresponds to F
 *   itask = 2 corresponds to G or J
 *   itask = 3 corresponds to FG or FJ
 *   itask = 4 corresponds to XS
 *   itask = 5 corresponds to XL
 *   itask = 6 corresponds to XU
 * @param_type
 *   param = 1: usual problem
 *   param = 2: least-squares problem
 */
void task_int2str(int itask, char *stask, int prob_type)
{
  switch (itask) {
  case 1:
    strcpy(stask, "F");
    break;
  case 2:
    switch (prob_type) {
    case 2:
      strcpy(stask, "J");
    default:
      strcpy(stask, "G");
    }
    break;
  case 3:
    switch (prob_type) {
    case 2:
      strcpy(stask, "FJ");
    default:
      strcpy(stask, "FG");
    }
    break;
  case 4:
    strcpy(stask, "XS");
    break;
  case 5:
    strcpy(stask, "XL");
    break;
  case 6:
    strcpy(stask, "XU");
    break;
  }
}

/**
... follow dminfg
 */
double minpack2_wrapper(int n, int nx, int ny, double *x,
			double *f, double *g, int itask,
			int prob, double par, double *w)
{
  char stask[3];
  double v = 0.0;

  switch(prob) {
  case DAER:
    
    wdaerfj_(&ny, &nx, x, f, g+nx, &ny, &itask);
    /* compute f as the mean of f[1] to f[nx] */
    if (itask == 1 || itask == 3) {
      //      nt_vect_print(stdout, "f", f, ny);
      v = st_mean(f, 1, 0, ny);
    }
    /* compute the mean gradient */
    if (itask == 2 || itask == 3) {
      nlp_ls_grad(ny, nx, g, ny, f, 1, 1);
      //      nt_vect_print(stdout, "g", g, nx*(ny+1));
    }
    if (itask == 1 || itask == 3) return v;
    break;
  case DEPT:
    break;
  case DGL1:
    break;
  case DGL2:
    break;
  case DLJ2:
    break;
  case DLJ3:
    break;
  case DMSA:
    task_int2str(itask, stask, 1);
    dmsabc_(&nx, &ny, w, &(w[nx+2]), &(w[2*nx+4]), &(w[2*nx+ny+6]));
    dmsafg_(&nx, &ny, x, f, g, stask, w,
	    &(w[nx+2]), &(w[2*nx+4]), &(w[2*nx+ny+6]));
    break;
  case DMSO:
    break;
  case DODC:
    break;
  case DPJB:
    break;
  case DSSC:
    break;
  default:
    c_error("unsupported problem");
  }

  return 0.0;
}

int minpack2_starting_point(Minpack2Problem *mp, double *x)
{
  int nx, ny;
  double f;
  double work[100000];

  int task;

  nx = mp->nx;
  ny = mp->ny;

  switch(mp->id) {
  case DAER:
    task = 4;
    wdaerfj_(&ny, &nx, x, work, work, &ny, &task);
    break;
  case DMSA:
    dmsabc_(&nx, &ny,
	    mp->w, &(mp->w[nx+2]), &(mp->w[2*nx+4]), &(mp->w[2*nx+ny+6]));
    dmsafg_(&nx, &ny, x, &f, NULL, "XS",
	    mp->w, &(mp->w[nx+2]), &(mp->w[2*nx+4]), &(mp->w[2*nx+ny+6]));
    break;
  default:
    c_error("unsupported problem");
  }

  return f;
}

inline double minpack2_wrapper_f(double *x, int n, void *problem)
{
  double f = 0.0;
  Minpack2Problem *mp = (Minpack2Problem *)problem;

  minpack2_wrapper(n, mp->nx, mp->ny, x, &f, NULL, 1, mp->id, mp->par,
		   mp->w);

  return f;
}

inline double minpack2_wrapper_f_values(double *x, int n, void *problem,
					double *values, int m)
{
  double f = 0.0;
  Minpack2Problem *mp = (Minpack2Problem *)problem;

  /* comment obtenir values??? */
  f = minpack2_wrapper(n, n, m, x, values, NULL, 1, mp->id,
		       mp->par, mp->w);

  return f;
}

inline double minpack2_wrapper_fg_values(double *x, int n, void *problem,
					 double *values, int m, double *grad)
{
  double f = 0.0;
  Minpack2Problem *mp = (Minpack2Problem *)problem;

  /* comment obtenir values??? */
  f = minpack2_wrapper(n, n, m, x, values, grad, 3,
		       mp->id, mp->par, mp->w);

  return f;
}

inline void minpack2_wrapper_g(double *x, int n, void *problem,
			       double *grad)
{
  double f = 0.0;

  Minpack2Problem *mp = (Minpack2Problem *)problem;

  minpack2_wrapper(n, mp->nx, mp->ny, x, &f, grad, 2, mp->id, mp->par,
		   mp->w);
}
