/**
 * algebra_test.c
 *
 * ORATIO Library
 * Copyright (C) 1999-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "ntools_test.h"

double fct(double *x, int n, void *data)
{
  /* x^2+y^2+2x+xy */
  return x[0]*x[0]+x[1]*x[1]+2*x[0]+x[0]*x[1];
}

void grad_eval(double *x, int n, void *data, double *grad)
{
  grad[0] = 2*x[0]+x[1]+2;
  grad[1] = 2*x[1]+x[0];
}

int derivatives()
{
  double x[2];
  double grad[2];
  double **H;
  double h[2];
  double tol = 0.0;
  double scaling[2];

  /* x^2+y^2+2x+xy */

  x[0] = x[1] = 1.0;
  nt_vect_ones(2, scaling, 1);

  H = nt_matrix_new(2, 2);

  nt_derive_hess_fd(grad_eval, x, grad, H, h, tol, scaling, 2,
		    NULL, NULL, NULL);

  nt_vect_print(stdout, "x", x, 2);
  nt_vect_print(stdout, "g", grad, 2);

  nt_matrix_print(stdout, "H", H, 2, 2);

  nt_derive_hess_cd(fct, x, H, h, tol, scaling, 2, NULL, NULL, NULL);

  nt_matrix_print(stdout, "H", H, 2, 2);

  nt_matrix_free(H);

  return 0;
}

int nt_algebra_test()
{
  derivatives();

  return 0;
}
