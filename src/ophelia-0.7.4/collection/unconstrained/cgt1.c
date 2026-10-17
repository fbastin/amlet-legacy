/**
 * cgt1.c
 *
 * OPHELIA Library
 * Copyright (C) 2004-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
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
 */
double cgt_1(double *x, int n, void *data)
{
  double x_square = x[0]*x[0];
  double y_square = x[1]*x[1];

  return -10*x_square + 10*y_square + 4*sin(x[0]*x[1]) - 2*x[0]
    + x_square*x_square;
}

/**
 * Problem CGT-1 gradient.
 *
 * @param data is unused.
 */
void cgt_1_deriv1(double *x, int n, void *data, double *grad)
{
  double x_square = x[0]*x[0];

  grad[0] = -20*x[0]+4*x[1]*cos(x[0]*x[1])-2+4*x[0]*x_square;
  grad[1] = 20*x[1]+4*x[0]*cos(x[0]*x[1]);
}

void cgt_1_deriv2(double *x, int n, void *data, double *H)
{
  double p = x[0]*x[1];
  int i;

  for (i = 0; i < n*n; i++) H[i] = 0.0;

  H[0] = 4*(3*x[0]*x[0]-x[1]*x[1]*sin(x[0]*x[1])-5);
  H[1] = H[2] = 4*(cos(p)-p*sin(p));
  H[3] = 20-4*x[0]*x[0]*sin(x[0]*x[1]);
}
