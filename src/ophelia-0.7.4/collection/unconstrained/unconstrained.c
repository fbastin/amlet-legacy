/**
 * unconstrained.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2004-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <math.h>

#include <oratio/ntools.h>

/**
 * Home-made problems.
 */

/**
 * @param x is the vector where the function is evaluated.
 * @param n is the problem dimension.
 *
 * @return the value of the function unc1.
 */
double unc1(double *x, int n, void *data)
{
  return x[0]*(x[0]-1)+2*x[1]*(x[1]-2)+x[2]*(x[2]+1)+x[3]*(x[3]-4)
    +cos(x[0]*x[1])+sin(x[1]*x[2])+cos(x[2]*x[3]);
}

/**
 * @param x is the vector where the function is evaluated.
 * @param n is the problem dimension.
 */
void unc1_grad(double *x, int n, void *data, double *grad)
{
  grad[0] = 2*x[0]-1-x[1]*sin(x[0]*x[1]);
  grad[1] = 4*x[1]-4-x[0]*sin(x[0]*x[1])+x[2]*cos(x[1]*x[2]);
  grad[2] = 2*x[2]+1+x[1]*cos(x[1]*x[2])-x[3]*sin(x[2]*x[3]);
  grad[3] = 2*x[3]-4-x[2]*sin(x[2]*x[3]);
}

/**
 * @param x is the vector where the function is evaluated.
 * @param n is the problem dimension.
 *
 * @return the value of the function unc1.
 */
double unc2(double *x, int n, void *data)
{
  return (x[0]-1)*(x[0]-1) +
    x[1]*(x[1]-1)+2*x[2]*(x[2]-2)+x[3]*(x[3]+1)+x[4]*(x[4]-4)
    + cos(x[1]*x[2])+sin(x[2]*x[3])+cos(x[3]*x[4])+x[5]*x[5]
    + x[7]*x[7]+x[6]*x[6]+cos(x[6]+1);
}

/**
 * @param x is the vector where the function is evaluated.
 * @param n is the problem dimension.
 */
void unc2_grad(double *x, int n, void *data, double *grad)
{
  grad[0] = 2*x[0]-2;
  grad[1] = 2*x[1]-1-x[2]*sin(x[1]*x[2]);
  grad[2] = 4*x[2]-4-x[1]*sin(x[1]*x[2])+x[3]*cos(x[2]*x[3]);
  grad[3] = 2*x[3]+1+x[2]*cos(x[2]*x[3])-x[4]*sin(x[3]*x[4]);
  grad[4] = 2*x[4]-4-x[3]*sin(x[3]*x[4]);
  grad[5] = 2*x[5];
  grad[6] = 2*x[6]-sin(x[6]+1);
  grad[7] = 2*x[7];
}
