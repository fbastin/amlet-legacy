#include <oratio/constants.h>

/**
 * Chichinadze test function.
 *
 * Reference: V. K. Chichinadze, "The $\psy$-transform for Solving Linear and Nonlinear Programming Problems", Automatica 5(3), pp. 347--355, 1969.
 *
 * @param x is the vector (of dimension 2) where the function should be evaluated.
 */
inline double chichinadze(double *x)
{
  return x[0]*(x[0]-12)+11+10*cos(M_PI_2*x[0])+8*sin(M_PI*x[0])-M_SQRT1_2/M_SQRTPI*exp(-(x[1]-0.5)*(x[1]-0.5)/2);
}
