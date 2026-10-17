/**
 * Rosenbrock function.
 */

double rosenbrock(double *x, int n, void *data)
{
  double t1, t2;

  t1 = (x[1]-x[0]*x[0]);
  t2 = (1-x[0]);

  return 100*t1*t1+t2*t2;
}

void rosenbrock_deriv1(double *x, int n, void *data, double *grad)
{
  double t = (x[1]-x[0]*x[0]);

  grad[0] = -400*x[0]*t-2*(1-x[0]);;
  grad[1] = 200*t;
}
