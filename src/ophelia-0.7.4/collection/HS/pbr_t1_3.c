/**
 *
 * Willi Hock and Klaus Schittkowski, "Test Examples for Nonlinear
 * Programming Codes", Springer-Verlag, 1981.
 *
 */

/**
 * PBR_T1_3
 */

double pbr_t1_3(double *x, int n, void *data)
{
  double t;

  t = x[0]+1;

  t = t*t*t/3.0+x[1];

  return t;
}

void pbr_t1_3_set_bounds(double *b)
{
  b[0] = 1;
  b[2] = 0;
}

/**
 * Set the lower and upper bounds.
 *
 * @param xl is a pointer to the lower bound.
 * @param xu is a pointer to the upper bound.
 */
void pbr_t1_3_set_bounds_st(double *xl, double *xu)
{
  xl[0] = 1;
  xl[1] = 0;
  xu[0] = xu[1] = 100;
}
