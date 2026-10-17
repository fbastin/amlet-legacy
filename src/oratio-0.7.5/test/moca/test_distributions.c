/**
 * test_distributions.c
 *
 * ORATIO Library
 * Copyright (C) 2006 <fabian.bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "caffeine.h"

#include <oratio/oratio.h>

#include <stdio.h>

int moca_test_distributions(double **sample)
{
  FILE *f;
  double u, inc, val = 0.0;

  double knots[] = {0, 1.0/3, 2.0/3, 1.0};
  //double knots[] = {0, 0.25, 0.5, 0.75, 1.0};
  double *P;
  double a[200000], b[200000];
  int dim;
  int i;

  int error_code = 0;

  NT_BSpline *bs;

  //  double params[] = {0.0, 0.5};
  // double params[] = {0.0, 2.0};
  /* Estimation of a spline with a normal. */
  //  double params[] = {-0.09899688, 3.08547};
  /* Estimation of a spline with a normal - panel. */
  //  double params[] = {-0.045597712, 3.0789612};

  /* Estimation of a spline with a normal. */
  // double params[] = {2.378881, 1.1954589};
  /* Estimation of a spline with a normal - panel. */
  //  double params[] = {2.2624609, 1.0609617};

  /* Estimation of a positive spline by estimation with a lognormal. */
  // double params[] = {0.78285792, 0.84213603};
  /* Estimation of a postive spline by estimation with a lognormal - panel. */
  // double params[] = {0.76211241, 0.78678539};
  /* Estimation of a standard by estimation with a lognormal. */
  // double params[] = {-221.40941, 195.50534};
  /* Estimation of a standard by estimation with a lognormal - panel. */
  // double params[] = {-19.933735, 23.647358};

  /* Estimation of real data - normal. */
  // double params[] = {-0.070651772, 0.22180886};
  /* Estimation of real data - lognormal. */
  double params[] = {-4.9094585, 4.5704215};

  P = _mallocs(6, double);

  /* Standard spline. */
  P[0] = -15.0;
  P[1] = -3;
  P[2] = -5e-1;
  P[3] = 5e-1;
  P[4] = 3;
  P[5] = 15.0;

  /* Positive spline. */
  /*  P[0] = 0.0;
  P[1] = 1.0;
  P[2] = 1.5;
  P[3] = 3.0;
  P[4] = 3.0;
  P[5] = 15.0; */

  /* Estimation of a Lognormal L(0.0, 5.0). */
  /*P[0] = 0.56648517;
  P[1] = 0.69978206;
  P[2] = 0.69978206;
  P[3] = 0.69978206;
  P[4] = 2.3439802;
  P[5] = 2.3439802; */

  /* Estimation of a Normal N(0,4). */
  /*  P[0] = -3.6771307;
  P[1] = -3.6771307;
  P[2] = -0.074608489;
  P[3] = 0.36500275;
  P[4] = 2.4455607;
  P[5] = 2.4455607; */

  /* Estimation of a Normal N(0,4) - panel. */
  /*  P[0] = -3.0638226;
  P[1] = -2.8132104;
  P[2] = -0.76157993;
  P[3] = 0.88705575;
  P[4] = 1.7126744;
  P[5] = 3.9004939; */

  /* Estimation of a Lognormal of parameters (0.0, 0.5). */
  /*  P[0] = 0.56648517;
  P[1] = 0.69978206;
  P[2] = 0.69978206;
  P[3] = 0.69978206;
  P[4] = 2.3429802;
  P[5] = 2.3429802; */

  /* Estimation of a Lognormal of parameters (0.0, 0.5) - panel. */
  /*  P[0] = 0.61091044;
  P[1] = 0.61091044;
  P[2] = 0.61091044;
  P[3] = 1.2236644;
  P[4] = 1.990672;
  P[5] = 1.990672; */

  /* Estimation of a spline with another spline. */
  P[0] = -15.924551;
  P[1] = -2.4689725;
  P[2] = -1.4908432;
  P[3] =  1.5407503;
  P[4] =  1.5407503;
  P[5] =  12.024644;

  /* Estimation of a spline with another spline - panel. */
  /*  P[0] = -12.04302;
  P[1] = -3.5660613;
  P[2] = -0.58255938;
  P[3] =  0.73265751;
  P[4] =  2.2508159;
  P[5] = 17.890703; */

  /* Estimation of real data with a spline. */
  P[0] = -2.8900049;
  P[1] = -0.16892952;
  P[2] = -0.11920471;
  P[3] = 0.030428079;
  P[4] = 0.030428079;
  P[5] = 0.7530208;

  dim = 0;
  u = 1.0e-6;
  inc = 0.00001;
  while (u < 1.0) {
    val = moca_affect_value(LOGNORMAL, u, params, NULL);
    a[dim] = u;
    b[dim] = val;
    dim++;
    u += inc;
  }

  f = fopen("lognormal_pdf.dat", "w");

  for (i = 0; i < dim-1; i++) {
    fprintf(f, "%g %g\n", 0.5*(b[i]+b[i+1]), (a[i+1]-a[i])/(b[i+1]-b[i]));
  }

  fclose(f);

  f = fopen("lognormal_cdf.dat", "w");

  for (i = 0; i < dim-1; i++) {
    fprintf(f, "%g %g\n", b[i], a[i]);
  }

  fclose(f);

  f = fopen("lognormal_icdf.dat", "w");

  for (i = 0; i < dim-1; i++) {
    fprintf(f, "%g %g\n", a[i], b[i]);
  }

  fclose(f);

  dim = 0;
  u = 1.0e-6;
  inc = 0.00001;
  while (u < 1.0) {
    val = moca_affect_value(NORMAL, u, params, NULL);
    a[dim] = u;
    b[dim] = val;
    dim++;
    u += inc;
  }

  f = fopen("normal_pdf.dat", "w");

  for (i = 0; i < dim-1; i++) {
    fprintf(f, "%g %g\n", b[i], st_gaussian_pdf(b[i]));
  }

  fclose(f);

  f = fopen("normal_cdf.dat", "w");

  for (i = 0; i < dim-1; i++) {
    fprintf(f, "%g %g\n", b[i], a[i]);
  }
  

  fclose(f);

  f = fopen("normal_icdf.dat", "w");

  for (i = 0; i < dim-1; i++) {
    fprintf(f, "%g %g\n", a[i], b[i]);
  }

  fclose(f);

  bs = nt_bsplines_new(3, 4, knots, 1);

  dim = 0;
  u = 1.0e-6;
  inc = 0.00001;
  while (u < 1.0) {
    /* val = ntools_bsplines_curve_point(bs->n, bs->p, u, bs->U, P, bs->N); */
    val = moca_affect_value(RAN_SP_SPLINE, u, P, bs);
    a[dim] = u;
    b[dim] = val;
    dim++;
    u += inc;
  }

  f = fopen("splines_cdf.dat", "w");

  for (i = 0; i < dim-1; i++) {
    fprintf(f, "%g %g\n", b[i], a[i]);
  }

  fclose(f);

  f = fopen("splines_icdf.dat", "w");

  for (i = 0; i < dim-1; i++) {
    fprintf(f, "%g %g\n", a[i], b[i]);
  }

  fclose(f);

  f = fopen("splines_pdf.dat", "w");

  u = -20.0;
  while (u < 0.5*(b[0]+b[1])) {
    fprintf(f, "%g %g\n", u, 0.0);
    u += 0.01;
  }
  for (i = 0; i < dim-1; i++) {
    fprintf(f, "%g %g\n", 0.5*(b[i]+b[i+1]), (a[i+1]-a[i])/(b[i+1]-b[i]));
  }
  u = 0.5*(b[dim-2]+b[dim-1])+0.01;
  while (u <= 20.0) {
    fprintf(f, "%g %g\n", u, 0.0);
    u += 0.1;
  }

  fclose(f);

  nt_bsplines_free(bs);
  free(P);

  return error_code;
}
