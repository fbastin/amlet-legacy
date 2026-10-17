/**
 * acceptance_rejection.c
 *
 * ORATIO Library
 * Copyright (C) 1999-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "test_ran.h"

#include <oratio/ran.h>

double f(double x, double *params)
{
  return 30*x*(x*(1-x*(2-x)));
}

double t(double x, double *params)
{
  return 15.0/8.0;
}

void test_acceptance_rejection()
{
  int n = 1000000;
  int i;
  double m = 0.0;

  Random *ran = ran_random();

  printf("*****************************************\n");
  printf("* Acceptance-rejection test\n");
  printf("*****************************************\n");

  for (i = 0; i < n; i++)
    m += (ran_acceptance_rejection(ran, ran, NULL, NULL, f, t)-m)/(i+1);

  printf("Mean: %g\n", m);

  ran_random_free(ran);
}
