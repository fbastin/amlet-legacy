/**
 * ntools_test.c
 *
 * ORATIO Library
 * Copyright (C) 1999-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "ntools_test.h"

#include <stdio.h>
#include <string.h>

int utilities()
{
  int code;
  clock_t start, time;
  double cpu_time;

  double x = 10.0;
  double y = 3;
  int i;

  printf("* Machine characteristics *\n");
  code = nt_machar_test();

  printf("* Natural numbers *\n");
  nt_natural_test();

  printf("* Power functions *\n");

  start = clock();

  for (i = 0; i < 50000000; i++) cpu_time = pow(x,y);

  time = clock();

  cpu_time = ((double)(time-start))/CLOCKS_PER_SEC;

  printf("%g^%g = %g; time spent: %g seconds\n",x, y, pow(x,y), cpu_time);

  start = time;
  for (i = 0; i < 50000000; i++) cpu_time = nt_ipow(x,y);

  time = clock();

  cpu_time = ((double)(time-start))/CLOCKS_PER_SEC;

  printf("%g^%g = %g; time spent: %g seconds\n", x, y, nt_ipow(x,y), cpu_time);

  return code ;
}

int main()
{
  int i;

  NTLog *log;
  int code[NT_TEST_END];

  log = nt_log_new(NULL);

  *code = 0;

  nt_log_section(log, "Algebra routines");
  code[NT_TEST_ALGEBRA] = nt_algebra_test();

  nt_log_section(log, "Modular arithmetic");
  code[NT_TEST_MODULAR] = nt_modular_check();

  nt_log_section(log, "Polynomial routines");
  code[NT_TEST_POLYNOMIAL] = nt_polynomial_test();

  nt_log_section(log, "Prime numbers routines");
  code[NT_TEST_PRIMES] = primes_test();

  nt_log_section(log, "Utility routines");
  code[NT_TEST_UTILITIES] = utilities();

  nt_log_section(log, "B-splines routines");
  code[NT_TEST_BSPLINES] = nt_bsplines_check();

  for (i = NT_TEST_START+1; i < NT_TEST_END; i++) {
    if (code[i]) {
      *code = 1;
      nt_log_section(log, "%d test failure.", i);
    }
  }
  if (!(*code))
    nt_log_section(log, "All tests have been successfully completed.");

  nt_log_free(log);

  return *code;
}
