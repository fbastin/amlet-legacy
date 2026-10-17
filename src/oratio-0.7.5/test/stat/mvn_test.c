/**
 * mvn_test.c
 *
 * ORATIO Library
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <oratio/stat.h>

int st_mvn_test()
{
  double x[] = {0.0, 0.0};

  printf("MVN test: %g\n", st_mvn_standard_pdf(2, x, 1));

  return 0;
}
