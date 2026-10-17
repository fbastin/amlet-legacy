/**
 * test_ctools.c
 *
 * ORATIO Library - CTools
 * Copyright (C) 2005-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "test_ctools.h"

#include <oratio/ctools.h>

int main()
{
  char *stamp = NULL;

  stamp = ctools_timestamp(NULL);
  printf("%s\n", stamp);

  ctools_file_tests();

  return 0;
}
