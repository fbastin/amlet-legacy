/**
 * file.c
 *
 * ORATIO Library - CTools
 * Copyright (C) 2005-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License version.
 */

#include "test_ctools.h"

/**
 * CTools file functions tests.
 * 
 * @return 0 if no error has been encountered.  
 */
int ctools_file_tests()
{
  char *buff;
  long length;

  printf("file.c exist: %d\n", c_file_check_exist("file.c"));

  buff = c_file_get_text("file.c", &length);

  free(buff);

  return 0;
}
