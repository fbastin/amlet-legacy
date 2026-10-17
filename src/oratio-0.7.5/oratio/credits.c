/**
 * credits.c
 *
 * ORATIO Library
 * Copyright (C) 2007-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "config.h"

#include <stdio.h>

/**
 * Print library credits on standard output.
 */
void oratio_credits()
{
  printf("ORATIO Library v.%s\n", PACKAGE_VERSION);
  printf("Please report any bug to %s\n", PACKAGE_BUGREPORT);
}
