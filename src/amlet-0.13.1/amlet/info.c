/**
 * info.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <stdio.h>

#include <config.h>

#include "info.h"
#include <amlet/amlet.h>
#include "io.h"

#include <oratio/ctools.h>

/**
 * @param info is the numerical identifier of the message to print.
 */
char *amlet_info_message(int info)
{
  return amlet_error_msgs[info].name;
}

/**
 * Information routines
 */

/**
 * Print a message describing the basic usage of the program.
 */
void amlet_help(char *s)
{
  if (s) printf("%s\n", s);
  printf("Usage %s [OPTIONS] FILE\n", COMMAND);
  printf("\n");
  printf("Supported options:\n");
  printf(" -o FILE, --output FILE: save log information into a file\n");
  printf(" -h, --help: display this help message\n");
  printf(" -V, --version: display the program version\n");
  printf(" -v, --verbose: verbose mode\n");
  printf("\nPlease signal any bugs to <%s>\n", EMAIL);
}
