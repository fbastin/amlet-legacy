/**
 * credits.c
 *
 * AMLET Library
 * Copyright (C) 2002-2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the MIT License
 */

#include "config.h"  /* produced by the configure command */
#include "info.h"

#include <stdio.h>

/**
 * Print library version.
 */
void amlet_version()
{
  printf("%s %s\n%s %s\n", PACKAGE, VERSION, COPYRIGHT, AUTHORS);
  printf("\nPlease signal any bugs to <%s>\n", EMAIL);
}

/**
 * Print library credits.
 *
 * @param f is the output file (assumed to be open).
 * @param start is the string to use at the beginning of each line.
 *   This variable will often be set to "".
 * @param end is the string to use at the end of each line.
 *   This variable will often be set to "".
 */
void amlet_credits(FILE *f, char *start, char *end)
{
  int i;

  fprintf(f, "%s", start);
  fprintf(f, "* %s %s *", PROGRAM, VERSION);
  fprintf(f, "%s\n", end);
  fprintf(f, "%s", start);
  fprintf(f, "%s\n", end);
  fprintf(f, "%s Fabian Bastin (et al.)", COPYRIGHT);
  fprintf(f, "%s\n", end);
  fprintf(f, "%s", start);
  fprintf(f, "Licensed under the MIT License");
  fprintf(f, "%s\n", end);
  fprintf(f, "%s", start);
  fprintf(f, "Department of Mathematics, University of Namur (2002-2004)");
  fprintf(f, "%s\n", end);
  fprintf(f, "%s", start);
  fprintf(f, "Imperial College, London (2004)");
  fprintf(f, "%s\n", end);
  fprintf(f, "%s", start);
  fprintf(f, "Cerfacs, Toulouse (2005-2006)");
  fprintf(f, "%s\n", end);
  fprintf(f, "%s", start);
  fprintf(f, "Department of Computing Science and Operational Research (2007-)");
  fprintf(f, "%s\n", end);
  fprintf(f, "%s", start);
  fprintf(f, "Initial development was supported by the Belgian");
  fprintf(f, "%s\n", end);
  fprintf(f, "%s", start);
  fprintf(f, "National Fund for Scientific Research.");
  fprintf(f, "%s\n", end);
  fprintf(f, "%s", start);
  for (i = 0; i < 70; i++) fprintf(f, "-");
  fprintf(f, "%s\n\n", end);
}
