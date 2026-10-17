/**
 * error.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <stdio.h>
#include <string.h>

#include <oratio/ctools.h>

#include "amlet.h"

/**
 * Error and warnings functions
 */

/**
 * @param f is the file where the warning will be issued.
 * @param err is the error code.
 */
int amlet_warning_va(FILE *f, int err, char *type,
		     const char *msg, va_list ap)
{
  char *buff;
  int n = 1024;

  /* Create the message */
  buff = _mallocs(n, char);
  while (vsnprintf(buff, n, msg, ap) == -1) {
    n <<= 1;
    buff = _reallocs(buff, n, char);
  }

  if (f != stderr && err) fprintf(stderr, "* AMLET %s *: %s\n", type, buff);
  fprintf(f, "* AMLET %s *: %s\n", type, buff);

  n = strlen(buff);

  free(buff);

  return n;
}

/**
 * @param f is the file where the warning will be issued.
 */
void amlet_warning(FILE *f, int verbose, char *msg, ...)
{
  va_list curseur;

  va_start(curseur, msg);

  amlet_warning_va(f, verbose, "WARNING", msg, curseur);

  va_end(curseur);
}
