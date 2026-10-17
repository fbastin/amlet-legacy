/**
 * error.c
 *
 * ORATIO Library
 * Copyright (C) 2000-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/memory.h"

#include <stdarg.h>
#include <stdio.h>
#include <unistd.h>

/**
 * Print a message on the standard error terminal.
 */
void c_print_msg(const char *type, const char *msg, va_list ap)
{
  char *buff;
  int n = 1024;

  /* Create the message */
  buff = _mallocs(n, char);
  while (vsnprintf(buff, n, msg, ap) == -1) {
    n *= 2;
    buff = _reallocs(buff, n, char);
  }

  fprintf(stderr, "** %s **: %s\n", type, buff);

  free(buff);
}

void c_error(char *msg, ...)
{
  va_list curseur;

  va_start(curseur, msg);

  c_print_msg("ERROR", msg, curseur);

  va_end(curseur);

  exit(1);
}

void c_warning(char *msg, ...)
{
  va_list curseur;

  va_start(curseur, msg);

  c_print_msg("WARNING", msg, curseur);

  va_end(curseur);
}
