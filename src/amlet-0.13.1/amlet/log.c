/**
 * log.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2004 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the MIT License
 */

#include <log.h>

#include <stdarg.h>

#include <io.h>

void amlet_log(char *msg, ...)
{
  va_list curseur;

  va_start(curseur, msg);

  amlet_print_msg(mclog->f, msg, curseur);

  va_end(curseur);
}

