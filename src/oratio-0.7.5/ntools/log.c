/**
 * log.c
 *
 * NTools Library
 * Copyright (C) 2002 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <ntools/log.h>

#include <stdarg.h>
#include <string.h>

#include <oratio/ctools.h>

/**
 * Create a new log object.
 *
 * @param name is the output file name, or NULL is stdout has to be used.
 *
 * @return a NT_Log object.
 */
NTLog *nt_log_new(char *name)
{
  NTLog *ntlog;
  FILE *f;

  ntlog = _malloc(NTLog);

  if (name) {
    f = fopen(name, "w");
    ntlog->f = f;
    if (!f) c_error("unable to open log file %s", name);
    ntlog->name = strdup(name);
  }
  else {
    ntlog->f = stdout;
    ntlog->name = NULL;
  }

  return ntlog;
}

/**
 * Set the output file the log objet 'log'.
 *
 * @param name must be NULL or assigned to a valid file name.
 */
void nt_log_set_file(NTLog *log, char *name)
{
  FILE *f;

  if (log->name) {
    free(log->name);
    fclose(log->f);
  }

  if (name) {
    f = fopen(name, "w");
    if (!f) c_error("unable to open log file %s", name);
    log->f = f;
    log->name = strdup(name);
  }
  else {
    log->f = stdout;
    log->name = NULL;
  }
}

void nt_log_free(NTLog *ntlog)
{
  if (ntlog) {
    if (ntlog->name) fclose(ntlog->f);
    free(ntlog->name);
    free(ntlog);
  }
}

/**
 * Input functions
 */

void nt_log_section(NTLog *ntlog, const char *msg, ...)
{
  va_list curseur;

  fprintf(ntlog->f,
	  "\n*************************************************************\n");
  va_start(curseur, msg);
  vfprintf(ntlog->f, msg, curseur);
  va_end(curseur);
  fprintf(ntlog->f,
	  "\n*************************************************************\n");
}

void nt_log_subsection(NTLog *ntlog, const char *msg, ...)
{
  va_list curseur;

  fprintf(ntlog->f,
	  "\n* ");
  va_start(curseur, msg);
  vfprintf(ntlog->f, msg, curseur);
  va_end(curseur);
  fprintf(ntlog->f,
	  " *\n");
}

void nt_log_write(NTLog *ntlog, const char *msg, ...)
{
  va_list curseur;

  va_start(curseur, msg);
  vfprintf(ntlog->f, msg, curseur);
  va_end(curseur);
}
