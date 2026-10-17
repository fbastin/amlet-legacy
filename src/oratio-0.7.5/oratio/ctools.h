/**
 * ctools.h
 *
 * ORATIO Library
 * Copyright (C) 2000-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_CTOOLS_H__
#define __ORATIO_CTOOLS_H__

#include <errno.h>

#include "values.h"

#include <sys/types.h>

#include <stdio.h>
#include <stdlib.h>

#if TIME_WITH_SYS_TIME
# include <time.h>
# include <sys/time.h>
# include <limits.h>
#else
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
#  include <limits.h>
# else
#  include <time.h>
# endif
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /* Generic function */
  typedef void (*C_GENERIC)(void);

  typedef struct _C_CSV C_CSV;

  struct _C_CSV {
    char *blanks;   /* blank characters */
    char separator; /* items separator */
    char *comments; /* character indicating a comment line */
  };

  C_CSV *c_csv_new(char separator);

  void c_csv_free(C_CSV *csv);

  char **c_csv_read_line(C_CSV *csv, FILE *input, int *length);

  void c_csv_write_line(C_CSV *csv, FILE *output, char **items, int nb_items);

  /**
   * Errors management.
   */

  void c_warning(char *msg, ...);

  void c_error(char *msg, ...);

#include "ctime.h"
#include "ctimer.h"
#include "exception.h"
#include "cfile.h"
#include "cstring.h"
#include "ctools_macros.h"
#include "memory.h"
#include "symbol.h"
#include "csystem.h"
#include "values.h"

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_CTOOLS_H__ */
