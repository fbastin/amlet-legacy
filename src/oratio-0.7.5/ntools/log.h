/**
 * log.h
 *
 * NTools Library
 * Copyright (C) 2002 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __NTOOLS_LOG_H__
#define __NTOOLS_LOG_H__

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  typedef struct {
    FILE *f;
    char *name;
  } NTLog;

  /**
   * Create a new NTLog file.
   *
   * @param name must be NULL or assigned to a valid file name.
   *
   * @return a newly allocated NTLog structure.
   */
  NTLog *nt_log_new(char *name);

  /**
   * Set the output file the log objet 'log'.
   *
   * @param name must be NULL or assigned to a valid file name.
   */
  void nt_log_set_file(NTLog *log, char *name);

  void nt_log_free(NTLog *ntlog);

  /**
   * Output functions
   */

  void nt_log_write(NTLog *ntlog, const char *msg, ...);

  void nt_log_section(NTLog *ntlog, const char *msg, ...);

  void nt_log_subsection(NTLog *ntlog, const char *msg, ...);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __NTOOLS_LOG_H__ */
