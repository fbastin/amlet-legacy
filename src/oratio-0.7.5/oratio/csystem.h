/**
 * system.h
 *
 * ORATIO Library
 * Copyright (C) 2000-2002 Fabian Bastin <slashbin@namurlug.org>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_CTOOLS_SYSTEM_H__
#define __ORATIO_CTOOLS_SYSTEM_H__

/**
 * Routines for system calls
 */

extern char **environ;

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Execute an external Unix command.
   *
   * @param command is the command to execute.
   *
   * @return the status of the command; in particural, return -1 if the command
   *   cannot be executed, and 1 if command is NULL.
   */
  int c_system(char *command);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_CTOOLS_SYSTEM_H__ */
