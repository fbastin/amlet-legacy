/**
 * system.c
 *
 * ORATIO Library
 * Copyright (C) 2000-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <oratio/ctools.h>

/**
 * Execute an external Unix command.
 *
 * @param command is the command to execute.
 *
 * @return the status of the command; in particural, return -1 if the command
 *   cannot be executed, and 1 if command is NULL.
 */
int c_system(char *command)
{
#ifdef _WIN32
  fprintf(stderr, "c_system is not supported under Windows\n");
  return -1;
#else
  /**
   * This code is based on the system manpage of Linux
   */
  int pid, status;

  if (command == 0) return 1;

  pid = fork();
  if (pid == -1) return -1;
  if (pid == 0) {
    char *argv[4];
    argv[0] = "sh";
    argv[1] = "-c";
    argv[2] = command;
    argv[3] = 0;
    /* environ designs the environment variables */
    execve("/bin/sh", argv, environ);
    exit(127);
  }
  do {
    if (waitpid(pid, &status, 0) == -1) {
      if (errno != EINTR)
	return -1;
    }
    else
      return status;
  } while(1);
#endif
}
