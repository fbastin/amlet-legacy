/**
 * cfile.h
 *
 * C Tools
 * Copyright (C) 2000-2010 Fabian Bastin <slashbin@namurlug.org>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __CTOOLS_FILE_H__
#define __CTOOLS_FILE_H__

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Error codes
   */
  enum {
    C_FILE_NOT_EXISTS = 1,
    C_FILE_CANNOT_OPEN
  };

  /**
   * Compute the size of file "name".
   *
   * @param file is a read opened file.
   *
   * @return the file size.
   */
  long c_file_size(FILE *file);

   /**
   * Copy the content of a text file into a character array.
   *  
   * @param name is the name of the file to open.
   * @param length contains on return the number of characters composing the
   *   text, plus one, if length is not NULL on input. If the file could not be
   *   open, and length is not NULL, *length will contain 0 on output.
   * 
   * @return a newly allocated pointer to a character array that contains a
   *   copy of all the text, that composes the text, file of name given by the
   *   pointer name, without any modification, other than the addition of the
   *   special character EOF at the end of the text. If the file could not be
   *   open, return NULL.   
   */
  char *c_file_get_text(char *name, long *length);

  /**
   * Check if the file given by filename exists.
   *
   * @return 1 if the file exists, 0 else.
   */
  int c_file_check_exist(char *filename);

  /**
   * Check if the file given by filename and residing in the directory
   * directory exists.
   *
   * @return NULL is the file doesn't exist, the complete filename 
   *   otherwise. In the last case, the returned string is newly allocated.
   */
  char *c_file_check_exists(char *directory, char *filename);

  /**
   * Eliminate the comments of the line
   */
  void c_file_line_cut_comments(const char *line);

  /**
   * Read the current line of input
   *    
   * @param input is a file open for reading
   * @return a char pointer to a string newly allocated, formed by the
   *    content of the line
   */
  char *c_file_get_line(FILE *input);

  /**
   * Read the current line of input
   *
   * @param buff is a pointer to a string, of size *size; the string
   *    will be reallocated if the preliminary defined size is unsufficient
   * @param input is a file open for reading
   * @param size is a pointer to an integer containing the size of
   *    *buff on input, and set to the size of *buff on output
   * @return a char pointer to a string made of the content of
   *    the line
   */
  char *c_file_get_line_b(char **buff, FILE *input, int *size);

  /**
   * Try to open a file, in basic mode.
   */
  int c_file_open(char *filename, int *fd);

  /**
   * Open a text file.
   *
   * @param name is the name of the file to open.
   */
  FILE *c_fopen(char *name, char *mode);

  /**
   * Close an open text file.
   *
   * @param f is a pointer to an open file, that will be closed on output.
   */
  void c_fclose(FILE *f);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __CTOOLS_FILE_H__ */
