/**
 * file.c
 *
 * C Tools
 * Copyright (C) 2000-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/cfile.h"
#include "oratio/memory.h"

#include <fcntl.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#if defined(WINDOWS)
#ifndef C_DIR_SEPARATOR_S
#define C_DIR_SEPARATOR_S "\\"
#endif
#else
#ifndef C_DIR_SEPARATOR_S
#define C_DIR_SEPARATOR_S "/"
#endif
#endif

#ifndef C_COMMENTS
#define C_COMMENTS '#'
#endif

/**
 * Compute the size of a file.
 *
 * @param file is a read opened file.
 *
 * @return the file size.
 */
long c_file_size(FILE *file)
{
  long pos, size;

  pos = ftell(file);
  fseek(file, 0L, SEEK_END);
  size = ftell(file);
  fseek(file, pos, SEEK_SET);

  return size;
}

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
char *c_file_get_text(char *name, long *length)
{
  FILE *file;
  long size;
  char *buf = NULL;

  if ((file = fopen(name, "r"))) {
    fseek(file, 0L, SEEK_END);
    size = ftell(file);
    fseek(file, 0L, SEEK_SET);

    buf = _mallocs(size+1, char);
    fread(buf, size, sizeof(char), file);
    buf[size] = EOF;

    if (length) *length = size;

    fclose(file);
  }
  else
    if (length) *length = 0;

  return buf;
}

/**
 * Eliminate the comments of the line
 */
void c_file_line_cut_comments(const char *line)
{
  char *s;

  s = strchr(line, '#');
  if (s) *s = '\0';
}

/**
 * Read the current line of input
 *
 * @param input is a read-open file.
 *
 * @return a char pointer to a string newly allocated, formed by the
 *    content of the line
 */
char *c_file_get_line(FILE *input)
{
  char *buff;
  int i = 512, j=0, c='1';

  if (feof(input)) return NULL;

  buff = _mallocs(i, char);

  do {
    c = getc(input);
    if (c == EOF) c = '\0';
    if (j < i-1)
      buff[j] = c;
    else {
      /* need more memory */
      i <<= 1;
      buff = _reallocs(buff, i, char);
      buff[j] = c;
    }
    j++;
  } while (c != '\n' && c != '\0');
  buff[j] = '\0';

  return buff;
}

/**
 * Read the current line of input
 *
 * @param buff is a pointer to a string, of size *size; the string
 *    will be reallocated if the preliminary defined size is unsufficient
 * @param input is a file open for reading
 * @param size is a pointer to an integer containing the size of
 *    *buff on input, and set to the size of *buff on output
 *
 * @return a char pointer to a string made of the content of
 *    the line
 */
char *c_file_get_line_b(char **buff, FILE *input, int *size)
{
  int j=0, c='1';

  if (!(*size)) *size = 2;

  if (feof(input)) {
    free(*buff);
    *buff = NULL;
  }
  else {
    do {
      c = getc(input);
      if (c == EOF) c = '\0';
      if (j < *size-1)
	(*buff)[j] = c;
      else {
	/* need more memory; we double *size */
	*size <<= 1;
	*buff = _reallocs(*buff, *size, char);
	(*buff)[j] = c;
      }
      j++;
    } while (c != '\n' && c != '\0');
    (*buff)[j] = '\0';
  }

  return *buff;
}

/**
 * Check if the file given by filename exists.
 * 
 * @param filename is the name of the file.
 *  
 * @return 1 if the file exists, 0 else.
 */
int c_file_check_exist(char *filename)
{
  struct stat s;

#ifdef S_ISLNK
  int status = lstat(filename, &s);
  if (!status && (S_ISREG(s.st_mode) || S_ISLNK(s.st_mode)))
    return 1;
#else
  if (!_stat(filename, &s)) return 1;
#endif

  return 0;
}

/**
 * Check if the file given by filename and residing in the directory
 * directory exists.
 *
 * @param filename is the name of the file.
 *
 * @return NULL is the file doesn't exist, the complete filename 
 *   otherwise. In the last case, the returned string is newly allocated.
 */
char *c_file_check_exists(char *directory, char *filename)
{
  char *full_filename;

  /* TODO: use dirname and basename defined in libgen.h? */

  full_filename = _mallocs(strlen(directory)+strlen(filename)+2, char);
  strcpy(full_filename, directory);
  strcat(full_filename, C_DIR_SEPARATOR_S);
  strcat(full_filename, filename);

  if (c_file_check_exist(full_filename))
    return full_filename;
  else
    return NULL;
}

/**
 * Try to open a file, in basic mode.
 *
 * @param filename is the name of the file.
 */
int c_file_open(char *filename, int *fd)
{
  *fd = open(filename, O_RDONLY);
  if (*fd < 0) return C_FILE_CANNOT_OPEN;
  else return 0;
}

/**
 * Open a text file.
 *
 * @param name is the name of the file to open.
 */
FILE *c_fopen(char *name, char *mode)
{
  FILE *f = NULL;
  char *buff;

  if (!(f = fopen(name, mode))) {
    if ((buff = c_file_check_exists(getenv("PWD"), name))) {
      f = fopen(buff, mode);
      free(buff);
    }
  }

  return f;
}

/**
 * Close an open text file.
 *
 * @param f is a pointer to an open file, that will be closed on output.
 */
void c_fclose(FILE *f)
{
  if (f) fclose(f);
}
