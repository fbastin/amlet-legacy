/**
 * cstring.h
 *
 * ORATIO library
 * Copyright (C) 2000-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef  __ORATIO_CTOOLS_CSTRING_H__
#define  __ORATIO_CTOOLS_CSTRING_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * String management functions.
   */

  /**
   * Compute the length of the characters string str.
   *
   * @return return the length of str, or 0 if str is NULL.
   */
#define _strlen(str) ((str) ? strlen(str) : 0)

  /**
   * Reverse a string.
   *
   * @param s is the string to reverse.
   * @param len is the len of the string.
   */
  void c_string_reverse(char *s, int len);

  /**
   * Convert an integer into a (newly memory-allocated) string.
   *
   * @param s is the string that will containt the representation of the
   *   integer; s can be set to NULL on output, and will have the
   *   exactly needed size on output.
   * @param i is the integer to convert.
   * 
   * @return a newly allocated string containing a text representation of
   *   the integer.
   */
  char *c_string_itoa(char *s, int i);

  /**
   * Skip initial spaces in string s.
   *
   * @param s is the string under consideration.
   */
  char *c_string_skip_spaces(char *s);

  /**
   * Prepend the string pre to the string dest.
   *
   * @param dest is the original string.
   * @param pre is to string to add at the end of dest.
   *
   * @return a pointer to the newly formed film.
   */
  char *c_string_prepend(char *dest, char *pre);

  /**
   * Append the string post to the string dest.
   *
   * @param dest is a pointer to the original string. dest can be set to
   *   NULL.
   * @param post is a pointer to the string to append.
   *
   * @return a pointer to the newly formed string.
   */
  char *c_string_append(char *dest, char *post);

#ifdef __cplusplus
}
#endif

#endif /*  __ORATIO_CTOOLS_CSTRING_H__ */
