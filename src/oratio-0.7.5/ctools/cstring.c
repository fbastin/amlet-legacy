/**
 * cstring.c
 *
 * ORATIO library
 * Copyright (C) 2000-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include <oratio/ctools.h>

/**
 * Reverse a string.
 *
 * @param s is the string to reverse.
 * @param len is the len of the string.
 */
void c_string_reverse(char *s, int len)
{
  int i, j;
  char c;
 
  for (i = 0, j = len-1; i<j; i++, j--) {
    c = s[i];
    s[i] = s[j];
    s[j] = c;
  }
}

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
char *c_string_itoa(char *s, int i)
{
  /* We just follow K&R implementation */
  int sign;
  int len, size;
  int k;
 
  size = 4;
  s = _reallocs(s, size, char);

  if ((sign = i) < 0)  /* record sign */
    i *= -1;           /* make n positive */

  len = 0;
  do {       /* generate digits in reverse order */
    s[len++] = i % 10 + '0';   /* get next digit */
    if (len == size) {
      size *= 2;
      s = _reallocs(s, size, char);
    }
  } while ((i /= 10) > 0);     /* delete it */

  /* compute how many characters have to be added and possibly adjust size */
  k = (sign < 0) ? 2 : 1;
  if (len+k >= size) {
    size += k;
    s = _reallocs(s, size, char);    
  }

  if (sign < 0)
    s[len++] = '-';

  s[len] = '\0';
  s = _reallocs(s, len, char);

  c_string_reverse(s, len);

  return s;
}

/**
 * Skip initial spaces in string s.
 *
 * @param s is the string under consideration.
 */
char *c_string_skip_spaces(char *s)
{
  if (s) while (isblank(*s)) s++;

  return s;
}

/**
 * Prepend the string pre to the string dest.
 *
 * @param dest is the original string.
 * @param pre is to string to add at the end of dest.
 *
 * @return a pointer to the newly formed film.
 */
char *c_string_prepend(char *dest, char *pre)
{
  int n = _strlen(dest);

  dest = _reallocs(dest, _strlen(dest)+_strlen(pre)+1, char);
  memmove(dest+_strlen(pre), dest, n*sizeof(char));
  memcpy(dest, pre, _strlen(pre)*sizeof(char));

  return dest;
}

/**
 * Append the string post to the string dest.
 *
 * @param dest is a pointer to the original string. dest can be set to
 *   NULL.
 * @param post is a pointer to the string to append.
 *
 * @return a pointer to the newly formed string.
 */
char *c_string_append(char *dest, char *post)
{
  dest = _reallocs(dest, _strlen(dest)+_strlen(post)+1, char);
  strcat(dest, post);

  return dest;
}
