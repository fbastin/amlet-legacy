/**
 * symbol.h
 *
 * ORATIO Library
 * Copyright (C) 2002-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_CTOOLS_SYMBOL_H__
#define __ORATIO_CTOOLS_SYMBOL_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * c_symbol structure.
   */
  typedef struct {
    char *name;
    int token;
  } c_symbol;

  /**
   * Get the name of a token.
   *
   * @param token is a token identifier.
   * @param symbols is the set of tokens where we search the token.
   * @param n is the number of tokens.
   *
   * @return the name associated to the symbol defined by token, or NULL
   *   is the token has not been found.
   */
  char *c_symbol_get_name(int token, c_symbol *symbols, int n);

  /**
   * Get a token associated to a given name.
   *
   * @param name is some token name (case sensitive).
   * @param symbols is the set of tokens where we search the token.
   * @param n is the number of tokens.
   *
   * @return the token index or INT_MIN (a quite improbable - but not
   *   impossible - value) if no token corresponding to name has been
   *   found.
   */
  int c_symbol_get_token(char *name, const c_symbol *symbols, int n);

  /**
   * Check if a specified string matches one in some set of strings.
   *
   * @param s is the string we want to check.
   * @param symbols is a set of c_symbol's.
   * @param n is the number of elements in symbols.
   * @param pos is a pointer to an integer that will contain on output
   *  the position of s in the symbols array, or -1 if it has not been found.
   *
   * @return the corresponding token index or INT_MIN (a quite improbable -
   *   but not impossible - value) if no token corresponding to s has been
   *   found.
   */
  int c_symbol_search(char *s, const c_symbol symbols[], int n, int *pos);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_CTOOLS_SYMBOL_H__ */
