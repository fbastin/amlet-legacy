/**
 * symbol.h
 *
 * C Tools
 * Copyright (C) 2002-2006 Fabian Bastin <bastin@cerfacs.fr>
 *
 * Licensed under the MIT License
 */

#ifndef __CTOOLS_SYMBOL_H__
#define __CTOOLS_SYMBOL_H__

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

#ifdef __cplusplus
}
#endif

#endif  /* __CTOOLS_SYMBOL_H__ */
