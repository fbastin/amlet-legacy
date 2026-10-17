/**
 * credits.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2015 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __AMLET_CREDITS_H__
#define __AMLET_CREDITS_H__

#include <oratio/ctools.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Print library version.
   */
  void amlet_version();

  /**
   * Print library credits.
   *
   * @param f is the output file (assumed to be open).
   * @param start is the string to use at the beginning of each line.
   *   This variable will often be set to "".
   * @param end is the string to use at the end of each line.
   *   This variable will often be set to "".
   */
  void amlet_credits(FILE *f, char *start, char *end);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_CREDITS_H__ */
