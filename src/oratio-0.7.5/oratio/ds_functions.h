/**
 * ds_functions.h
 *
 * ORATIO Library
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_DSTRUCT_FUNCTIONS_H__
#define __ORATIO_DSTRUCT_FUNCTIONS_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * General utilities
   */

  /**
   * Swap variables a and b; in the calling function, the variable
   * swap must be declared, and a, b, and swap must share the same type.
   */
#define SWAP(a,b) { swap=(a); (a)=(b); (b)=swap; }

  /**
   * Functions type definition.
   */

  typedef void (*Dfct)(void *data);

  /**
   * Comparison functions type.
   *
   * @return 0 if a == b, a striclty negative number if a < b,
   *   a strictly positive number if a > b.
   */
  typedef int (*DSCompareFunc)(const void *a, const void *b);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_DSTRUCT_FUNCTIONS_H__ */
