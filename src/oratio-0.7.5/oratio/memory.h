/**
 * memory.h
 *
 * C Tools
 * Copyright (C) 2000-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_CTOOLS_MEMORY_H__
#define __ORATIO_CTOOLS_MEMORY_H__

#include <stdlib.h>
#include "ctools_macros.h"

/* Macros destinées à faciliter la gestion de la mémoire */
#define _malloc(type) ((type *)c_malloc(sizeof(type)))
#define _mallocs(nb, type) ((type *)c_malloc((nb)*sizeof(type)))
#define _realloc(ptr, type) ((type *)c_realloc(ptr, sizeof(type)))
#define _reallocs(ptr, nb, type) ((type *)c_realloc(ptr, (nb)*sizeof(type)))
#define _free(var) (if (var) free(var));

#define _sizeenum(enum) (sizeof(enum)/sizeof(enum[0]))

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Allocation memory function. This function checks if the memory has been
   * correctly assigned. If not, the program is stopped.
   *
   * @return a pointer to the newly allocated memory space.
   */
  void *c_malloc(size_t size);

  /**
   * Reallocation memory. This function checks if the memory has been
   * correctly assigned. If not, the program is stopped.
   *
   * @return a pointer to the reallocated memory space.
   */
  void *c_realloc(void *ptr, size_t size);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_CTOOLS_MEMORY_H__ */
