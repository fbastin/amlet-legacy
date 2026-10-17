/**
 * memory.c
 *
 * C Tools
 * Copyright (C) 2000-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/ctools.h"

/**
 * Allocation memory function. This function checks if the memory has been
 * correctly assigned. If not, the program is stopped.
 *
 * @return a pointer to the newly allocated memory space.
 */
void *c_malloc(size_t size)
{
  void *p;

  p = malloc(size);
  if (!p) c_error("sorry, not enough memory");

  return p;
}

/**
 * Reallocation memory. This function checks if the memory has been
 * correctly assigned. If not, the program is stopped.
 *
 * @return a pointer to the reallocated memory space.
 */
void *c_realloc(void *ptr, size_t size)
{
  void *p;

  p = realloc(ptr, size);
  if (!p && !size) c_warning("unable to do memory reallocation");

  return p;
}
