/**
 * dssarray.h
 *
 * ORATIO Library
 * Copyright (C) 2000-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_DSTRUCT_DSSARRAY_H__
#define __ORATIO_DSTRUCT_DSSARRAY_H__

#include <unistd.h>

/**
 * Sorted array routines
 */

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Create a new array structure.
   *
   * @param size_elem is the size of one element in the array.
   * @param size is a pointer to an integer designing the number of
   *   elements to be allocated; is 0 is given, a default value is used
   *   and stored in *size.
   *
   * @return a pointer to the newly allocated array.
   */
  void *ds_sarray_new(size_t size_elem, int *size);

  /**
   * Free the memory associated to a generic array.
   *
   * @param n is the number of elements in array.
   * @param size_elem is the size (in bytes) of an element of the array.
   * @param ff is the function to apply to the an element of the array
   *   (can be set to NULL).
   */
  void ds_sarray_free(int n, void *array, size_t size_elem, Dfct ff);

  /**
   * Insert a new element into a generic array.
   *
   * @param n is the number of elements in the array.
   *        size_elem: size of one element in the array
   * #POST: n incremented by 1
   */
  void *ds_sarray_insert(int *n,
			 void *array, void *data, size_t size_elem,
			 DSCompareFunc f, int *size_max);

  /**
   * @param n is the number of elements in array.
   */
  int ds_sarray_search(int n, void *array, size_t size_elem, void
		       *data, DSCompareFunc f);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_DSTRUCT_DSSARRAY_H__ */
