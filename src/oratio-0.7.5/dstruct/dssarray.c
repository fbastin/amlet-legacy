/**
 * dssarray.c
 *
 * ORATIO Library
 * Copyright (C) 2000-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/ctools.h"
#include "oratio/dstruct.h"

#include <string.h>

#define BASE_MEM 32

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
void *ds_sarray_new(size_t size_elem, int *size)
{
  void *array;

  if (*size <= 1) *size = BASE_MEM;  
  array = c_malloc((*size)*size_elem);

  return array;
}

/**
 * Free the memory associated to a generic array.
 *
 * @param n is the number of elements in array.
 * @param array is a pointer to the sorted array.
 * @param size_elem is the size (in bytes) of an element of the array.
 * @param ff is the function to apply to the an element of the array
 *   (can be set to NULL).
 */
void ds_sarray_free(int n, void *array, size_t size_elem, Dfct ff)
{
  int i;

  if (array) {
    if (ff) for (i = 0; i < n; i++) ff(array+i*size_elem);
    free(array);
  }

  array = NULL;
}

/**
 * Insert a new element into a generic sorted array. If the element is
 * already in the array, no insertion is performed.
 *
 * @param n is a pointer to the number of elements in the array; the value *n
 *   is incremented if the element has been properly inserted.
 * @param array is a pointer to the sorted array.
 * @param size_elem is the memory size of one arrey element.
 * @param f is a comparison function between two elements a, b, considered in
 *   that order, that returns a strictly negative value if a < b, a strictly
 *   positive value if a > b, and 0 if a = b.
 * @param size_max if a pointer to the maximum number of array elements.
 *
 * @return a pointer to the sorted array, possibly different different to the
 *   original pointer address if the memory allocated to the array has been
 *   increased.
 */
void *ds_sarray_insert(int *n, void *array, void *data, size_t size_elem,
		       DSCompareFunc f, int *size_max)
{
  int pos, i;

  /* Search the position for the element */
  pos = ds_find_inter(array, *n, 1, size_elem, data, f);

  if (pos > -1) {
    /* Check if the element is already present in the array */
    if (f(array+pos*size_elem, data) == 0) return array;
  }

  (*n)++;
  /* Check if the memory limit has to be increased */
  if (*n > *size_max) {
    if (!(*size_max)) *size_max = 1;
    *size_max = 2*(*size_max);
    array = realloc(array, *(size_max)*size_elem);
  }

  /* Shift the elements greater than data one position to the right. */
  for (i = *n; i > pos+1; i--)
    memcpy (array+i*size_elem, array+(i-1)*size_elem, size_elem);
  /* Copy the element */
  memcpy(array+(pos+1)*size_elem, data, size_elem);

  return array;
}

int ds_sarray_search(int n, void *array, size_t size_elem, void
		     *data, DSCompareFunc f)
{
  /* Dichotomic search */
  return ds_find_dicho(array, n, 1, size_elem, data, f);
}
