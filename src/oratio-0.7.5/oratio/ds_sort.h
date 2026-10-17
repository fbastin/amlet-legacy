/**
 * dssort.h
 *
 * D Struct
 * Copyright (C) 2005-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_DSTRUCT_DSSORT_H__
#define __ORATIO_DSTRUCT_DSSORT_H__

#include "dstruct.h"

#include <unistd.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * SORTING FUNCTIONS
   */

  /**
   * Heapsort (double array).
   *
   * Sort the vector array[0..n-1] in numerically ascendent order
   * using the heaposrt algorithm.
   *
   * @param n is the dimension of array.
   * @param array is the vector to sort, and on output, it is the
   *   sorted vector.
   */
  void ds_sort_dheapsort(int n, double *array);

  /**
   * Heapsort (integer array).
   *
   * Sort the vector array[0..n-1] in numerically ascendent order
   * using the heaposrt algorithm.
   *
   * @param n is the dimension of array.
   * @param array is the vector to sort, and on output, it is the
   *   sorted vector.
   */
  void ds_sort_iheapsort(int n, int *array);

  /**
   * Sort a generic array using the heapsort algorithm.
   *
   * @param base is a pointer to the start of the array to sort on
   *   input, and to the sorted array on output.
   * @param nmemb is the number of elements in the array pointed by base.
   * @param size is the size (in bytes) of an element of the array.
   * @param compar is a comparison function that must return an integer
   *   less  than, equal  to,  or  greater than zero if the first argument
   *   is considered to be respectively  less  than,  equal  to, or
   *   greater than the second.
   */
  void ds_sort_gheapsort(void *base, size_t nmemb, size_t size,
			 int(*compar)(const void *, const void *));

  /**
   * Heapsort with permutation vector (double version).
   *
   * The use of this function is similar to a call to
   * ds_sort_heapsort, but it also computes the permutation vector
   * corresponding to the sorting.
   *
   * @param n is the dimension of array.
   * @param array is on input the array to sort, and on output the
   *   sorted vector, such that array_orig[perm[i]] = array[i], where
   *   array_orig[i] design the unsorted vector.
   */
  void ds_sort_dheapsort_with_perm(int n, double *array, int *perm);

  /**
   * Heapsort with permutation vector (integer version).
   *
   * The use of this function is similar to a call to
   * ds_sort_heapsort, but it also computes the permutation vector
   * corresponding to the sorting.
   *
   * @param n is the dimension of array.
   * @param array is on input the array to sort, and on output the
   *   sorted vector, such that array_orig[perm[i]] = array[i], where
   *   array_orig[i] design the unsorted vector.
   */
  void ds_sort_iheapsort_with_perm(int n, int *array, int *perm);

  /**
   * Sort a generic array using the heapsort algorithm, and generates
   * the corresponding permutation vector.
   *
   * @param base is a pointer to the start of the array to sort on
   *   input, and to the sorted array on output.
   * @param nmemb is the number of elements in the array pointed by base.
   * @param size is the size (in bytes) of an element of the array.
   * @param compar is a comparison function that must return an integer
   *   less  than, equal  to,  or  greater than zero if the first argument
   *   is considered to be respectively  less  than,  equal  to, or
   *   greater than the second.
   */
  void ds_sort_gheapsort_with_perm(void *base, size_t nmemb, size_t size,
				   int(*compar)(const void *, const void *),
				   int *perm);

  /**
   * PERMUTATION FUNCTIONS.
   */

  /**
   * Array permutation.
   *
   * Apply the permutation perm to the vector array, and produce a newly
   * allocated vector.
   *
   * @param n is the dimension of array.
   * @param array is the vector to permute. On output, array is
   *   unchanged.
   * @param size_elem is the size (in byte) of one element of array.
   * @param perm is the vector of permutation indexes. perm[i] indicates
   *   the position of the element originally at i in the permuted
   *   vector, i.e. if v is the permuted vector, v[i] = array[p[i]].
   *
   * @return a newly allocated vector, corresponding to permuted array vector.
   */
  void *ds_perm(int n, void *array, size_t size_elem, int *perm);

  /**
   * Generate the inverse permutation corresponding to a given permutation.
   *
   * @param n is the number of elements in p and pinv.
   * @param p is a permutation vector.
   * @param pinv is a vector whose content is undefined on input, and
   *   that will contains the inverse permutation on output.
   */
  void ds_perm_inv(int n, int *p, int *pinv);

  /**
   * In place array permutation.
   *
   * Apply the permutation perm to the vector array, which is directly
   * modified.
   *
   * @param n is the dimension of array.
   * @param array is the vector to permute. On output, array is
   *   the permuted vector.
   * @param size_elem is the size (in byte) of one element of array.
   * @param perm is the vector of permutation indexes. perm[i] indicates
   *   the position of the element originally at i in the permuted
   *   vector, i.e. if v is the permuted vector, v[i] = array[p[i]].
   *
   * @return 0 if an error has occured (memory allocation), 0 else.
   */
  int ds_perm_inp(int n, void *array, size_t size_elem, int *perm);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_DSTRUCT_DSSORT_H__ */
