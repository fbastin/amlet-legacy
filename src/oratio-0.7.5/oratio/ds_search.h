/**
 * dssearch.h
 *
 * ORATIO Library
 * Copyright (C) 2000-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_DSTRUCT_DSSEARCH_H__
#define __ORATIO_DSTRUCT_DSSEARCH_H__

#include "dstruct.h"

#include <unistd.h>

/**
 * Sorted array routines
 */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Compare the double values contained at adresses a and b.
   */
  int ds_search_double(const void *a, const void *b);

  /**
   * Compare the integer values contained at adresses a and b.
   */
  int ds_search_int(const void *a, const void *b);

  int ds_search_int_opposite(const void *a, const void *b);

  /**
   * Dichotomy search.
   *
   * @param base is a pointer to the start of the array on which the
   *   search has to be performed. The array must be sorted in
   *   ascending order.
   * @param n is the number of elements if the array.
   * @param inc is the array increment.
   * @param size is the size in bytes of one array element.
   * @param search is a comparison function between two elements, say a and b,
   *   that returns 0 if a == b, a striclty negative number if a < b, a
   *   strictly positive number if a > b. The type of the elements must
   *   correspond to the array type.
   *
   * @return the element index corresponding to value, or -1 if it has not
   *   been found.
   */
  int ds_find_dicho(void *base, int n, int inc, size_t size,
		    void *value, DSCompareFunc search);

  /**
   * Interval search
   *
   * @param base is a pointer to the start of the array on which the
   *   search has to be performed. The array must be sorted in
   *   ascending order.
   * @param n is the number of elements if the array.
   * @param inc is the array increment.
   * @param size is the size in bytes of one array element.
   * @param search is a comparison function between two elements, say a and b,
   *   that returns 0 if a == b, a striclty negative number if a < b, a
   *   strictly positive number if a > b. The type of the elements must
   *   correspond to the array type.
   *
   * @return an integer i greater or equal to 0, such that the value is in 
   *         [base+i*size, base+(i+1)*size[, or [base+(nmemb-1)*size, \infty[
   *         or -1 if no interval can be determined.
   */
  int ds_find_inter(void *base, int n, int inc, size_t size,
		    void *value, DSCompareFunc search);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_DSTRUCT_DSSEARCH_H__ */
