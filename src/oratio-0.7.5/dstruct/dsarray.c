/**
 * dsarray.c
 *
 * ORATIO Library
 * Copyright (C) 2000-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/dstruct.h"

/**
 * Compute the cumulated frequency of the array elements that are
 * inferior to some criterion.
 *
 * @param base is a pointer towards a sorted array of nmemb elements,
 *    each of size "size".
 * @param search is a comparison function between two elements, say a and b,
 *   that returns 0 if a == b, a striclty negative number if a < b, a
 *   strictly positive number if a > b. The type of the elements must
 *   correspond to the array type.
 */
double ds_array_cum(void *base, size_t nmemb, size_t size, void *criterion,
		    DSCompareFunc search)
{
  int idx;
  double cum;

  idx = ds_find_dicho(base, nmemb, 1, size, criterion, search);
  if (idx == -1)
    cum = (search(base, criterion) < 0) ? 1.0 : 0.0;
  else {
    int done = 0;
    while (!done) {
      if (idx <= 0)
	done = 1;
      else {
	if (search((char *)base+(idx-1)*size, criterion) == 0)
	  idx--;
	else
	  done = 1;
      }
    }
    cum = (double)idx/(double)nmemb;
  }

  return cum;
}
