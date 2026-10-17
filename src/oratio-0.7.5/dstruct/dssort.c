/**
 * dssort.c
 *
 * ORATIO Library
 * Copyright (C) 2000-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/oratio.h"

#include <string.h>

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
void *ds_perm(int n, void *array, size_t size_elem, int *perm)
{
  int i;
  void *arperm = NULL;

  if (n > 0) {
    arperm = malloc(n*size_elem);
    if (!arperm) return NULL;
    for (i = 0; i < n; i++)
      memcpy(arperm+i*size_elem, array+perm[i]*size_elem, size_elem);
  }

  return arperm;
}

/**
 * Generate the inverse permutation corresponding to a given permutation.
 *
 * @param n is the number of elements in p and pinv.
 * @param p is a permutation vector.
 * @param pinv is a vector whose content is undefined on input, and
 *   that will contains the inverse permutation on output.
 */
void ds_perm_inv(int n, int *p, int *pinv)
{
  int i;

  for (i = 0; i < n; i++) pinv[p[i]] = i;
}

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
int ds_perm_inp(int n, void *array, size_t size_elem, int *perm)
{
  int i, pos, next;
  void *t;

  if (!(t = malloc(size_elem))) return 0;

      /*  if (n > 0) {
    arperm = malloc(n*size_elem);
    if (!arperm) return NULL;
    for (i = 0; i < n; i++)
      memcpy(arperm+i*size_elem, array+perm[i]*size_elem, size_elem);
      }*/

  for (i = 0; i < n; i++) {
    pos = i;

    /* We walk through one cycle. */
    while (perm[pos] >= 0) {
      next = perm[pos];
      perm[pos] = -perm[pos]-1;
      if (perm[next] >= 0) {
	memcpy(t, array+pos*size_elem, size_elem);
	memcpy(array+pos*size_elem, array+next*size_elem, size_elem);
	memcpy(array+next*size_elem, t, size_elem);
	pos = next;
      }
    }
  }

  /* Reset permutation vector */
  for (i = 0; i < n; i++) perm[i] = -perm[i]-1;

  free(t);

  return 1;
}

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
void ds_sort_dheapsort(int n, double *array)
{
  /**
   * The implementation is the one proposed on Wikipedia (2006/08/17).
   */
  int i = n/2, parent, child;
  double t;
 
  for (;;) { /* Loops until arr is sorted */
    if (i > 0) { /* First stage - Sorting the heap */
      i--;           /* Save its index to i */
      t = array[i];    /* Save parent value to t */
    } else {     /* Second stage - Extracting elements in-place */
      n--;           /* Make the new heap smaller */
      if (n == 0) return; /* When the heap is empty, we are done */
      t = array[n];    /* Save last value (it will be overwritten) */
      array[n] = array[0]; /* Save largest value at the end of arr */
    }
 
    parent = i; /* We will start pushing down t from parent */
    child = i*2 + 1; /* parent's left child */
 
    /* Shift operation - pushing the value of t down the heap */
    while (child < n) {
      if (child + 1 < n  &&  array[child + 1] > array[child]) {
	child++; /* Choose the largest child */
      }
      if (array[child] > t) { /* If any child is bigger than the parent */
	array[parent] = array[child]; /* Move the largest child up */
	parent = child; /* Move parent pointer to this child */
	child = parent*2 + 1; /* Find the next child */
      } else
	break; /* t's place is found */
    }
    array[parent] = t; /* We save t in the heap */
  }
}

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
void ds_sort_iheapsort(int n, int *array)
{
  /**
   * The implementation is the one proposed on Wikipedia (2006/08/17).
   */

  int i = n/2, parent, child;
  int t;

  for (;;) { /* Loops until arr is sorted */
    if (i > 0) { /* First stage - Sorting the heap */
      i--;           /* Save its index to i */
      t = array[i];    /* Save parent value to t */
    } else {     /* Second stage - Extracting elements in-place */
      n--;           /* Make the new heap smaller */
      if (n == 0) return; /* When the heap is empty, we are done */
      t = array[n];    /* Save last value (it will be overwritten) */
      array[n] = array[0]; /* Save largest value at the end of arr */
    }
 
    parent = i; /* We will start pushing down t from parent */
    child = i*2 + 1; /* parent's left child */
 
    /* Shift operation - pushing the value of t down the heap */
    while (child < n) {
      if (child + 1 < n  &&  array[child + 1] > array[child]) {
	child++; /* Choose the largest child */
      }
      if (array[child] > t) { /* If any child is bigger than the parent */
	array[parent] = array[child]; /* Move the largest child up */
	parent = child; /* Move parent pointer to this child */
	child = parent*2 + 1; /* Find the next child */
      } else
	break; /* t's place is found */
    }
    array[parent] = t; /* We save t in the heap */
  }
}

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
		       int(*compar)(const void *, const void *))
{
  int i = nmemb/2, parent, child;
  void *t;

  t = malloc(size);
 
  for (;;) {
    if (i > 0) {
      i--;
      memcpy(t, base+i*size, size);
    } else {
      nmemb--;
      if (nmemb == 0) return;
      memcpy(t, base+nmemb*size, size);
      memcpy(base+nmemb*size, base, size);
    }
 
    parent = i;
    child = i*2 + 1;
 
    while (child < nmemb) {
      if (child + 1 < nmemb  &&
	  (compar(base+(child+1)*size, base+child*size) > 0)) {
	child++;
      }
      if (compar(base+child*size, t) > 0) {
	memcpy(base+parent*size, base+child*size, size);
	parent = child;
	child = parent*2 + 1;
      } else
	break;
    }
    memcpy(base+parent*size, t, size);
  }

  free(t);
}

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
void ds_sort_dheapsort_with_perm(int n, double *array, int *perm)
{
  /**
   * The implementation follows the one proposed on Wikipedia (2006/08/17).
   */
  int i = n/2, parent, child, idx;
  double t;
 
  /* Initialize the permutation vector. */
  for (i = 0; i < n; i++) perm[i] = i;

  for (;;) { /* Loops until arr is sorted */
    if (i > 0) { /* First stage - Sorting the heap */
      i--;           /* Save its index to i */
      t = array[i];    /* Save parent value to t */
      idx = perm[i];
    } else {     /* Second stage - Extracting elements in-place */
      n--;           /* Make the new heap smaller */
      if (n == 0) return; /* When the heap is empty, we are done */
      t = array[n];    /* Save last value (it will be overwritten) */
      array[n] = array[0]; /* Save largest value at the end of arr */
      idx = perm[n];
      perm[n] = perm[0];
     }
 
    parent = i; /* We will start pushing down t from parent */
    child = i*2 + 1; /* parent's left child */
 
    /* Shift operation - pushing the value of t down the heap */
    while (child < n) {
      if (child + 1 < n  &&  array[child + 1] > array[child]) {
	child++; /* Choose the largest child */
      }
      if (array[child] > t) { /* If any child is bigger than the parent */
	array[parent] = array[child]; /* Move the largest child up */
	perm[parent] = perm[child];
	parent = child; /* Move parent pointer to this child */
	child = parent*2 + 1; /* Find the next child */
      } else
	break; /* t's place is found */
    }
    array[parent] = t; /* We save t in the heap */
    perm[parent] = idx;
  }
}

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
void ds_sort_iheapsort_with_perm(int n, int *array, int *perm)
{
  /**
   * The implementation follows the one proposed on Wikipedia (2006/08/17).
   */
  int i = n/2, parent, child, idx;
  int t;
 
  /* Initialize the permutation vector. */
  for (i = 0; i < n; i++) perm[i] = i;

  for (;;) { /* Loops until arr is sorted */
    if (i > 0) { /* First stage - Sorting the heap */
      i--;           /* Save its index to i */
      t = array[i];    /* Save parent value to t */
      idx = perm[i];
    } else {     /* Second stage - Extracting elements in-place */
      n--;           /* Make the new heap smaller */
      if (n == 0) return; /* When the heap is empty, we are done */
      t = array[n];    /* Save last value (it will be overwritten) */
      array[n] = array[0]; /* Save largest value at the end of arr */
      idx = perm[n];
      perm[n] = perm[0];
     }
 
    parent = i; /* We will start pushing down t from parent */
    child = i*2 + 1; /* parent's left child */
 
    /* Shift operation - pushing the value of t down the heap */
    while (child < n) {
      if (child + 1 < n  &&  array[child + 1] > array[child]) {
	child++; /* Choose the largest child */
      }
      if (array[child] > t) { /* If any child is bigger than the parent */
	array[parent] = array[child]; /* Move the largest child up */
	perm[parent] = perm[child];
	parent = child;       /* Move parent pointer to this child */
	child = parent*2 + 1; /* Find the next child */
      } else
	break; /* t's place is found */
    }
    array[parent] = t; /* We save t in the heap */
    perm[parent] = idx;
  }
}

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
				 int *perm)
{
  int i = nmemb/2, parent, child, idx;
  void *t;

  t = malloc(size);
  
  /* Initialize the permutation vector. */
  for (i = 0; i < nmemb; i++) perm[i] = i;

  for (;;) {
    if (i > 0) {
      i--;
      memcpy(t, base+i*size, size);
      idx = perm[i];
    } else {
      nmemb--;
      if (nmemb == 0) return;
      memcpy(t, base+nmemb*size, size);
      memcpy(base+nmemb*size, base, size);
      idx = perm[nmemb];
      perm[nmemb] = perm[0];
    }
 
    parent = i;
    child = i*2 + 1;
 
    while (child < nmemb) {
      if (child + 1 < nmemb  &&
	  (compar(base+(child+1)*size, base+child*size) > 0)) {
	child++;
      }
      if (compar(base+child*size, t) > 0) {
	memcpy(base+parent*size, base+child*size, size);
	perm[parent] = perm[child];
	parent = child;
	child = parent*2 + 1;
      } else
	break;
    }
    memcpy(base+parent*size, t, size);
    perm[parent] = idx;
  }

  free(t);
}
