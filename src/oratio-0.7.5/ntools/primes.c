/**
 * primes.c
 *
 * ORATIO Library
 * Copyright (C) 2004-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/ctools.h"

#include <math.h>

/**
 * Compute the n first prime numbers using successive divisions.
 *
 * @param n is a non-null natural.
 * @param primes is an allocated array of at least n elements. On
 *   output, primes contains the n first prime numbers, by ascending
 *   order, in the first n positions.
 */
void nt_primes_simple(unsigned int n, unsigned int *primes)
{
  int i, j, k, isprime;

  primes[0] = 2;
  i = 2;
  j = 3;

  while (i < n) {
    k = 0;
    isprime = 1;

    while (primes[k]*primes[k] <= j && isprime)
      if (!(j % primes[k++])) isprime = 0;

    if (isprime) primes[i++] = j;

    j += 2;
  }
}

/**
 * Compute the n first prime numbers using the Eratoshtene sieve.
 *
 * @param n is a non-null natural.
 * @param primes is an allocated array of at least n elements. On
 *   output, primes contains the n first prime numbers, by ascending
 *   order, in the first n positions.
 */
void nt_primes_eratosthene(unsigned int n, unsigned int *primes)
{
  /* The implementation could probably be refined in order to obtain
     better performances, both in terms of speed and of memory
     consumption. Any suggestion/modification is welcome. */

  /* In order to save space, we should record the array isprime
     as an array of bits */
  /* We could also use linked lists, in order to free the first part
     of the sequence when possible, and to manage to size of primes. */

  unsigned int max;     /* upper bound of the n-th prime number */
  char *isprime;
  int i = 0;            /* number of found prime numbers */
  int j;                /* array loop index */
  int p;

  /* We use the upper bounds given by Massias and Robin:
     J. Massias and G. Robin, "Bornes effectives pour certaines
     fonctions concernant les nombres premiers," J. Théorie Nombres
     Bordeaux, 8 (1996) 215-242. */
  if (n <= 12)
    /* The 12th prime number is 37 */
    max = 37;
  else if (n < 15985)
    max = (n*(log(n)+log(log(n-1))+1.8*(log(log(n))/log(n))));
  else
    max = (n*(log(n)+log(log(n))-0.9427));

  isprime = _mallocs(max+1, char);

  /* we initialize the array isprime as if each number (>= 2) was prime */
  for (j = 2; j <= max; j++) isprime[j] = 1;

  p = 2;
  for (i = 0; i < n; i++) {
    while (!isprime[p]) p++;
    /* record the prime number */
    primes[i] = p;
    /* mark the multiples of p as being non prime */
    j = p << 1;
    for (j = p << 1 ; j <= max; j += p) isprime[j] = 0;
    /* go to the next integer */
    p++;
  }

  free(isprime);
}
