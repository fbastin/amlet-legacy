/**
 * primes_test.c
 *
 * ORATIO Library
 * Copyright (C) 1999-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "ntools_test.h"

#include <stdio.h>

enum {
  ERROR = 1,
  ERROR_PRIMES
};

#define PRIME        10000
#define PRIME_VALUE  541      /* 100th prime number */

int primes_test()
{
  unsigned int *primes;

  printf("Generate the %d first prime numbers\n", PRIME);

  primes = _mallocs(PRIME, unsigned int);

  nt_primes_eratosthene(PRIME, primes);

  printf("10000th prime: %d\n", primes[PRIME-1]);

  if (primes[99] != PRIME_VALUE) {
    free(primes);
    return ERROR_PRIMES;
  }

  free(primes);

  return 0;
}
