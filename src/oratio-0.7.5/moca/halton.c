/**
 * halton.c
 *
 * ORATIO Library - Monte Carlo and Quasi-Monte Carlo module
 * Copyright (C) 2009-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/moca.h"
#include "oratio/oratio.h"

#define HALTON_BW_MAXDIM       16
#define HALTON_BW_PRIME_MAXDIM 53

/**
 * Braaten Weller permutations.
 *
 * These permutations are quite old, but still popular in transportation
 * research.
 */
const int braaten_weller_permutations
  [HALTON_BW_MAXDIM][HALTON_BW_PRIME_MAXDIM] = {
  {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0},
  {0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0},
  {0, 3, 1, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0},
  {0, 4, 2, 6, 1, 5, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0},
  {0, 5, 8, 2, 10, 3, 6, 1, 9, 7, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0},
  {0, 6, 10, 2, 8, 4, 12, 1, 9, 5, 11, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0},
  {0, 8, 13, 3, 11, 5, 16, 1, 10, 7, 14, 4, 12, 2, 15, 6, 9, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0},
  {0, 9, 14, 3, 17, 6, 11, 1, 15, 7, 12, 4, 18, 8, 2, 16, 10, 5, 13, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0},
  {0, 11, 17, 4, 20, 7, 13, 2, 22, 9, 15, 5, 18, 1, 14, 10, 21, 6, 16, 3,
   19, 8, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 15, 7, 24, 11, 20, 2, 27, 9, 18, 4, 22, 13, 26, 5, 16, 10, 23, 1, 19,
   28, 6, 14, 17, 3, 25, 12, 8, 21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 15, 23, 5, 27, 9, 18, 2, 29, 12, 20, 7, 25, 11, 17, 3, 30, 14, 22,
   1, 21, 8, 26, 10, 16, 28, 4, 19, 6, 24, 13, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 18, 28, 6, 23, 11, 34, 3, 25, 14, 31, 8, 20, 36, 1, 16, 27, 10, 22,
   13, 32, 4, 29, 17, 7, 35, 19, 2, 26, 12, 30, 9, 24, 15, 33, 5, 21, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 20, 31, 7, 26, 12, 38, 3, 23, 34, 14, 17, 29, 5, 40, 10, 24, 1, 35,
   18, 28, 9, 33, 15, 21, 4, 37, 13, 30, 8, 39, 22, 2, 27, 16, 32, 11, 25,
   6, 36, 19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 21, 32, 7, 38, 13, 25, 3, 35, 17, 28, 10, 41, 5, 23, 30, 15, 37, 1,
   19, 33, 11, 26, 42, 8, 18, 29, 4, 39, 14, 22, 34, 6, 24, 12, 40, 2, 31,
   20, 27, 9, 36, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 24, 12, 39, 6, 33, 20, 44, 3, 29, 16, 36, 10, 42, 22, 8, 31, 26, 14,
   46, 1, 35, 18, 28, 5, 40, 19, 37, 11, 25, 43, 4, 30, 15, 34, 9, 45, 21,
   2, 32, 17, 41, 13, 27, 7, 38, 23, 0, 0, 0, 0, 0, 0},
  {0, 26, 40, 9, 33, 16, 49, 4, 36, 21, 45, 12, 29, 6, 51, 23, 38, 14, 43,
   1, 30, 19, 47, 10, 34, 24, 42, 3, 27, 52, 15, 18, 39, 7, 46, 31, 11, 35,
   20, 48, 2, 28, 41, 8, 22, 50, 13, 32, 17, 44, 5, 37, 25}
};

int moca_halton_std_permutations_max_dim_get()
{
  return HALTON_BW_MAXDIM;
}

int moca_halton_std_permutations_prime_max_dim_get()
{
  return HALTON_BW_PRIME_MAXDIM;
}

int *moca_halton_std_permutations()
{
  return (int *)braaten_weller_permutations;
}

/**
 * Reset a Halton draws generator.
 *
 * @param h is an initialized Halton generator.
 */
void moca_halton_reset(MOCA_Halton *h)
{
  h->in = 0;
}

/**
 * @param prime is the prime used to construct the Halton draw.
 * @param n is the integer to transform in a Halton draw.
 */
double moca_halton(int prime, int s, int n, int *scrambling, int lds)
{
  int digit;
  double half, h;

  h = 0.0;
  half = 1.0/prime;

  while (n > 0) {
    digit = (n % prime);
    if (scrambling)
      /* do the permutation */
      digit = scrambling[s*lds+digit];
    h += digit*half;
    n /= prime;
    half /= prime;
  }

  return h;
}

/**
 * Generate a Halton draws sample.
 *
 * @param h is a Halton draws generator.
 * @param n is the numbers of vectors forming the sample.
 * @param S contains on output the halton sample.
 * @param ldS is the leading dimension of S; ldS must be grater or equal to
 *   sample dimension, and is typically equal to this dimesion.
 */
void moca_halton_sample(MOCA_Halton *h,
			int n, double *S, int ldS)
{
  int i;
  
  moca_halton_reset(h);
  for (i = 0; i < n; i++) moca_halton_next(h, S+i*ldS);
}

/**
 * Generate a sample based on Halton sequences
 *
 * @param ran is random number generator.
 * @param loose is equal to 1 if for some errors can be managed by
 *   only raising a warning, 0 else.
 */
void moca_halton_generate_sample(Random *ran, void *sample,
                                 unsigned int s, int n,
                                 int randomization, int loose)
{
  int i, j, k;
  unsigned int *primes;   /* s first prime numbers */
  /* array of draws, for one individual */
  double **draws;
  double *shifts;

  int *scrambling = NULL;

  int scrambled;

  primes = _mallocs(s, unsigned int);
  shifts = _mallocs(s, double);
  draws = nt_matrix_new(s, n);

  /* We generate the primes associated to the s dimensions. */
  nt_primes_eratosthene(s, primes);

  if (randomization == MOCA_HALTON_RANDOMIZATION_SHIFTING ||
      randomization == MOCA_HALTON_RANDOMIZATION_SHIFTING_SCRAMBLING ||
      randomization == MOCA_HALTON_RANDOMIZATION_SHIFTING_SHUFFLING)
    for (i = 0; i < s; i++) shifts[i] = ran_random_get_val(ran);

  scrambled = (randomization == MOCA_HALTON_RANDOMIZATION_SCRAMBLING ||
	       randomization == MOCA_HALTON_RANDOMIZATION_SHIFTING_SCRAMBLING)
    ? 1 : 0;
  if (scrambled) {
    if (s > HALTON_BW_MAXDIM) {
      if (loose)
	c_warning("scrambled Halton sequences are only available up to dimension 16; no scrambling will be made");
      else
	c_error("scrambled Halton sequences only available up to dimension 16; stop");
    }
    scrambling = moca_halton_std_permutations();
  }

  /* For each dimension, generate the corresponding Halton sequence */
  for (j = 0; j < s; j++) {
    for (k = 0; k < n; k++) {
      draws[j][k] = moca_halton(primes[j], j, k+1,
				scrambling, HALTON_BW_PRIME_MAXDIM);
      if (randomization == MOCA_HALTON_RANDOMIZATION_SHIFTING ||
	  randomization == MOCA_HALTON_RANDOMIZATION_SHIFTING_SCRAMBLING ||
	  randomization == MOCA_HALTON_RANDOMIZATION_SHIFTING_SHUFFLING)
	draws[j][k] = fmod(draws[j][k]+shifts[j],1);
    }
    if (randomization == MOCA_HALTON_RANDOMIZATION_SHUFFLING ||
	randomization == MOCA_HALTON_RANDOMIZATION_SHIFTING_SHUFFLING)
      ran_shuffle_double(ran, n, draws[j], 1);
    if (sample) {
      for (k = 0; k < n; k++)
	  ((double **)sample)[k][j] = draws[j][k];
    }
  }

  free(primes);
  free(shifts);
  nt_matrix_free(draws);
}

/**
 * Generate the next point of a Halton sequence.
 *
 * @param h is a MOCA_Halton initialized generator.
 * @param x is a double vector that contains on output the next point
 *   of the Halton sequence.
 */
void moca_halton_next(MOCA_Halton *h, double *x)
{
  int i;

  h->in++;
  for (i = 0; i < h->s; i++)
    x[i] = moca_halton(h->primes[i], i, h->in,
		       h->scrambling, h->primes[h->lds-1]);
}

/**
 * @param h is an initialized Halton structure.
 */
void moca_halton_set_scrambling(MOCA_Halton *h, int sid)
{
  switch (sid) {
  case MOCA_HALTON_NO_SCRAMBLING:
    h->maxdim = -1;
    h->scrambling = NULL;
    break;
  case MOCA_HALTON_BW_SCRAMBLING:
  default:
    h->maxdim = HALTON_BW_MAXDIM;
    h->scrambling = (int *)braaten_weller_permutations;
    break;
  }
}

/**
 * Create a new Halton sequence generator.
 *
 * @param s is the dimension of the vectors to generate.
 * @param sid an identifier of the scrambling method to use, 0 meaning
 *   that no scrambling has to be performed.
 *
 * @return a newly allocate MOCA_Halton structure.
 */
MOCA_Halton *moca_halton_new(unsigned int s, int sid)
{
  MOCA_Halton *h;

  h = _malloc(MOCA_Halton);

  h->s = s;
  h->in = 0;

  moca_halton_set_scrambling(h, sid);
  h->lds = s;

  h->primes = _mallocs(h->lds, unsigned int);
  nt_primes_eratosthene(h->lds, h->primes);

  return h;
}

/**
 * Free the memory associated to a MOCA_Halton structure.
 *
 * @param h is an initialized MOCA_Halton structure.
 */
void moca_halton_free(MOCA_Halton *h)
{
  if (h) {
    free(h->primes);
    free(h);
  }
}
