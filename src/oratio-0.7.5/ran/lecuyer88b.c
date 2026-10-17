/**
 * lecuyer88b.c
 *
 * Random Numbers Library
 * Copyright (C) 1999-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/ran.h"
#include "oratio/memory.h"
#include <lecuyer88b.h>

#include <string.h>

/** TAKE THE REFERENCE FROM NR **/

/**
 * L'Ecuyer algorithm constants
 */
#define IM1  2147483563  /* near INT_MAX */
#define IM2  2147483399
#define IMM1 2147483562
#define IA1  40014
#define IA2  40692
#define IQ1  53668
#define IQ2  52774
#define IR1  12211
#define IR2  3791
#define ndiv (1+(int) floor(IMM1/32))

/**
 * L'ECUYER ALGORITHM IMPLEMENTATION
 */

/**
 * Return a pseudo U(0,1)-random number following the L'Ecuyer's
 * 1988 algorithm and a postponed output technique.
 *
 * @param r is a point to a LEcuyer88b structure. It is declared as
 *   null to conform to the common definition of rng_get_val.
 */
double ran_lecuyer88b_get_val(void *r)
{
  int k;
  int i;

  unsigned long s1, s2;
  int iy;

  LEcuyer88b *ecu = (LEcuyer88b *)r;

  s1 = ecu->s[0];
  s2 = ecu->s[1];
  iy = ecu->iy;

  /* Génération d'un nombre sur base d'un algorithme construit sur la
   * méthode congruentielle, et utilisant l'algorithme de Schrage */
  k = s1/IQ1;
  s1 = IA1*(ecu->s[0]-k*IQ1)-k*IR1;
  /* S1=40014*(S1-IQ1*Floor(S1/IQ1))-Floor(S1/IQ1)*12211 */
  if (s1<0) s1 += IM1;

  /* Génération de la valeur de S2 */
  k = s2/IQ2;
  s2 = IA2*(s2-k*IQ2)-k*IR2;
  if (s2<0) s2 += IM2;

  /* Choice of the element to take out. */
  i = iy/ndiv;
  iy = ecu->iv[i] - s2;
  if (iy<1) iy += IMM1;
  /* insertion du nouvel élément */
  ecu->iv[i] = s1;

  ecu->s[0] = s1;
  ecu->s[1] = s2;
  ecu->iy = iy;

  return (double)iy/IM1;
}

/**
 */
unsigned long *ran_lecuyer88b_get_state(void *rng)
{
  return ((LEcuyer88b *)rng)->s;
}

/**
 * It is currently assumed that state represents only 1 long integer,
 * while it should represent two!.
 */
void ran_lecuyer88b_set_state(LEcuyer88b *rng, unsigned long *state)
{
  int i;

  memcpy(rng->s, state, 2*sizeof(double));

  /* Generation of an array of 32 numbers in order to allow postponed output */
  for (i = ntab+7; i >= 0; i--) {
    int k;

    k = rng->s[0]/IQ1;
    rng->s[0] = IA1*(rng->s[0]-k*IQ1)-k*IR1;
    if (rng->s[0] < 0) rng->s[0] +=IM1;
    if (i < ntab) rng->iv[i] = rng->s[0];
  }

  rng->iy = rng->iv[0];
}

void ran_lecuyer88b_init(LEcuyer88b *rng, unsigned long *seed, int n)
{
  if (n == 1)
    rng->s[0] = rng->s[1] = *seed;
  else /* n == 2 */
    ran_lecuyer88b_set_state(rng, seed);
}

LEcuyer88b *ran_lecuyer88b_new(unsigned long *seed, int n)
{
  LEcuyer88b *rng;

  if (n < 1 || n > 2) return NULL;

  rng = _malloc(LEcuyer88b);

  if (seed) {
    int i;

    for (i = 0; i < n; i++) if (seed[i] == 0) seed[i] = ran_init_random();
    ran_lecuyer88b_init(rng, seed, n);
  }
  else {
    unsigned long s;

    s = ran_init_random();
    ran_lecuyer88b_init(rng, &s, 1);
  }

  return rng;
}

void ran_lecuyer88b_free(LEcuyer88b *rng)
{
  if (rng) free((LEcuyer88b *)rng);
}
