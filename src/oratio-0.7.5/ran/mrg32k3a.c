/**
 * mrg32k3a.c
 *
 * ORATIO Library
 * Copyright (C) 2007-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/ran.h"

#include <oratio/ntools.h>

#include "mrg32k3a.h"

#include <string.h>

#define norm    2.328306549295727688e-10
#define m1      4294967087.0
#define m2      4294944443.0
#define a12     1403580.0
#define a13n    810728.0
#define a21     527612.0
#define a23n    1370589.0

#define fact    5.9604644775390625e-8    /* 1 / 2^24 */

/* Default initial seed of the package. Will be updated to become
   the seed of the next created stream. */
static double nextSeed[6] = { 12345, 12345, 12345, 12345, 12345, 12345 };

/* The following are the transition matrices of the two MRG components */
/* (in matrix form), raised to the powers -1, 1, 2^76, and 2^127, resp.*/
static double InvA1[9] = {          /* Inverse of A1p0 */
  184888585.0,   0.0,  1945170933.0,
          1.0,   0.0,           0.0,
          0.0,   1.0,           0.0
};

static double InvA2[9] = {          /* Inverse of A2p0 */
  0.0,  360363334.0,  4225571728.0,
  1.0,          0.0,           0.0,
  0.0,          1.0,           0.0
};

static double A1p0[9] = {
        0.0,        1.0,       0.0 ,
        0.0,        0.0,       1.0 ,
  -810728.0,  1403580.0,       0.0
};

static double A2p0[9] = {
         0.0,        1.0,       0.0,
         0.0,        0.0,       1.0,
  -1370589.0,        0.0,  527612.0
};

static double A1p76[9] = {
    82758667.0, 1871391091.0, 4127413238.0, 
  3672831523.0,   69195019.0, 1871391091.0, 
  3672091415.0, 3528743235.0,   69195019.0
};

static double A2p76[9] = {
  1511326704.0, 3759209742.0, 1610795712.0,
  4292754251.0, 1511326704.0, 3889917532.0,
  3859662829.0, 4292754251.0, 3708466080.0
};

static double A1p127[9] = {
  2427906178.0, 3580155704.0,  949770784.0,
   226153695.0, 1230515664.0, 3580155704.0,
  1988835001.0,  986791581.0, 1230515664.0
};

static double A2p127[9] = {
  1464411153.0,  277697599.0, 1610723613.0,
    32183930.0, 1464411153.0, 1022607788.0,
  2824425944.0,   32183930.0, 2093834863.0
};

void ran_mrg32k3a_resetstartstream (MRG32K3a *mrg)
{
  int i;

  for (i = 0; i < 6; ++i) mrg->s[i] = mrg->bs[i] = mrg->is[i];
}

void ran_mrg32k3a_resetnextsubstream(MRG32K3a *mrg)
{
  int i;

  nt_modular_matvec(3, A1p76, 3, mrg->bs, mrg->bs, m1);
  nt_modular_matvec(3, A2p76, 3, mrg->bs+3, mrg->bs+3, m2);

  for (i = 0; i < 6; ++i) mrg->s[i] = mrg->bs[i];
}

void ran_mrg32k3a_resetstartsubstream(MRG32K3a *mrg)
{
  int i;

  for (i = 0; i < 6; ++i) mrg->s[i] = mrg->bs[i];
}

void ran_mrg32k3a_advancestate(MRG32K3a *mrg, long e, long c)
{
  double B1[9], C1[9], B2[9], C2[9];

  if (e > 0) {
    nt_modular_mattwopow(3, A1p0, 3, B1, 3, m1, e);
    nt_modular_mattwopow(3, A2p0, 3, B2, 3, m2, e);
  }
  else if (e < 0) {
    nt_modular_mattwopow(3, InvA1, 3, B1, 3, m1, -e);
    nt_modular_mattwopow(3, InvA2, 3, B2, 3, m2, -e);
  }

  if (c >= 0) {
    nt_modular_matpow(3, A1p0, 3, C1, 3, m1, c);
    nt_modular_matpow(3, A2p0, 3, C2, 3, m2, c);
  }
  else {
    nt_modular_matpow(3, InvA1, 3, C1, 3, m1, -c);
    nt_modular_matpow(3, InvA2, 3, C2, 3, m2, -c);
  }

  if (e) {
    nt_modular_matmat(3, B1, 3, C1, 3, C1, 3, m1);
    nt_modular_matmat(3, B2, 3, C2, 3, C2, 3, m2);
  }

  nt_modular_matvec(3, C1, 3, mrg->s, mrg->s, m1);
  nt_modular_matvec(3, C2, 3, mrg->s+3, mrg->s+3, m2);
}

/**
 * Draw a pseudo U(0,1)-random number following the MRG32KA algorithm.
 *
 * @param mrg is an initialized MRG32K3a random numbers generator.
 */
double ran_mrg32k3a_get_val(void *mrg)
{
  long   k;
  double p1, p2;

  MRG32K3a *rng = (MRG32K3a *)mrg;

  /* Component 1 */
  p1 = a12*rng->s[2] - a13n*rng->s[0];
  k = p1/m1;   p1 -= k*m1;
  if (p1 < 0.0) p1 += m1;
  /* Update the state. */
  rng->s[0] = rng->s[1];
  rng->s[1] = rng->s[2];
  rng->s[2] = p1;

  /* Component 2 */
  p2 = a21*rng->s[5] - a23n*rng->s[3];
  k  = p2/m2;  p2 -= k*m2;
  if (p2 < 0.0) p2 += m2;
  /* Update the state. */
  rng->s[3] = rng->s[4];
  rng->s[4] = rng->s[5];
  rng->s[5] = p2;

  /* Combination */
  return (p1 <= p2) ? ((p1 - p2 + m1) * norm) : ((p1 - p2) * norm);
}

/**
 * Get the state of the random number generator.
 *
 * @param mrg is an initialized MRG32K3a random number generator.
 *
 * @return the current of the generator mrg.
 */
double *ran_mrg32k3a_get_state(MRG32K3a *mrg)
{
  return mrg->s;
}

/**
 * This also reinitializes the stream.
 */
void ran_mrg32k3a_set_state(MRG32K3a *mrg, void *state)
{
  memcpy(mrg->s, state, 6*sizeof(double));
  memcpy(mrg->bs, state, 6*sizeof(double));
  memcpy(mrg->is, state, 6*sizeof(double));

  nt_modular_matvec(3, A1p127, 3, nextSeed, nextSeed, m1);
  nt_modular_matvec(3, A2p127, 3, nextSeed+3, nextSeed+3, m2);
}

void ran_mrg32k3a_init(MRG32K3a *mrg, unsigned long *seed, int n)
{
  double s[6];
  int i;

  if (n <= 1)
    for (i = 0; i < 6; i++) s[i] = *seed;
  else
    /* n == 6 */
    for (i = 0; i < 6; i++) s[i] = seed[i];

  ran_mrg32k3a_set_state(mrg, s);
}

/**
 * Initialize the random generator using a default, pre-implemented,
 * initial state.
 */
void ran_mrg32k3a_init_default(MRG32K3a *mrg)
{
  ran_mrg32k3a_set_state(mrg, nextSeed);
}

/**
 * Create a new mrg based on the MRG32K3a method.
 *
 * @param seed is the seed of the generater. It can be set to NULL, 
 * @param n is the size of the vector pointed by seed. n must be equal
 *   to 0, 1 or 6. If n = 0, a pre-implemented initial state is used.
 *
 * @return a pointer to an initialized MRG32K3a structure if n is
 *   equal to 0, 1 or 6, NULL otherwise.
 */
MRG32K3a *ran_mrg32k3a_new(unsigned long *seed, int n)
{
  MRG32K3a *mrg;

  if (!(n == 0 || n == 1 || n == 6)) return NULL;

  mrg = _malloc(MRG32K3a);

  if (!seed)
    ran_mrg32k3a_init_default(mrg);
  else if (n != 1 && n != 6)
    ran_mrg32k3a_init_default(mrg);
  else {
    if (*seed == 0 && n <= 1) *seed = ran_init_random();
    else if (n == 6) {
      if (!(seed[0] || seed[1] || seed[2])) seed[2] = ran_init_random();
      if (!(seed[3] || seed[4] || seed[5])) seed[5] = ran_init_random();
    }
    ran_mrg32k3a_init(mrg, seed, n);
  }
  return mrg;
}

/**
 * Free memory associated to the MRG32K3a random numbers generator.
 *
 * @param mrg must be an initialized pointer to a MRG32K3a structure.
 */
void ran_mrg32k3a_free(void *mrg)
{
  if (mrg) free((MRG32K3a *)mrg);
}
