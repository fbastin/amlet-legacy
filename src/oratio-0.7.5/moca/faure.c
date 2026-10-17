/**
 * faure.c
 *
 * Monte Carlo and Quasi-Monte Carlo Library
 * Copyright (C) 2005 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/moca.h"
#include "error.h"

#include <math.h>

#include <oratio/oratio.h>

/**
 * Reference:
 *   Algorithm 647 from ACM.
 *   Algorithm appeared in ACM-Trans. Math. Software, Vol.12(4),
 *   p. 362, 1986.
 */

/**
 * Initialize the Faure quasi-random generator.
 *
 * @param s is the vectors dimension.
 *
 * @return MOCA_ERROR_DIMENSION is the dimension is not strictly
 *   between 1 and 41 and MOCA_ERROR_DIGITS if too many digits are
 *   needed to represent the numbers.
 */
int moca_faure_init(MOCAFaure *mcf, unsigned int s, int atmost)
{
  unsigned int *primes;

  double r;

  int i, j;

  /* Then it calculates an upper summation limit "hisum" based on
     "s" and the user-supplied number "atmost" of quasirandom
     vectors required. flag(2)=.true. if atmost is ok. */

  primes = _mallocs(s-1, unsigned int);

  /* Dimension check */
  mcf->s = s;
  if (s < 2 || s > 40) return MOCA_ERROR_DIMENSION;

  /* Compute the prime used for the sequence */
  nt_primes_eratosthene(s-1, primes);
  mcf->prime = primes[s-2];

  /* Computing 4th power */
  mcf->testn = mcf->prime;
  mcf->testn *= mcf->testn;
  mcf->testn *= mcf->testn;

  /* Compute log(atmost+testn) in base qs using a ratio of natural
     logs to get an upper bound on (the number of digits in the base
     qs representation of atmost+testn) minus one. */

  r = log((double)(atmost+mcf->testn))/log((double) mcf->prime);
  /* we round r to the nearest integer */
  mcf->hisum = (int)floor(r+0.5);
  if (mcf->hisum >= 20) return MOCA_ERROR_DIGITS;

  /* Now find binomial coefficients mod qs in a lower-triangular
     matrix "coef" using recursion
     binom(i,j)=binom(i-1,j)+binom(i-1,j-1) and a=b+c implies
     mod(a,d)= mod(mod(b,d)+mod(c,d),d) */

  mcf->coef[0] = 1;
  for (j = 1; j <= mcf->hisum; ++j) {
    mcf->coef[j] = 1;
    mcf->coef[j + j * 20] = 1;
  }
  for (j = 1; j <= mcf->hisum; ++j) {
    for (i = j + 1; i <= mcf->hisum; ++i) {
      mcf->coef[i + j * 20] = (mcf->coef[i - 1 + j * 20] + 
			       mcf->coef[i - 1 + (j - 1) * 20]) % mcf->prime;
    }
  }

  /* Calculating these coefficients mod qs avoids possible overflow
     problems with raw binomial coefficients */

  /* Now complete initialization as described in section 2.
     nextn has 4 digits in base qs, so hisum equals 3. */

  mcf->nextn = mcf->testn - 1;
  mcf->hisum = 3;
  mcf->invprime = 1./(double)mcf->prime;

  free(primes);

  return MOCA_OK;
}

/**
 * Generate a new Faure quasirandom vector.
 *
 * Ref.: H. Faure, Acta Arithmetica, Vol. XLI, pp. 337-351, 1982.
 *
 * The user must call moca_faure_init before calling moca_faure_next.
 *
 * @return a new quasirandom vector
 */
int moca_faure_next(MOCAFaure *mcf, double *x)
{
  int i, j, k;
  double r;
  int ktemp, ltemp, mtemp, ytemp[20], ztemp;

  ktemp = mcf->testn;
  ltemp = mcf->nextn;

  /* We now compute the ytemp(j)'s */
  for (i = mcf->hisum; i >= 0; i--) {
    ktemp /= mcf->prime;
    mtemp = ltemp % ktemp;
    ytemp[i] = (ltemp - mtemp) / ktemp;
    ltemp = mtemp;
  }

  /* x(k) has the form sum over j from zero to hisum of
     ytemp(j)*(qs**(-(j+1))) 
     Ready to compute x(1) using nested multiplication */

  /* Compute x[0] */
  r = (double) ytemp[mcf->hisum];
  for (i = mcf->hisum - 1; i >= 0; --i) {
    r = ytemp[i] + mcf->invprime * r;
  }
  x[0] = r * mcf->invprime;

  /* Find the other s-1 components of the quasi-random vector */
  for (k = 1; k < mcf->s; k++) {
    x[k] = 0.;
    r = mcf->invprime;
    for (j = 0; j <= mcf->hisum; ++j) {
      ztemp = 0;
      for (i = j; i <= mcf->hisum; ++i) {
	ztemp += mcf->coef[i + j * 20] * ytemp[i];
      }

      ytemp[j] = ztemp % mcf->prime;
      x[k] += ytemp[j] * r;
      r *= mcf->invprime;
    }
  }

  /* Update nextn and, if needed, testn and hisum */
  mcf->nextn++;
  if (mcf->nextn == mcf->testn) {
    mcf->testn *= mcf->prime;
    mcf->hisum++;
  }

  return 0;
}

MOCAFaure *moca_faure_new()
{
  MOCAFaure *mcf;

  mcf = _malloc(MOCAFaure);

  return mcf;
}

void moca_faure_free(MOCAFaure *mcf)
{
  if (mcf)
    free(mcf);
}
