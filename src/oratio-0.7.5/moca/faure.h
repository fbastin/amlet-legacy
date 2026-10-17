/**
 * faure.h
 *
 * Monte Carlo and Quasi-Monte Carlo Library
 * Copyright (C) 2005-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __MOCA_FAURE_H__
#define __MOCA_FAURE_H__

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
    int s;                /* dimension */
    int prime;            /* smallest prime greater or equal to s */
    double invprime;      /* 1.0/prime */
    int testn, nextn, hisum;
    int coef[400];        /* table of binomial coefficients */
  } MOCAFaure;

  /**
   * Initialize the Faure quasi-random generator.
   *
   * @param s is the vectors dimension.
   *
   * @return MOCA_ERROR_DIMENSION is the dimension is not strictly
   *   between 1 and 41 and MOCA_ERROR_DIGITS if too many digits are
   *   needed to represent the numbers.
   */
  int moca_faure_init(MOCAFaure *mcf, unsigned int s, int atmost);

  MOCAFaure *moca_faure_new();

  void moca_faure_free(MOCAFaure *mcf);

  int moca_faure_next(MOCAFaure *mcf, double *x);

#ifdef __cplusplus
}
#endif

#endif  /* __MOCA_FAURE_H__ */
