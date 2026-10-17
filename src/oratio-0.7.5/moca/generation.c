/**
 * generation.c
 *
 * ORATIO Library
 * Copyright (C) 2005-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/moca.h"
#include "oratio/values.h"

#include <stdio.h>

/**
 * Return a pointer to the set of moca methods symbols.
 *
 * @return a pointer to the methods symbols array.
 */
c_symbol *moca_get_methods()
{
  return (c_symbol *)moca_methods;
}

/**
 * Give the number of implemented methods.
 *
 * @return the number of implemented methods.
 */
int moca_methods_number()
{
  return MOCA_METHODS_END;
}

/**
 * Return the name of the sampling method.
 *
 * @param method is a sampling method identifier.
 *
 * @param method the name of the sampling method associated to the
 *   given identifier, or NULL if this identifier is invalid.
 */
char *moca_get_method_name(int method)
{
  return c_symbol_get_name(method, moca_get_methods(), moca_methods_number());
}

void moca_convert_vector(int s, void *x, int *type)
{
  int i;

  for (i = 0; i < s; i++)
    ((double *)x)[i] = moca_convert_random_value(type[i], ((double *)x)[i]);
}

/**
 * Generate a Monte Carlo sample.
 *
 * @param ran is an initialized random numbers generator.
 * @param n is the number of vectors to generate.
 * @param s is the dimension of each vector.
 * @param S is a n by s matrix that will contain on output the Monte
 *   Carlo sample.
 */
void moca_mc_generate_sample(Random *ran, int n, int s,
				    double **S)
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < s; j++)
      S[i][j] = ran_random_get_val(ran);
}
