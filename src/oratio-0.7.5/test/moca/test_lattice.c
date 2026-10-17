/**
 * test_lattice.c
 *
 * ORATIO Library
 * Copyright (C) 2005-2007 <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "caffeine.h"

/**
 * Lattice tests
 *
 * @param ran is an initialized random numbers generator.
 *
 * @return 0 if the test has been successful, an error code otherwise.
 */
int moca_lattice_check(Random *ran)
{
  FILE *f;

  int n;

  double *quasi;
  int dimen = 2;
  int a[] = {1,149,115,87,55};
  double shift[5];

  int i;

  MOCA_Rank1Lattice *l;

  n = 512;
  quasi = _mallocs(dimen, double);

  l = moca_rank1_lattice_new(n, dimen, a);
  moca_shift_init32(ran, 5, shift);

  f = fopen("LatticeDraws", "w");

  for (i = 0; i < n; i++) {
    moca_rank1_lattice_next(l, quasi, 1);
    moca_shift_apply(dimen, quasi, shift);
    moca_baker_point(dimen, quasi, 1);
    ran_write_vector(f, quasi, dimen);
  }

  if (f) fclose(f);

  free(quasi);
  moca_rank1_lattice_free(l);

  return 0;
}
