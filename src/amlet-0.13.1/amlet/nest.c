/**
 * nest.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <nest.h>

#include <oratio/ctools.h>

nest *amlet_nest_new(int n, int nalt)
{
  nest *ne;
  int i;

  ne = _malloc(nest);
  ne->n = n;
  ne->alt = _mallocs(nalt, int);
  for (i = 0; i < nalt; i++) ne->alt[i] = 0;
  ne->mult = _mallocs(n, double);
  for (i = 0; i < n; i++) ne->mult[i] = 1.0;

  return ne;
}

void amlet_nest_free(nest *n)
{
  if (n) {
    free(n->alt);
    free(n->mult);
    free(n);
  }
}
