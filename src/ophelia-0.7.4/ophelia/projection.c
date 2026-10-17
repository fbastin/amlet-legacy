/**
 * projection.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 1999-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "projection.h"

#include "btr_struct.h"
#include "isotonic.h"

#include <oratio/ctools.h>
#include <oratio/ntools.h>

#include <float.h>
#include <string.h>

int nlp_proj_identity(int n, double *x, double *w,
		      double *proj, void *data)
{
  memcpy(proj, x, n*sizeof(double));

  return 0;
}

/**
 * Perform some prefined projection within the basic trust-region algorithm.
 *
 * @param b is an initialized BTR structure.
 * @param y is the point to project. y is replaced on output by its projection.
 *
 * @return 0 is no problem has been encountered, otherwise some error code.
 */
int nlp_btr_proj_standard(BTR *b, double *y)
{
  DS_List *l = b->constraints;
  NLP_Constraint *c;

  /* We loop on the constraint to apply classical projections. */
  while (l) {
    c = (NLP_Constraint *)(l->data);
    switch (c->type) {
    case NLP_BOUNDS:
      /* TODO! */
      break;
    case NLP_MONOTONICITY:
      nlp_isotonics_wrapper(b->n, y, NULL, y, c->nidx, c->cidx);
      break;
    default:
      /* no projection is performed. */
      break;
    }
    l = l->next;
  }

  return 0;
}
