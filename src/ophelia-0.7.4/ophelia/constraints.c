/**
 * constraints.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "constraints.h"

#include <oratio/ctools.h>

NLP_Constraint *nlp_constraint_new(int type, int nidx,
					  int *cidx, int id)
{
  NLP_Constraint *c = _malloc(NLP_Constraint);

  c->type = type;
  c->nidx = nidx;
  c->cidx = cidx;
  c->id = id;

  return c;
}

void nlp_constraint_free(NLP_Constraint *c)
{
  if (c) {
    free(c->cidx);
    free(c);
  }
}

