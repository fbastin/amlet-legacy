/**
 * utilities.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 1999-2006 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the MIT License
 */

int nlp_feasible_default(double *x, int n, void *data)
{
  return 1;
}
