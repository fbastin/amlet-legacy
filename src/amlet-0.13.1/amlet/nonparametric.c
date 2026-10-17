/**
 * nonparametric.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/oratio.h>

/**
 * Determine if the type is considered as nonparamtric.
 *
 * @param type is an integer identifier.
 *
 * @return 1 is the type is considered as nonparametric, 0 else.
 */
int amlet_is_nonparametric(int type)
{
  return (type < MOCA_END_DISTRIBUTIONS) ? 0 : 1;
}

void *amlet_nonparametric_set(int type)
{
  void *object = NULL;

  switch (type) {
  case RAN_SP_SPLINE: {
    /* The current gestion is too much static. */
    NT_BSpline *bs;
    /*    int n = 4;
	  double knots[] = {0.0, 1.0/3.0, 2.0/3.0, 1.0}; */
    /*    int n = 5;
	  double knots[] = {0.0, 0.25, 0.5, 0.75, 1.0}; */
    int n = 6;
    double knots[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};

    bs = nt_bsplines_new(3, n, knots, 1);

    return bs;
  }
  default:
    c_error("invalid nonparametric choice");
  }

  return object;
}

void amlet_nonparametric_free(void *np_object, int type)
{
  switch (type) {
  case RAN_SP_SPLINE:
    nt_bsplines_free((NT_BSpline *)np_object);
    break;
  default:
    free(np_object);
  }
}
