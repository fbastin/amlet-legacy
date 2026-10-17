/**
 * alternative.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2004-2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the MIT License
 */

#ifndef __AMLET_ALTERNATIVE_H__
#define __AMLET_ALTERNATIVE_H__

#include "amlet_driver.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct _amlet_alt AMLET_Alt;

  struct _amlet_alt {
    int id;

    char *name;

    double *prop;
    int nprop;
  };

  int amlet_alternatives_sampling(AMLET_Driver *dr, double ***obs,
				  int **availability, int *choices,
				  int np, int nobs, int nalt, int nvar,
				  int *rep, int n);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_ALTERNATIVE_H__ */
