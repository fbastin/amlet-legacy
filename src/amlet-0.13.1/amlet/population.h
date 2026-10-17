/**
 * population.h
 *
 * AMLET Library
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __AMLET_POPULATION_H__
#define __AMLET_POPULATION_H__

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct _AMLET_Population {
    double ***obs;
    int *rep;
    int **availability;
    int *choices;
  } AMLET_Population;

#ifndef __cplusplus
}
#endif

#endif  /* __AMLET_POPULATION_H__ */
