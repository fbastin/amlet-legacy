/**
 * nest.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __AMLET_NEST_H__
#define __AMLET_NEST_H__

#include <amlet/data.h>
#include <amlet/dcl.h>

#ifdef __cplusplus
extern "C" {
#endif

  nest *amlet_nest_new(int n, int nalt);

  void amlet_nest_free(nest *n);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_NEST_H__ */
