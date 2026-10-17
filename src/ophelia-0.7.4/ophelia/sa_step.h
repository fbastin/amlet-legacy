/**
 * sa_step.h
 *
 * ATHOS Library
 * Copyright (C) 1999-2004 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ATHOS_SA_STEP_H__
#define __ATHOS_SA_STEP_H__

#include "btr_struct.h"

#include <sys/types.h>

#include <oratio/ran.h>

#ifdef __cplusplus
extern "C" {
#endif

  void nlp_sa_step_bounds(Random *ran, double *x, int n, void *data,
			  double step_size);

  void nlp_sa_step_default(Random *ran, double *x, int n, void *data,
			   double step_size);

#ifdef __cplusplus
}
#endif

#endif  /* __ATHOS_SA_STEP_H__ */
