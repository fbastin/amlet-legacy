/**
 * amlet_optimization.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2004 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __AMLET_OPTIMIZATION_H__
#define __AMLET_OPTIMIZATION_H__

#include <ophelia/nlp.h>

#include <amlet_driver.h>

#ifdef __cplusplus
extern "C" {
#endif

#define AMLET_DEFAULT_TRUST_REGION_RADIUS 1.0

  void amlet_init_btr(BTR *b, AMLET_Driver *dr);

  void amlet_init_bfgs(BFGS *bf, AMLET_Driver *dr);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_OPTIMIZATION_H__ */
