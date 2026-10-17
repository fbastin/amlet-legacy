/**
 * amlet_stopping.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2005 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the MIT License
 */

#ifndef __AMLET_STOPPING_H__
#define __AMLET_STOPPING_H__

#include <amlet/amlet.h>

#include <ophelia/nlp.h>
#include <amlet/dcl.h>

#ifdef __cplusplus
extern "C" {
#endif

  int amlet_bfgs_stopping_test(BFGS *bf, int info);

  int amlet_btr_stopping_test(BTR *btr, double rho);

  int amlet_btrda_stopping_test(BTR *btr, double *step, int nlast,
				dcl_data *md, double accuracy,
				double rfct, double rho);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_STOPPING_H__ */
