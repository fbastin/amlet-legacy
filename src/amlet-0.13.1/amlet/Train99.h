/**
 * Train99.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2004 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __AMLET_TRAIN_1999_H__
#define __AMLET_TRAIN_1999_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Read observations written using the Train's format (1999).
   *
   * @param times: times[i] (i > 0) represents the number of observations for
   *   individual i-1.
   */
  int import_train_1999(const char *obsname, double ***obs,
			int **availability, int *choices,
			int np, int nobs, int nalt, int nvar,
			int idcensor, int iddep, int *times);

  /**
   * Save observations written using the Train's format (1999).
   */
  int export_train_1999(const char *obsname, double ***obs,
			int **availability, const int *choices,
			int np, int nobs, int nalt, int nvar,
			int idcensor, int *times);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_TRAIN_1999_H__ */
