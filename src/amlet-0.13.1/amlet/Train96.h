/**
 * Train96.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2004 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the MIT License
 */

#ifndef __AMLET_TRAIN_1996_H__
#define __AMLET_TRAIN_1996_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Train's format file importation (1996 Gauss code).
   *
   * @param idcensor must be set to 0 if there is no line giving the
   *   availibility of the alternatives
   */
  int import_train_1996(char *obsname, double ***obs,
			int **availability, int *choices,
			int np, int nobs, int nalt, int nvar,
			int idcensor, int iddep);

  /**
   * Save observations at the Gauss format.
   *
   * @param idcensor indexes the line where the alternatives accessibility
   *   is expressed.
   */
  int export_train_1996(char *obsname, double ***obs,
			int **availability, const int *choices,
			int np, int nobs, int nalt, int nvar,
			int idcensor);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_TRAIN_1996_H__ */
