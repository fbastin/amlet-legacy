/**
 * amlet_data.h
 *
 * AMLET Library
 * Copyright (C) 2002-2015 Fabian Bastin <fabian.bastin@gmail.com>
 */

#ifndef __AMLET_DATA_H__
#define __AMLET_DATA_H__

/**
 * Data import and export
 */

#include <oratio/ctools.h>

enum {
  TRAIN96 = 1,
  TRAIN99,
  AMLET,
  END_DATA_FORMAT
};

static const c_symbol formats[] = {
  {"Train96", TRAIN96},
  {"Train99", TRAIN99},
  {"amlet", AMLET}
};

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Check that expressed choices correspond to available alternatives
   *
   * @return 0 if no error has been found, AMLET_ERROR_CHOICES otherwise.
   */
  int amlet_error_choice(FILE *f, int nobs, int *choices, int **availability);

  /**
   * Check that the number of observations is coherent, i.e. that the sum of the
   * components of the repetition vectors is equal to the number of
   * observations.
   *
   * @param np is number of individuals.
   * @param nobs is the number of observations.
   * @param rep is the vector of repetitions; rep[i] should give the number of
   *   observations delivered by individual i.
   *
   * @return 0 is no error has been found, AMLET_ERROR_REPETITIONS otherwise.
   */
  int amlet_check_repetitions(int np, int nobs, int *rep);

  /**
   * Read observations from an external file.
   *
   * @param format is a valid data file format identifier.
   */
  int import(int format, char *obsname, double ***obs,
	     int **availability, int *choices,
	     int np, int nind, int nalt, int nobs,
	     int idcensor, int iddep, int *times);

  /**
   * Export observation in an external file.
   *
   * @param format is a valid data file format identifier.
   */
  int export(int format, char *obsname, double ***obs,
	     int **availability, const int *choices,
	     int np, int nind, int nalt, int nobs,
	     int idcensor, int *times);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_DATA_H__ */
