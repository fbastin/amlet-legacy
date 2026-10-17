/**
 * error.h
 *
 * NTools Library
 * Copyright (C) 2002 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_NTOOLS_ERROR_H__
#define __ORATIO_NTOOLS_ERROR_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /* Use also the routines in ctools */

  /* Définition de constantes d'erreur */
  enum {
    NT_EFAILED = 1,
    NT_EBADTOL,
    NT_EROUND,
    NT_EMAXITER,
    NT_ESING,
    NT_MATRIX_BAD_ORIENTATION,
    NT_MATRIX_BAD_DIM1,
    NT_MATRIX_BAD_DIM2,
    NT_VECTOR_BAD_DIM
  };

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_NTOOLS_ERROR_H__ */
