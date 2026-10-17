/**
 * nonparametric.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2006 Fabian Bastin <bastin@cerfacs.fr>
 *
 * Licensed under the MIT License
 */

#ifndef __OPHELIA_NONPARAMETRIC_H__
#define __OPHELIA_NONPARAMETRIC_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Determine if the type is considered as nonparamtric.
   *
   * @param type is an integer identifier.
   *
   * @return 1 is the type is considered as nonparametric, 0 else.
   */
  int amlet_is_nonparametric(int type);

  void *amlet_nonparametric_set(int type);

  void amlet_nonparametric_free(void *np_object, int type);

#ifdef __cplusplus
}
#endif

#endif  /* __OPHELIA_NONPARAMETRIC_H__ */
