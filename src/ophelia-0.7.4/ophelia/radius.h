/**
 * radius.h
 *
 * ATHOS Library
 * Copyright (C) 1999-2006 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ATHOS_BTR_RADIUS_H__
#define __ATHOS_BTR_RADIUS_H__

#include "btr_struct.h"

#ifdef __cplusplus
extern "C" {
#endif


  /**
   * Accept or reject the iterate candidate.
   *
   * @param b is an initialized BTR structure.
   * @param rho is the ratio between actuel and predicted reductions.
   *
   * @return 1 if the iterate candidate is accepted, 0 else.
   */
  int nlp_btr_accept_candidate(BTR *b, double rho);

  /**
   * Radius initialization procedures
   */

  /**
   * Compute the initial trust-region radius on basis of the gradient value.
   *
   * Ref.: Lin and More, "Newton's Method For Large Bound-Constrained
   * Optimization Problems", Siam Journal on Optimization 9(4), 1999.
   *
   * @param b is an initialized BTR structure.
   */
  double btr_radius_init_grad(BTR *b, void *data);

  /**
   * Radius update procedures
   */

  void btr_radius_update_basic(BTR *b, double rho);

  void btr_radius_update_usual(BTR *b, double rho);

#ifdef __cplusplus
}
#endif

#endif  /* __ATHOS_RADIUS_H__ */
