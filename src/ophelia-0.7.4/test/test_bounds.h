/**
 * test_bounds.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2004 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __NLP2_TESTING_BOUNDS_H__
#define __NLP2_TESTING_BOUNDS_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * PBR_T1_1
   */

  double pbr_t1_1(double *x, int n, void *data);
  void pbr_t1_1_set_bounds(double *b);

  /**
   * PBR_T1_3
   */

  double pbr_t1_3(double *x, int n, void *data);
  void pbr_t1_3_set_bounds(double *b);

#ifdef __cplusplus
}
#endif

#endif  /* __NLP2_TESTING_BOUNDS_H__ */
