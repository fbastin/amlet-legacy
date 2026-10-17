/**
 * projection.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 1999-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __NLP_PROJECTION_H__
#define __NLP_PROJECTION_H__

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct _NlpProj {
    void *data;   /* date used for the projection */
    int id;       /* constraints identifiers */
  } NlpProj;

  /**
   * Projection function type.
   *
   * @param n is the problem dimension.
   * @param x is the point to project.
   * @param w is a weighting vector.
   * @param proj is, on output, the projected vector.
   * @param data represent problem-specific data.
   *
   * @return an error code
   */
  typedef int (*NlpProjection)(int n, double *x, double *w,
			       double *proj, void *data);

  int nlp_proj_identity(int n, double *x, double *w,
			double *proj, void *data);

#ifdef __cplusplus
}
#endif

#endif  /* __NLP_PROJECTION_H__ */
