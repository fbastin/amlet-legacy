/**
 * constraints.h
 *
 * OPHELIA Library
 * Copyright (C) 2006-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __NLP_CONSTRAINTS_H__
#define __NLP_CONSTRAINTS_H__

#ifdef __cplusplus
extern "C" {
#endif

  enum {
    NLP_NONE,
    NLP_BOUNDS,
    NLP_MONOTONICITY
  };

  typedef struct _NLP_Constraint NLP_Constraint;

  struct _NLP_Constraint {
    int type;   /* constraint type; see previous enumeration */
    int nidx;   /* number of constraints sets. */
    int *cidx;  /* arrays of constraints sets; each set i is
		   identified as starting at cidx[2*i] (inclusive) and
		   ending at cidx[2*i+1] (exclusive). */
    int id;     /* identifier */
  };

  NLP_Constraint *nlp_constraint_new(int type, int nidx, int *cidx, int id);

  void nlp_constraint_free(NLP_Constraint *c);

#ifdef __cplusplus
}
#endif

#endif  /* __NLP_CONSTRAINTS_H__ */
