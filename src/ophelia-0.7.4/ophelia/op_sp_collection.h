/**
 * op_sp_collection.h
 *
 * OPHELIA Library
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <ophelia/problem.h>
#include <oratio/ran.h>

#ifndef __OP_SP_COLLECTION_H__
#define __OP_SP_COLLECTION_H__

#ifdef __cplusplus
extern "C" {
#endif

  OP_Problem *op_probl_sp_quad1_new(Random *r, int n);

#ifdef __cplusplus
}
#endif

#endif  /* __OP_SP_COLLECTION_H__ */
