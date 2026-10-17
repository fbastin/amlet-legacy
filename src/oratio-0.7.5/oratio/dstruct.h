/**
 * dstruct.h
 *
 * Data Structures Library
 * Copyright (C) 2000-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_DSTRUCT_H__
#define __ORATIO_DSTRUCT_H__

#include <oratio/ctools.h>

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NULL
#define NULL 0
#endif

#include "ds_functions.h"
#include "dslist.h"
#include "dsslist.h"
#include "ds_sarray.h"
#include "ds_search.h"
#include "ds_sort.h"


#include "dsqueue.h"

#include "dsarray.h"

  /**
   * Tree structure definition
   */

  typedef struct _DS_Tree DS_Tree;

  struct _DS_Tree {
    void *root;
    DS_Tree *parent;
    DS_List *children;
  };

  DS_Tree *ds_tree_new();

  DS_Tree *ds_tree_new_with_root(void *obj);

  void ds_tree_free(DS_Tree *tree, Dfct ff);

  DS_Tree *ds_tree_set_root(DS_Tree *tree, void *obj);

  DS_Tree *ds_tree_insert(DS_Tree *parent, void *obj);

  void ds_tree_foreach(DS_Tree *tree, Dfct ff);

  DS_Tree *ds_tree_search(DS_Tree *tree, void *obj, DSCompareFunc fct);

#include "dsbtree.h"

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_DSTRUCT_H__ */
