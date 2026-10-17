/**
 * dstree.c
 *
 * ORATIO Library
 * Copyright (C) 2000-2005 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/memory.h"
#include "oratio/dstruct.h"

DS_Tree *ds_tree_new()
{
  return NULL;
}

DS_Tree *ds_tree_new_with_root(void *obj)
{
  DS_Tree *tree = _malloc(DS_Tree);

  tree->root = obj;
  tree->children = ds_list_new();
  tree->parent = NULL;

  return tree;
}

void ds_tree_free(DS_Tree *tree, Dfct ff)
{
  DS_List *children;

  if (tree) {
    if (ff) ff(tree->root);
    children = tree->children;
    while (children) {
      ds_tree_free(children->data, ff);
      children = children->next;
    }
    ds_list_free(tree->children, NULL);
    free(tree);
  }
}

/*****************************************************************************/

DS_Tree *ds_tree_set_root(DS_Tree *tree, void *obj)
{
  if (tree)
    tree->root = obj;
  else
    tree = ds_tree_new_with_root(obj);

  return tree;
}

DS_Tree *ds_tree_insert(DS_Tree *parent, void *obj)
{
  DS_Tree *child;

  child = ds_tree_new_with_root(obj);
  if (!parent) return child;
  child->parent = parent;
  parent->children = ds_list_append(parent->children, child);

  return parent;
}

/*****************************************************************************/

void ds_tree_foreach(DS_Tree *tree, Dfct ff)
{
  DS_List *children;
  if (tree) {
    ff(tree->root);
    children = tree->children;
    while (children) {
      ds_tree_foreach(children->data, ff);
      children = children->next;
    }
  }
}

/*****************************************************************************/

DS_Tree *ds_tree_search(DS_Tree *tree, void *obj, DSCompareFunc fct)
{
  if (tree) {
    DS_Tree *found;
    DS_List *children = tree->children;

    if (!fct(tree->root, obj)) return tree;

    while (children) {
      if ((found = ds_tree_search(children->data, obj, fct)))
	  return found;
      children = children->next;
    }
  }

  return NULL;
}
