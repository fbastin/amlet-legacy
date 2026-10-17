/**
 * dsbtree.c
 *
 * ORATIO
 * Copyright (C) 2000-2002 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/memory.h"
#include "oratio/dstruct.h"

DS_BTree *ds_btree_new()
{
  return NULL;
}

DS_BTree *ds_btree_new_with_root(void *obj)
{
  DS_BTree *tree = _malloc(DS_BTree);

  tree->root = obj;
  tree->left = NULL;
  tree->right = NULL;

  return tree;
}

void ds_btree_free(DS_BTree *tree, Dfct ff)
{
  if (tree) {
    if (ff) ff(tree->root);
    ds_btree_free(tree->left, ff);
    ds_btree_free(tree->right, ff);

    free(tree);
  }
}

DS_BTree *ds_btree_set_root(DS_BTree *tree, void *obj)
{
  if (tree)
    tree->root = obj;
  else
    tree = ds_btree_new_with_root(obj);
  return tree;
}

DS_BTree *ds_btree_insert_left(DS_BTree *tree, void *obj)
     /* FIXME: dans l'implémentation actuelle, s'il existe un sous-arbre,
      * celui-ci est recopié à gauche de l'objet inséré. Dès lors, des appels
      * successifs de cette fonction risquent de désiquilibrer l'arbre
      */
{
  if (tree) {
    DS_BTree *tmp_btree = tree->left;
    tree->left = ds_btree_new_with_root(obj);
    tree->left->left = tmp_btree;
  }
  else
    tree = ds_btree_new_with_root(obj);
  return tree;
}

DS_BTree *ds_btree_insert_right(DS_BTree *tree, void *obj)
{
  if (tree) {
    DS_BTree *tmp_btree = tree->right;
    tree->right = ds_btree_new_with_root(obj);
    tree->right->right = tmp_btree;
  }
  else
    tree = ds_btree_new_with_root(obj);
  return tree;
}

DS_BTree *ds_btree_insert(DS_BTree *tree, void *obj, int place)
{
  if (place == LEFT) return ds_btree_insert_left(tree, obj);
  if (place == RIGHT) return ds_btree_insert_right(tree, obj);
  return NULL;
}

void ds_btree_foreach(DS_BTree *tree, Dfct fct)
{
  /* Traitement de l'arbre en parcours préfixe (?) */
  if (tree) {
    fct(tree->root);
    ds_btree_foreach(tree->left, fct);
    ds_btree_foreach(tree->right, fct);
  }
}

int ds_btree_count_leaves(DS_BTree *tree)
{
  int leaves = 0;

  if (tree) {
    if (!tree->left && !tree->right) leaves++;
    if (tree->left) leaves += ds_btree_count_leaves(tree->left);
    if (tree->right) leaves += ds_btree_count_leaves(tree->right);
  }

  return leaves;
}

void ds_btree_expand_leaf(DS_SList ***branches, int *nb_branches, void *obj)
{
  *nb_branches = 1;
  *branches = _malloc(DS_SList *);
  (*branches)[0] = NULL;
  (*branches)[0] = ds_slist_prepend((*branches)[0], obj);
}

void ds_btree_expand(DS_BTree *tree, DS_SList ***branches, int *nb_branches)
{
  if(tree) {
    int i;
    DS_SList **branches_left, **branches_right;
    int nb_branches_left, nb_branches_right;

    ds_btree_expand(tree->left, &branches_left, &nb_branches_left);
    ds_btree_expand(tree->right, &branches_right, &nb_branches_right);

    if (!nb_branches_left && !nb_branches_right)
      /* on a atteint une feuille, i.e. on termine une branche */
      ds_btree_expand_leaf(branches, nb_branches, tree->root);
    else {
      /* Ajout de la racine aux branches obtenues */
      for (i = 0; i < nb_branches_left; i++)
	branches_left[i] = ds_slist_prepend(branches_left[i], tree->root);
      for (i = 0; i < nb_branches_right; i++)
	branches_left[i] = ds_slist_prepend(branches_right[i], tree->root);
      /* Fusion des deux tableaux obtenus */
      branches_left = (DS_SList **)realloc(branches_left,
					 (nb_branches_left+nb_branches_right)*
					 sizeof(DS_SList *));
      for (i = 0; i < nb_branches_right; i++)
	branches_left[i+nb_branches_left] = branches_right[i];
      free(branches_right);
      *branches = branches_left;
      *nb_branches = nb_branches_left + nb_branches_right;
    }
  }
  else {
    *branches = NULL;
    *nb_branches = 0;
  }
}

/* void ds_btree_boolean_expand_core(DS_BTree *tree, DS_SList ***branches,
				 int *nb_branches, int level)
{
  if(tree) {
    int i;
    DS_SList **branches_left, **branches_right;
    int nb_branches_left, nb_branches_right;

    ds_btree_boolean_expand_core(tree->left, &branches_left,
				&nb_branches_left, level+1);
    ds_btree_boolean_expand_core(tree->right, &branches_right,
				&nb_branches_right, level+1);

    if (!nb_branches_left && !nb_branches_right) {
      if (tree->root == (void *)1)
	ds_btree_expand_leaf(branches, nb_branches, (void *)level);
      else {
	*nb_branches = 0;
	*branches = NULL;
      }
    }
    else {
      // Ajout de la racine aux branches obtenues si elle vaut 1
    if (tree->root == (void *)1) {
	if (nb_branches_left || nb_branches_right) {
	  for (i = 0; i < nb_branches_left; i++)
	    branches_left[i] = ds_slist_prepend(branches_left[i],
					       (void *)level);
	  for (i = 0; i < nb_branches_right; i++)
	    branches_right[i] = ds_slist_prepend(branches_right[i],
						(void *)level);
	}
	else
	  ds_btree_expand_leaf(&branches_left, &nb_branches_left,
			       (void *)level);
    }
    // Fusion des deux tableaux obtenus
    branches_left = _reallocs(branches_left,
			      nb_branches_left+nb_branches_right,
			      DS_SList *);
    for (i = 0; i < nb_branches_right; i++)
      branches_left[i+nb_branches_left] = branches_right[i];
    free(branches_right);
    *branches = branches_left;
    *nb_branches = nb_branches_left + nb_branches_right;
    }
  }
  else {
    *branches = NULL;
    *nb_branches = 0;
  }
  }

void ds_btree_boolean_expand(DS_BTree *tree, DS_SList ***branches,
			    int *nb_branches)
{
  int root_level = -1;
  ds_btree_boolean_expand_core(tree, branches, nb_branches, root_level);
}
*/
