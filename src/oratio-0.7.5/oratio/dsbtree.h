/**
 * dsbtree.h
 *
 * D Struct
 * Copyright (C) 2000-2002 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_DSTRUCT_DSBTREE_H__
#define __ORATIO_DSTRUCT_DSBTREE_H__

#include "dstruct.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Binary tree
   */

  enum {
    LEFT = 1,
    RIGHT
  };

  typedef struct _DS_BTree DS_BTree;

  struct _DS_BTree {
    void *root;
    DS_BTree *left;
    DS_BTree *right;
  };

  /**
   * Create a new tree
   */
  DS_BTree *ds_btree_new();

  /**
   * Création d'un nouvel arbre en mettant obj comme racine
   */
  DS_BTree *ds_btree_new_with_root(void *obj);

  /**
   * Suppression d'un arbre.
   * @param ff est une fonction à appliquer à un élément de l'arbre pour
   * libérer la mémoire qui lui est associée. ff peut être mis à NULL si
   * aucune fonction ne doit être utilisée.
   */
  void ds_btree_free(DS_BTree *tree, Dfct ff);

  /**
   * Modification de la valeur de la racine.
   * Au besoin, la routine crée un nouvel arbre.
   * Renvoie l'arbre modifié/créé.
   */
  DS_BTree *ds_btree_set_root(DS_BTree *tree, void *obj);

  /**
   * Insertion d'un élément à gauche de la racine
   * Au besoin, la routine crée un nouvel arbre.
   * Renvoie l'arbre modifié/créé.
   */
  DS_BTree *ds_btree_insert_left(DS_BTree *tree, void *obj);

  /**
   * Insertion d'un élément à droite de la racine
   * Au besoin, la routine crée un nouvel arbre.
   * Renvoie l'arbre modifié/créé.
   */
  DS_BTree *ds_btree_insert_right(DS_BTree *tree, void *obj);

  /**
   * Insertion d'un élément, comme fils de la racine
   * @param place indique où on veut insérer l'élément: LEFT si on souhaite
   * le placer à gauche de la racine; RIGHT pour le mettre à droite
   * Au besoin, la routine crée un nouvel arbre.
   * Renvoie l'arbre modifié/créé.
   */
  DS_BTree *ds_btree_insert(DS_BTree *tree, void *obj, int place);

  /**
   * Applique la fonction fct à chaque élément de l'arbre
   */
  void ds_btree_foreach(DS_BTree *tree, Dfct fct);

  /**
   * Compte le nombre de feuilles de l'arbre
   * @return le nombre de feuilles de l'arbre
   */
  int ds_btree_count_leaves(DS_BTree *tree);

  /**
   * Explicite chaque branche de l'arbre.
   */
  void ds_btree_expand(DS_BTree *tree, DS_SList ***branches, int *nb_branches);

  /**
   * Explicite chaque branche d'un arbre de décision. Chaque élément de l'arbre
   * vaut 0 ou 1. Seul les éléments marqués 1 sont pris en considération,
   * et est repris dans la branche formés en notant le numéro de son niveau,
   * en mettant à -1 la racine de l'arbre.
   */
  void ds_btree_boolean_expand(DS_BTree *tree, DS_SList ***branches,
			      int *nb_branches);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_DSTRUCT_DSBTREE_H__ */
