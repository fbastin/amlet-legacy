/**
 * dslist.h
 *
 * D Struct
 * Copyright (C) 2000-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_DSTRUCT_DSLIST_H__
#define __ORATIO_DSTRUCT_DSLIST_H__

#include <sys/types.h>

#include "ds_functions.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Doubled linked list
   */

  typedef struct _DS_List DS_List;

  struct _DS_List {
    void *data;
    DS_List *next;
    DS_List *prev;
  };

  DS_List *ds_list_new();

  /**
   * Free the memory associated to list.
   *
   * The function frees memory associated to list and the next
   * elements, but not the previous ones. ds_list_first can be used to
   * first set the pointer to the list beginning.
   *
   * @param ff is a function that must be applied to each list
   *   element, in order to automize the memory freeing process. ff
   *   can be set to NULL if no function has to be called.
   */
  void ds_list_free(DS_List *list, Dfct ff);

  /**
   * Return a pointer to the first list element.
   */
  DS_List *ds_list_first(DS_List *list);

  /**
   * Return a pointer to the last list element.
   */
  DS_List *ds_list_last(DS_List *list);

  /**
   * Return the list length.
   *
   * @param list is an initialized simply-linked list.
   *
   * @return the list length.
   */
  int ds_list_length(DS_List *list);

  /**
   * Return 1 if the list is empty, 0 otherwise
   */
#define ds_list_empty(list) ((list) ? 0 : 1)

  /**
   * Modifie l'élément de donnée du premier maillon de la liste list par obj.
   * #PRE: list différent de NULL
   */
  void ds_list_set_data(DS_List *list, void *obj);

  /**
   * Return the content of the first element of a double linked list.
   *
   * @return the content of the first element of list, or NULL is the
   *   list is empty.
   */
  void *ds_list_get_data(DS_List *list);

  /**
   * Return the n-th element of the list, or NULL if this element does
   * not exist.
   */
  DS_List *ds_list_nth(DS_List *list, int n);

  /**
   * Ajoute un maillon à liste list, en tête de celle-ci. La maillon
   * nouvellement créé reçoit obj comme donnée. list peut être vide, auquel
   * cas le pointeur doit posséder la valeur NULL.
   */
  DS_List *ds_list_prepend(DS_List *list, void *obj);

  /**
   * Add an element at the end of a double linked list.
   *
   * The addition of an element to the end of the list is far less
   * efficient than the addition in the beginning of this list.
   *
   * @param list is a pointer to a double linked list, empty (and then
   *   set to NULL), or not.
   * @param obj is the object that must be pointed by the new element list.
   */
  DS_List *ds_list_append(DS_List *list, void *obj);

  /**
   * Remove the first element of list.
   *
   * @param list is an initialized doubly-linked list.
   *
   * @return the first element of list is list is not empty, NULL otherwise.
   */
  DS_List *ds_list_pop(DS_List *list);

  /**
   * Crée une copie de la liste. Dans l'implémentation actuelle, les données
   * des deux listes sont partagés, ce qui peut conduire à certains problèmes.
   */
  DS_List *ds_list_copy(DS_List *list);

  /**
   * Create a list copy using a clone function.
   */
  DS_List *ds_list_clone_with_function(DS_List *list, size_t size_data,
				       void *(*f)(void *));

  /**
   * Crée une copie de la liste, dont les données sont indépendantes de celles
   * de la liste originale
   *
   * #PRE:  size_data indique le nombre de bytes occupés par les données d'un
   *        noeud de la liste; si size_data est nul, cela signifie que
   *        le champ data du noeud n'est pas un pointeur, mais un entier.
   */
  DS_List *ds_list_clone(DS_List *list, size_t size_data);

  /**
   * Concatène les listes list1 et list2; plus explicitement, la fin de la
   * liste list1 est concaténée avec le début de la liste contenant list2.
   */
  DS_List *ds_list_concat(DS_List *list1, DS_List *list2);

  /**
   * Invert the order of a double linked list.
   *
   * @param list is a pointer to a double linked list, empty (and then
   *   set to NULL), or not. If list has element before the one pointed
   *   by list, the order of these previous elements is not
   *   modified. Only the list part beginning at the current element is
   *   reversed.
   *
   * @return a pointer to the first element of the reversed list.
   */
  DS_List *ds_list_reverse(DS_List *list);

  /**
   * Renvoie le premier maillon de la liste list ayant comme donnée obj.
   */
  DS_List *ds_list_find(DS_List *list, void *obj);

  /**
   * Remove the first element having obj as data from the list,
   * starting from list (forward removal).
   *
   * If the suppressed element corresponds to list, list is set to the
   * next element if it exists, and to the previous one else.
   */
  DS_List *ds_list_remove(DS_List *list, void *obj);

  /**
   * Remove the node element form the list, starting from list
   * (forward removal).
   *
   * @param list is a pointer to the list under consideration; on output,
   *   the element 'node' has been removed from this list.
   * @param node is the node to suppress.
   *
   * @return the original pointer list, except if the suppressed element
   *   corresponds to list; in this last case, the returned pointer is the
   *   next element if it exists, and the previous one else.
   */
  DS_List *ds_list_remove_node(DS_List *list, DS_List *node);

  /**
   * Remove the current element pointed by a list.
   *
   * @param list is a pointer to the list under consideration; on output,
   *   the element 'node' has been removed from this list.
   *
   * @return a pointer to the next element if it exists, and to the previous
   *   one else.
   */
  DS_List *ds_list_drop(DS_List *list);

  /**
   * Convert the list into an array.
   *
   * Each array element will be a pointer to the data structure
   * previously referenced by the corresponding list element.
   */
  void **ds_list_to_array(DS_List *list);

  /**
   * Convert a list into an array.
   *
   * @param list is an initialized double-linked list.
   * @param size is the number of bytes that are copied from a list
   *   element to a cell of the array.
   *
   * @return a pointer to the array beginning.
   */
  void *ds_list_to_array_with_size(DS_List *list, size_t size);

  /**
   * Apply the function ff to each element of the list
   *
   * @param list is an initialized double-linked list.
   * @param ff is a function to be applied to each list member.
   */
  void ds_list_foreach(DS_List *list, Dfct ff);

  /**
   * Insert an element in a sorted list.
   *
   * @param list is a sorted double-linked list.
   */
  DS_List *ds_list_insert_sorted(DS_List *list, void *data,
				 DSCompareFunc f);

  /**
   * Insert an element in a sorted list. In case two elements have the same
   * order, the new one is put after the previous one.
   *
   * @param list is a sorted double-linked list.
   */
  DS_List *ds_list_insert_sorted_last(DS_List *list, void *data,
				      DSCompareFunc f);

  /**
   * Renvoie le maillon suivant, ou NULL s'il la liste n'existe pas
   */
#define ds_list_next(list)	((list) ? (((DS_List *)(list))->next) : NULL)
#define ds_list_prev(list)      ((list) ? (((DS_List *)(list))->prev) : NULL)

#define _ds_list_get_data(list) ((list) ? (((DS_List *)(list))->data) : NULL)

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_DSTRUCT_DSLIST_H__ */
