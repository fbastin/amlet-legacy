/**
 * dsslist.h
 *
 * D Struct
 * Copyright (C) 2000-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_DSTRUCT_DSSLIST_H__
#define __ORATIO_DSTRUCT_DSSLIST_H__

#include "ds_functions.h"

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**
   * Singly linked chain
   */

  typedef struct _DS_SList DS_SList;

  struct _DS_SList {
    void *data;
    DS_SList *next;
  };

  /**
   * Create a new singly linked list.
   *
   * @return return a newly constructed singly linked list.        
   */
  DS_SList *ds_slist_new();

  /**
   * Free the memory associated to a singly linked list.
   *
   * @param list is the singly list to free.
   * @param ff is a function that will be applied to each list element, in
   *   order to free associated memory. ff can be set to NULL if no function
   *   has to be called.               
   */
  void ds_slist_free(DS_SList *list, Dfct ff);

  /**
   * Return a pointer to the last list position.
   *
   * @param list is the singly linked list under consideration.   
   *
   * @return a pointer to the last list position.
   */
  DS_SList *ds_slist_last(DS_SList *list);

  /**
   * Check if is the list is empty.
   *
   * @param list is the singly linked list under consideration.
   *
   * @return 1 if the list is empty, 0 otherwise.
   */
#define ds_slist_empty(list) ((list) ? 0 : 1)

  /**
   * Length of a list.
   *    
   * @param list is the singly linked list under consideration.
   *
   * @return the list length. 
   */
  int ds_slist_length(DS_SList *list);

  /**
   * Replace the first list element pointer.
   *     
   * @param list is the singly linked list under consideration.
   * @param obj a pointer to the object to become the first list element
   *   pointer.
   */
  void ds_slist_set_data(DS_SList *list, void *obj);

  /**
   * Return the first list element pointer.
   *
   * @param list is the singly linked list under consideration.
   *         
   * @return a pointer to the first list element pointer, or NULL if the list
   *   is empty.   
   */
  void *ds_slist_get_data(DS_SList *list);

  /**
   * Return a pointer to the nth position in a singly linked list.  
   *
   * @param list is the singly linked list under consideration.
   * @param n is the position (>=0) to look for.
   *
   * @return a singly linked list pointer corresponding the the nth position
   *   of the original list, or NULL if the list has strictly less than n
   *   elements.      
   */
  DS_SList *ds_slist_nth(DS_SList *list, int n);

  /**
   * Return a pointer to the nth element of a singly linked list.
   *    
   * @param list is the singly linked list under consideration.
   * @param n is the position (>=0) to look for.
   *
   * @return a pointer to the element stored at position n.   
   */
  void *ds_slist_nth_data(DS_SList *list, int n);

  /**
   * Search a node index in a singly linked list.
   *
   * @param list is the singly linked list under consideration.
   * @param node is the node list to search for.
   *
   * @return position of node in the list, if found, or -1 the search has not
   *   produced has been unsuccessful.                 
   */
  int ds_slist_index(DS_SList *list, DS_SList *node);

  /**
   * Search the position of some element in a singly linked list.
   *      .  
   * @param list is the singly linked list under consideration.
   * @param obj is a pointer to the element to search for.
   *   
   * @return position of element obj in the list, if found, or -1 the search has
   *   has produced has been unsuccessful.
   */                    
  int ds_slist_index_data(DS_SList *list, void *obj);

  /**
   * Add an element in front of a singly linked list.
   *
   * @param list is the singly linked list under consideration (possibly empty).
   * @param obj is the element pointer to add to the list.
   *
   * @return a pointer to the new list.
   */
  DS_SList *ds_slist_prepend(DS_SList *list, void *obj);

  /**
   * Add an element in front of a singly linked list (macro version).
   *
   * @param list is the singly linked list under consideration (possibly empty).
   * @param obj is the element pointer to add to at the beginning of list.
   *
   * @return a pointer to the new list, with the element obj added at the
   *   beginning of the list.
   */
#define ds_slist_push(list, obj) ds_slist_prepend(list, obj);

  /**
   * Add an element at the end of a singly linked list.
   *
   * @param list is the singly linked list under consideration (possibly empty).
   * @param obj is the element pointer to add at the of the list.
   *
   * @return a pointer to the new list, with the element obj added at the
   *   end of the list.
   */
  DS_SList *ds_slist_append(DS_SList *list, void *obj);

  /**
   * Remove the first element of a list.
   *
   * @param list is the singly linked list under consideration.
   *
   * @return a pointer to the list with the first element removed (no action
   *   is taken to ensure that associated memory has been freed).
   */
  DS_SList *ds_slist_pop(DS_SList *list);

  /**
   * Copy a list, sharing elements pointers.
   *
   * @param list is the singly linked list under consideration.
   *
   * @return a pointer to the list copy.
   */
  DS_SList *ds_slist_copy(DS_SList *list);

  /**
   * Clone a list, using a cloning function for the elements.
   *
   * @param list is the singly linked list under consideration.
   * @param size_data is the size of the list elements.
   * @param clone is a cloning function.
   *
   * @return a pointer to the newly-created, clone list.
   */
  DS_SList *ds_slist_clone_with_function(DS_SList *list, size_t size_data,
					 void *(*clone)(void *));

  /**
   * Clone a list.
   *
   * @param list is the singly linked list under consideration.
   * @param size_data is the size of the list elements.
   *
   * @return a pointer to the newly-created, clone list.
   */
  DS_SList *ds_slist_clone(DS_SList *list, size_t size_data);

  /**
   * Concatenate two lists.
   */
  DS_SList *ds_slist_concat(DS_SList *list1, DS_SList *list2);

  /**
   * Inverse the list order.
   */
  DS_SList *ds_slist_reverse(DS_SList *list);

  /**
   * Renvoie le premier maillon de la liste list ayant comme donnée obj.
   * @param list is the singly linked list under consideration.
   */
  DS_SList *ds_slist_find(DS_SList *list, void *obj);

  /**
   * Enlève de la liste list le premier maillon ayant comme donnée obj.
   * @param list is the singly linked list under consideration.
   */
  DS_SList *ds_slist_remove(DS_SList *list, void *obj);

  /**
   * Enlève de la liste list le maillon node.
   * @param list is the singly linked list under consideration.
   */
  DS_SList *ds_slist_remove_node(DS_SList *list, DS_SList *node);

  /**
   * Convertit la liste list en tableau. Chaque élément du tableau sera
   * un pointeur vers la structure de données pointée par le maillon de
   * la liste correspondant.
   * @param list is the singly linked list under consideration.
   */
  void **ds_slist_to_array(DS_SList *list);

  /**
   * Convertit la liste list en tableau, en copiant size bytes depuis
   * le maillon de la liste dans la case du tableau correspondante.
   *
   * The data are duplicated.
   *
   * @param list is the singly linked list under consideration.
   * @param size is the size in bytes of one element.
   */
  void *ds_slist_to_array_with_size(DS_SList *list, size_t size);

  /**
   * Apply the function ff to each list element.
   *
   * @param list is the singly linked list under consideration.
   */
  void ds_slist_foreach(DS_SList *list, Dfct ff);

  /**
   * Insère un élément en respectant l'ordre établi par f
   *
   * @param list is the singly linked list under consideration.
   */
  DS_SList *ds_slist_insert_sorted(DS_SList *list, void *data,
				   DSCompareFunc f);

  /**
   * Return the next link, or NULL if the list doesn't exist.
   *
   * @param list is the singly linked list under consideration.
   */
#define ds_slist_next(slist)	((slist) ? (((DS_SList *)(slist))->next) : NULL)

#define _ds_slist_get_data(slist) ((slist) ? (((DS_SList *)(slist))->data) : NULL)

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ORATIO_DSTRUCT_DSSLIST_H__ */
