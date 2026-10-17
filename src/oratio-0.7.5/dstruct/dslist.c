/**
 * dslist.c
 *
 * D Struct - ORATIO Library
 * Copyright (C) 2000-2005 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/ctools.h"
#include "oratio/dstruct.h"

#include <string.h>

DS_List *ds_list_new()
{
  return NULL;
}

DS_List *ds_list_new_with_obj(void *obj)
{
  DS_List *list;

  list = _malloc(DS_List);
  list->data = obj;
  list->next = NULL;
  list->prev = NULL;

  return list;
}

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
void ds_list_free(DS_List *list, Dfct ff)
{
  if (list) {
    DS_List *root;
    if (list->prev) list->prev->next = NULL;

    while (list) {
      root = list;
      list = root->next;
      if (ff) ff(root->data);
      free(root);
    }
  }
}

DS_List *ds_list_first(DS_List *list)
{
  if (list)
    while (list->prev) list = list->prev;

  return list;
}

DS_List *ds_list_last(DS_List *list)
{
  if (list)
    while (list->next) list = list->next;

  return list;
}

/**
 * Compute the list length.
 *
 * @param list is an initialized simply-linked list.
 *
 * @return the list length.
 */
int ds_list_length(DS_List *list)
{
  int length = 0;

  while (list) {
    length++;
    list = list->next;
  }

  return length;
}

void ds_list_set_data(DS_List *list, void *obj)
{
  if (list)
    list->data = obj;
  else
    c_warning("list undefined in routine ds_list_set_data");
}

/**
 * Return the content of the first element of a double linked list.
 *
 * @return the content of the first element of list, or NULL is the
 *   list is empty.
 */
void *ds_list_get_data(DS_List *list)
{
  if (list)
    return list->data;
  else
    return NULL;
}

DS_List *ds_list_nth(DS_List *list, int n)
{
  while ((n-- > 0) && list)
    list = list->next;
  return list;
}

void *ds_list_nth_data(DS_List *list, int n)
{
  while ((n-- > 0) && list)
    list = list->next;
  return list ? list->data : NULL;
}

long int ds_list_index(DS_List *list, DS_List *node)
{
  long int i = 0;

  while (list) {
    if (list == node) return i;
    i++;
    list = list->next;
  }

  return -1;
}

long int ds_list_index_data(DS_List *list, void *obj)
{
  long int i = 0;

  while (list) {
    if (list->data == obj) return i;
    i++;
    list = list->next;
  }

  return -1;
}

DS_List *ds_list_prepend(DS_List *list, void *obj)
{
  DS_List *new_list;

  new_list = ds_list_new_with_obj(obj);

  if (list) {
    if (list->prev) {
      list->prev->next = new_list;
      new_list->prev = list->prev;
    }
    list->prev = new_list;
    new_list->next = list;
  }

  return new_list;
}

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
DS_List *ds_list_append(DS_List *list, void *obj)
{
  DS_List *new_list;
  DS_List *last;

  new_list = ds_list_new_with_obj(obj);

  if (list) {
    last = ds_list_last(list);
    last->next = new_list;
    new_list->prev = last;
    return list;
  }
  else
    return new_list;
}

/**
 * Remove the first element of list.
 *
 * @param list is an initialized doubly-linked list.
 *
 * @return the first element of list is list is not empty, NULL otherwise.
 */
DS_List *ds_list_pop(DS_List *list)
{
  if (list) {
    DS_List *root;

    root = list;
    if (list->prev) list->prev->next = list->next;
    if (list->next) list->next->prev = list->prev;
    list = list->next;

    free(root);
  }
  return list;
}

DS_List *ds_list_copy(DS_List *list)
{
  DS_List *new_list = NULL;

  if (list) {
    DS_List *last;

    new_list = ds_list_new_with_obj(list->data);
    last = new_list;
    list = list->next;
    while(list) {
      last->next = ds_list_new_with_obj(list->data);
      last->next->prev = last;
      last = last->next;
      list = list->next;
    }
  }

  return new_list;
}

DS_List *ds_list_clone_item(DS_List *list, size_t size_data)
{
  DS_List *item = NULL;

  if (list) {
    item = _malloc(DS_List);
    item->next = NULL;
    item->prev = NULL;
    if (size_data) {
      item->data = malloc(size_data);
      memcpy(item->data, list->data, size_data);
    }
    else
      item->data = list->data;
  }

  return item;
}

DS_List *ds_list_clone_item_with_function(DS_List *list, size_t size_data,
					  void *(*f)(void *))
{
  DS_List *item = NULL;

  if (list) {
    item = _malloc(DS_List);
    item->next = NULL;
    item->prev = NULL;
    item->data = f(list->data);
  }

  return item;
}

DS_List *ds_list_clone(DS_List *list, size_t size_data)
{
  DS_List *new_list = NULL;

  if (list) {
    DS_List *last;

    new_list = ds_list_clone_item(list, size_data);
    last = new_list;
    list = list->next;
    while(list) {
      last->next = ds_list_clone_item(list, size_data);
      last->next->prev = last;
      last = last->next;
      list = list->next;
    }
  }

  return new_list;
}

/**
 * Create a list copy using a clone function.
 */
DS_List *ds_list_clone_with_function(DS_List *list, size_t size_data,
				     void *(*f)(void *))
{
  DS_List *new_list = NULL;

  if (list) {
    DS_List *last;

    new_list = ds_list_clone_item_with_function(list, size_data, f);
    last = new_list;
    list = list->next;
    while (list) {
      last->next = ds_list_clone_with_function(list, size_data, f);
      last->next->prev = last;
      last = last->next;
      list = list->next;
    }
  }

  return new_list;
}

DS_List *ds_list_concat(DS_List *list1, DS_List *list2)
{
  if (list2) {
    if (list1) {
      DS_List *first, *last;
      first = ds_list_first(list2);
      last = ds_list_last(list1);
      last->next = first;
      first->prev=last;
    }
    else
      list1 = list2;
  }
  return list1;
}

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
DS_List *ds_list_reverse(DS_List *list)
{
  DS_List *tmp;
  DS_List *prev = NULL;
  DS_List *last = NULL;
  DS_List *before;

  if (list)
    before = list->prev;
  else
    before = NULL;

  while (list) {
    last = list;
    tmp = list->next;
    list->next = prev;
    list->prev = tmp;

    prev = list;
    list = tmp;
  }

  if (before) {
    before->next = last;
    last->prev = before;
  }

  return last;
}

DS_List *ds_list_find(DS_List *list, void *obj)
{
  while (list) {
    if (list->data == obj) break;
    list = list->next;
  }

  return list;
}

/**
 * Remove the first element having obj as data from the list,
 * starting from list (forward removal).
 *
 * If the suppressed element corresponds to list, list is set to the
 * next element if it exists, and to the previous one else.
 */
DS_List *ds_list_remove(DS_List *list, void *obj)
{
  if (list) {
    DS_List *tmp;

    tmp = list;

    /* We ensure that list is still valid (possibly by setting it to NULL) */
    if (list->data == obj) {
      list = (list->next) ? list->next : list->prev;
      return list;
    }

    while (tmp) {
      if (tmp->data == obj) {
	if (tmp->prev)
	  tmp->prev->next = tmp->next;
	else
	  list = list->next;
	if (tmp->next)
	  tmp->next->prev = tmp->prev;
	free(tmp);
	break;
      }
      tmp = tmp->next;
    }
  }

  return list;
}

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
DS_List *ds_list_remove_node(DS_List *list, DS_List *node)
{
  DS_List *tmp;

  tmp = list;

  /* We ensure that list is still valid (possibly by setting it to NULL) */
  if (list && list == node) {
    list = (list->next) ? list->next : list->prev;
    return list;
  }

  while (tmp) {
    if (tmp == node) {
      if (tmp->prev)
	tmp->prev->next = tmp->next;
 
      if (tmp->next)
	tmp->next->prev = tmp->prev;

      free(tmp);
      break;
    }
    tmp = tmp->next;
  }

  return list;
}

/**
 * Remove the current element pointed by a list.
 *
 * @param list is a pointer to the list under consideration; on output,
 *   the element 'node' has been removed from this list.
 *
 * @return a pointer to the next element if it exists, and to the previous
 *   one else.
 */
DS_List *ds_list_drop(DS_List *list)
{
  return ds_list_remove_node(list, list);
}

/**
 * Convert the list into an array.
 *
 * Each array element will be a pointer to the data structure
 * previously referenced by the corresponding list element.
 */
void **ds_list_to_array(DS_List *list)
{
  void **array = NULL;

  if (list) {
    int i = 0;
    DS_List *node = list;
    array = _mallocs(ds_list_length(list), void *);
    while (node) {
      array[i++] = node->data;
      node = node->next;
    }
  }

  return array;
}

void *ds_list_to_array_with_size(DS_List *list, size_t size)
{
  void *array = NULL;

  if (list) {
    int i = 0;
    DS_List *node = list;
    array = malloc(ds_list_length(list)*size);
    while (node) {
      memcpy(array+(i++)*size, node->data, size);
      node = node->next;
    }
  }

  return array;
}

void ds_list_foreach(DS_List *list, Dfct ff)
{
  while (list) {
    if (ff) ff(list->data);
    list = list->next;
  }
}

DS_List *ds_list_insert_sorted(DS_List *list, void *data,
			       DSCompareFunc f)
{
  DS_List *tmp_list = list;
  DS_List *new_list, *prev;

  new_list = ds_list_new_with_obj(data);

  if (!list) return new_list;
  
  prev = list->prev;

  while (f(data, tmp_list->data) > 0) {
    prev = tmp_list;
    tmp_list = tmp_list->next;
    if (!tmp_list) break;
  }

  new_list->next = tmp_list;
  if (tmp_list) tmp_list->prev = new_list;
  if (prev) {
    prev->next = new_list;
    new_list->prev = prev;
    return list;
  }
  else
    return new_list;
}

DS_List *ds_list_insert_sorted_last(DS_List *list, void *data,
				    DSCompareFunc f)
{
  DS_List *tmp_list = list;
  DS_List *new_list, *prev;

  new_list = ds_list_new_with_obj(data);

  if (!list) return new_list;
  
  prev = list->prev;

  while (f(data, tmp_list->data) >= 0) {
    prev = tmp_list;
    tmp_list = tmp_list->next;
    if (!tmp_list) break;
  }

  new_list->next = tmp_list;
  if (tmp_list) tmp_list->prev = new_list;
  if (prev) {
    prev->next = new_list;
    new_list->prev = prev;
    return list;
  }
  else
    return new_list;
}
