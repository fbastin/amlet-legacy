/**
 * dsslist.c
 *
 * ORATIO Library
 * Copyright (C) 2000-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/ctools.h>
#include "oratio/dstruct.h"
#include "oratio/dsslist.h"

#include <string.h>

/**
 * Create a new simply linked list.
 *
 * @return a pointed to a newly created simply linked-list.      
 */
DS_SList *ds_slist_new()
{
  return NULL;
}

DS_SList *ds_slist_new_with_obj(void *obj)
{
  DS_SList *l;

  l = _malloc(DS_SList);

  l->data = obj;
  l->next = NULL;

  return l;
}

void ds_slist_free(DS_SList *list, Dfct ff)
{
  DS_SList *root;

  while (list) {
    root = list;
    list = root->next;
    if (ff) ff(root->data);
    free(root);
  }
}

DS_SList *ds_slist_last(DS_SList *list)
{
  if (list) {
    while (list->next)
      list = list->next;
  }
  return list;
}

int ds_slist_length(DS_SList *list)
{
  int length;

  length=0;
  while (list) {
    length++;
    list = list->next;
  }

  return length;
}

void ds_slist_set_data(DS_SList *list, void *obj)
{
  if (list)
    list->data = obj;
  else
    c_warning("list undefined in routine ds_slist_set_data");
}

void *ds_slist_get_data(DS_SList *list)
{
  if (list)
    return list->data;
  else
    return NULL;
}

DS_SList *ds_slist_nth(DS_SList *list, int n)
{
  while ((n-- > 0) && list)
    list = list->next;
  return list;
}

void *ds_slist_nth_data(DS_SList *list, int n)
{
  while ((n-- > 0) && list)
    list = list->next;
  return list ? list->data : NULL;
}

int ds_slist_index(DS_SList *list, DS_SList *node)
{
  int i = 0;

  while (list) {
    if (list == node) return i;
    i++;
    list = list->next;
  }

  return -1;
}

int ds_slist_index_data(DS_SList *list, void *obj)
{
  int i = 0;

  while (list) {
    if (list->data == obj) return i;
    i++;
    list = list->next;
  }

  return -1;
}

DS_SList *ds_slist_prepend(DS_SList *list, void *obj)
{
  DS_SList *new_list;

  new_list = ds_slist_new_with_obj(obj);
  new_list->next = list;

  return new_list;
}

DS_SList *ds_slist_append(DS_SList *list, void *obj)
{
  DS_SList *new_list;
  DS_SList *last;

  new_list = ds_slist_new_with_obj(obj);

  if (list) {
    last = ds_slist_last(list);
    last->next = new_list;
    return list;
  }
  else
    return new_list;
}

DS_SList *ds_slist_pop(DS_SList *list)
{
  if (list) {
    DS_SList *root;
    root = list;
    list = list->next;
    free(root);
  }
  return list;
}

DS_SList *ds_slist_copy(DS_SList *list)
{
  DS_SList *new_list = NULL;

  if (list) {
    DS_SList *last;

    new_list = ds_slist_new_with_obj(list->data);
    last = new_list;
    list = list->next;
    while(list) {
      last->next = ds_slist_new_with_obj(list->data);
      last = last->next;
      list = list->next;
    }
  }

  return new_list;
}

DS_SList *ds_slist_clone_item(DS_SList *list, size_t size_data,
			      void *(*f)(void *))
{
  DS_SList *item = NULL;

  if (list) {
    item = _malloc(DS_SList);
    item->next = NULL;
    if (f)
      item->data = f(list->data);
    else {
      if (size_data) {
	item->data = malloc(size_data);
	memcpy(item->data, list->data, size_data);
      }
      else
	item->data = list->data;
    }
  }

  return item;
}

DS_SList *ds_slist_clone(DS_SList *list, size_t size_data)
{
  DS_SList *new_list = NULL;

  if (list) {
    DS_SList *last;

    new_list = ds_slist_clone_item(list, size_data, NULL);
    last = new_list;
    list = list->next;
    while(list) {
      last->next = ds_slist_clone_item(list, size_data, NULL);
      last = last->next;
      list = list->next;
    }
  }

  return new_list;
}

DS_SList *ds_slist_clone_with_function(DS_SList *list, size_t size_data,
				       void *(*f)(void *))
{
  DS_SList *new_list = NULL;

  if (list) {
    DS_SList *last;

    new_list = ds_slist_clone_item(list, size_data, f);
    last = new_list;
    list = list->next;
    while (list) {
      last->next = ds_slist_clone_with_function(list, size_data, f);
      last = last->next;
      list = list->next;
    }
    last->next = NULL;
  }

  return new_list;
}

DS_SList *ds_slist_concat(DS_SList *list1, DS_SList *list2)
{
  if (list2) {
    if (list1)
      ds_slist_last(list1)->next = list2;
    else
      list1 = list2;
  }
  return list1;
}

DS_SList *ds_slist_reverse(DS_SList *list)
{
  DS_SList *tmp;
  DS_SList *prev = NULL;
  DS_SList *last = NULL;

  while (list) {
    last = list;
    tmp = list->next;
    list->next = prev;

    prev = list;
    list = tmp;
  }

  return last;
}

DS_SList *ds_slist_find(DS_SList *list, void *obj)
{
  while (list) {
    if (list->data == obj) break;
    list = list->next;
  }

  return list;
}

DS_SList *ds_slist_remove(DS_SList *list, void *obj)
{
  DS_SList *prev, *tmp;

  prev = NULL;
  tmp = list;

  while (tmp) {
    if (tmp->data == obj) {
      if (prev)
	prev->next = tmp->next;
      else
	list = list->next;
      tmp->next = NULL;
      free(tmp);
      break;
    }
    prev = tmp;
    tmp = tmp->next;
  }

  return list;
}

DS_SList *ds_slist_remove_node(DS_SList *list, DS_SList *node)
{
  DS_SList *prev, *tmp;

  prev = NULL;
  tmp = list;

  while (tmp) {
    if (tmp == node) {
      if (prev)
	prev->next = tmp->next;
      else
	list = list->next;
      tmp->next = NULL;
      break;
    }
    prev = tmp;
    tmp = tmp->next;
  }

  return list;
}

void **ds_slist_to_array(DS_SList *list)
{
  void **array = NULL;

  if (list) {
    int i = 0;
    DS_SList *node = list;
    array = _mallocs(ds_slist_length(list), void *);
    while (node) {
      array[i++] = node->data;
      node = ds_slist_next(node);
    }
  }

  return array;
}

/**
 * Convertit la liste list en tableau, en copiant size bytes depuis
 * le maillon de la liste dans la case du tableau correspondante.
 *
 * The data are duplicated.
 *
 * @param size is the size in bytes of one element.
 */
void *ds_slist_to_array_with_size(DS_SList *list, size_t size)
{
  void *array = NULL;

  if (list) {
    int i = 0;
    DS_SList *node = list;
    array = malloc(ds_slist_length(list)*size);
    while (node) {
      memcpy(array+(i++)*size, node->data, size);
      node = ds_slist_next(node);
    }
  }

  return array;
}

void ds_slist_foreach(DS_SList *list, Dfct ff)
{
  while (list) {
    if (ff) ff(list->data);
    list = list->next;
  }
}

DS_SList *ds_slist_insert_sorted(DS_SList *list, void *data,
				 DSCompareFunc f)
{
  DS_SList *tmp_list = list;
  DS_SList *new_list, *prev = NULL;

  new_list = ds_slist_new_with_obj(data);

  if (!list)
    return new_list;

  while (f(data, tmp_list->data) > 0) {
    prev = tmp_list;
    tmp_list = tmp_list->next;
    if (!tmp_list) break;
  }

  new_list->next = tmp_list;
  if (prev) {
    prev->next = new_list;
    return list;
  }
  else
    return new_list;
}
