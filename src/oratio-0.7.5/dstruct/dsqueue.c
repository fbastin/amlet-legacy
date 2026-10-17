/**
 * dsqueue.c
 *
 * D Struct
 * Copyright (C) 2000-2005 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <oratio/ctools.h>
#include "oratio/dstruct.h"

DS_Queue *ds_queue_new()
{
  DS_Queue *queue;

  queue = _malloc(DS_Queue);
  queue->front = queue->rear = NULL;

  return queue;
}

void ds_queue_free(DS_Queue *queue, Dfct ff)
{
  if (queue) {
    ds_slist_free(queue->front, ff);
    free(queue);
  }
}

void ds_queue_clear(DS_Queue *queue, Dfct ff)
{
  if (queue) {
    ds_slist_free(queue->front, ff);
    queue->front = queue->rear = NULL;
  }
}

DS_SList *ds_queue_front(DS_Queue *queue)
{
  return queue->front;
}

int ds_queue_is_empty(DS_Queue *queue)
{
  return (queue->front == NULL);
}

int ds_queue_is_full(DS_Queue *queue)
{
  return 0;
}

void *ds_queue_dequeue(DS_Queue *queue)
{
  void *data;

  if (ds_queue_is_empty(queue))
    return NULL;
  else {
    data = queue->front->data;
    queue->front = ds_slist_pop(queue->front);
    if (!queue->front) queue->rear = NULL;
    return data;
  }
}

void ds_queue_enqueue(DS_Queue *queue, void *data)
{
  queue->rear = ds_slist_append(queue->rear, data);
  if (queue->front)
    queue->rear = ds_slist_next(queue->rear);
  else
    queue->front = queue->rear;
}

void *ds_queue_first(DS_Queue *queue)
{
  return (queue->front) ? queue->front->data : NULL;
}

void *ds_queue_end(DS_Queue *queue)
{
  return (queue->rear) ? queue->rear->data : NULL;
}

/*****************************************************************************/

void ds_queue_foreach(DS_Queue *queue, Dfct ff)
{
  DS_SList *list;

  list = queue->front;
  ds_slist_foreach(list, ff);
}

/*****************************************************************************/

/* untested function */
DS_Queue *queue_insert_sorted(DS_Queue *queue, void *data, DSCompareFunc f)
{
  if (!queue) {
    queue = ds_queue_new();
    ds_queue_enqueue(queue, data);
  }
  else {
    queue->front = ds_slist_insert_sorted(queue->front, data, f);
    if (queue->rear->next) queue->rear = ds_slist_next(queue->rear);
  }

  return queue;
}
