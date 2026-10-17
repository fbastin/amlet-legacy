/**
 * dsqueue.h
 *
 * D Struct
 * Copyright (C) 2000-2002 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_DSTRUCT_QUEUE_H__
#define __ORATIO_DSTRUCT_QUEUE_H__

#include "dstruct.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  typedef struct _DS_Queue DS_Queue;

  struct _DS_Queue {
    DS_SList *front, *rear;
  };

  DS_Queue *ds_queue_new();

  void ds_queue_free(DS_Queue *queue, Dfct ff);

  void ds_queue_clear(DS_Queue *queue, Dfct ff);

  DS_SList *ds_queue_front(DS_Queue *queue);

  int ds_queue_is_empty(DS_Queue *queue);

  int ds_queue_is_full(DS_Queue *queue);

  void *ds_queue_dequeue(DS_Queue *queue);

  void ds_queue_enqueue(DS_Queue *queue, void *data);

  /**
   * Return the first element of a queue, or NULL if the queue is empty.
   */
  void *ds_queue_first(DS_Queue *queue);

  /**
   * Renvoie le dernier élément d'une queue, ou NULL si la queue est vide
   */
  void *ds_queue_end(DS_Queue *queue);

  /**
   * Applique la fonction ff à chaque élément de la queue queue.
   */
  void ds_queue_foreach(DS_Queue *queue, Dfct ff);

  /**
   * Insert an element inside the queue, while respecting some order
   * in the queue.
   */
  DS_Queue *queue_insert_sorted(DS_Queue *queue, void *data, DSCompareFunc f);
  /* !!! untested function */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __ORATIO_DSTRUCT_QUEUE_H__ */
