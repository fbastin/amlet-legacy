/**
 * dsconvert.h
 *
 * D Struct
 * Copyright (C) 2000-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <oratio/dstruct.h>

/**
 * Modify a queue so that it contains the same elements as some
 * simply linked chained list.
 *
 * @param list is initialized list
 * @param queue is an initialized queue
 *
 * #POST: the queue 'queue' contains the same elements as the list 'list';
 *        these elements are shared between both the structures.
 *        The beginning of the queue corresponds to the beginning of
 *        the list, who is unchanged.
 */
void ds_slist_to_queue(DS_SList *list, DS_Queue *queue)
{
  if (!ds_queue_is_empty(queue)) ds_queue_clear(queue, 0);

  queue->front = list;
  if (list) while (list->next) list = list->next;

  queue->rear = list;
}
