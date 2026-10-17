/**
 * dsconvert.h
 * 
 * ORATIO Library
 * Copyright (C) 2000-2002 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __DSTRUCT_CONVERT_H__
#define __DSTRUCT_CONVERT_H__

#include <dsslist.h>
#include <dsqueue.h>

#ifdef __cplusplus
extern "C" {
#endif  /* __cplusplus */

  /**
   * Modify a queue so that it contains the same elements as some
   * simply linked chained list.
   *
   * #PRE:  initialized list and queue
   *
   * #POST: the queue 'queue' contains the same elements as the list 'list';
   *        these elements are shared between both the structures.
   *        The beginning of the queue corresponds to the beginning of
   *        the list, who is unchanged.
   */
  void ds_slist_to_queue(DS_SList *list, DS_Queue *queue);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __DSTRUCT_CONVERT_H__ */
