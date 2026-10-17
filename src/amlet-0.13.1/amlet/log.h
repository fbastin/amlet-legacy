/**
 * log.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __AMLET_LOG_H__
#define __AMLET_LOG_H__

#include <oratio/ntools.h>

#ifdef __cplusplus
extern "C" {
#endif

  extern NTLog *mclog;

  void amlet_log(char *msg, ...);

#ifdef __cplusplus
}
#endif

#endif /* __AMLET_LOG_H__ */
