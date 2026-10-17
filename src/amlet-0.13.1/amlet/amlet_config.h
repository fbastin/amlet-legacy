/**
 * amlet_config.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __AMLET_CONFIG_H__
#define __AMLET_CONFIG_H__

#include <oratio/oratio.h>

#include <amlet/nest.h>
#include <amlet/dcl.h>

#ifdef __cplusplus
extern "C" {
#endif

  extern int loose;       /* 1 if we are not strict on the syntax of driver files,
		      0 else */

  extern NTLog *mclog;
  extern int nfctev;

  extern dcl_data data;

  extern unsigned int twopowers[31];

  typedef struct {
    char *logname;    /* name of the log file */
  } amlet_config;

#ifdef __cplusplus
}
#endif

#endif
