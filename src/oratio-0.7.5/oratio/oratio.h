/**
 * oratio.h
 *
 * "Operational Research, A Toolbox in Operation"
 * Copyright (C) 2008-2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_H__
#define __ORATIO_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Print library credits on standard output.
   */
  void oratio_credits();

  /**
   * Headers of the various components.
   */
#include "ctools.h"
#include "dstruct.h"
#include "moca.h"
#include "ntools.h"
#include "ran.h"
#include "stat.h"

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_H__ */
