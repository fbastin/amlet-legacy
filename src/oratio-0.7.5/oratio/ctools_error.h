/**
 * ctools_error.h
 *
 * C Tools
 * Copyright (C) 2000-2002 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ATHOR_CTOOLS_ERROR_H__
#define __ATHOR_CTOOLS_ERROR_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  void c_warning(char *msg, ...);
  void c_error(char *msg, ...);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __ATHOR_CTOOLS_ERROR_H__ */
