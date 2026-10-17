/**
 * ctime.h
 *
 * ORATIO Library - CTools module
 * Copyright (C) 2005 Fabian Bastin <slashbin@slashbin.be>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_CTOOLS_TIME_H__
#define __ORATIO_CTOOLS_TIME_H__

#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Read a date written in standard form: YYYY/MM/DD HH:MM:SS.
   *
   * @param d is a string containing a date written in standard form.
   */
  time_t ctools_read_european_date_time(char *d);

  /**
   * Read a date written in standard form: YYYY/MM/DD HH:MM:SS.
   */
  time_t ctools_read_standard_date_time(char *d);

  /**
   * Read a time written in standard form: HH:MM:SS.
   *
   * @param d is a string containing a date written in standard form.
   */
  time_t ctools_read_standard_time(char *d, struct tm *t);

  /**
   * Read a date written in american form: MM/DD/YYYY HH:MM:SS.
   */
  time_t ctools_read_american_date_time(char *d);

  /**
   * Create a stamp based on the current time.
   */
  char *ctools_timestamp(char *buff);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_CTOOLS_TIME_H__ */

