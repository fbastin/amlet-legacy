/**
 * time.c
 *
 * ORATIO Library
 * Copyright (C) 2005-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/memory.h"

#include <stdio.h>
#include <time.h>

#define TIME_SIZE 15

/**
 * Read a date written in standard form: YYYY/MM/DD HH:MM:SS.
 *
 * @param d is a string containing a date written in standard form.
 *
 * @return the corresponding time as a time_t structure.
 */
time_t ctools_read_standard_date_time(char *d)
{
  time_t tt;
  struct tm t;

  t.tm_sec = 0;
  t.tm_wday = t.tm_yday = 0;
  t.tm_isdst = -1;

  sscanf(d, "%d/%d/%d %d:%d:%d", &(t.tm_year), &(t.tm_mon),
	 &(t.tm_mday), &(t.tm_hour), &(t.tm_min), &(t.tm_sec));

  t.tm_year -= 1900;

  tt = mktime(&t);

  return tt;
}

/**
 * Read a time written in standard form: HH:MM:SS.
 *
 * @param d is a string containing a date written in standard form.
 */
time_t ctools_read_standard_time(char *d, struct tm *t)
{
  time_t tt;

  t->tm_wday = t->tm_yday = 0;
  t->tm_isdst = -1;
  t->tm_year = 80;
  t->tm_mon = 1;
  t->tm_mday = 1;

  sscanf(d, "%d:%d:%d", &(t->tm_hour), &(t->tm_min), &(t->tm_sec));

  tt = mktime(t);

  return tt;
}

/**
 * Read a date written in european form: DD/MM/YYYY HH:MM:SS.
 *
 * @param d is a string containing a date written as DD/MM/YYYY HH:MM:SS.
 */
time_t ctools_read_european_date_time(char *d)
{
  time_t tt;
  struct tm t;

  t.tm_sec = 0;
  t.tm_wday = t.tm_yday = 0;
  t.tm_isdst = -1;

  sscanf(d, "%d/%d/%d %d:%d:%d", &(t.tm_mday), &(t.tm_mon),
	 &(t.tm_year), &(t.tm_hour), &(t.tm_min), &(t.tm_sec));

  t.tm_year -= 1900;

  tt = mktime(&t);

  return tt;
}

/**
 * Read a date written in american form: MM/DD/YYYY HH:MM:SS.
 *
 * @param d is a string containing a date written as MM/DD/YYYY HH:MM:SS.
 */
time_t ctools_read_american_date_time(char *d)
{
  time_t tt;
  struct tm t;

  t.tm_sec = 0;
  t.tm_wday = t.tm_yday = 0;
  t.tm_isdst = -1;

  sscanf(d, "%d/%d/%d %d:%d:%d", &(t.tm_mon), &(t.tm_mday),
	 &(t.tm_year), &(t.tm_hour), &(t.tm_min), &(t.tm_sec));

  t.tm_year -= 1900;

  tt = mktime(&t);
 
  return tt;
}

/**
 * Create a stamp based on the current time.
 */
char *ctools_timestamp(char *buff)
{
  struct tm *tm;
  time_t now;

  if (!buff) buff = _mallocs(TIME_SIZE, char);

  now = time(NULL);
  tm = localtime(&now);

  strftime(buff, TIME_SIZE, "%Y%m%d%H%M%S", tm);

  return buff;
}
