/**
 * ctimer.c
 *
 * ORATIO Library - CTools module
 * Copyright (C) 2002-2015 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <config.h>

#include <errno.h>
#include <string.h>

#include "oratio/ctools.h"

#if defined(__CYGWIN__) || defined(_WIN32) || defined(__MINGW32__)
#define MAX_TIMERVAL  100000
#define ITIMER_AMLET  ITIMER_REAL  /* windows has no other itimer :( */
  clock_t start, end;
  double cpu_time;
#else

  struct itimerval old_t, new_t;
#ifndef timersub
#define timersub(a, b, result)						      \
  do {									      \
    (result)->tv_sec = (a)->tv_sec - (b)->tv_sec;			      \
    (result)->tv_usec = (a)->tv_usec - (b)->tv_usec;			      \
    if ((result)->tv_usec < 0) {					      \
      --(result)->tv_sec;						      \
      (result)->tv_usec += 1000000;					      \
    }									      \
  } while (0)
#endif 

#endif

/**
 * Initialize the internal timer.
 */
void ctools_timer_init()
{
#if defined(__CYGWIN__) || defined(_WIN32) || defined(__MINGW32__)
  /* We can easily get some overflow with the clock function due to
     the time needed by the program */
  start = clock();
#else
  /* We record the clock in order to get the ellapsed computation time */
  int r;

  /* Set up interval timer; we start the timer now */
  old_t.it_value.tv_sec = MAX_TIMERVAL;
  old_t.it_value.tv_usec = 0;
  old_t.it_interval.tv_sec = old_t.it_interval.tv_usec = 0;

  r = setitimer(ITIMER_AMLET, &old_t, NULL);

  /* adjust the value of old_t */
  getitimer(ITIMER_AMLET, &new_t);
  old_t.it_value.tv_sec = new_t.it_value.tv_sec;
  old_t.it_value.tv_usec = new_t.it_value.tv_usec;

  if (r == EINVAL || r == EFAULT)
    c_warning("Unable to set interval timer.\n");
#endif
}

/**
 * Return the time used (in seconds) since the timer initialization.
 *
 * @param output is a pointer to the file where the time has to recorded,
 *    or NULL if it does not to be recorded.
 * @param sec contains on output the number of seconds.
 * @param ms  contains on output the number of millseconds.
 *
 * @return the elapsed time in seconds.
 */
long ctools_timer(FILE *output, long *sec, long *ms)
{
#if defined(__CYGWIN__) || defined(_WIN32) || defined(__MINGW32__)
  end = clock();
  cpu_time = ((double)(end-start))/CLOCKS_PER_SEC;

  *sec = cpu_time;
  *ms = 0.0;

  if (output)
    fprintf(output, "Time spent: %ld seconds\n", *sec);

  return *sec;
#else
  struct itimerval tmp;

  memset(&tmp, 0, sizeof(tmp));
  getitimer(ITIMER_AMLET, &new_t);
  timersub(&(old_t.it_value), &(new_t.it_value), &(tmp.it_value));

  if (output)
    fprintf(output, "Time spent: %ld.%ld seconds\n",
	    tmp.it_value.tv_sec, (tmp.it_value.tv_usec/1000));

  *sec = tmp.it_value.tv_sec;
  *ms = tmp.it_value.tv_usec/1000;

  return *sec + ((tmp.it_value.tv_usec >= 500000) ? 1 : 0);
#endif
}
