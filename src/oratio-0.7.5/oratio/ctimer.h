/**
 * ctimer.h
 *
 * ORATIO Library - CTools module
 * Copyright (C) 2005-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __CTOOLS_TIMER_H__
#define __CTOOLS_TIMER_H__

#include <time.h>
#include <sys/time.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct _C_Timer C_Timer;

#define MAX_TIMERVAL  100000
#define ITIMER_AMLET  ITIMER_REAL

  struct _C_Timer {
#if defined(__CYGWIN__) || defined(_WIN32) || defined(__MINGW32__)
    clock_t start, end;
    double cpu_time;
#else
    
    struct itimerval old_t, new_t;
#ifndef timersub
#define timersub(a, b, result)						\
    do {								\
      (result)->tv_sec = (a)->tv_sec - (b)->tv_sec;			\
      (result)->tv_usec = (a)->tv_usec - (b)->tv_usec;			\
      if ((result)->tv_usec < 0) {					\
	--(result)->tv_sec;						\
	(result)->tv_usec += 1000000;					\
      }									\
    } while (0)
#endif 
    
#endif
  };

  /**
   * Initialize the internal timer.
   */
  void ctools_timer_init();

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
  long ctools_timer(FILE *output, long *sec, long *ms);

#ifdef __cplusplus
}
#endif

#endif  /* __CTOOLS_TIMER_H__ */
