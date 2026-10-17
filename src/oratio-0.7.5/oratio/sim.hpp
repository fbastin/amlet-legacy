/**
 * sim.hpp
 *
 * ORATIO Library
 * Copyright (C) 2010-2013 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License
 */

#ifndef __ORATIO_SIM_HPP__
#define __ORATIO_SIM_HPP__

#include <oratio/dstruct.h>
  
/**
 * Discrete events simulation
 */
class ESim {
public:
  void *system;    /* system under consideration */
  double clock;    /* simulation clock */
  DS_List *events; /* events list */
  int status;      /* indicate the simulation status */

  ESim();
  ~ESim();

  double time() {return clock;};

  void init();
  void start();
  void stop();
};

enum {
  ENDED = 0,
  PAUSED,
  RUNNING
};

/**
 * Create a new discrete events simulation instance, and initialize it.
 *
 * @return the newly created discrete events simulation instance.        
 */
ESim *esim_new();
  
/**
 * Reset a discrete events simulation instance.
 *
 * @param esim is the pointer to a discrete events simulation instance.  
 */
void esim_reset(ESim *sim);
  
/**
 * Event prototype.
 *
 * @param system is a pointer to the system under consideration.
 *
 * @return 0 if the event has been correctly handled, some error code
 *   otherwise.         
 */
struct ESim_Event {
  ESim *sim;
  double time;   /* time when the event has to occur. */
  int priority;  /* priority of the event, default 0. */
  void *data;

public:
  int id;

  ESim_Event(ESim *sim);
  virtual ~ESim_Event();

  virtual void actions() {};

  void schedule(double time, int priority = 0);
};

ESim_Event *esim_get_current_event(ESim *sim);

/**
 * Free memory allocated to an event.
 *
 * @param ev is a pointer to a memory-allocated event. 
 */
void esim_event_free(void *ev);

/**
 * Get the data associated to the current event.
 *
 * @param esim is the pointer to a discrete events simulation instance.
 *
 * @return a pointer to the data of the current event.
 */
void *esim_get_current_event_data(ESim *sim);

#endif  /* __ORATIO_SIM_HPP__ */
