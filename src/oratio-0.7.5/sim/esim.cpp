/**
 * esim.c
 *
 * ORATIO Library - Simulation module
 * Copyright (C) 2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License
 */

#include "oratio/sim.hpp"
#include "oratio/memory.h"

ESim::ESim()
{
  init();

  return;
}

ESim::~ESim()
{
  ds_list_free(events, esim_event_free);
}

/**
 * Initialize a discrete events simulation instance.
 */
void ESim::init()
{
  clock = 0.0;
  events = ds_list_new();
  status = PAUSED;
}

void ESim::stop()
{
  ds_list_free(events, esim_event_free);
  events = NULL;
  status = ENDED;
}

/**
 * Launch a discrete events simulation.
 *     
 * @param esim is the pointer to a discrete events simulation instance.  
 */  
void ESim::start()
{
  ESim_Event *ev;

  status = RUNNING;
  while (events && status == RUNNING) {
    /* Generate the current event. */
    ev = esim_get_current_event(this);
    clock = ev->time;
    ev->actions();
    /* Since the event has been treated, we free it */
    if (status == RUNNING) {
      delete (ev);
      /* Go to the next event. */ /* Chercher l'ID */
      events = ds_list_pop(events);
    }
  }
}

/** The code below is kept for compatibility issues */

void esim_init(ESim *sim)
{
  sim->clock = 0.0;
  sim->events = ds_list_new();
  sim->status = PAUSED;
}

/**
 * Create a new discrete events simulation instance, and initialize it.
 *
 * @return the newly created discrete events simulation instance.        
 */
ESim *esim_new()
{
  ESim *sim = _malloc(ESim);
  
  esim_init(sim);
  
  return sim;
}

/**
 * Reset a discrete events simulation instance.
 *
 * @param esim is the pointer to a discrete events simulation instance.  
 */
 void esim_reset(ESim *sim)
 {
   ds_list_free(sim->events, esim_event_free);
   sim->clock = 0.0;
 }

/**
 * @param esim is the pointer to a discrete events simulation instance.  
 */
ESim_Event *esim_get_current_event(ESim *sim)
{
  return (ESim_Event *)(sim->events->data);
}

/**
 * Get the data associated to the current event.
 *
 * @param esim is the pointer to a discrete events simulation instance.
 *
 * @return a pointer to the data of the current event.
 */
void *esim_get_current_event_data(ESim *sim)
{
  return ((esim_get_current_event(sim))->data);
}
