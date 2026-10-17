/**
 * event.cpp
 *
 * ORATIO Library - Simulation module
 * Copyright (C) 2010-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License
 */

#include "oratio/sim.hpp"

static int event_id = 0;

int ds_compare_events(const void *a, const void *b)
{
  const double time_a = ((ESim_Event *)a)->time;
  const double time_b = ((ESim_Event *)b)->time;

  const int priority_a = ((ESim_Event *)a)->priority;
  const int priority_b = ((ESim_Event *)b)->priority;

  int res = ds_search_double(&time_a, &time_b);

  // if ds_search_double returns 0 (time_a = time_b), this function returns
  // the opposite of the difference between priorities a & b.
  if (!res)
    return ds_search_int_opposite(&priority_a, &priority_b);

  return res;
}

ESim_Event::ESim_Event(ESim *sim)
{
  this->sim = sim;
  time = 0.0;
  priority = 0;
  data = NULL;
  this->id = event_id++;
}

ESim_Event::~ESim_Event()
{}

void ESim_Event::schedule(double t, int prior)
{
  time = sim->clock+t;
  priority = prior;
  sim->events = ds_list_insert_sorted_last(sim->events, this,
					   ds_compare_events);  
}

/*
reschedule

cancel

we simply use the real order to sort the event, as the simulation clock can be described by reals.
*/

/**
 * Free memory allocated to an event.
 *
 * @param ev is a pointer to a memory-allocated event. 
 */
void esim_event_free(void *ev)
{
  if (ev) free(ev);
}
