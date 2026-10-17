/**
 * accumulate.cpp
 *
 * ORATIO Library - Simulation module
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License
 */

#include <string.h>

#include "oratio/sim_accumulate.hpp"

#include "oratio/ctools.h"
#include "oratio/ntools.h"

/**
 * Initializes the statistical collector and puts the current value of
 * the corresponding variable to 0.
 */
void ESim_Accumulate::init()
{
  max_value = -nt_infinity();;
  min_value = nt_infinity();
  last_value = 0.0;
  sum_value = sum_compensation = 0.0;
  init_time = last_time = sim->time();
}

void ESim_Accumulate::report(FILE *output)
{
  update(last_value);
  fprintf(output, "Report on Accumulate collector");
  if (name) fprintf(output, " %s", name);
  fprintf(output, "\n");
  fprintf(output, "Simulation time: %g to %g\n", init_time, last_time);
  fprintf(output, "Minimum value: %g\tMaximum value: %g\n",
	  min_value, max_value);
  fprintf(output, "Average: %g\n", average());
}

/**
 * Update the information stored by the accumulate collector.
 */
void ESim_Accumulate::update(double x)
{
  double y, t;

  double time = sim->time();
  if (x < min_value) min_value = x;
  if (x > max_value) max_value = x;

  /* We use Kahan's method to update the sum */
  y = last_value*(time-last_time) - sum_compensation;
  t = sum_value + y;
  sum_compensation = (t-sum_value)-y;
  sum_value = t;

  last_value = x;
  last_time = time;
}

/**
 * Set the name of the accumulate.
 *
 * @param name is the name to associated to the accumulate; a copy of the
 *   string will be created.
 */
void ESim_Accumulate::set_name(char *name)
{
  if (this->name) free(this->name);
  this->name = _mallocs(strlen(name)+1, char);
  memcpy(this->name, name, strlen(name)*sizeof(char));
}
