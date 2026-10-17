/**
 * esim_accumulate.hpp
 *
 * ORATIO Library - Simulation module
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License
 */

#include <stdio.h>

#include <oratio/ntools.h>
#include <oratio/sim.hpp>

class ESim_Accumulate
{ 
  double init_time;    // Initialization time.
  double last_time;    // Last update time.
  double last_value;   // Value since last update.
		   
  double max_value;
  double min_value;
  double sum_value;
  double sum_compensation;  // compensation of low order bits during update  

  ESim *sim;

  char *name;

public:

  ESim_Accumulate() {}

  ESim_Accumulate(ESim *sim) {
    this->sim = sim;
    this->name = NULL;
    init();
  }

  ~ESim_Accumulate() {if (name) free(name);};

  /**
   * Initializes the statistical collector and puts the current value of
   * the corresponding variable to 0.
   */
  void init();

  /**
   * Update the information stored by the accumulate collector.
   */
  void update(double x);

  /**
   * Returns the time-average since the last initialization.
   */
  double average()  {
    double period = last_time - init_time;
    update(last_value);
    return (period > 0.0) ? sum_value/period : 0.0;
  }

  void report(FILE *output);

  /**
   * Get the name of the accumulate collector.
   *
   * @return un pointer to the accumulate name.
   */
  char *get_name() {return name;}

  /**
   * Set the name of the accumulate.
   *
   * @param name is the name to associated to the accumulate; a copy of the
   *   string will be created.
   */
  void set_name(char *name);
};
