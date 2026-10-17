/**
 * test_sim.cpp
 *
 * ORATIO Library
 * Copyright (C) 2009 <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License.
 */

#include <oratio/sim.hpp>

#include <math.h>

class basic: public ESim_Event {
public:
  basic(ESim *sim) : ESim_Event(sim) {printf("Hallo\n");};

  void action() {
    printf("Clock: %g\n", sim->clock);
  };
};

int main()
{
  ESim *sim;
  ESim_Event *ev;

  sim = new ESim();

  printf("Salut 1\n");
  ev = new basic(sim);
  printf("Salut 2\n");
  ev->schedule(2.0);
  printf("Salut 3\n");
  ev = new basic(sim);
  printf("Salut 4\n");
  ev->schedule(1.0);
  printf("Salut 100\n");

  sim->start();

  delete(sim);

  return 0;
}
