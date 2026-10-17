/**
 * mixedlogit.c
 *
 * AMLET Library
 * Copyright (C) 2015 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the BSD License.
 **/

#include <math.h>

#include "mixedlogit.h"

/**
 * Compute the simulation bias.
 *
 * Ref.: F. Bastin, C. Cirillo, and Ph. L. Toint, An adaptive Monte Carlo
 * algorithm for computing mixed logit estimators, Computational Management
 * Science, 3(1):55--79 (2006).
 * Special thanks to Adrien Bonneu, INSA, Toulouse. First bias correction.
 *
 * @param popsize is the population size
 * @param accuracy is the estimated precision.
 * @param alpha is the significance level used for the accuracy estimation.
 *
 * @return the bias of simulation
 */
double simulation_bias(int popsize, double alpha, double accuracy)
{
  double t = accuracy/alpha;
  if (isnan(t)) t = 0.0;
  return -0.5*popsize*t*t;
}
