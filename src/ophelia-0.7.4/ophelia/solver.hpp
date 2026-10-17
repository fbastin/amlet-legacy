/**
 * Solver.hpp
 *
 * OPHELIA
 * Copyright (C) 2015 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the MIT License
 */

#include <string>
using namespace std;

#ifndef __OPHELIA_SOLVER_HPP__
#define __OPHELIA_SOLVER_HPP__

class OPSolver {
  bool verbose = false;

  string *name;

  double *sol;
  double value;

public:
  // Solve the optimization problem
  virtual int solve();

  OPSolver();
  ~OPSolver();
};

#endif  /* __OPHELIA_SOLVER_HPP__ */
