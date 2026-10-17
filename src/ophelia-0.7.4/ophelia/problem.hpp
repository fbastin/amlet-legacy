/**
 * problem.hpp
 *
 * OPHELIA Library
 * Copyright (C) 2005-2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the MIT License
 */

#include <string>
using namespace std;

#ifndef __OPHELIA_PROBLEM_HPP__
#define __OPHELIA_PROBLEM_HPP__

class OPProblem {
  string *name;

  int n;       // problem dimension
  double *lo;  // lower bounds
  double *up;  // upper bounds

public:
  double evaluate();
  double gradient();
  double hessian();

  bool is_feasible();

  OPProblem();
  OPProblem(string *name);
  ~OPProblem();
};

#endif  /* __OPHELIA_PROBLEM_HPP__ */
