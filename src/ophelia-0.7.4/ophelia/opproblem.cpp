#include "problem.hpp"

OPProblem::OPProblem()
{
  name = NULL;

  lo = up = NULL;
}

OPProblem::OPProblem(string *name)
{
  this->name = name;

  lo = up = NULL;
}

OPProblem::~OPProblem()
{
}
