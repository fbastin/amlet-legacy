#include <stdio.h>

#include "btr.hpp";

double OP_BTR::ratio()
{
  return 10.0;
}

void OP_BTR::update_radius()
{
  printf("Hello %g\n", ratio());
}
