/**
 * amlet_check.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2005 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the MIT License
 */

#include <stdio.h>
#include <stdlib.h>

#include <amlet/amlet.h>

#include "src/help.h"

int main()
{
  printf("AMLET tests set.\n");

  amlet_version();
  /*
  system("amlet -v -s driver_australia_mnl");
  system("amlet -v -s driver_australia_mmnl");

  system("amlet -v -s driver_simulated_panel");*/

  system("amlet -v driver_simulated_errc");

  return 0;
}
