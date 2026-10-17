/**
 * individual.c
 *
 * AMLET Library
 * Copyright (C) 2004-2015 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the BSD License.
 **/

#include "individual.h"

#include <oratio/ctools.h>

AMLET_Ind *amlet_ind_new()
{
  AMLET_Ind *i = _malloc(AMLET_Ind);

  //  i->prob_choice = i->grad_prob_choice = i->alternatives = NULL;

  return i;
}

void amlet_ind_free(AMLET_Ind *i)
{
  if (i) {
    /*    if (i->prob_choice) free(i->prob_choice);
	  if (i->grad_prob_choice) free(i->grad_prob_choice); */
    //    if (i->alternatives) free(i->alternatives);
    free(i);
  }
}
