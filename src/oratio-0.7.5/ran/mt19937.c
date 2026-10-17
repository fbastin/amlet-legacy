/**
 * mt19937.c
 *
 * ORATIO Library
 * Copyright (C) 2004-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "mt19937.h"
#include "oratio/ran.h"

#include <oratio/ctools.h>

#include <string.h>

double ran_mt19937_get_val(void *mt)
{
  /* TO IMPROVE... */
  double u;

  do {
    u = ran_mt19937_genrand_real3(mt);
  } while (u == 1.0);

  return u;
}

/**
 * Initialize the Mersenne twister generator with a specified seed vector.
 *
 * @param seed is the seed vector.
 * @param n is the dimension of seed.
 */
void ran_mt19937_init(MT19937 *mt, unsigned long *seed, int n)
{
  n == 1 ? ran_mt19937_init_genrand(mt, *seed) :
    ran_mt19937_init_by_array(mt, seed, n);
}

/**
 * UNTESTED FUNCTION!
 */
MT19937 *ran_mt19937_new_with_state(unsigned long *state)
{
  MT19937 *mt;

  mt = _malloc(MT19937);

  memcpy(mt->mt, state, MT19937_STATE_LENGTH*sizeof(unsigned long));

  return mt;
}

/**
 * Retrieve the generator state.
 */
unsigned long *ran_mt19937_get_state(void *mt)
{
  return ((MT19937 *)mt)->mt;
}

/**
 * UNTESTED FUNCTION!
 */
void ran_mt19937_set_state(MT19937 *mt, void *state)
{
  memcpy(mt->mt, state, MT19937_STATE_LENGTH*sizeof(unsigned long));
} 

/**
 * @param seed is the seed vector.
 * @param n is the dimension of seed.
 */
MT19937 *ran_mt19937_new(unsigned long *seed, int n)
{ 
  MT19937 *mt;

  mt = _malloc(MT19937);

  if (seed)
    ran_mt19937_init(mt, seed, n);
  else {
    unsigned long s;

    s = 12345;
    ran_mt19937_init(mt, &s, 1);
  }

  return mt;
}

void ran_mt19937_free(void *mt)
{
  free((MT19937 *)mt);
}
