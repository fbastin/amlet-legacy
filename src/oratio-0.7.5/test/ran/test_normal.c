/**
 * test_normal.c
 *
 * Random Numbers Library
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License.
 */

#include "test_ran.h"

#include <stdio.h>

/**
 * Test of normal distribution generation routines.
 */
void test_normal()
{
  FILE *f;
  Random *normal;
  int classes[52];
  double d, size;
  int i,j;
  unsigned long state = 0;

  int nclasses = 52;

  f = fopen("normal.txt", "w");

  for (i = 0; i < nclasses; i++) classes[i] = 0.0;
  normal = ran_normal_new(&state, 1, -5.0, 2);

  size = 20.0/50.0;

  for (i = 0; i < 500000; i++) {
    d = ran_normal_get_val(normal);
    if (d < -15.0) classes[0]++;
    else if (d >= 5.0) classes[51]++;
    else {
      for (j = 0; j < 50; j++) {
	if (d >= (-15.0 + j*size) && d < (-15.0 +(j+1)*size)) {
	  classes[j+1]++;
	  break;
	}
      }
    }
  }

  for (i = 0; i < nclasses; i++) fprintf(f, "%d %d\n", i, classes[i]);

  fclose(f);

  f = fopen("normal_bm.txt", "w");

  for (i = 0; i < nclasses; i++) classes[i] = 0.0;

  size = 20.0/50.0;

  for (i = 0; i < 500000; i++) {
    d = ran_gaussian_box_muller(normal);
    if (d < -10.0) classes[0]++;
    else if (d >= 10.0) classes[51]++;
    else {
      for (j = 0; j < 50; j++) {
	if (d >= (-10.0 + j*size) && d < (-10.0 +(j+1)*size)) {
	  classes[j+1]++;
	  break;
	}
      }
    }
  }

  for (i = 0; i < nclasses; i++) fprintf(f, "%d %d\n", i, classes[i]);

  fclose(f);

  ran_random_free(normal);
}

