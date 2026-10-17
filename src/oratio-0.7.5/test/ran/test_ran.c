/**
 * test_ran.c
 *
 * Random Numbers Library
 * Copyright (C) 1999-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <oratio/oratio.h>

#include "cauchy_test.h"
#include "test_ran.h"
#include "test_m_unit_sphere.h"
#include "test_ran_sp.h"
#include "uniforms_test.h"

#include <stdio.h>
#include <stdlib.h>

#define SIZE 20000

unsigned long state;

void Bernouilli_test()
{
  Random *b;
  int prop[2];
  int i;

  prop[0] = prop[1] = 0;

  state = 0;
  b = ran_bernouilli_new(&state, 1, 0.6);

  for (i = 0; i < SIZE; i++)
    prop[ran_bernouilli_get_val(b)]++;

  printf("*****************************************\n");
  printf("* Bernouilli test\n");
  printf("*****************************************\n");
  printf("Parameter: %g\n", ran_bernouilli_get_param(b));
  printf("Proportions: %g %g\n",
	 prop[0]/((double)SIZE), prop[1]/((double)SIZE));

  ran_random_free(b);
}

void Gumbel_test()
{
  Random *g;
  double v, w, x;
  int prop[2];
  int i;

  prop[0] = prop[1] = 0;

  state=0;
  g = ran_gumbel_new(&state, 1, -M_EULER, 1.0);

  w = st_gumbel_median(-M_EULER, 1.0);

  printf("*****************************************\n");
  printf("* Gumbel test\n");
  printf("*****************************************\n");

  for (i = 0; i < SIZE; i++) {
    v = ran_gumbel_get_val(g);
    if (v <= w) prop[0]++;
    else prop[1]++;
  }

  printf("Proportions: %g %g\n",
	 prop[0]/((double)SIZE), prop[1]/((double)SIZE));

  ran_random_free(g);
}

void triangular_test()
{
  FILE *f = NULL;
  Random *r;
  int classes[52];
  int i, j;
  unsigned long state = 0;
  double d, size;

  f = fopen("triangular.txt", "w");

  for (i = 0; i < 52; i++) classes[i] = 0.0;
  r = ran_triangular_new(&state, 1, 0., 1., 0.25);

  size = 1.0/50.0;

  for (i = 0; i < 500000; i++) {
    d = ran_triangular_get_val(r);
    if (d < 0.0) classes[0]++;
    else if (d >= 1.0) classes[51]++;
    else {
      for (j = 0; j < 50; j++) {
	if (d >= (0.0 + j*size) && d < (0.0 +(j+1)*size)) {
	  classes[j+1]++;
	  break;
	}
      }
    }
  }

  for (i = 0; i < 52; i++) fprintf(f, "%d %d\n", i, classes[i]);

  fclose(f);

  ran_random_free(r);
}

int main()
{
  Random *r;
  FILE *f = NULL;

  r = ran_random_new(0, &state, 1);
  uniforms_check(f, r);
  test_m_unit_sphere(r);

  Bernouilli_test();
  Gumbel_test();
  
  test_exponential();
  
  test_normal();

  triangular_test();
  test_truncated_normal();

  test_ran_sp();

  test_mvn();

  test_acceptance_rejection();

  ran_random_free(r);

  return 0;
}
