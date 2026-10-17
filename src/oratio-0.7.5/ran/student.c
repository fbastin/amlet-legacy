/**
 * student.c
 *
 * Random Numbers Library
 * Copyright (C) 2007-2008 Adrien Bonneu <abonneu@etud.insa-toulouse.fr>
 *                         Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Random number generation for student distributions
 */

#include "oratio/ran.h"
#include "oratio/ntools.h"

#include <math.h>
#include <string.h>

/**
 * Student Distribution
 */

/**
 * Create a new generator object dedicated to Student distribution.
 *
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 * @param n is the dimension of the seed vector.
 * @param nu is the number of degrees of freedom of the Student distribution.
 *
 * @return the newly allocated generator object.
 */
Random *ran_student_new(unsigned long *seed, int n,
			double nu)
{
  Random *ran=ran_random_new_with_size(seed, n, sizeof(double));

  memcpy(ran->parameters, &nu, sizeof(double));

  ran->type = RAN_STUDENT;

  return ran;
}

/**
 * Generate a pseudo-random from a Student distribution with the parameter nu.
 *
 * @param nu is the number of degrees of freedom of the Student distribution.
 *
 * @return a number draw from a Student with the number of degrees of
 *   freedom equal to nu.
 */
double ran_student(const Random *random, double nu)
{
  double x, z, y = 0;
  int i;

  for (i=0;i<nu;i++) {
    x = ran_gaussian_icdf(random);
    y += (x*x-y)/(i+1);
  }
  z=sqrt(y);

  x = ran_gaussian_icdf(random);

  return x/z;
}

/**
 * Generate a pseudo-random from a Student distribution using a
 * specialized generator.
 *
 * @return a number draw from a Student with the number of degrees of
 *   freedom equal to nu.
 */
double ran_student_get_val(const Random *random)
{
  return ran_student(random, *((double *)(random->parameters)));
}


/**
 * Get the number of degrees of freedoms of a Student random variate
 * generator.
 *
 * @param random is an initialized random numbers generator, dedicated
 *   to the generation from a Student random variate with degrees of
 *   freedoms equal to nu.
 * @param nu is on output the number of degrees of freedom of the
 *   Student distribution.
 */
void ran_student_get_parameters(const Random *random, double *nu)
{
  *nu = *((double *)(random->parameters));
}

/**
 * Set the number of degrees of freedom of a Student random variate
 * generator.
 *
 * @param random is an initialized random numbers generator.
 * @param nu is the number of degrees of freedom of the Student distribution
 */
void ran_student_set_parameters(Random *random, double nu)
{
  ran_random_set_parameters(random, &nu, sizeof(double));
}
