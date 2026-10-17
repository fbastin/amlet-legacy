/**
 * m_unit_sphere.c
 *
 * ORATIO Library
 * Copyright (C) 2007 Fabian Bastin <bastin@iro.umontreal>
 *
 * Licensed under the Open Software License version 2.1
 */

#include "oratio/ran.h"
#include "oratio/cblas.h"

/**
 * Create a new uniform generator on the unit hypersphere with some
 * initial state.
 *
 * @param state is an array of long integers defining the initial
 *   state. If defined to NULL, a new array will be allocated and will
 *   contain an arbitrary initial state.
 * @param m is the size of the random vector.
 *
 * @return the newly allocated random generator object.
 */
Random *ran_m_unit_sphere_new(unsigned long *seed, int n,
				     int m)
{
  Random *random = ran_random_new_with_size(seed, n, 0);

  return random;
}

/**
 * Draw a vector uniformly distributed on the unit hypersphere.
 *
 * @param ran is an initialized random numbers generator.
 * @param m is the dimension of the sphere.
 * @param x is on output a realization of the random vector defined on
 *   the unit sphere.
 *
 * @return x, a realization of the random vector defined uniformly on
 *   the unit sphere and the state of ran is updated.
 */
void ran_m_unit_sphere(Random *ran, int m, double *x)
{
  /**
   * Ref.: Luc Devroye, p.230.
   */

  double norm;
  int i;

  for (i = 0; i < m; i++) x[i] = ran_gaussian_icdf(ran);
  nt_vect_print(stdout, "X", x, m);
  norm = cblas_dnrm2(m, x, 1);
  cblas_dscal(m, 1.0/norm, x, 1);
} 

/**
 * Draw a vector uniformly distributed on the unit hypersphere.
 *
 * @param ran is an initialized random numbers generator.
 * @param m is the dimension of the sphere.
 * @param x is on output a realization of the random vector defined on
 *   the unit sphere.
 *
 * @return x, a realization of the random vector defined uniformly on
 *   the unit sphere and the state of ran is updated.
 */
void ran_m_unit_sphere_get_val(Random *ran, int m, double *x)
{
  return ran_m_unit_sphere(ran, m, x);
}
