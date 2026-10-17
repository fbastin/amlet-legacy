/**
 * lattice.h
 *
 * OPHELIA Library
 * Copyright (C) 2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *                    David Munger  <mungerd@gmail.com>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __AMLET_LATTICE_H__
#define __AMLET_LATTICE_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Return the default lattice.
   *
   * @return a pointer to the default lattice, described a single vecteur, of
   *   (AMLET_DEFAULT_LATTICE_MAX_DIMENSION+1)x(AMLET_DEFAULT_LATTICE_NUMBER).
   *   The number of points of lattices are given in position
   *   k(AMLET_DEFAULT_LATTICE_MAX_DIMENSION+1),
   *   k=0,...,(AMLET_DEFAULT_LATTICE_NUMBER)-1, and followed by the correspond
   *   lattice parameters for AMLET_DEFAULT_LATTICE_MAX_DIMENSION dimensions.
   */
  int *amlet_get_default_lattice();

  /**
   * Return the maximum dimension for which one can use one predefined
   * rank-1 lattice.
   *
   * @return the maximum dimension of predefined rank-1 lattices.
   */
  int amlet_get_default_lattice_max_dimension();

  int amlet_get_default_lattice_number();

  /**
   * Read the lattices dimension.
   *
   * @param f is the file containing the lattices parameters.
   *
   * @return the dimension of the lattice to consider in the file, starting
   *   from the line below the current one.
   */
  int amlet_read_lattices_dimension(FILE *f);

  /**
   * Read lattices parameters.
   *
   * @param f is the file containing the lattices parameters.
   * @param s is the lattices dimension.
   * @param lattices is a pointer to an unallocated integer array on input,
   *   and to the lattices description on output.
   * @param nlattices is a pointer to an integer containing on output the
   *   number of parameters sets (each one defined for a certain number of
   *   points).
   *
   * @return the number of defined lattices.
   */
  int amlet_read_lattices(FILE *f, int s, int **lattices, int *nlattices);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_LATTICE_H__ */
