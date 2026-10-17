/**
 * ntools.h
 *
 * ORATIO Library
 * Copyright (C) 2002-2015 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_NTOOLS_H__
#define __ORATIO_NTOOLS_H__

#include "bsplines.h"
#include "cblas.h"
#include "constants.h"
#include "nt_coordinates.h"
#include "nt_derive.h"
#include "error.h"
#include "functions.h"
#include "nt_io.h"
#include "matrix.h"
#include "natural.h"
#include "norm.h"
#include "nt_gamma.h"
#include "nt_polynomial.h"
#include "nt_simplex.h"
#include "nt_utilities.h"
#include "vector.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Graph algorithms.
   */

  /**
   * Compute the shortest path using the label-correcting
   * algorithm. The graph may present negative arcs, but no negative
   * cycles.
   *
   * @param n is the number of nodes in the graph.
   * @param s is the origin.
   * @param t is the end.
   * @param adj is the collection of adjency lists, i.e. adj[i] is the
   *   adjency list associated list associated to node i. The end of the
   *   list must be represented by -1.
   * @param v_adj is the values of the arcs corresponding to adj.
   * @param u is an array such that u[i] is an underestimate of the
   *   shortest past from i to t.
   * @param parent is an array that will contain on output the parents
   *   of the nodes on the shortest path.
   * @param w are the corresponding weigths on the shortest path.
   *
   * @return 1 if the shortest path has been found, 0 else.
   */
  int nt_graph_labelcorrecting_adj(int n, int s, int t,
				   int **adj, double **v_adj, double *u,
				   int *parent, double *w);

  /**
   * Log file gestion.
   */

  typedef struct {
    FILE *f;
    char *name;
  } NTLog;

  /**
   * Create a new NTLog file.
   *
   * @param name must be NULL or assigned to a valid file name.
   *
   * @return a newly allocated NTLog structure.
   */
  NTLog *nt_log_new(char *name);

  /**
   * Set the output file the log objet 'log'.
   *
   * @param name must be NULL or assigned to a valid file name.
   */
  void nt_log_set_file(NTLog *log, char *name);

  void nt_log_free(NTLog *ntlog);

  /**
   * Output functions
   */

  void nt_log_write(NTLog *ntlog, const char *msg, ...);

  void nt_log_section(NTLog *ntlog, const char *msg, ...);

  void nt_log_subsection(NTLog *ntlog, const char *msg, ...);

  /**
   * CONVERSION ROUTINES
   */

  /**
   * Convertit la matrice A en vecteur x, en écrivant A ligne par ligne si
   * orientation vaut NT_ROW, ou colonne par colonne si orientation vaut
   * NT_COLUMN
   */
  void nt_matrix_2_vector(double **A, double x[],
			  int n, int m, int orientation);

  void nt_vect_2_matrix(double x[], double **A,
			int n, int m, int orientation);


  /**
   * Modular arithmetic.
   */

  /**
   * Compute ax+y mod m.
   *
   * @param a is an integer satisfying 0 <= a*m < 2^53.
   * @param x is an integer.
   * @param y is an integer.
   * @param m is a strictly positive integer satisfying m < 2^35.
   *
   * @return ax+x mod m.
   */
  double nt_modular_daxpy(double a, double x, double y, double m);

  /**
   * Compute ax mod m using Schrage's method.
   *
   * @param a is an integer such that a*a<m or a = floor(m/i), where i*i<m.
   * @param x is a positive integer between 0 and m-1 (included).
   * @param m is a positive integer.
   *
   * @return ax mod m.
   */
  long nt_modular_schrage(long a, long x, long m);

  /**
   * Compute (a'x + y) mod m.
   *
   * @param n is the dimension of the vectors a and x.
   * @param a is a vector of dimension n.
   * @param x is a vector of dimension n.
   * @param y is a double with some initialized value on input.
   * @param m is a positive integer such that -m < x[i] < m.
   *
   * @return the value (a'x + y) mod m.
   */
  double nt_modular_vecvec(int n, double *a, double *x,
			   double y, double m);

  /**
   * Returns v = A*x mod m.
   *
   * @param A is a matrix of leading dimension incA, stored by rows. In
   *   other terms incA is the number of columns of A.
   * @param x is a vector of dimension n. Its content can be replaced on
   *   output if v is a pointer to the same vector.
   * @param v is a vector whose value on output will be equal to Ax mod m.
   *   On input is content is undefined. v can point to the same pointer
   *   than x.
   * @param m is an integer such that -m < x[i] < m.
   *
   * @return v = A*x mod m.
   */
  void nt_modular_matvec(int n, double *A, int incA, double *x,
			 double *v, double m);

  /**
   * Compute C = A*B mod m.
   *
   * @param is the effective dimension of A, B, C, i.e. they are viewed
   *   as nxn matrices.
   * @param A is a matrix of dimension at least nxn.
   * @param incA is the leading dimension of A.
   * @param B is a matrix of dimension at least nxn.
   * @param incB is the leading dimension of A.
   * @param C is a matrix of dimension at least nxn.
   *   On output, C = A*B mod m. Work even if A = C or B = C or A = B = C.
   * @param incC is the leading dimension of A.
   * @param m is a strictly positive integer.
   */
  void nt_modular_matmat(int n, double *A, int incA, double *B, int incB,
			 double *C, int incC, double m);

  /**
   * Compute matrix B = A^p mod m.
   *
   * @param n is the dimensions of matrices A and B, which are assumed
   *   to be square.
   * @param A is a matrix of dimension nxn.
   * @param incA is the leading dimension of A.
   * @param B is a matrix of dimension nxn. B may point to A.
   * @param incB is the leading dimension of B.
   * @param m is a strictly positive integer.
   * @param p is a positive integer.
   *
   * @returns B = A^p mod m.
   */
  void nt_modular_matpow(int n, double *A, int incA,
			 double *B, int incB, double m, long p);

  /**
   * Compute matrix B = (A^(2^e) % m).
   *
   * @param n is the dimensions of matrices A and B, which are assumed
   *   to be square.
   * @param A is a matrix of dimension nxn.
   * @param incA is the leading dimension of A.
   * @param B is a matrix of dimension nxn. B may point to A.
   * @param incB is the leading dimension of B.
   * @param m is a strictly positive integer.
   * @param p is a positive integer.
   *
   * @returns B = A^(2^e) mod m.
   */
  void nt_modular_mattwopow(int n, double *A, int incA,
			    double *B, int incB, double m, long p);

  /**
   * Primes numbers.
   */

  /**
   * Compute the n first prime numbers using successive divisions.
   *
   * @param n is a non-null natural.
   * @param primes is an allocated array of at least n elements. On
   *   output, primes contains the n first prime numbers, by ascending
   *   order, in the first n positions.
   */
  void nt_primes_simple(unsigned int n, unsigned int *primes);

  /**
   * Compute the n first prime numbers using the Eratoshtene sieve.
   *
   * @param n is a non-null natural.
   * @param primes is an allocated array of at least n elements. On
   *   output, primes contains the n first prime numbers, by ascending
   *   order, in the first n positions.
   */
  void nt_primes_eratosthene(unsigned int n, unsigned int *primes);

  /**
   * Special values.
   */

  /**
   * Return the infinity representation.
   *
   * @return infinity;
   */
  double nt_infinity();

  /**
   * Return the Not a Number representation.
   *
   * This representation is architecture dependent.
   * Be careful: this will produce a floating exception on some
   * architectures, in particular when no direct double representation
   * is available.
   */
  double nt_nan();

  /**
   * Misc utilities.
   */

  /**
   * Compute x^n, where x is a real and n is a positive integer.
   *
   * @param x is a real number.
   * @param n is a positive (unsigned) integer.
   *
   * @return x^n.
   */
  double nt_ipow(double x, unsigned int n);

  /**
   * Compute ln(n!). This can be useful if n is high, since we obtain
   * then a better precision.
   *
   * @param n is a positive integer.
   *
   * @return len(n!) as a double.
   */
  double nt_factorial_ln(int n);

  /**
   * Compute consecutive powers of 2, starting from 2^0 = 1.
   *
   * @param n is the number of powers to compute (including 0).
   * @param p is an allocated integer array. On output, this array
   *   contains the successive powers of 2.
   */
  void nt_two_powers_array(int n, unsigned int *p);

  /**
   * Compute the closest power of two to some integer x, but less or equal to x.
   *
   * @param x is the considered integer.
   *
   * @return i, where 2^i <= x, and for all j > i, 2^j > x.
   */
  int nt_two_power_floor(unsigned int x);

  /**
   * Compute the closest power of two to some integer x, but greater or
   * equal to x.
   *
   * @param x is the considered integer.
   *
   * @return i, where 2^i >= x, and for all j < i, 2^j < x.
   */
  int nt_two_power_ceil(int x);

  /**
   * Compute the radical inverse of integer i in base b.
   *
   * @param b is the base.
   * @param i is the integer to develop.
   *
   * @return the radical inverse of i in base b.
   */
  double nt_radical_inverse(int b, int i);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_NTOOLS_H__ */
