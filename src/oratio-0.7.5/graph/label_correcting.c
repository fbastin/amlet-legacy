/**
 * label_correcting.c
 *
 * ORATIO Library
 * Copyright (C) 2009-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License
 */

#include "oratio/ntools.h"
#include "oratio/stat.h"

/**
 * Computes the shortest path using the label-correcting
 * algorithm. The graph may present negative arcs, but no negative
 * cycle. The offset of the graph has to be equal to 0.
 *
 * @param n is the number of nodes in the graph.
 * @param origin is the internal identifier of the origin node.
 * @param destination is the internal identifier of the destination node.
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
int gr_graph_labelcorrecting_adj(int n, int origin, int destination,
				 int **adj, double **v_adj, double *u,
				 int *parent, double *w)
{
  int open[n];
  int last = 0;

  int present[n];
  double d[n];

  double upper;

  int i, j, c;

  /* Initialization. */
  open[0] = origin;

  upper = nt_infinity();

  for (i = 0; i < n; i++) {
    d[i] = nt_infinity();
    present[i] = 0;
    parent[i] = i;
    w[i] = 0;
  }
  d[origin] = 0.0;
  present[origin] = 1;

  while (last >= 0) {
    /* Remove a node i from open. */
    i = open[last--];
    present[i] = 0;

    /* Go along the adjacency list of i. */
    j = 0;
    c = adj[i][0];
    while (c != -1) {
      /*      printf("Node i: %d, Next: %d\n", i+1, c+1);
	      printf("d_i: %g, weight: %g, d_next: %g, upper: %g, u_next: %g\n", d[i], v_adj[i][j], d[c], upper, u[c]); */
      if (d[i] + v_adj[i][j] < nt_min(d[c], upper-u[c])) {
	d[c] = d[i] + v_adj[i][j];
	parent[c] = i;
	w[c] = v_adj[i][j];
	if (c == destination)
	  upper = d[c];
	else {
	  if (!present[c]) {
	    present[c] = 1;
	    open[++last] = c;
	  }
	}
      }
      c = adj[i][++j];
    }
  }

  return (upper == nt_infinity()) ? 0 : 1;
}

/**
 * Computes the shortest path using the label-correcting
 * algorithm. The graph may present negative arcs, but no negative
 * cycle. The offset of the graph has to be equal to 0.
 *
 * @param n is the number of nodes in the graph.
 * @param origin is the internal identifier of the origin node.
 * @param destination is the internal identifier of the destination node.
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
int gr_graph_label_adj_slflll(int n, int origin, int destination,
			      int **adj, double **v_adj, double *u,
			      int *parent, double *w)
{
  int open[n];
  int open_start, open_last;
  int size;

  int present[n];
  double d[n];

  double upper, average;

  int i, j, c;

  /* Initialization. */
  open[0] = origin;
  open_start = open_last = 0;
  size = 1;
  average = 0.0;

  upper = nt_infinity();

  for (i = 0; i < n; i++) {
    d[i] = nt_infinity();
    present[i] = 0;
    parent[i] = i;
    w[i] = 0;
  }
  d[origin] = 0.0;
  present[origin] = 1;

  while (size) {
    /* Remove a node i from open. We take the first one in the list. */
    i = open[open_start++];
    if (open_start == n) open_start = 0;
    present[i] = 0;

    /* Go along the adjacency list of i. */
    j = 0;
    c = adj[i][0];
    while (c != -1) {
      if (d[i] + v_adj[i][j] < nt_min(d[c], upper-u[c])) {
	d[c] = d[i] + v_adj[i][j];
	parent[c] = i;
	w[c] = v_adj[i][j];
	if (c == destination)
	  upper = d[c];
	else {
	  if (!present[c]) {
	    present[c] = 1;
	    average = st_mean_update(size++, average, d[c]);
	    /* Apply the SLF technique */
	    if (d[c] < d[open[open_start]]) {
	      open_start--;
	      if (open_start == -1) open_start = n-1;
	      open[open_start] = c;
	    }
	    else {
	      open_last++;
	      if (open_last == n) open_last = 0;
	      open[open_last] = c;
	    }
	  }
	}
      }
      c = adj[i][++j];
    }

    average = st_mean_remove(size--, average, d[i]);
    /* LLL */
    if (size) {
      while (d[open[open_start]] > average) {
	open_last++;
	if (open_last == n) open_last = 0;
	open[open_last] = open[open_start];
	open_start++;
	if (open_start == n) open_start = 0;
      }
    }
  }

  return (upper == nt_infinity()) ? 0 : 1;
}
