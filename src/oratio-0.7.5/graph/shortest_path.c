/**
 * shortest_path.c
 *
 * ORATIO Library
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License.
 */

#include "oratio/graph.h"
#include "oratio/ntools.h"

/**
 * Compute a shortest path in a topologically ordered graph, using
 * backward linking algorithm.
 *
 * @param graph is the graph under consideration.
 * @param origin is the internal identifier of the origin node.
 * @param destination is the internal identifier of the destination node.
 * @param J is a double array of size at least equal to n; on output,
 *   J[i] contains the distance from node i to the destination, using
 *   an offset equal to 0.
 * @param next is a int array of size at least equal to n; on output,
 *   next[i] contains the successor of node i on the shortest path
 *   from i to the destination, using an offset equal to 0.
 */
void gr_graph_backward_linking(GR_Graph *graph,
			       int origin, int destination,
			       double *J, int *next)
{
  int i, j;

  double Jcand;
  DS_List *neighbor;

  int n = graph->n; /* number of nodes in the graph. */
  GR_Node *nodes = graph->nodes; /* array of nodes in the graph. */

  origin -= graph->offset;
  destination -= graph->offset;

  J[n-1] = 0;
  next[n-1] = -1;

  for (i = n-2; i >= 0; i--) {
    J[i] = nt_infinity();
    next[i] = -1;

    neighbor = nodes[i].out;
    while (neighbor) {
      j = ((GR_Edge *)(neighbor->data))->destination->internal_id;
      Jcand = ((GR_Edge *)(neighbor->data))->weight + J[j];
      if (Jcand < J[i]) {
	J[i] = Jcand;
	next[i] = j;
      }
      neighbor = neighbor->next;
    }
  }
}

/**
 * Compute a shortest path in a topologically ordered graph, using
 * forward linking algorithm.
 *
 * @param graph is the graph under consideration.
 * @param origin is the internal identifier of the origin node.
 * @param destination is the internal identifier of the destination node.
 * @param D is a double array of size at least equal to n; on output,
 *   D[i] contains the distance from the origin to node i, using an offset
 *   equal to 0.
 * @param previous is a int array of size at least equal to n; on output,
 *   previous[i] contains the predecessor of node i on the shortest path
 *   from the origin to i, using an offset equal to 0.
 */
void gr_graph_forward_linking(GR_Graph *graph,
			      int origin, int destination,
			      double *D, int *previous)
{
  int i, j;

  double Dcand;
  DS_List *neighbor;

  int n = graph->n; /* number of nodes in the graph. */
  GR_Node *nodes = graph->nodes; /* array of nodes in the graph. */

  origin -= graph->offset;
  destination -= graph->offset;

  D[0] = 0;
  previous[0] = -1;

  for (j = 1; j < n; j++) {
    D[j] = nt_infinity();
    previous[j] = -1;

    neighbor = nodes[j].in;
    while (neighbor) {
      i = ((GR_Edge *)(neighbor->data))->origin->internal_id;
      Dcand = ((GR_Edge *)(neighbor->data))->weight + D[i];
      if (Dcand < D[j]) {
	D[j] = Dcand;
	previous[j] = i;
      }
      neighbor = neighbor->next;
    }
  }
}

/**
 * Compute a shortest path in a topologically ordered graph, using
 * the accession method.
 *
 * @param graph is the graph under consideration.
 * @param origin is the internal identifier of the origin node.
 * @param destination is the internal identifier of the destination node.
 * @param D is a double array of size at least equal to n; on output,
 *   D[i] contains the distance from the origin to node i.
 * @param previous is a int array of size at least equal to n; on output,
 *   previous[i] contains the predecessor of node i on the shortest path
 *   from the origin to i.
 */
void gr_graph_accession(GR_Graph *graph,
			int origin, int destination,
			double *D, int *previous)
{
  int i, j;

  double Dcand;
  DS_List *neighbor;

  int n = graph->n; /* number of nodes in the graph. */
  GR_Node *nodes = graph->nodes; /* array of nodes in the graph. */

  origin -= graph->offset;
  destination -= graph->offset;

  D[0] = 0;
  previous[0] = -1;

  for (j = 1; j < n; j++) {
    D[j] = nt_infinity();
    previous[j] = -1;
  }

  for (i = 0; i < n-1; i++) {
    if (D[i] < D[n-1]) {
      neighbor = nodes[i].out;
      while (neighbor) {
	j = ((GR_Edge *)(neighbor->data))->destination->internal_id;
	Dcand = ((GR_Edge *)(neighbor->data))->weight + D[i];
	if (Dcand < D[j]) {
	  D[j] = Dcand;
	  previous[j] = i;
	}
	neighbor = neighbor->next;
      }
    }
  }
}
