/**
 * dijkstra.c
 *
 * ORATIO Library
 * Copyright (C) 2009-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License.
 */

#include "oratio/graph.h"
#include "oratio/ntools.h"

/**
 * Dijkstra algorithm.
 *
 * @param graph is the graph under consideration.
 * @param origin is the internal identifier of the origin node.
 * @param destination is the internal identifier of the destination node.
 * @param dist is a double array of size at least equal to n; on output,
 *   dist[i] contains the distance from the origin to node i+offset.
 * @param previous is a int array of size at least equal to n; on output,
 *   previous[i] contains the predecessor of i on the shortest path is i
 *   is present on this path, -1 else, using internal identifiers.
 */
void gr_graph_dijkstra(GR_Graph *graph,
		       int origin, int destination,
		       double *dist, int *previous)
{
  int i, j;

  int u, v;  /* node indexes */

  double min_dist, alt;

  DS_List *neighbor;

  int n = graph->n; /* number of nodes in the graph. */
  GR_Node *nodes = graph->nodes; /* array of nodes in the graph. */

  int present[n];

  /* Initialization */
  for (i = 0; i < n; i++) {
    dist[i] = nt_infinity();  /* unknown distance function from source to i */
    previous[i] = -1;       /* previous node in optimal path from source */
    present[i] = 0;
  }
  dist[origin] = 0.0;  /* distance from source to source */

  for (i = 0; i < n; i++) {
    u = -1;
    min_dist = nt_infinity();

    /* We search the next node to consider. This is not efficient for
       large networks - there is place for improvment here. */
    for (j = 0; j < n; j++) {
      if (dist[j] < min_dist && !present[j]) {
	u = j;
	min_dist = dist[j];
      }
    }
    if (u == -1)
      return;  /* all remaining nodes are inaccessible from origin */

    present[u] = 1;
    /* We consider the edges going out the current node. */
    neighbor = nodes[u].out;

    while (neighbor) {
      alt = dist[u] + ((GR_Edge *)(neighbor->data))->weight;
      v = ((GR_Edge *)(neighbor->data))->destination->internal_id;
      if (alt < dist[v]) {
	dist[v] = alt;
	previous[v] = u;
      }
      neighbor = neighbor->next;
    }
  }
}

/**
 * Dijkstra algorithm (edges version).
 *
 * @param graph is the graph under consideration.
 * @param origin is the internal identifier of origin node.
 * @param destination is the internal identifier of destination node.
 * @param dist is a double array of size at least equal to n; on output,
 *   dist[i] contains the distance from the origin to node i, where i is
 *   the internal identifier.
 * @param previous is an array of pointer to edges, of size at least
 *   equal to n; on output previous [i] is a pointer to the edge leading to
 *   i in the shortest path, if i is present on this path, else it is set
 *   to NULL.
 * @param weight is a function giving the weight to consider for an edge.
 */
void gr_graph_dijkstra_edges(GR_Graph *graph,
			     int origin, int destination,
			     double *dist, GR_Edge **previous,
			     double (*weight)(GR_Edge *edge))
{
  int i, j;

  int u, v;  /* node indexes */

  double min_dist, alt;

  DS_List *neighbor;

  int n = graph->n; /* number of nodes in the graph. */
  GR_Node *nodes = graph->nodes; /* array of nodes in the graph. */

  int present[n];

  /* Initialization */
  for (i = 0; i < n; i++) {
    dist[i] = nt_infinity();  /* unknown distance function from source to i */
    previous[i] = NULL;
    present[i] = 0;
  }
  dist[origin] = 0.0;  /* distance from source to source */

  for (i = 0; i < n; i++) {
    u = -1;
    min_dist = nt_infinity();

    /* We search the next node to consider. This is not efficient for
       large networks - there is place for improvment here. */
    for (j = 0; j < n; j++) {
      if (dist[j] < min_dist && !present[j]) {
	u = j;
	min_dist = dist[j];
      }
    }
    if (u == -1)
      return;  /* all remaining nodes are inaccessible from origin */

    present[u] = 1;
    /* We consider the edges going out the current node. */
    neighbor = nodes[u].out;

    while (neighbor) {
      alt = dist[u] + weight((GR_Edge *)(neighbor->data));
      v = ((GR_Edge *)(neighbor->data))->destination->internal_id;
      if (alt < dist[v]) {
	dist[v] = alt;
	previous[v] = (GR_Edge *)(neighbor->data);
      }
      neighbor = neighbor->next;
    }
  }
}
