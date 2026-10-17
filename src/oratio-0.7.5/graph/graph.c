/**
 * graph.c
 *
 * ORATIO Library
 * Copyright (C) 2009-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License.
 */

#include "oratio/graph.h"
#include "oratio/ctools.h"
#include "oratio/ntools.h"

#include <limits.h>

/**
 * Return the internal identifier of a node with some specific id. If two
 × nodes have the same identifier, the function returns the identifier of
 × the first node found.
 *
 * @param graph is the graph where we are looking for the node.
 × @param id is the identifier of the searched node.
 *
 × @return the internal identifier of the node with identifier id, or -1
 *   if the node has not been found.
 */
int gr_graph_get_node_internal_id(GR_Graph *graph, int id)
{
  int i;

  for (i = 0; i < graph->n; i++) {
    if (graph->nodes[i].id == id) return i;
  }

  return -1;
}

/**
 * Correct the nodes identifiers in a graph with respect to the offset.
 */
void gr_graph_correct_offset(GR_Graph *graph)
{
  int i;

  for (i = 0; i < graph->n; i++) graph->nodes[i].id -= graph->offset;
}

/**
 * Compute the offset inside a graph.
 */
int gr_graph_compute_offset(GR_Graph *graph)
{
  int i;
  int offset = INT_MAX;

  for (i = 0; i < graph->n; i++)
    offset = nt_min(offset, graph->nodes[i].id);
  graph->offset = offset;

  return offset;
}

void gr_graph_finalize(GR_Graph *graph)
{
  graph->nodes = _reallocs(graph->nodes, graph->n, GR_Node);
  graph->nmax = graph->n;
  gr_graph_compute_offset(graph);
  gr_graph_correct_offset(graph);
}

void gr_graph_add_node_id(GR_Graph *graph, int id)
{
  int i;

  i = graph->n++;

  if (graph->n > graph->nmax) {
    graph->nmax = 2*graph->n;
    graph->nodes = _reallocs(graph->nodes, graph->nmax, GR_Node);
  }
  graph->nodes[i].in = ds_list_new();
  graph->nodes[i].out = ds_list_new();

  graph->nodes[i].label = NULL;
  graph->nodes[i].value = NULL;

  graph->nodes[i].internal_id = i;

  graph->nodes[i].id = id;
}

void gr_graph_add_node(GR_Graph *graph)
{
  gr_graph_add_node_id(graph, graph->n);
}

/**
 * Allocate memory for a new graph and initialize the nodes.
 *
 * @param n is the number of nodes.
 * @param offset it the number of the first number, assuming that nodes
 *   numbers are consecutive.
 *
 * @return a pointer to the newly created graph.
 */
GR_Graph *gr_graph_new_offset(int n, int offset)
{
  int i;

  GR_Graph *graph = _malloc(GR_Graph);

  graph->n = graph->nmax = n;
  graph->offset = offset;

  if (n) {
    graph->nodes = _mallocs(n, GR_Node);
    for (i = 0; i < n; i++) {
      graph->nodes[i].in = ds_list_new();
      graph->nodes[i].out = ds_list_new();

      graph->nodes[i].label = NULL;
      graph->nodes[i].value = NULL;
      
      graph->nodes[i].id = graph->nodes[i].internal_id = i;
    }
  }
  else
    graph->nodes = NULL;

  return graph;
}

/**
 * Allocate memory for a new graph and initialize the nodes.
 *
 * @param n is the number of nodes.
 *
 * @return a pointer to the newly created graph.
 */
GR_Graph *gr_graph_new(int n)
{
  return gr_graph_new_offset(n, 0);
}

/**
 * Free the memory associated to a graph.
 *
 * @param graph is an initialized graph.
 */
void gr_graph_free(GR_Graph *graph)
{
  int i;

  if (graph) {
    for (i = 0; i < graph->n; i++) gr_free_node_content(graph->nodes+i);

    free(graph);
  }
}

/**
 * Create an edge between two nodes in a graph.
 *
 * @param origin is the identifier of the origin node.
 * @param destination is the identifier of the destination node.
 * @param weight is the weight of the link.
 *
 * @return the newly created edge.
 */
GR_Edge *gr_graph_connect(GR_Graph *graph, int origin, int destination,
			  double weight)
{
  int a, b;
  GR_Edge *edge;

  a = gr_graph_get_node_internal_id(graph, origin);
  b = gr_graph_get_node_internal_id(graph, destination);

  edge = gr_edge_new(graph->nodes+a, graph->nodes+b);
  edge->weight = weight;

  return edge;
}

/**
 * Create an edge in a graph.
 *
 * @param origin is the internal identifier of the origin node.
 * @param destination is the internal identifier of the destination node.
 * @param weight is the weight of the link.
 *
 * @return the newly created edge.
 */
GR_Edge *gr_graph_create_edge(GR_Graph *graph,
			      int origin, int destination,
			      double weight)
{
  GR_Edge *edge = gr_edge_new(graph->nodes+origin, graph->nodes+destination);
  edge->weight = weight;

  return edge;
}
