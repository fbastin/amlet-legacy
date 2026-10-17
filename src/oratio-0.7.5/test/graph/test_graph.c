/**
 * test_graph.c
 *
 * ORATIO Library
 * Copyright (C) 2009-2010 <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License.
 */

#include <oratio/graph.h>

#include <math.h>

#define N 9

void likelihood_edge(double theta, GR_Node *node)
{
  GR_Edge *e;
  GR_Node *dest;
  DS_List *edges = node->out;

  while (edges) {
    e = (GR_Edge *)(edges->data);
    dest = e->destination;
    if (node->value[0] < dest->value[0] && node->value[1] > dest->value[1])
      e->value[0] = exp(theta*(dest->value[0] - node->value[0] - e->weight));
    else
      e->value[0] = 0.0;
    edges = edges->next;
  }
}

void link_weight(int origin, GR_Node *node)
{
  DS_List *edge_out = node->out;
  DS_List *edge_in;
  GR_Edge *e;

  while (edge_out) {
    e = (GR_Edge *)(edge_out->data);
    if (node->id != origin) {
      e->value[1] = 0.0;
      edge_in = node->in;
      while (edge_in) {
	e->value[1] += ((GR_Edge *)(edge_in->data))->value[1];
	edge_in = edge_in->next;
      }
      e->value[1] *= e->value[0];
    }
    else
      e->value[1] = e->value[0];      
    edge_out = edge_out->next;
  }
}

void link_flow(int destination, double flow, GR_Node *node)
{
  DS_List *edge_out = node->out;
  DS_List *edge_in = node->in;
  GR_Edge *e;

  double sum_in, sum_out;

  sum_in = 0.0;
  while (edge_in) {
    sum_in += ((GR_Edge *)(edge_in->data))->value[1];
    edge_in = edge_in->next;
  }

  if (node->id == destination)
    sum_out = flow;
  else {
    sum_out = 0.0;
    while (edge_out) {
      sum_out += ((GR_Edge *)(edge_out->data))->value[2];
      edge_out = edge_out->next;
    }
  }

  edge_in = node->in;
  while (edge_in) {
    e = (GR_Edge *)(edge_in->data);
    e->value[2] = (e->value[1]*sum_out);
    if (e->value[2] != 0.0) e->value[2] /= sum_in;
    edge_in = edge_in->next;
  }
}

int main()
{
  double dist[N];
  int previous[N];
  GR_Edge *previous_edges[N];

  int perm[N];

  int i;

  GR_Graph *graph = gr_graph_new(N);
  GR_Node *nodes = graph->nodes;
  GR_Node *temp_node;
  GR_Edge *e;

  DS_List *edges = ds_list_new(), *temp_edge, *start;

  for (i = 0; i < N; i++) graph->nodes[i].value = _mallocs(2, double);

  /* Connect the nodes */
  edges = ds_list_prepend(edges, gr_graph_connect(graph, 0, 1, 2.0));
  edges = ds_list_prepend(edges, gr_graph_connect(graph, 0, 3, 2.0));
  edges = ds_list_prepend(edges, gr_graph_connect(graph, 1, 2, 2.0));
  edges = ds_list_prepend(edges, gr_graph_connect(graph, 1, 4, 2.0));
  edges = ds_list_prepend(edges, gr_graph_connect(graph, 2, 5, 2.0));
  edges = ds_list_prepend(edges, gr_graph_connect(graph, 3, 4, 1.0));
  edges = ds_list_prepend(edges, gr_graph_connect(graph, 3, 6, 2.0));
  edges = ds_list_prepend(edges, gr_graph_connect(graph, 4, 5, 1.0));
  edges = ds_list_prepend(edges, gr_graph_connect(graph, 4, 7, 2.0));
  edges = ds_list_prepend(edges, gr_graph_connect(graph, 5, 8, 2.0));
  edges = ds_list_prepend(edges, gr_graph_connect(graph, 6, 7, 2.0));
  edges = ds_list_prepend(edges, gr_graph_connect(graph, 7, 8, 2.0));

  start = edges;
  while (edges) {
    ((GR_Edge *)(edges->data))->value = _mallocs(3, double);
    edges = edges->next;
  }

  /*  gr_graph_dijkstra(N, graph->nodes, 0, 8, dist, previous); */
  gr_graph_dijkstra_edges(graph, 0, 8, dist, previous_edges,
			  gr_edge_get_weight);

  /* Label the nodes */
  for (i = 0; i < N; i++) nodes[i].value[0] = dist[i];

  /* Reverse the connections. */
  for (i = 0; i < N; i++) {
    temp_edge = nodes[i].in;
    nodes[i].in = nodes[i].out;
    nodes[i].out = temp_edge;
  }

  edges = start;
  while (edges) {
    e = (GR_Edge *)(edges->data);
    temp_node = e->origin;
    e->origin = e->destination;
    e->destination = temp_node;
    edges = edges->next;
  }

  gr_graph_dijkstra(graph, 8, 0, dist, previous);

  /* Print the shortest path */
  printf("Shortest path: ");
  i = 0;
  do {
    printf("%d ", i);
    i = previous[i];
  } while (i >= 0);
  printf("\n");

  for (i = 0; i < N; i++) nodes[i].value[1] = dist[i];

  /* Reverse the connections. We return to the original graph */
  edges = start;

  for (i = 0; i < N; i++) {
    temp_edge = nodes[i].in;
    nodes[i].in = nodes[i].out;
    nodes[i].out = temp_edge;
  }
  while (edges) {
    e = (GR_Edge *)(edges->data);
    temp_node = e->origin;
    e->origin = e->destination;
    e->destination = temp_node;
    edges = edges->next;
  }

  /* Compute the initial link likelihoods */
  for (i = 0; i < N; i++) likelihood_edge(1.0, nodes+i);

  /* Step 1*/

  /* Forward pass */
  for (i = 0; i < N; i++) dist[i] = nodes[i].value[0];
  ds_sort_dheapsort_with_perm(N, dist, perm);
  for (i = 0; i < N; i++) link_weight(0, nodes+perm[i]);

  /* Backward pass */
  for (i = 0; i < N; i++) dist[i] = nodes[i].value[1];
  ds_sort_dheapsort_with_perm(N, dist, perm);
  for (i = 0; i < N; i++) link_flow(8, 1000, nodes+perm[i]);

  /* Print the flows */
  edges = start;
  while (edges) {
    e = (GR_Edge *)(edges->data);
    printf("Flow over (%d,%d): %g\n",
	   e->origin->id+1, e->destination->id+1, e->value[2]);
    edges = edges->next;
  }

  gr_graph_free(graph);

  return 0;
}
