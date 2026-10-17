/**
 * edge.c
 *
 * ORATIO Library
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License.
 */

#include "oratio/graph.h"

/**
 * Free the memory associated to edge values.
 *
 * @param edge is the edge under consideration.
 */
void gr_free_edge_content(GR_Edge *edge)
{
  if (edge->value) free(edge->value);
}

/**
 * Create a new edge.
 *
 * @param origin is the origin node, with an offset equal to 0.
 * @param destination is the destination node, with an offset equal to 0.
 *
 * @return the newly memory-allocated edge, with a weight initialized to 0.
 */
GR_Edge *gr_edge_new(GR_Node *origin, GR_Node *destination)
{
  GR_Edge *edge = _malloc(GR_Edge);

  edge->origin = origin;
  edge->destination = destination;

  edge->label = NULL;

  edge->weight = 0.0;
  edge->value = edge->data = NULL;

  /* Connect the origin to the destination. */
  origin->out = ds_list_prepend(origin->out, edge);
  destination->in = ds_list_prepend(destination->in, edge);

  return edge;
}

/**
 * Suppress an edge.
 *
 * @param edge is a memory-allocated edge.
 */
void gr_edge_free(GR_Edge *edge)
{
  if (edge) {
    if (edge->value) free(edge->value);

    edge->origin->out = ds_list_remove(edge->origin->out, edge);
    edge->destination->in = ds_list_remove(edge->destination->in, edge);

  free(edge);
  }
}

/**
 * Return the weight of an edge.
 *
 * @param edge is a memory-allocated edge.
 *
 * @return the weight of the edge.
 */
double gr_edge_get_weight(GR_Edge *edge)
{
  return edge->weight;
}
