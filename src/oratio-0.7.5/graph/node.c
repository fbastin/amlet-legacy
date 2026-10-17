/**
 * node.c
 *
 * ORATIO Library
 * Copyright (C) 2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License.
 */

#include "oratio/graph.h"

/**
 * Create a new node.
 *
 * @return a pointer to the newly created node.
 */
GR_Node *gr_node_new()
{
  GR_Node *node = _malloc(GR_Node);

  node->in = ds_list_new();
  node->out = ds_list_new();

  node->label = NULL;

  node->id = 1;

  node->value = node->data = NULL;

  return node;
}

/**
 * Free the memory associated to the content of a node.
 *
 * @param node is an initialized node.
 */
void gr_free_node_content(GR_Node *node)
{
  if (node->label) free(node->label);
  if (node->value) free(node->value);

  ds_list_free(node->in, NULL);
  ds_list_free(node->out, (Dfct)gr_free_edge_content);
}

/**
 * Free the memory associated to a node.
 *
 * @param node is an initialized node.
 */
void gr_node_free(GR_Node *node)
{
  if (node) {
    gr_free_node_content(node);

    free(node);
  }
}

/**
 * Return the node identifier.
 *
 * @param node is an initialized node.
 *
 * @return the node identifier.
 */
int gr_node_get_id(GR_Node *node)
{
  return node->id;
}

/**
 * Return the node internal identifier.
 *
 * @param node is an initialized node.
 *
 * @return the node internal identifier.
 */
int gr_node_get_internal_id(GR_Node *node)
{
  return node->internal_id;
}
