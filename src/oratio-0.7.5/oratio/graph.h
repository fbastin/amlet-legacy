/**
 * graph.h
 *
 * ORATIO Library
 * Copyright (C) 2009-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License
 */

#ifndef __ORATIO_GRAPH_H__
#define __ORATIO_GRAPH_H__

#include <oratio/dstruct.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Data structures.
   */

  typedef struct _GR_Node GR_Node;

  typedef struct _GR_Edge GR_Edge;

  /**
   * Node description.
   */
  struct _GR_Node {
    int id;          /* node identifier */
    int internal_id; /* internal identifier */

    DS_List *in;     /* input edges */
    DS_List *out;    /* output edges */

    char *label;

    double *value;

    void *data;
  };

  /**
   * Edge description.
   */
  struct _GR_Edge {
    GR_Node *origin;
    GR_Node *destination;

    char *label;

    double weight;

    double *value;

    void *data;
  };

  /**
   * Graph structure.
   */
  typedef struct _GR_Graph {
    int n;           /* number of nodes */
    int nmax;        /* maximum number of nodes */
    int offset;      /* offset: the node numbering starts at offset. */

    GR_Node *nodes;  /* set of nodes */
  } GR_Graph;

  /**
   * Create a new node.
   *
   * @return a pointer to the newly created node.
   */
  GR_Node *gr_node_new();

  /**
   * Free the memory associated to a node.
   *
   * @param node is a memory-allocated node.
   */
  void gr_node_free(GR_Node *node);

  /**
   * Free the memory associated to the content of a node.
   *
   * @param node is an initialized node.
   */
  void gr_free_node_content(GR_Node *node);

  /**
   * Return the node identifier.
   *
   * @param node is an initialized node.
   *
   * @return the node identifier.
   */
  int gr_node_get_id(GR_Node *node);

  /**
   * Return the node internal identifier.
   *
   * @param node is an initialized node.
   *
   * @return the node internal identifier.
   */
  int gr_node_get_internal_id(GR_Node *node);

  /**
   * Return the internal identifier of a node with some specific id. If two
   × nodes have the same identifier, the function returns the identifier of
   × the first node found.
   *
   * @param graph is the graph where we are looking for the node.
   × @param id is the identifier of the searched node.
   *
   × @return the internal identifier of the node with identifier id.
  */
  int gr_graph_get_node_internal_id(GR_Graph *graph, int id);

  /**
   * Create a new edge.
   *
   * @param origin is the origin node.
   * @param destination is the destination node.
   *
   * @return the newly memory-allocated edge.
   */
  GR_Edge *gr_edge_new(GR_Node *origin, GR_Node *destination);

  /**
   * Suppress an edge.
   *
   * @param edge is a memory-allocated edge.
   */
  void gr_edge_free(GR_Edge *edge);

  /**
   * Free the memory associated to edge values.
   *
   * @param edge is the edge under consideration.
   */
  void gr_free_edge_content(GR_Edge *edge);

  /**
   * Return the weight of an edge.
   *
   * @param edge is a memory-allocated edge.
   *
   * @return the weight of the edge.
   */
  double gr_edge_get_weight(GR_Edge *edge);

  /**
   * Allocate memory for a new graph and initialize the nodes.
   *
   * @param n is the number of nodes.
   * @param offset it the number of the first number, assuming that nodes
   *   numbers are consecutive.
   *
   * @return a pointer to the newly created graph.
   */
  GR_Graph *gr_graph_new_offset(int n, int offset);

  /**
   * Allocate memory for a new graph and initialize the nodes.
   *
   * @param n is the number of nodes.
   *
   * @return a pointer to the newly created graph.
   */
  GR_Graph *gr_graph_new(int n);

  /**
   * Free the memory associated to a graph.
   *
   * @param graph is an initialized graph.
   */
  void gr_graph_free(GR_Graph *graph);

  /**
   * Create an edge between two nodes in a graph.
   *
   * @param origin is the identifier of the origin node.
   * @param destination is the identifier of the destination node.
   * @param weight is the weight of the link.
   *
   * @return the newly created edge.
   */
  GR_Edge *gr_graph_connect(GR_Graph *graph,
			    int origin, int destination,
			    double weight);

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
				double weight);

  /**
   * Correct the nodes identifiers in a graph with respect to the offset.
   */
  void gr_graph_correct_offset(GR_Graph *graph);

  /**
   * Compute the offset inside a graph.
   */
  int gr_graph_compute_offset(GR_Graph *graph);

  void gr_graph_finalize(GR_Graph *graph);

  void gr_graph_add_node(GR_Graph *graph);

  void gr_graph_add_node_id(GR_Graph *graph, int id);

  /**
   * Shortest paths algorithms.
   */

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
				   int *parent, double *w);

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
				int *parent, double *w);

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
			 double *dist, int *previous);

  /**
   * Dijkstra algorithm (edges version).
   *
   * @param graph is the graph under consideration.
   * @param origin is the internal identifier of the origin node.
   * @param destination is the internal identifier of the destination node.
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
			       double (*weight)(GR_Edge *edge));

  /**
   * Compute a shortest path in a topologically ordered graph, using
   * backward linking algorithm.
   *
   * @param graph is the graph under consideration.
   * @param origin is the internal identifier of the origin node.
   * @param destination is the internal identifier of the destination node.
   * @param J is a double array of size at least equal to n; on output,
   *   J[i] contains the distance from node i to the destination.
   * @param next is a int array of size at least equal to n; on output,
   *   next[i] contains the successor of node i on the shortest path
   *   from i to the destination.
   */
  void gr_graph_backward_linking(GR_Graph *graph,
				 int origin, int destination,
				 double *J, int *next);

  /**
   * Compute a shortest path in a topologically ordered graph, using
   * forward linking algorithm.
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
  void gr_graph_forward_linking(GR_Graph *graph,
				int origin, int destination,
				double *D, int *previous);

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
			  double *D, int *previous);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_GRAPH_H__ */

