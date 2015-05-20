#ifndef GRAPH_H
#define GRAPH_H


#include <fstream>
#include <iostream>
#include <vector>

#include <bits/succinctvector.h>
#include <misc/definitions.h>


namespace CSA
{

//--------------------------------------------------------------------------

struct GraphNode
{
  uint label; // We assume that CHAR_BIT bits is enough for the label.
  uint value; // Used for sampling in GCSA.

  GraphNode() {}
  GraphNode(uint _label, uint _value) : label(_label), value(_value) {}
};

struct GraphEdge
{
  uint from, to;

  GraphEdge() {}
  GraphEdge(uint _from, uint _to) : from(_from), to(_to) {}
  GraphEdge(const GraphEdge& original, uint adjustment) :
    from(original.from + adjustment),
    to(original.to + adjustment)
  {
  }
};

//--------------------------------------------------------------------------

struct PathNode
{
  uint      from, to;
  pair_type key;

  void setSorted()      { this->to = this->from; }
  bool isSorted() const { return (this->to == this->from); }

  // These can be used only when PathGraph.status == ready.
  // Backbone assumes that outdegree does not need the highest order bit in an usint.
  // When backbone is set, outdegree cannot be used anymore.
  uint  value() const         { return this->to; }
  usint outdegree() const     { return this->key.first; }
  bool  isBackbone() const    { return (this->key.first & (((usint)1) << (WORD_BITS - 1))); }
  void  setBackbone()         { this->key.first |= ((usint)1) << (WORD_BITS - 1); }
};

struct PathEdge
{
  usint from, rank; // rank equals to.
  uint label;       // Of from.

  PathEdge(usint _from, usint _rank, uint _label) : from(_from), rank(_rank), label(_label) {}
};

//--------------------------------------------------------------------------

class Graph
{
  public:
    Graph(std::vector<GraphNode>& node_vector, std::vector<GraphEdge>& edge_vector);
    Graph(const std::string& base_name, bool is_full_name = false);
    ~Graph();

    void write(const std::string& base_name, bool is_full_name = false);

    void printInfo();

    void sortEdges(bool from);
    pair_type getEdges(uint node, bool from);  // Use sortEdges() first.
    pair_type getNextEdgeRange(pair_type range, bool from); // Use sortEdges() first.

    // Marks predecessor labels in GraphEdges.
    void markPredecessors();

    // Lower case labels indicate that the node is not a part of the backbone.
    void createBackbone();  // Create the backbone structure; convert labels to upper case.
    void removeBackbone();  // Convert labels to upper case.
    void encodeBackbone();  // Encode backbone information into labels.

    // Make the graph reverse deterministic. Updates backbone, if present.
    // The value of a merged node is the maximum value of the original below threshold,
    // or the maximum of all original nodes if the values are all >= threshold.
    // Threshold is the number of backbone nodes, if backbone is present, or the value of the 
    // current node otherwise.
    void determinize();
    bool isDeterministic();

    uint maxLabel();

    // Updates the labels when there are multiple automata to be indexed.
    // Final nodes get labels 0..automata-1.
    // The labels of regular nodes are incremented by automata.
    // Initial nodes get labels max_label+automata..max_label+2automata-1.
    void changeLabels(uint automata, uint max_label);
    void restoreLabels(uint automata, uint max_label);

    // Adds the nodes and edges from the given graph to current graph.
    // If both graphs contain a backbone, the value fields of the nodes are updated accordingly.
    // The backbone of this graph will be deleted, as it is no longer valid.
    // To avoid losing backbone information, call encodeBackbone() on both graphs before addGraph().
    void addGraph(const Graph& another);

    GraphNode* nodes;
    GraphEdge* edges;
    uint       node_count, edge_count;

    SuccinctVector* backbone;

    bool       ok;  // Also tells that edges are sorted by from.

    const static usint BACKBONE_BLOCK_SIZE = 32;
};

//--------------------------------------------------------------------------

class PathGraph
{
  public:
    PathGraph(Graph& parent);       // Creates a generation 0 PathGraph. Updates the labels of parent.
    PathGraph(PathGraph& previous); // Creates a next generation PathGraph. Invalidates previous.

    ~PathGraph();

    void printInfo();

    typedef std::vector<PathNode> nvector;
    typedef std::vector<PathEdge> evector;

    // status must be sorted. Calling this invalidates parent and
    // sets status to ready.
    // Writes outdegree to PathNode.key.second, value to PathNode.to, and
    // predecessor labels to PathNode.key.first.
    // Restores the labels of parent.
    bool generateEdges(Graph& parent);

    // status must be ready.
    std::vector<pair_type>* getSamples(usint sample_rate, usint& max_sample, const Graph& parent);

    nvector nodes;
    evector edges;
    usint   node_count, edge_count, ranks, automata;
    uint    max_label;  // Node label of initial nodes.
    usint   temp_nodes; // Total number of nodes created before sorting.

    usint   generation; // Sorted by paths of length 2^generation.

    enum status_t { error, ok, sorted, ready, edges_sorted } status;
    bool has_stabilized;  // The number of nodes will probably not explode in subsequent doublings.
    
    bool debug;

    // Can create an index by using key.second in PathNodes.
    // If the graph is not ready, its status becomes error.
    // Sorting edges by from actually sorts them by (from, to).
    void      sortEdges(bool by_from, bool create_index);
    pair_type getEdges(usint node, bool by_from); // Create index first.

  private:
    void      createPathNode(const PathNode& left, const PathNode& right);
    void      sort();
    pair_type nextMaximalSet(pair_type node_range);

    void      sortByKey();
    void      sortEdges();  // By (from.label, to.rank)


    // Can create an index by using key.second in PathNodes.
    void      sortByFrom(bool create_index);
    pair_type getNodesFrom(uint node);        // Use sortByFrom(true) first.
    pair_type getNextRange(pair_type range);  // Use sortByFrom() first.

    // As in Graph.
    void      restoreLabels();
};

//--------------------------------------------------------------------------


} // namespace CSA


#endif // GRAPH_H

