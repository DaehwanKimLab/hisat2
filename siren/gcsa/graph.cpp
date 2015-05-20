#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <sstream>
#include <deque>
#include <map>
#include <stack>
#include <vector>

#include "graph.h"
#include "gcsa.h"

#include <misc/utils.h>
#include <bits/bitbuffer.h>


namespace CSA
{

//--------------------------------------------------------------------------

Graph::Graph(std::vector<GraphNode>& node_vector, std::vector<GraphEdge>& edge_vector) :
  nodes(0), edges(0),
  node_count(node_vector.size()), edge_count(edge_vector.size()),
	backbone(0),
  ok(false)
{
  GraphNode* node = this->nodes = new GraphNode[this->node_count];
  for(std::vector<GraphNode>::iterator iter = node_vector.begin();
      iter != node_vector.end(); ++iter, ++node)
  {
    *node = *iter;
  }

  GraphEdge* edge = this->edges = new GraphEdge[this->edge_count];
  for(std::vector<GraphEdge>::iterator iter = edge_vector.begin();
      iter != edge_vector.end(); ++iter, ++edge)
  {
    *edge = *iter;
  }

  this->ok = true;
}

Graph::Graph(const std::string& base_name, bool is_full_name) :
  nodes(0), edges(0),
  node_count(0), edge_count(0),
  backbone(0),
  ok(false)
{
  std::string automaton_name = (is_full_name ? base_name : base_name + AUTOMATON_EXTENSION);
  std::ifstream infile(automaton_name.c_str(), std::ios_base::binary);
  if(!infile)
  {
    std::cerr << "Graph: Cannot open input file (" << automaton_name << ")!" << std::endl;
    return;
  }

  infile.read((char*)&(this->node_count), sizeof(this->node_count));
  infile.read((char*)&(this->edge_count), sizeof(this->edge_count));
  this->nodes = new GraphNode[this->node_count];
  infile.read((char*)(this->nodes), this->node_count * sizeof(GraphNode));
  this->edges = new GraphEdge[this->edge_count];
  infile.read((char*)(this->edges), this->edge_count * sizeof(GraphEdge));

  infile.close();
  this->sortEdges(true);  // FIXME unnecessary, as the edges are already sorted?
  this->ok = true;
}


Graph::~Graph()
{
  delete[] this->nodes; this->nodes = 0;
  delete[] this->edges; this->edges = 0;
  delete this->backbone; this->backbone = 0;
}

void
Graph::write(const std::string& base_name, bool is_full_name)
{
  if(!(this->ok)) { return; }

  std::string automaton_name = (is_full_name ? base_name : base_name + AUTOMATON_EXTENSION);
  std::ofstream outfile(automaton_name.c_str(), std::ios_base::binary);
  if(!outfile)
  {
    std::cerr << "Graph: Cannot open output file (" << automaton_name << ")!" << std::endl;
    return;
  }

  outfile.write((char*)&(this->node_count), sizeof(this->node_count));
  outfile.write((char*)&(this->edge_count), sizeof(this->edge_count));
  largeWrite(outfile, (char*)(this->nodes), this->node_count * sizeof(GraphNode), 1024);
  largeWrite(outfile, (char*)(this->edges), this->edge_count * sizeof(GraphEdge), 1024);

  outfile.close();
}

void
Graph::printInfo()
{
  std::cout << "Initial automaton (" << this->node_count << " nodes, " << this->edge_count << " edges)" << std::endl;
}

struct GraphEdgeFromComparator
{
  bool operator() (const GraphEdge& a, const GraphEdge& b) const
  {
    return (a.from < b.from);
  }
} ge_from_comparator;

struct GraphEdgeToComparator
{
  bool operator() (const GraphEdge& a, const GraphEdge& b) const
  {
    return (a.to < b.to);
  }
} ge_to_comparator;

//--------------------------------------------------------------------------

void
Graph::sortEdges(bool from)
{
  if(from)
  {
    parallelSort(this->edges, this->edges + this->edge_count, ge_from_comparator);
  }
  else
  {
    parallelSort(this->edges, this->edges + this->edge_count, ge_to_comparator);
  }
}

pair_type
Graph::getEdges(uint node, bool from)
{
  pair_type range;
  uint      temp;

  // Lower bound for the range.
  sint low = 0, high = (usint)(this->edge_count) - 1;
  while(low < high)
  {
    sint mid = low + (high - low) / 2;
    temp = (from ? this->edges[mid].from : this->edges[mid].to);
    if(node == temp)     { high = mid; }
    else if(node < temp) { high = mid - 1; }
    else                 { low = mid + 1; }
  }
  temp = (from ? this->edges[low].from : this->edges[low].to);
  if(node == temp) { range.first = low; }
  else { return EMPTY_PAIR; }

  // Upper bound for the range.
  high = (usint)(this->edge_count) - 1;
  while(low < high)
  {
    sint mid = low + (high - low + 1) / 2;
    temp = (from ? this->edges[mid].from : this->edges[mid].to);
    if(node == temp) { low = mid; }
    else             { high = mid - 1; }
  }
  range.second = high;

  return range;
}

pair_type
Graph::getNextEdgeRange(pair_type range, bool from)
{
  if(range.second >= this->edge_count - 1) { return EMPTY_PAIR; }

  if(isEmpty(range)) { range = pair_type(0, 0); } // Create the first range.
  else { range.second++; range.first = range.second; }

  if(from)
  {
    while(range.second < this->edge_count - 1 && this->edges[range.second + 1].from == this->edges[range.first].from)
    {
      range.second++;
    }
  }
  else
  {
    while(range.second < this->edge_count - 1 && this->edges[range.second + 1].to == this->edges[range.first].to)
    {
      range.second++;
    }
  }

  return range;
}

//--------------------------------------------------------------------------

void
Graph::createBackbone()
{
  if(!(this->ok) || this->backbone != 0) { return; }

  SuccinctEncoder encoder(BACKBONE_BLOCK_SIZE);
  for(usint i = 0; i < this->node_count; i++)
  {
    if(islower(this->nodes[i].label)) { this->nodes[i].label = toupper(this->nodes[i].label); }
    else { encoder.setBit(i); }
  }
  this->backbone = new SuccinctVector(encoder, this->node_count);
}

void
Graph::removeBackbone()
{
  if(!(this->ok)) { return; }

  for(usint i = 0; i < this->node_count; i++)
  {
    this->nodes[i].label = toupper(this->nodes[i].label);
  }
}

void
Graph::encodeBackbone()
{
  if(!(this->ok) || this->backbone == 0) { return; }

  SuccinctVector::Iterator iter(*(this->backbone));
  for(usint i = 0; i < this->node_count; i++)
  {
    if(!(iter.isSet(i))) { this->nodes[i].label = tolower(this->nodes[i].label); }
  }
}

//--------------------------------------------------------------------------

struct TempNode
{
  std::basic_string<uint> nodes;
  uint id, label, value;
  bool backbone, potential;

  TempNode(const GraphNode& original, uint original_id, bool _backbone) :
    id(0), label(original.label), value(original.value), backbone(_backbone), potential(false)
  {
    this->nodes.push_back(original_id);
  }

  GraphNode getNode() const
  {
    return GraphNode(label, value);
  }
};

struct TempEdge
{
  TempNode* from;
  TempNode* to;

  TempEdge(TempNode* _from, TempNode* _to) : from(_from), to(_to) {}

  GraphEdge getEdge() const
  {
    return GraphEdge(from->id, to->id);
  }
};

struct NodeLabelIdComparator
{
  bool operator() (const GraphNode* a, const GraphNode* b) const
  {
    return ((a->label < b->label) || (a->label == b->label && a < b));
  }
} nli_comparator;

struct TempEdgeFromComparator
{
  bool operator() (const TempEdge& a, const TempEdge& b) const
  {
    return (a.from < b.from);
  }
} tef_comparator;

struct TempEdgeToComparator
{
  bool operator() (const TempEdge& a, const TempEdge& b) const
  {
    return (a.to < b.to);
  }
} tet_comparator;

void
Graph::determinize()
{
  if(!(this->ok)) { return; }

  bool create_backbone = (this->backbone != 0);
  uint backbone_nodes = (create_backbone ? this->backbone->getNumberOfItems() : 0);

  // Find final node(s) and put them into active queue.
  uint* outdegrees = new uint[this->node_count];
  for(uint i = 0; i < this->node_count; i++) { outdegrees[i] = 0; }
  for(uint i = 0; i < this->edge_count; i++) { outdegrees[this->edges[i].from]++; }
  std::deque<TempNode*> active;
  std::map<std::basic_string<uint>, TempNode*> new_nodes;
  std::vector<TempEdge> new_edges;
  for(uint i = 0; i < this->node_count; i++)
  {
    if(outdegrees[i] == 0)
    {
      TempNode* temp = new TempNode(this->nodes[i], i, create_backbone);
      active.push_back(temp);
      new_nodes[temp->nodes] = temp;
    }
  }
  delete[] outdegrees; outdegrees = 0;


  // Actual work.
  this->sortEdges(false); // Sort edges by destination node.
  while(!(active.empty()))
  {
    // Find predecessors of all original nodes contained in the next active node.
    TempNode* curr = active.front(); active.pop_front();
    std::vector<GraphNode*> predecessors;
    for(std::basic_string<uint>::iterator iter = curr->nodes.begin(); iter != curr->nodes.end(); ++iter)
    {
      pair_type edge_range = this->getEdges(*iter, false);
      for(uint i = edge_range.first; i <= edge_range.second; i++)
      {
        predecessors.push_back(this->nodes + this->edges[i].from);
      }
    }
    removeDuplicates(predecessors, false);

    // Build new tempNodes from the predecessors.
    sequentialSort(predecessors.begin(), predecessors.end(), nli_comparator);
    SuccinctVector::Iterator* bb_iter = 0;
    if(create_backbone) { bb_iter = new SuccinctVector::Iterator(*(this->backbone)); }
    for(std::vector<GraphNode*>::iterator iter = predecessors.begin(); iter != predecessors.end(); )
    {
      GraphNode* first = *iter; ++iter;
      uint threshold = (create_backbone ? backbone_nodes : curr->value);
      TempNode* temp = new TempNode(*first, first - this->nodes, false);
      if(create_backbone) { temp->potential = bb_iter->isSet(first - this->nodes); }
      while(iter != predecessors.end() && (*iter)->label == first->label) // Other original predecessors.
      {
        uint new_value = temp->value;
        if((temp->value < threshold) != ((*iter)->value < threshold))
        {
          new_value = std::min(temp->value, (*iter)->value);
        }
        else { new_value = std::max(temp->value, (*iter)->value); }
        if(create_backbone)
        {
          temp->potential = bb_iter->isSet(*iter - this->nodes) |
            (new_value == temp->value ? temp->potential : 0);
        }
        temp->value = new_value;
        temp->nodes.push_back(*iter - this->nodes);
        ++iter;
      }
      std::map<std::basic_string<uint>, TempNode*>::iterator existing = new_nodes.find(temp->nodes);
      if(existing == new_nodes.end()) // Create new node and make it active.
      {
        new_nodes[temp->nodes] = temp;
        active.push_back(temp);
        new_edges.push_back(TempEdge(temp, curr));
      }
      else
      {
        delete temp;
        new_edges.push_back(TempEdge((*existing).second, curr));
      }
      curr->id++; // Increase indegree.
    }
    delete bb_iter; bb_iter = 0;
  }

  // Create the new backbone.
  if(create_backbone)
  {
    for(std::map<std::basic_string<uint>, TempNode*>::iterator iter = new_nodes.begin(); iter != new_nodes.end(); ++iter)
    {
      TempNode* node = (*iter).second;
      if(node->backbone) { active.push_back(node); }
    }
    parallelSort(new_edges.begin(), new_edges.end(), tet_comparator);
    while(!(active.empty()))
    {
      TempNode* curr = active.front(); active.pop_front();
      TempEdge dummy(curr, curr);
      std::pair<std::vector<TempEdge>::iterator, std::vector<TempEdge>::iterator> edge_range =
        equal_range(new_edges.begin(), new_edges.end(), dummy, tet_comparator);
      while(edge_range.first != edge_range.second)
      {
        TempNode* node = edge_range.first->from;
        if(node->potential && curr->value == node->value + 1)
        {
          node->backbone = true; node->potential = false;
          active.push_back(node);
        }
        ++(edge_range.first);
      }
    }
  }

  // Topological sort. Create new nodes in sorted order.
  delete[] this->nodes; this->nodes = new GraphNode[new_nodes.size()]; this->node_count = 0;
  SuccinctEncoder* encoder = 0;
  if(create_backbone) { encoder = new SuccinctEncoder(BACKBONE_BLOCK_SIZE); }
  parallelSort(new_edges.begin(), new_edges.end(), tef_comparator);
  for(std::map<std::basic_string<uint>, TempNode*>::iterator iter = new_nodes.begin(); iter != new_nodes.end(); ++iter)
  {
    TempNode* node = (*iter).second;
    if(node->id == 0)
    {
      active.push_back(node);
      node->id = this->node_count;
      this->nodes[this->node_count] = node->getNode();
      if(create_backbone && node->backbone) { encoder->setBit(this->node_count); }
      this->node_count++;
    }
  }
  while(!(active.empty()))
  {
    TempNode* curr = active.front(); active.pop_front();
    TempEdge dummy(curr, curr);
    std::pair<std::vector<TempEdge>::iterator, std::vector<TempEdge>::iterator> edge_range =
      equal_range(new_edges.begin(), new_edges.end(), dummy, tef_comparator);
    while(edge_range.first != edge_range.second)
    {
      TempNode* node = edge_range.first->to;
      node->id--;
      if(node->id == 0)
      {
        active.push_back(node);
        node->id = this->node_count;
        this->nodes[this->node_count] = node->getNode();
        if(create_backbone && node->backbone) { encoder->setBit(this->node_count); }
        this->node_count++;
      }
      ++(edge_range.first);
    }
  }
  if(create_backbone)
  {
    delete this->backbone;
    this->backbone = new SuccinctVector(*encoder, this->node_count);
    usint new_bb = this->backbone->getNumberOfItems();
    if(backbone_nodes != new_bb)
    {
      std::cerr << "Error: Backbone nodes: " << backbone_nodes << " -> " << new_bb << std::endl;
    }
    delete encoder; encoder = 0;
  }

  // Create new edges.
  delete[] this->edges; this->edges = new GraphEdge[new_edges.size()]; this->edge_count = 0;
  for(std::vector<TempEdge>::iterator iter = new_edges.begin(); iter != new_edges.end(); ++iter)
  {
    this->edges[this->edge_count] = iter->getEdge(); this->edge_count++;
  }
  this->sortEdges(true);  // Edges must be sorted by from.

  // Delete TempNodes.
  for(std::map<std::basic_string<uint>, TempNode*>::iterator iter = new_nodes.begin(); iter != new_nodes.end(); ++iter)
  {
    delete (*iter).second;
  }
}

bool
Graph::isDeterministic()
{
  if(!(this->ok)) { return false; }

  this->sortEdges(false); // Sort edges by destination node.
  pair_type range = EMPTY_PAIR;
  for(range = this->getNextEdgeRange(range, false); !isEmpty(range); range = this->getNextEdgeRange(range, false))
  {
    WriteBuffer buf(BITS_TO_WORDS(CHARS));
    for(usint i = range.first; i <= range.second; i++)
    {
      uint c = this->nodes[this->edges[i].from].label;
      if(buf.isSet(c)) { return false; }
      buf.setBit(c);
    }
  }

  return true;
}

//--------------------------------------------------------------------------

uint
Graph::maxLabel()
{
  uint max_label = 0;
  for(usint i = 0; i < this->node_count; i++) { max_label = std::max(max_label, this->nodes[i].label); }
  return max_label;
}


void
Graph::changeLabels(uint automata, uint max_label)
{
  uint inits = 0, finals = 0;
  max_label += automata;
  for(uint i = 0; i < this->node_count; i++)
  {
    if(this->nodes[i].label == 0) { this->nodes[i].label += finals; finals++; }
    else
    {
      this->nodes[i].label += automata;
      if(this->nodes[i].label == max_label) { this->nodes[i].label += inits; inits++; }
    }
  }
}

void
Graph::restoreLabels(uint automata, uint max_label)
{
  for(usint i = 0; i < this->node_count; i++)
  {
    if(this->nodes[i].label < automata) { this->nodes[i].label = 0; }
    else
    {
      this->nodes[i].label -= automata;
      if(this->nodes[i].label >= max_label) { this->nodes[i].label = max_label; }
    }
  }
}

//--------------------------------------------------------------------------

void
Graph::addGraph(const Graph& another)
{
  if(!(this->ok) || !(another.ok)) { return; }

  bool update_values = (this->backbone != 0 && another.backbone != 0);
  GraphNode* new_nodes = new GraphNode[this->node_count + another.node_count];

  uint threshold = (update_values ? this->backbone->getNumberOfItems() : 0);
  uint aligned_adjustment = 0;
  uint nonaligned_adjustment = (update_values ? another.backbone->getNumberOfItems() : 0);
  for(usint i = 0; i < this->node_count; i++)
  {
    new_nodes[i] = this->nodes[i];
    new_nodes[i].value += (new_nodes[i].value >= threshold ? nonaligned_adjustment : aligned_adjustment);
  }

  threshold = (update_values ? another.backbone->getNumberOfItems() : 0);
  aligned_adjustment = (update_values ? this->backbone->getNumberOfItems() : 0);
  nonaligned_adjustment = (update_values ? this->backbone->getNumberOfItems() : 0);
  for(usint i = 0, j = this->node_count; i < another.node_count; i++, j++)
  {
    new_nodes[j] = another.nodes[i];
    new_nodes[j].value += (new_nodes[j].value >= threshold ? nonaligned_adjustment : aligned_adjustment);
  }
  delete[] this->nodes; this->nodes = new_nodes;

  GraphEdge* new_edges = new GraphEdge[this->edge_count + another.edge_count];
  for(usint i = 0; i < this->edge_count; i++) { new_edges[i] = this->edges[i]; }
  for(usint i = 0; i < another.edge_count; i++)
  {
    new_edges[i + this->edge_count] = GraphEdge(another.edges[i], this->node_count);
  }
  delete[] this->edges; this->edges = new_edges;

  delete this->backbone; this->backbone = 0;
  this->node_count += another.node_count;
  this->edge_count += another.edge_count;
}

//--------------------------------------------------------------------------

PathGraph::PathGraph(Graph& parent) :
  node_count(0), edge_count(0), ranks(0), automata(0), max_label(0),
  temp_nodes(0),
  generation(0),
  status(error), has_stabilized(false)
{
  if(!(parent.ok)) { return; }

  std::vector<uint> end_markers;
  for(uint i = 0; i < parent.node_count; i++)
  {
    if(parent.nodes[i].label == 0) { this->automata++; end_markers.push_back(i); }
    this->max_label = std::max(this->max_label, parent.nodes[i].label);
  }
  if(this->automata == 0)
  {
    std::cerr << "Error: The graph has no final nodes!" << std::endl;
    return;
  }
  parent.changeLabels(this->automata, this->max_label);

  // One PathNode per edge and one extra for each of the end markers.
  // First this->automata key values are for end markers.
  // Initial nodes also have distinct keys.
  this->temp_nodes = this->node_count = parent.edge_count + this->automata;
  this->nodes.reserve(this->node_count);
  for(GraphEdge* edge = parent.edges; edge < parent.edges + parent.edge_count; ++edge)
  {
    PathNode pn;
    pn.from = edge->from; pn.to = edge->to;
    pn.key = pair_type(parent.nodes[edge->from].label, 0);
    this->nodes.push_back(pn);
  }
  for(uint i = 0; i < this->automata; i++)  // Final nodes.
  {
    PathNode pn;
    pn.from = pn.to = end_markers[i];
    pn.key = pair_type(i, 0);
    this->nodes.push_back(pn);
  }
    
    // daehwan - for debugging purposes
    debug = parent.node_count <= 20;
    
    if(debug) {
        std::cerr << "Nodes:" << std::endl;
        for(size_t i = 0; i < parent.node_count; i++) {
            const GraphNode& node = parent.nodes[i];
            std::cerr << "\t" << i << "\t" << (char)node.label << "\t" << node.value << std::endl;
        }
        std::cerr << std::endl;
        std::cerr << "Edges: " << std::endl;
        for(size_t i = 0; i < parent.edge_count; i++) {
            const GraphEdge& edge = parent.edges[i];
            std::cerr << "\t" << i << "\t" << edge.from << " --> " << edge.to << std::endl;
        }
        std::cerr << std::endl;
    }

  parent.restoreLabels(this->automata, this->max_label);
  this->status = ok;
  this->sort();
    
}

PathGraph::PathGraph(PathGraph& previous) :
  node_count(0), edge_count(0), ranks(0), automata(previous.automata), max_label(previous.max_label),
  temp_nodes(0),
  generation(previous.generation + 1),
  status(error), has_stabilized(false)
{
  if(previous.status != ok) { return; }
    
    // daehwan - for debugging purposes
    debug = previous.debug;

  previous.sortByFrom(true);

  // A heuristic to determine, whether the number of new nodes should be counted.
  usint new_nodes = previous.node_count + previous.node_count / 8;
  if(previous.ranks >= previous.node_count / 2 && !(previous.has_stabilized))
  {
    new_nodes = 0;
    for(nvector::iterator left = previous.nodes.begin(); left != previous.nodes.end(); ++left)
    {
      if(left->isSorted()) { new_nodes++; continue; }
      pair_type pn_range = previous.getNodesFrom(left->to);
      new_nodes += length(pn_range);
    }
    std::cout << "Trying to allocate space for " << new_nodes << " nodes." << std::endl;
  }
  this->nodes.reserve(new_nodes);

  for(nvector::iterator left = previous.nodes.begin(); left != previous.nodes.end(); ++left)
  {
    if(left->isSorted()) { this->nodes.push_back(*left); continue; }

    pair_type pn_range = previous.getNodesFrom(left->to);
    for(usint pn = pn_range.first; pn <= pn_range.second; pn++)
    {
      this->createPathNode(*left, previous.nodes[pn]);
    }
  }
  this->temp_nodes = this->node_count = this->nodes.size();

  this->status = ok;
  this->sort();

  if(previous.has_stabilized || (this->generation >= 11 && this->ranks >= 0.8 * new_nodes))
  {
    this->has_stabilized = true;
  }
}

PathGraph::~PathGraph()
{
  this->status = error;
}

void
PathGraph::printInfo()
{
  std::cout << "Generation " << this->generation << " (" << this->temp_nodes << " -> " << this->node_count << " nodes, " << this->ranks << " ranks)" << std::endl;
}

//--------------------------------------------------------------------------

bool
PathGraph::generateEdges(Graph& parent)
{
  if(this->status != sorted) { return false; }

  this->sortByFrom(false);  // Sort nodes by from, do not create index.
    
  parent.sortEdges(false);  // Sort parent edges by to.
  parent.changeLabels(this->automata, this->max_label);

  // Create edges (from, to) as (from.from, to = rank(to)).
  pair_type pn_range = this->getNextRange(EMPTY_PAIR);
  pair_type ge_range = parent.getNextEdgeRange(EMPTY_PAIR, false);
  this->edges.reserve(this->node_count + this->node_count / 4);
  while(!isEmpty(pn_range) && !isEmpty(ge_range))
  {
    if(this->nodes[pn_range.first].from == parent.edges[ge_range.first].to)
    {
      for(usint node = pn_range.first; node <= pn_range.second; node++)
      {
        for(usint edge = ge_range.first; edge <= ge_range.second; edge++)
        {
          uint from = parent.edges[edge].from;
          this->edges.push_back(PathEdge(from, this->nodes[node].key.first, parent.nodes[from].label));
        }
      }
      pn_range = this->getNextRange(pn_range);
      ge_range = parent.getNextEdgeRange(ge_range, false);
    }
    else if(this->nodes[pn_range.first].from < parent.edges[ge_range.first].to)
    {
      pn_range = this->getNextRange(pn_range);
    }
    else
    {
      ge_range = parent.getNextEdgeRange(ge_range, false);
    }
  }
  this->edge_count = this->edges.size();

  this->sortByKey(); // Restore correct node order.
  this->sortEdges(); // Sort edges by (from.label, to.rank).

  // Sets PathNode.to = GraphNode.value and PathNode.key.first to outdegree.
  // Replaces (from.from, to) with (from, to)
  nvector::iterator node = this->nodes.begin(); node->key.first = 0;
  evector::iterator edge = this->edges.begin();
  while(node != this->nodes.end() && edge != this->edges.end())
  {
    if(edge->from == node->from)
    {
      edge->from = node - this->nodes.begin(); ++edge;
      node->key.first++;
    }
    else
    {
      node->to = parent.nodes[node->from].value;
      ++node; node->key.first = 0;
    }
  }
  if(node != this->nodes.end())
  {
    node->to = parent.nodes[node->from].value;
  }
    
    if(debug) {
        // this->restoreLabels();
        std::cerr << "Path edges" << std::endl;
        for(size_t i = 0; i < edge_count; i++) {
            const PathEdge& edge = edges[i];
            std::cerr << "\t" << i << "\tfrom: " << edge.from << "\tranking: " << edge.rank << "\t" << (char)edge.label << std::endl;
        }
    }
    
    if(debug) {
        std::cerr << "Path nodes" << std::endl;
        for(size_t i = 0; i < node_count; i++) {
            const PathNode& node = nodes[i];
            std::cerr << "\t" << i << "\t(" << node.key.first << ", " << node.key.second << ")\t"
            << node.from << " --> " << node.to << std::endl;
        }
    }

  this->sortEdges(false, true);
    
  parent.restoreLabels(this->automata, this->max_label);
  this->restoreLabels();
  this->status = ready;
    
    
  return true;
}

//--------------------------------------------------------------------------

std::vector<pair_type>*
PathGraph::getSamples(usint sample_rate, usint& max_sample, const Graph& parent)
{
  if(this->status != ready && this->status != edges_sorted) { return 0; }

  this->sortEdges(true, true);
  WriteBuffer* found_nodes = new WriteBuffer(this->node_count, 1);
  ReadBuffer* found = found_nodes->getReadBuffer();
  WriteBuffer* sampled_nodes = new WriteBuffer(this->node_count, 1);
  ReadBuffer* sampled = sampled_nodes->getReadBuffer();
  std::vector<pair_type>* sample_pairs = new std::vector<pair_type>;

  std::stack<usint> active; // Push the initial nodes into stack.
  for(usint i = 1; i <= this->automata; i++) { active.push(this->node_count - i); }
  usint unsampled = 0, current = 0;
  max_sample = 0;
  usint total_found = 0;
  while(!active.empty())
  {
    current = active.top(); active.pop();
    unsampled = 1;

    while(true)
    {
      if(found->readItem(current))  // Nodes with multiple predecessors must be sampled.
      {
        if(!sampled->readItem(current))
        {
          sample_pairs->push_back(pair_type(current, this->nodes[current].value()));
          sampled_nodes->goToItem(current); sampled_nodes->writeItem(1);
          max_sample = std::max(max_sample, (usint)(this->nodes[current].value()));
        }
        break;
      }
      found_nodes->goToItem(current); found_nodes->writeItem(1); total_found++;

      pair_type successors = this->getEdges(current, true);
      // No nearby samples.
      // Multiple or no outgoing edges.
      // Discontinuity in values.
      if(unsampled >= sample_rate || length(successors) != 1 ||
         this->nodes[this->edges[successors.first].rank].value() != this->nodes[current].value() + 1)
      {
        sample_pairs->push_back(pair_type(current, this->nodes[current].value()));
        sampled_nodes->goToItem(current); sampled_nodes->writeItem(1);
        max_sample = std::max(max_sample, (usint)(this->nodes[current].value()));
        for(usint suc = successors.first + 1; suc <= successors.second; suc++)
        {
          active.push(this->edges[suc].rank);
        }
        unsampled = 0;
      }
      if(isEmpty(successors)) { break; }

      current = this->edges[successors.first].rank;
      unsampled++;
    }
  }
  delete found_nodes; delete found;
  delete sampled_nodes; delete sampled;

//  std::cout << "Found " << total_found << " nodes" << std::endl;
  return sample_pairs;
}

//--------------------------------------------------------------------------

void
PathGraph::createPathNode(const PathNode& left, const PathNode& right)
{
  PathNode pn;
  pn.from = left.from; pn.to = right.to;
  pn.key = pair_type(left.key.first, right.key.first);
  this->nodes.push_back(pn);
}

void
PathGraph::sort()
{
  this->sortByKey();

  // Update ranks.
  usint rank = 0;
  pair_type key = this->nodes[0].key;
  for(nvector::iterator iter = this->nodes.begin(); iter != this->nodes.end(); ++iter)
  {
    if(iter->key != key) { rank++; key = iter->key; }
    iter->key = pair_type(rank, 0);
  }
    
    // daehwan - for debugging purposes
    if(debug) {
        std::cerr << "Path nodes (" << generation << "-generation) before merge" << std::endl;
        for(size_t i = 0; i < node_count; i++) {
            const PathNode& node = nodes[i];
            std::cerr << "\t" << i << "\t(" << node.key.first << ", " << node.key.second << ")\t"
                      << node.from << " --> " << node.to << std::endl;
        }
    }

  // Merge equivalent nodes.
  usint top = 0;
  pair_type node_range = EMPTY_PAIR;
  while(true)
  {
    node_range = this->nextMaximalSet(node_range);
    if(isEmpty(node_range)) { break; }
    this->nodes[top] = this->nodes[node_range.first]; top++;
  }
  this->node_count = top; this->nodes.resize(this->node_count);

  // Check if sorted.
  bool is_sorted = true;
  PathNode* candidate = &(this->nodes[0]);
  key = candidate->key; this->ranks = 1;
  for(usint i = 1; i < this->node_count; i++)
  {
    if(this->nodes[i].key != key)
    {
      if(candidate != 0) { candidate->setSorted(); }
      candidate = &(this->nodes[i]);
      key = candidate->key; this->ranks++;
    }
    else { candidate = 0; is_sorted = false; }
  }
  if(candidate != 0) { candidate->setSorted(); }

  // Replace the ranks of a sorted graph so that rank(i) = i.
  // Merges may otherwise leave gaps in the ranks.
  if(is_sorted)
  {
    for(usint i = 0; i < this->node_count; i++)
    {
      this->nodes[i].key.first = i;
    }
    this->status = sorted;
  }
    
    if(debug) {
        std::cerr << "Path nodes (" << generation << "-generation) after merge" << std::endl;
        for(size_t i = 0; i < node_count; i++) {
            const PathNode& node = nodes[i];
            std::cerr << "\t" << i << "\t(" << node.key.first << ", " << node.key.second << ")\t"
                      << node.from << " --> " << node.to << (node.isSorted() ? "\tsorted" : "") << std::endl;
        }
    }
}

// Returns the next maximal mergeable set of PathNodes.
// A set of PathNodes sharing adjacent keys is mergeable, if each of the
// PathNodes begins in the same GraphNode, and no other PathNode shares
// the key. If the maximal set is empty, returns the next PathNode.
pair_type
PathGraph::nextMaximalSet(pair_type range)
{
  if(range.second >= this->node_count - 1) { return EMPTY_PAIR; }

  if(isEmpty(range)) { range = pair_type(0, 0); }
  else
  {
    range.first = range.second + 1; range.second = range.first;
    if(this->nodes[range.first - 1].key == this->nodes[range.first].key) { return range; }
  }

  for(usint i = range.second + 1; i < this->node_count; i++)
  {
    if(this->nodes[i - 1].key != this->nodes[i].key) { range.second = i - 1; }
    if(this->nodes[i].from != this->nodes[range.first].from) { return range; }
  }

  range.second = this->node_count - 1;
  return range;
}

//--------------------------------------------------------------------------

struct PathNodeComparator
{
  bool operator() (const PathNode& a, const PathNode& b) const
  {
    return (a.key < b.key);
  }
} pn_comparator;

void
PathGraph::sortByKey()
{
  parallelSort(this->nodes.begin(), this->nodes.end(), pn_comparator);
}

struct PathEdgeComparator
{
  bool operator() (const PathEdge& a, const PathEdge& b) const
  {
    return ((a.label < b.label) ||
            ((a.label == b.label) && (a.rank < b.rank)));
  }
} pe_comparator;

void
PathGraph::sortEdges()
{
  parallelSort(this->edges.begin(), this->edges.end(), pe_comparator);
}

//--------------------------------------------------------------------------

struct PathEdgeFromComparator
{
  bool operator() (const PathEdge& a, const PathEdge& b) const
  {
    return ((a.from < b.from) || (a.from == b.from && a.rank < b.rank));
  }
} pe_from_comparator;

struct PathEdgeToComparator
{
  bool operator() (const PathEdge& a, const PathEdge& b) const
  {
    return (a.rank < b.rank);
  }
} pe_to_comparator;

void
PathGraph::sortEdges(bool by_from, bool create_index)
{
  if(by_from)
  {
    if(this->status != edges_sorted)
    {
      parallelSort(this->edges.begin(), this->edges.end(), pe_from_comparator);
      if(this->status == ready) { this->status = edges_sorted; }
      else                      { this->status = error; }
    }
  }
  else
  {
    parallelSort(this->edges.begin(), this->edges.end(), pe_to_comparator);
    this->status = error;
  }

  if(create_index)
  {
    for(nvector::iterator iter = this->nodes.begin(); iter != this->nodes.end(); ++iter)
    {
      iter->key.second = 0;
    }
    if(by_from)
    {
      for(evector::iterator iter = this->edges.begin(); iter != this->edges.end(); ++iter)
      {
        this->nodes[iter->from].key.second++;
      }
    }
    else
    {
      for(evector::iterator iter = this->edges.begin(); iter != this->edges.end(); ++iter)
      {
        this->nodes[iter->rank].key.second++;
      }
    }
    for(usint i = 1; i < this->node_count; i++)
    {
      this->nodes[i].key.second += this->nodes[i - 1].key.second;
    }
  }
}

pair_type
PathGraph::getEdges(usint node, bool by_from)
{
  if(node >= this->node_count)
  {
    std::cout << "Error: Trying to get edges " << (by_from ? "from " : "to ") << node << std::endl;
  }
  if(this->nodes[node].key.second == 0) { return EMPTY_PAIR; }
  if(node == 0) { return pair_type(0, this->nodes[node].key.second - 1); }
  return pair_type(this->nodes[node - 1].key.second, this->nodes[node].key.second - 1);
}

//--------------------------------------------------------------------------

struct PathNodeFromComparator
{
  bool operator() (const PathNode& a, const PathNode& b) const
  {
    return (a.from < b.from);
  }
} pn_from_comparator;

void
PathGraph::sortByFrom(bool create_index)
{
  parallelSort(this->nodes.begin(), this->nodes.end(), pn_from_comparator);
  this->status = error;

  if(create_index)
  {
    usint current = 0;
    this->nodes[0].key.second = 0;
    for(nvector::iterator iter = this->nodes.begin(); iter != this->nodes.end(); ++iter)
    {
      while(current < iter->from)
      {
        current++;
        this->nodes[current].key.second = iter - this->nodes.begin();
      }
    }
    for(current++; current < this->node_count; current++)
    {
      this->nodes[current].key.second = this->nodes.size();
    }
  }
}

pair_type
PathGraph::getNodesFrom(uint node)
{
  if(node + 1 >= this->node_count) { return pair_type(this->nodes[node].key.second, this->node_count - 1); }
  if(this->nodes[node].key.second == this->nodes[node + 1].key.second) { return EMPTY_PAIR; }
  return pair_type(this->nodes[node].key.second, this->nodes[node + 1].key.second - 1);
}

pair_type
PathGraph::getNextRange(pair_type range)
{
  if(range.second >= this->node_count - 1) { return EMPTY_PAIR; }

  if(isEmpty(range)) { range = pair_type(0, 0); } // Create the first range.
  else { range.second++; range.first = range.second; }

  while(range.second < this->node_count - 1 && this->nodes[range.second + 1].from == this->nodes[range.first].from)
  {
    range.second++;
  }

  return range;
}

//--------------------------------------------------------------------------

void
PathGraph::restoreLabels()
{
  for(usint i = 0; i < this->edge_count; i++)
  {
    if(this->edges[i].label < this->automata) { this->edges[i].label = 0; }
    else
    {
      this->edges[i].label -= this->automata;
      if(this->edges[i].label >= this->max_label) { this->edges[i].label = this->max_label; }
    }
  }
}

//--------------------------------------------------------------------------

} // namespace CSA
