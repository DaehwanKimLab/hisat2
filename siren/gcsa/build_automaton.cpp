#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include <misc/utils.h>

#include "graph.h"


using namespace CSA;


/*
  Builds an automaton from aligned DNA sequences. Accepted characters are the
  bases ('A', 'C', 'G', 'T'), 'N' denoting any base, and '-' denoting a gap.

  Input format: n lines of m characters each, separated by '\n'. Use
  clean_alignment first to reformat the alignment into supported format.

  The program considers a context of k bases after each position.
  Two bases are merged to a common state, if they

    1) are identical,
    2) occur in the same position in the alignment, and
    3) share a context

  The reason for using post-context is that the automaton must be
  deterministic in reverse direction.
*/

const std::string alphabet = "ACGNTZ";
const usint alphabet_size = 7;

struct Node
{
  char  label;  // Character value for the node.
  usint key;    // Label + context.
  usint context_filter;

  usint  seq, pos;     // Current coordinates.
  usint  context_end;  // Last pos of the context.
  char*  data;
  usint  sequences;

  usint* positions;

  void init(char* _data, usint _sequences, usint _seq, usint context_length, usint* _positions)
  {
    this->context_filter = alphabet_size;
    this->seq = _seq; this->pos = 0;
    this->context_end = 0;
    this->data = _data; this->sequences = _sequences;
    this->positions = _positions;

    usint ptr = this->pos * (this->sequences + 1) + this->seq;
    this->label = this->data[ptr];
    
    for(usint i = 0; i < context_length; i++) { this->context_filter *= alphabet_size; }
    this->key = this->rank(this->label);
    for(usint i = 0; i < context_length; i++) { this->updateKey(); }
  }

  void print()
  {
    std::cout << "(" << this->label << ", " << this->key << ", " << this->context_filter << ", " << this->seq << ", " << this->pos << ", " << this->context_end << ")" << std::endl;
  }

  usint getPtr()     { return this->pos * (this->sequences + 1) + this->seq; }
  usint getContext() { return this->context_end * (this->sequences + 1) + this->seq; }
  bool atEnd()       { return (this->label == '\0'); }
  bool isGap()       { return (this->label == '-'); }

  GraphNode getNode() { return GraphNode(this->label, this->positions[this->pos]); }

  usint rank(char c)
  {
    usint temp = alphabet.find(c);
    return ((temp != std::string::npos) ? temp + 1 : 0);
  }

  void advance()
  {
    if(this->atEnd()) { return; }
    this->pos++;
    this->label = this->data[this->getPtr()];
    if(this->isGap()) { return; }

    this->updateKey();
  }

  void updateKey()
  {
    usint ptr = this->getContext();
    this->key = (alphabet_size * this->key) % this->context_filter;
    while(this->data[ptr] != '\0')
    {
      ptr += this->sequences + 1; this->context_end++;
      if(this->data[ptr] != '-') { break; }
    }
    this->key += this->rank(this->data[ptr]);
  }
};


struct NodeComparator
{
  bool operator() (const Node* a, const Node* b)
  {
    return ((a->key < b->key) || (a->key == b->key && a->seq < b->seq));
  }
} node_context_comparator;


void
buildAutomaton(std::ifstream& input, const std::string& base_name, usint sequences, usint context_length, bool backbone)
{
  usint  lines = fileSize(input) / (sequences + 1) + 2;
  usint  size = lines * (sequences + 1);
  char* data = new char[size];
  char* real_data = data + sequences + 1;
  char* data_end = data + size;
  char* real_end = data_end - (sequences + 1);

  // Read the input file.
  input.clear();
  input.seekg(0, std::ios::beg);
  for(usint i = 0; i < sequences; i++) { data[i] = 'Z'; }
  input.read(real_data, real_end - real_data);
  for(char* ptr = real_end; ptr < data_end; ++ptr) { *ptr = '\0'; }
  input.close();

  // Map positions in the multiple alignment to positions in the first sequence (backbone).
  // Initial node 0, sequence 1..n, final node n+1.
  // Alignment position that correspond to no position in the backbone are mapped to
  // separate positions starting from n+2 (or more).
  usint* positions = new usint[lines];
  if(backbone)
  {
    for(usint i = 0, pos = 0, extra = 0; i < lines; i++)
    {
      if(data[i * (sequences + 1)] == '-') { positions[i] = lines + extra; extra++; }
      else { positions[i] = pos; pos++; }
    }
  }
  else { for(usint i = 0; i < lines; i++) { positions[i] = i; } }

  Node nodes[sequences];
  usint previous[sequences]; // Previous node id for every sequence.
  for(usint i = 0; i < sequences; i++)
  {
    nodes[i].init(data, sequences, i, context_length, positions);
    previous[i] = 0;
  }
  std::vector<GraphNode> node_vector;
  std::vector<pair_type> edge_buffer;
  std::vector<GraphEdge> edge_vector;
  node_vector.push_back(GraphNode('Z', positions[0])); // Initial node.

  for(usint line = 1; line < lines; line++)
  {
    // Advance all nodes.
    Node* active[sequences];
    usint active_nodes = 0;
    for(usint i = 0; i < sequences; i++)
    {
      nodes[i].advance();
      if(!(nodes[i].isGap())) { active[active_nodes] = nodes + i; active_nodes++; }
    }

    // Create nodes and temporary edges.
    sequentialSort(active, active + active_nodes, node_context_comparator);
    edge_buffer.clear();
    for(usint i = 0; i < active_nodes; i++)
    {
      if(i == 0 || active[i]->key != active[i - 1]->key)
      {
        GraphNode n = active[i]->getNode();
        if(backbone && active[i]->seq > 0) { n.label = tolower(n.label); }
        node_vector.push_back(n);
      }
      edge_buffer.push_back(pair_type(previous[active[i]->seq], node_vector.size() - 1));
      previous[active[i]->seq] = node_vector.size() - 1;
    }

    // Create edges.
    sequentialSort(edge_buffer.begin(), edge_buffer.end());
    for(usint i = 0; i < edge_buffer.size(); i++)
    {
      if(i == 0 || edge_buffer[i] != edge_buffer[i - 1])
      {
        edge_vector.push_back(GraphEdge(edge_buffer[i].first, edge_buffer[i].second));
      }
    }
  }

  std::cout << "Lines: " << (lines - 2) << std::endl;
  std::cout << "Nodes: " << node_vector.size() << std::endl;
  std::cout << "Edges: " << edge_vector.size() << std::endl;
  std::cout << std::endl;

  delete[] data;
  delete[] positions;
  Graph graph(node_vector, edge_vector);
  graph.write(base_name);
}


int
main(int argc, char** argv)
{
  std::cout << "Automaton builder" << std::endl;
  std::cout << std::endl;
  if(argc < 3)
  {
    std::cout << "Usage: build_automaton [-b] alignment_file context_length" << std::endl;
		std::cout << "  -b  create backbone information" << std::endl;
    return 1;
  }

  bool backbone = false;
  usint name_arg = 1, context_arg = 2;
	if(argc >= 4 && argv[1][0] == '-' && argv[1][1] == 'b')
	{
    backbone = true; name_arg = 2; context_arg = 3;
	}

  std::cout << "Alignment file: " << argv[name_arg] << std::endl;
  std::string base_name = argv[name_arg];
  std::ifstream alignment_file(argv[name_arg], std::ios_base::binary);
  if(!alignment_file)
  {
    std::cerr << "Error opening alignment file!" << std::endl;
    return 2;
  }

  usint context_length = atoi(argv[context_arg]);
  std::cout << "Context length: " << context_length << std::endl;

  std::string line;
  std::getline(alignment_file, line);
  usint sequences = line.size();
  std::cout << "Number of sequences: " << sequences << std::endl;
  std::cout << std::endl;

  double start = readTimer();
  buildAutomaton(alignment_file, base_name, sequences, context_length, backbone);
  double time = readTimer() - start;
  std::cout << "Used " << time << " seconds." << std::endl;
  std::cout << "Memory: " << memoryUsage() << " kB" << std::endl;
  std::cout << std::endl;

  return 0;
}
