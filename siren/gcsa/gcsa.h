#ifndef GCSA_H
#define GCSA_H

#include <fstream>
#include <vector>

#include <bits/deltavector.h>
#include <bits/rlevector.h>
#include <alphabet.h>

#include "graph.h"


namespace CSA
{


const std::string AUTOMATON_EXTENSION = ".automaton";
const std::string GCSA_EXTENSION = ".gcsa";
const std::string BACKBONE_EXTENSION = ".backbone";


class Backbone;

class GCSA
{
  friend class Backbone;

  public:

//--------------------------------------------------------------------------
//  CONSTRUCTION
//--------------------------------------------------------------------------

    const static usint ARRAY_BLOCK_SIZE = 32;
    const static usint OUTGOING_BLOCK_SIZE = 32;
    const static usint SAMPLE_BLOCK_SIZE = 32;
    const static usint SAMPLE_RATE = 16;

    GCSA(const std::string& base_name);
    GCSA(PathGraph& graph, Graph& parent, bool print = false);
    ~GCSA();

    void writeTo(const std::string& base_name) const;

    bool isOk() const;

//--------------------------------------------------------------------------
//  REPORTING
//--------------------------------------------------------------------------

    // Returns the size of the data structure.
    usint reportSize(bool print = false) const;

    void printCounts() const;

    inline usint getSize() const { return this->node_count; }
    inline usint getNumberOfAutomata() const { return this->alphabet->countOf(0); }

//--------------------------------------------------------------------------
//  SEARCHING
//--------------------------------------------------------------------------

    // Returns the SA range for suffixes starting with the pattern.
    pair_type find(const std::string& pattern) const;

    // Returns the set of unique node values for the given SA range of nodes.
    // User must delete the returned vector.
    std::vector<usint>* locate(pair_type range) const;
    std::vector<usint>* locate(std::vector<pair_type>& ranges) const;

    bool isSampled(usint index) const;
    usint locate(usint index) const;

    // Returns the label of the given node.
    usint labelOf(usint index) const;

//--------------------------------------------------------------------------
//  LOW-LEVEL INTERFACE
//--------------------------------------------------------------------------

    pair_type getSARange() const;
    pair_type getBWTRange() const;
    pair_type getCharRange(usint c) const;

    void convertToBWTRange(pair_type& sa_range) const;
    void convertToSARange(pair_type& bwt_range) const;
    void convertToSARange(std::vector<pair_type>& bwt_ranges) const;

    // This is an unsafe function that does not check its parameters.
    // Note that LF will not work correctly, if there are multiple automata and
    // c is the label of initial nodes.
    pair_type LF(pair_type bwt_range, usint c) const;

    // User must free the returned vector.
    std::vector<usint>* locateRange(pair_type range) const;
    std::vector<usint>* locateRanges(std::vector<pair_type>& ranges) const;

    // The following are specific to GCSA.

    // This will not work correctly for an initial node, if there are multiple automata.
    std::vector<usint>* getSuccessors(usint index) const;

    DeltaVector::Iterator* getIterator(usint c) const;
    RLEVector::Iterator* getEdgeIterator() const;

//--------------------------------------------------------------------------
//  INTERNAL VARIABLES
//--------------------------------------------------------------------------

  private:
    bool ok, support_locate;

    usint node_count;

    DeltaVector* array[CHARS];  // BWT
    RLEVector* outgoing;        // M

    DeltaVector* sampled_positions;
    ReadBuffer* samples;

    Alphabet* alphabet;

    Backbone* backbone;                   // Only used during construction. FIXME: for now?

//--------------------------------------------------------------------------
//  INTERNAL STUFF
//--------------------------------------------------------------------------

    pair_type convertToNodeRange(pair_type edge_range) const;
    usint Psi(usint index) const; // Follows the first outgoing edge.
    usint LF(usint index, usint c) const;

    void locate(pair_type range, std::vector<usint>& vec) const;
    usint locateUnsafe(usint index) const;

    // These are not allowed.
    GCSA();
    GCSA(const GCSA&);
    GCSA& operator = (const GCSA&);
};

//--------------------------------------------------------------------------

class Backbone
{
  public:
    Backbone(const GCSA& _gcsa, PathGraph& graph, Graph& parent, bool print = false);
    Backbone(const std::string& base_name, const GCSA& _gcsa);
    ~Backbone();

    void writeTo(const std::string& base_name) const;
    usint getSize() const;
    bool isOk() const;

    // Returns the size of the data structure.
    usint reportSize(bool print = false) const;

    bool contains(usint index) const;

    bool originalContains(usint index) const;

    // The behavior of these functions is undefined if index is not in the backbone.
    usint next(usint index) const;
    usint previous(usint index) const;

    const static usint NODE_BLOCK_SIZE = 64;
    const static usint EDGE_BLOCK_SIZE = 16;

  private:
    const GCSA& gcsa;

    /*
      nodes[i] tells whether node i is in the backbone.
      edges.select(i) is the outgoing edge from node i to the next node in the backbone.
      If node i is not in the backbone, edges.select(i) is the first outgoing edge from the node.
    */
    SuccinctVector* nodes;
    RLEVector*   edges;

    SuccinctVector* original;

    usint size; // Number of backbone nodes.
    bool ok;

    // These are not allowed.
    Backbone();
    Backbone(const Backbone&);
    Backbone& operator = (const Backbone&);
};


} // namespace CSA


#endif // GCSA_H
