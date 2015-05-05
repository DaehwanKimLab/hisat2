/*
 * Copyright 2015, Daehwan Kim <infphilo@gmail.com>
 *
 * This file is part of HISAT 2.
 *
 * HISAT 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef GBWT_GRAPH_H_
#define GBWT_GRAPH_H_

#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

// Reference:
// Jouni Sirén, Niko Välimäki, and Veli Mäkinen: Indexing Graphs for Path Queries with Applications in Genome Research.
// IEEE/ACM Transactions on Computational Biology and Bioinformatics 11(2):375-388, 2014.
// http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6698337

//--------------------------------------------------------------------------

enum SNP_TYPE {
    SNP_SNS = 0, // single nucleotide substitution
    SNP_INS,     // small insertion wrt reference genome
    SNP_DEL,     // small deletion wrt reference genome
};

template <typename index_t>
struct SNP {
    index_t        pos;
    SNP_TYPE       type;
    EList<char, 4> diff;
};

template <typename index_t>
class RefGraph
{
public:
    struct Node {
        char    label; // ACGTN + #(=Z) + $(=0)
        index_t value; // location in a whole genome
        
        Node() {}
        Node(char label_, index_t value_) : label(label_), value(value_) {}
    };
    
    struct Edge {
        index_t from; // from Node
        index_t to;   // to Node
        
        Edge() {}
        Edge(index_t from_, index_t to_) : from(from_), to(to_) {}
    };
    
public:
    RefGraph(const string& ref_fname, const string& snp_fname, bool verbose);
    
#if 0
    RefGraph(std::vector<GraphNode>& node_vector, std::vector<GraphEdge>& edge_vector);
    RefGraph(const std::string& base_name);
    
    RefGraph();
    
    void write(const std::string& base_name);
    
    void printInfo();
    
    void sortEdges(bool from);
    pair_type getEdges(char node, bool from);  // Use sortEdges() first.
    pair_type getNextEdgeRange(pair_type range, bool from); // Use sortEdges() first.
    
    // Marks predecessor labels in GraphEdges.
    void markPredecessors();
    
    // Make the graph reverse deterministic. Updates backbone, if present.
    // The value of a merged node is the maximum value of the original below threshold,
    // or the maximum of all original nodes if the values are all >= threshold.
    // Threshold is the number of backbone nodes, if backbone is present, or the value of the
    // current node otherwise.
    void determinize();
    bool isDeterministic();
    
    char maxLabel();
    
    // Updates the labels when there are multiple automata to be indexed.
    // Final nodes get labels 0..automata-1.
    // The labels of regular nodes are incremented by automata.
    // Initial nodes get labels max_label+automata..max_label+2automata-1.
    void changeLabels(index_t automata, char max_label);
    void restoreLabels(index_t automata, char max_label);
    
    // Adds the nodes and edges from the given graph to current graph.
    // If both graphs contain a backbone, the value fields of the nodes are updated accordingly.
    // The backbone of this graph will be deleted, as it is no longer valid.
    // To avoid losing backbone information, call encodeBackbone() on both graphs before addGraph().
    void addGraph(const Graph<index_t>& another);
#endif
    
private:
    EList<Node> nodes;
    EList<Edge> edges;
    bool       ok;  // Also tells that edges are sorted by from.
};

/**
 * Load reference sequence file and snp information.
 * Construct a reference graph
 */
template <typename index_t>
RefGraph<index_t>::RefGraph(const string& ref_fname, const string& snp_fname, bool verbose) {
    EList<FileBuf*> is(MISC_CAT);
    int reverse = 0;
    bool nsToAs = false, bisulfite = false;
    RefReadInParams refparams(false, reverse, nsToAs, bisulfite);
    FILE *f = fopen(ref_fname.c_str(), "r");
    if (f == NULL) {
        cerr << "Error: could not open " << ref_fname.c_str() << endl;
        throw 1;
    }
    FileBuf *fb = new FileBuf(f);
    assert(fb != NULL);
    if(fb->peek() == -1 || fb->eof()) {
        cerr << "Warning: Empty fasta file: '" << ref_fname.c_str() << "'" << endl;
        throw 1;
    }
    assert(!fb->eof());
    assert(fb->get() == '>');
    ASSERT_ONLY(fb->reset());
    assert(!fb->eof());
    is.push_back(fb);
    if(is.empty()) {
        cerr << "Warning: All fasta inputs were empty" << endl;
        throw 1;
    }
    // Vector for the ordered list of "records" comprising the input
    // sequences.  A record represents a stretch of unambiguous
    // characters in one of the input sequences.
    EList<RefRecord> szs(MISC_CAT);
    bool bigEndian = false, sanityCheck = false;
    BitPairReference::szsFromFasta(is, string(), bigEndian, refparams, szs, sanityCheck);
    assert_gt(szs.size(), 0);
    
    EList<string> refnames;
    
    SString<char> s;
    index_t jlen = 0;
    for(size_t i = 0; i < szs.size(); i++) {
        jlen += (index_t)szs[i].len;
    }
    try {
        Timer _t(cerr, "  (1/5) Time reading reference sequence: ", verbose);
        
        s.resize(jlen);
        RefReadInParams rpcp = refparams;
        // For each filebuf
        assert_eq(is.size(), 1);
        FileBuf *fb = is[0];
        assert(!fb->eof());
        // For each *fragment* (not necessary an entire sequence) we
        // can pull out of istream l[i]...
        if(!fb->eof()) {
            // Push a new name onto our vector
            refnames.push_back("");
            TIndexOffU distoff = 0;
            fastaRefReadAppend(*fb, true, s, distoff, rpcp, &refnames.back());
        }
        fb->reset();
        assert(!fb->eof());
        // Joined reference sequence now in 's'
    } catch(bad_alloc& e) {
        // If we throw an allocation exception in the try block,
        // that means that the joined version of the reference
        // string itself is too larger to fit in memory.  The only
        // alternatives are to tell the user to give us more memory
        // or to try again with a packed representation of the
        // reference (if we haven't tried that already).
        cerr << "Could not allocate space for a joined string of " << jlen << " elements." << endl;
        // There's no point passing this exception on.  The fact
        // that we couldn't allocate the joined string means that
        // --bmax is irrelevant - the user should re-run with
        // ebwt-build-packed
        cerr << "Please try running hisat2-build on a computer with more memory." << endl;
        if(sizeof(void*) == 4) {
            cerr << "If this computer has more than 4 GB of memory, try using a 64-bit executable;" << endl
            << "this executable is 32-bit." << endl;
        }
        throw 1;
    }
    
    ifstream snp_file(snp_fname.c_str(), ios::in);
    if(!snp_file.is_open()) {
        cerr << "Error: could not open "<< snp_fname.c_str() << endl;
        throw 1;
    }
    
    EList<SNP<index_t> > snps;
    while(!snp_file.eof()) {
        string snp_id;
        snp_file >> snp_id;
        if(snp_id.empty() || snp_id[0] == '#') {
            string line;
            getline(snp_file, line);
            continue;
        }
        string chr, type, diff;
        index_t pos;
        snp_file >> chr >> pos >> type >> diff;
        
        snps.expand();
        SNP<index_t>& snp = snps.back();
        snp.pos = pos;
        snp.diff.clear();
        if(type == "SNS") {
            snp.type = SNP_SNS;
            if(diff.size() != 1) {
                cerr << "Error: SNS type takes only one base" << endl;
                throw 1;
            }
            snp.diff.push_back(diff[0]);
        } else if(type == "INS") {
            snp.type = SNP_INS;
            for(size_t i = 0; i < diff.size(); i++) {
                char base = toupper(diff[i]);
                if(base == 'A' || base == 'C' || base == 'G' || base == 'T') {
                    snp.diff.push_back(base);
                } else {
                    snps.pop_back();
                    continue;
                }
            }
        } else if(type == "DEL") {
            snp.type = SNP_DEL;
            for(size_t i = 0; i < diff.size(); i++) {
                assert_eq(diff[i], '-');
                snp.diff.push_back('-');
            }
        } else {
            cerr << "Error: unknown snp type " << type << endl;
            throw 1;
        }
    }
    snp_file.close();
    
    index_t num_predicted_nodes = jlen * 12 / 10;
    nodes.reserveExact(num_predicted_nodes);
    edges.reserveExact(num_predicted_nodes);
    
    // Created head node
    nodes.expand();
    nodes.back().label = 'Z';
    nodes.back().value = 0;
    // Create nodes and edges corresponding to reference genome
    for(size_t i = 0; i < s.length(); i++) {
        nodes.expand();
        nodes.back().label = "ACGT"[(int)s[i]];
        nodes.back().value = i;
        
        assert_geq(nodes.size(), 2);
        edges.expand();
        edges.back().from = nodes.size() - 2;
        edges.back().to = nodes.size() - 1;
    }
    
    // Create tail node
    nodes.expand();
    nodes.back().label = 0;
    nodes.back().value = s.length();
    edges.expand();
    edges.back().from = nodes.size() - 2;
    edges.back().to = nodes.size() - 1;
    
    // Create nodes and edges for SNPs
    for(size_t i = 0; i < snps.size(); i++) {
        const SNP<index_t>& snp = snps[i];
        if(snp.type == SNP_SNS) {
            assert_eq(snp.diff.size(), 1);
            nodes.expand();
            nodes.back().label = snp.diff[0];
            nodes.back().value = snp.pos;
            
            edges.expand();
            edges.back().from = snp.pos;
            edges.back().to = nodes.size() - 1;
            
            edges.expand();
            edges.back().from = nodes.size() - 1;
            edges.back().to = snp.pos + 2;
        } else if(snp.type == SNP_INS) {
            
        } else if(snp.type == SNP_DEL) {
            
        } else {
            assert(false);
        }
    }
    
    if(verbose) {
        cerr << "Nodes:" << endl;
        for(size_t i = 0; i < nodes.size(); i++) {
            const Node& node = nodes[i];
            cerr << "\t" << i << "\t" << node.label << "\t" << node.value << endl;
        }
        cerr << endl;
        cerr << "Edges: " << endl;
        for(size_t i = 0; i < edges.size(); i++) {
            const Edge& edge = edges[i];
            cerr << "\t" << i << "\t" << edge.from << " --> " << edge.to << endl;
        }
        cerr << endl;
    }
}

#if 0

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

#endif

#endif /*GBWT_GRAPH_H_*/
