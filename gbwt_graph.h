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
#include <map>
#include <deque>

using namespace std;

// Reference:
// Jouni Sirén, Niko Välimäki, and Veli Mäkinen: Indexing Graphs for Path Queries with Applications in Genome Research.
// IEEE/ACM Transactions on Computational Biology and Bioinformatics 11(2):375-388, 2014.
// http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6698337

//--------------------------------------------------------------------------

enum SNP_TYPE {
    SNP_SGL = 0, // single nucleotide substitution
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
class RefGraph {
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
    
    struct EdgeFromCmp {
        bool operator() (const Edge& a, const Edge& b) const {
            return a.from < b.from;
        }
    };

    struct EdgeToCmp {
        bool operator() (const Edge& a, const Edge& b) const {
            return a.to < b.to;
        };
    };
    
public:
    RefGraph(const string& ref_fname, const string& snp_fname, bool verbose);
    
#if 0
    void write(const std::string& base_name) {}
    void printInfo() {}
#endif
    
private:
    bool isReverseDeterministic();
    void reverseDeterminize();
    
    void sortEdgesFrom() {
        std::sort(edges.begin(), edges.end(), EdgeFromCmp());
    }
    void sortEdgesTo() {
        std::sort(edges.begin(), edges.end(), EdgeToCmp());
    }
    
    // Return edge ranges [begin, end)
    pair<index_t, index_t> findEdges(index_t node, bool from);
    pair<index_t, index_t> findEdgesFrom(index_t node) {
        return findEdges(node, true);
    }
    pair<index_t, index_t> findEdgesTo(index_t node) {
        return findEdges(node, false);
    }
    
#if 0
    pair_type getEdges(char node, bool from);  // Use sortEdges() first.
    pair_type getNextEdgeRange(pair_type range, bool from); // Use sortEdges() first.
    
    // Marks predecessor labels in GraphEdges.
    void markPredecessors();
    
    char maxLabel();
    
    // Updates the labels when there are multiple automata to be indexed.
    // Final nodes get labels 0..automata-1.
    // The labels of regular nodes are incremented by automata.
    // Initial nodes get labels max_label+automata..max_label+2automata-1.
    void changeLabels(index_t automata, char max_label);
    void restoreLabels(index_t automata, char max_label);
#endif
    
private:
    EList<Node> nodes;
    EList<Edge> edges;
    index_t     zNode;
    
private:
    // Following composite nodes and edges are used to reverse-determinize an automaton.
    struct CompositeNode {
        std::basic_string<index_t> nodes;
        index_t                    id;
        char                       label;
        index_t                    value;
        
        CompositeNode() : id(0), label(0), value(0) {}
        
        CompositeNode(char label_, index_t node_id) :
        id(0), label(label_)
        {
            nodes.push_back(node_id);
        }
        
        Node getNode() const {
            return Node(label, value);
        }
    };
    
    struct CompositeEdge {
        index_t from;
        index_t to;
        
        CompositeEdge() : from(0), to(0) {}
        CompositeEdge(index_t from_, index_t to_) : from(from_), to(to_) {}
        
        Edge getEdge(const EList<CompositeNode>& nodes) const
        {
            assert_lt(from, nodes.size());
            const CompositeNode& from_node = nodes[from];
            assert_lt(to, nodes.size());
            const CompositeNode& to_node = nodes[to];
            return Edge(from_node.id, to_node.id);
        }
        
        bool operator < (const CompositeEdge& o) const
        {
            if(from < o.from ||
               (from == o.from && to < o.to)) {
                return true;
            } else {
                return false;
            }
        }
    };
    
    struct TempNodeLabelCmp {
        TempNodeLabelCmp(const EList<Node>& nodes_) : nodes(nodes_) {}
        bool operator() (index_t a, index_t b) const {
            assert_lt(a, nodes.size());
            assert_lt(b, nodes.size());
            return nodes[a].label < nodes[b].label;
        }
        
        const EList<Node>& nodes;
    };
};

/**
 * Load reference sequence file and snp information.
 * Construct a reference graph
 */
template <typename index_t>
RefGraph<index_t>::RefGraph(const string& ref_fname, const string& snp_fname, bool verbose)
: zNode(0)
{
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
        // rs73387790	single	22:20000001-21000000	145	A
        string snp_id;
        snp_file >> snp_id;
        if(snp_id.empty() || snp_id[0] == '#') {
            string line;
            getline(snp_file, line);
            continue;
        }
        string type, chr, diff;
        index_t pos;
        snp_file >> type >> chr >> pos >> diff;
        
        snps.expand();
        SNP<index_t>& snp = snps.back();
        snp.pos = pos;
        snp.diff.clear();
        if(type == "single") {
            snp.type = SNP_SGL;
            if(diff.size() != 1) {
                cerr << "Error: SNS type takes only one base" << endl;
                throw 1;
            }
            snp.diff.push_back(diff[0]);
        } else if(type == "deletion") {
            snp.type = SNP_DEL;
            for(size_t i = 0; i < diff.size(); i++) {
                assert_eq(diff[i], '-');
                snp.diff.push_back('-');
            }
        } else if(type == "insertion") {
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
        } else {
            cerr << "Error: unknown snp type " << type << endl;
            throw 1;
        }
    }
    snp_file.close();
    
    index_t num_predicted_nodes = (index_t)(jlen * 1.2);
    nodes.reserveExact(num_predicted_nodes);
    edges.reserveExact(num_predicted_nodes);
    
    // Created head node
    nodes.expand();
    nodes.back().label = 'Z';
    nodes.back().value = 0;
    // Create nodes and edges corresponding to a reference genome
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
    nodes.back().label = '$';
    nodes.back().value = s.length();
    zNode = nodes.size() - 1;
    edges.expand();
    edges.back().from = nodes.size() - 2;
    edges.back().to = nodes.size() - 1;
    
    // Create nodes and edges for SNPs
    for(size_t i = 0; i < snps.size(); i++) {
        const SNP<index_t>& snp = snps[i];
        if(snp.type == SNP_SGL) {
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
        } else if(snp.type == SNP_DEL) {
            index_t deletionLen = snp.diff.size();
            assert_gt(deletionLen, 0);
            // assert_gt(snp.pos + deletionLen, *something*);
            edges.expand();
            edges.back().from = snp.pos;
            edges.back().to = snp.pos + deletionLen + 1;
        } else if(snp.type == SNP_INS) {
            assert_gt(snp.diff.size(), 0);
            for(size_t j = 0; j < snp.diff.size(); j++) {
                char inserted_bp = snp.diff[j];
                
                nodes.expand();
                nodes.back().label = inserted_bp;
                nodes.back().value = snp.pos;
                
                edges.expand();
                edges.back().from = (j == 0 ? snp.pos : nodes.size() - 2);
                edges.back().to = nodes.size() - 1;
            }
            edges.expand();
            edges.back().from = nodes.size() - 1;
            edges.back().to = snp.pos + 1;
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
    
    if(!isReverseDeterministic()) {
        cerr << "is not reverse-deterministic" << endl;
        reverseDeterminize();
    }
}

template <typename index_t>
pair<index_t, index_t> RefGraph<index_t>::findEdges(index_t node, bool from) {
    pair<index_t, index_t> range(0, 0);
    assert_gt(edges.size(), 0);
    
    // Find lower bound
    index_t low = 0, high = edges.size() - 1;
    index_t temp;
    while(low < high) {
        index_t mid = low + (high - low) / 2;
        temp = (from ? edges[mid].from : edges[mid].to);
        if(node == temp) {
            high = mid;
        } else if(node < temp) {
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }
    temp = (from ? edges[low].from : edges[low].to);
    if(node == temp) {
        range.first = low;
    } else {
        return range;
    }
    
    // Find upper bound
    high = edges.size() - 1;
    while(low < high)
    {
        index_t mid = low + (high - low + 1) / 2;
        temp = (from ? edges[mid].from : edges[mid].to);
        if(node == temp) {
            low = mid;
        } else {
            assert_lt(node, temp);
            high = mid - 1;
        }
    }
#ifndef NDEBUG
    temp = (from ? edges[high].from : edges[high].to);
    assert_eq(node, temp);
#endif
    range.second = high + 1;
    return range;
}

template <typename index_t>
bool RefGraph<index_t>::isReverseDeterministic()
{
    if(edges.size() <= 0) return true;
    
    // Sort edges by "to" nodes
    sortEdgesTo();
    
    index_t curr_to = edges.front().to;
    EList<char> seen;
    for(index_t i = 1; i < edges.size(); i++) {
        index_t from = edges[i].from;
        assert_lt(from, nodes.size());
        char nt = nodes[from].label;
        if(curr_to != edges[i].to) {
            curr_to = edges[i].to;
            seen.clear();
            seen.push_back(nt);
        } else {
            for(index_t j = 0; j < seen.size(); j++) {
                if(nt == seen[j]) {
                    return false;
                }
            }
            seen.push_back(nt);
        }
    }
    
    return true;
}


template <typename index_t>
void RefGraph<index_t>::reverseDeterminize()
{
    EList<CompositeNode> cnodes; cnodes.resize(nodes.size());
    map<basic_string<index_t>, index_t> cnode_map;
    deque<index_t> active_cnodes;
    EList<CompositeEdge> cedges; cedges.resize(edges.size());
    
    // Start from the final node ('$')
    assert_lt(zNode, nodes.size());
    const Node& z_node = nodes[zNode];
    cnodes.expand();
    cnodes.back().label = z_node.label;
    cnodes.back().nodes.push_back(zNode);
    active_cnodes.push_back(0);
    cnode_map[cnodes.back().nodes] = 0;
    
    sortEdgesTo();
    
    EList<index_t> predecessors;
    while(!active_cnodes.empty()) {
        index_t curr_cnode_id = active_cnodes.front(); active_cnodes.pop_front();
        assert_lt(curr_cnode_id, cnodes.size());
        CompositeNode& curr_cnode = cnodes[curr_cnode_id];
        
        // Find predecessors of this composite node
        predecessors.clear();
        for(size_t i = 0; i < curr_cnode.nodes.size(); i++) {
            index_t curr_node_id = curr_cnode.nodes[i];
            pair<index_t, index_t> edge_range = findEdgesTo(curr_node_id);
            for(index_t j = edge_range.first; j < edge_range.second; j++) {
                assert_lt(j, edges.size());
                predecessors.push_back(edges[j].from);
            }
        }
        
        predecessors.sort();
        index_t new_size = unique(predecessors.begin(), predecessors.end()) - predecessors.begin();
        predecessors.resize(new_size);
        
        // Create composite nodes by labels
        stable_sort(predecessors.begin(), predecessors.end(), TempNodeLabelCmp(nodes));
        for(size_t i = 0; i < predecessors.size(); i++) {
            index_t node_id = predecessors[i];
            assert_lt(node_id, nodes.size());
            const Node& node = nodes[node_id]; i++;
            cnodes.expand();
            cnodes.back().label = node.label;
            cnodes.back().nodes.push_back(node_id);
            
            while(i < predecessors.size()) {
                index_t next_node_id = predecessors[i];
                assert_lt(next_node_id, nodes.size());
                const Node& next_node = nodes[next_node_id];
                if(next_node.label != node.label) break;
                cnodes.back().nodes.push_back(next_node_id);
                i++;
            }
            
            // Create edges from this new composite node to current composite node
            typename map<basic_string<index_t>, index_t>::iterator existing = cnode_map.find(cnodes.back().nodes);
            if(existing == cnode_map.end()) {
                cnode_map[cnodes.back().nodes] = cnodes.size() - 1;
                active_cnodes.push_back(cnodes.size() - 1);
                cedges.push_back(CompositeEdge(cnodes.size() - 1, curr_cnode_id));
            } else {
                cnodes.pop_back();
                cedges.push_back(CompositeEdge((*existing).second, curr_cnode_id));
            }
            
            // Increment indegree
            curr_cnode.id++;
        }
    }
    
    // Create new nodes
    nodes.resizeExact(cnodes.size()); nodes.clear();
    sort(cedges.begin(), cedges.end());
    for(index_t i = 0; i < cnodes.size(); i++) {
        CompositeNode& node = cnodes[i];
        if(node.id == 0) {
            active_cnodes.push_back(i);
            node.id = nodes.size();
            nodes.expand();
            nodes.back() = node.getNode();
        }
    }
    
    while(!active_cnodes.empty()) {
        index_t curr_cnode_id = active_cnodes.front(); active_cnodes.pop_front();
        assert_lt(curr_cnode_id, cnodes.size());
        index_t i = cedges.bsearchLoBound(CompositeEdge(curr_cnode_id, 0));
        while(true) {
            assert_geq(cedges[i].from, curr_cnode_id);
            if(cedges[i].from != curr_cnode_id) break;
            index_t successor_cnode_id = cedges[i].to;
            assert_lt(successor_cnode_id, cnodes.size());
            CompositeNode& successor_cnode = cnodes[successor_cnode_id];
            assert_gt(successor_cnode.id, 0);
            successor_cnode.id--;
            if(successor_cnode.id == 0) {
                active_cnodes.push_back(successor_cnode_id);
                successor_cnode.id = nodes.size();
                nodes.expand();
                nodes.back() = successor_cnode.getNode();
            }
            i++;
        }
    }
    
    // Create new edges
    edges.resizeExact(cedges.size()); edges.clear();
    for(index_t i = 0; i < cedges.size(); i++) {
        const CompositeEdge& edge = cedges[i];
        edges.expand();
        edges.back() = edge.getEdge(cnodes);
    }
    
    sortEdgesFrom();
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
