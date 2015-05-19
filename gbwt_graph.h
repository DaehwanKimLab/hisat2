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

template <typename index_t> class PathGraph;

// Note: I wrote the following codes based on Siren's work, gcsa (see the reference above).
template <typename index_t>
class RefGraph {
    friend class PathGraph<index_t>;
public:
    struct Node {
        char    label; // ACGTN + #(=Z) + $(=0)
        index_t value; // location in a whole genome
        bool    backbone; // backbone node, which corresponds to a reference sequence
        
        Node() { reset(); }
        Node(char label_, index_t value_, bool backbone_ = false) : label(label_), value(value_), backbone(backbone_) {}
        void reset() { label = 0; value = 0; backbone = false; }
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
    
    bool repOk() { return true; }
    
    void write(const std::string& base_name) {}
    void printInfo() {}
    
private:
    bool isReverseDeterministic();
    void reverseDeterminize(bool debug = false);
    
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
    index_t     lastNode; // $
    
private:
    // Following composite nodes and edges are used to reverse-determinize an automaton.
    struct CompositeNode {
        std::basic_string<index_t> nodes;
        index_t                    id;
        char                       label;
        index_t                    value;
        bool                       backbone;
        bool                       potential;
        
        CompositeNode() { reset(); }
        
        CompositeNode(char label_, index_t node_id) :
        id(0), label(label_)
        {
            nodes.push_back(node_id);
        }
        
        Node getNode() const {
            return Node(label, value, backbone);
        }
        
        void reset() {
            nodes.clear();
            id = 0;
            label = 0;
            value = 0;
            backbone = potential = false;
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
            return from < o.from;
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
: lastNode(0)
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
    
    EList<pair<index_t, index_t> > chr_szs;
    index_t jlen = 0;
    for(index_t i = 0; i < szs.size(); i++) {
        if(szs[i].first) {
            chr_szs.expand();
            chr_szs.back().first = jlen;
            chr_szs.back().second = i;
        }
        jlen += (index_t)szs[i].len;
    }
    bool debug = (jlen <= 10);
    SString<char> s;
    EList<string> refnames;
    try {
        Timer _t(cerr, "  (1/5) Time reading reference sequence: ", verbose);
        
        s.resize(jlen);
        RefReadInParams rpcp = refparams;
        // For each filebuf
        assert_eq(is.size(), 1);
        FileBuf *fb = is[0];
        assert(!fb->eof());
        index_t szsi = 0;
        index_t distoff = 0;
        bool first = true;
        while(!fb->eof()) {
            if(szs[szsi].first) refnames.push_back("");
            RefRecord rec = fastaRefReadAppend(*fb, first, s, distoff, rpcp, &refnames.back());
            first = false;
            assert_eq(rec.off, szs[szsi].off);
            assert_eq(rec.len, szs[szsi].len);
            assert_eq(rec.first, szs[szsi].first);
            szsi++;
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
        
        index_t chr_idx = 0;
        for(; chr_idx < refnames.size(); chr_idx++) {
            if(chr == refnames[chr_idx])
                break;
        }
        if(chr_idx >= refnames.size()) continue;
        
        assert_eq(chr_szs.size(), refnames.size());
        assert_lt(chr_idx, chr_szs.size());
        pair<index_t, index_t> tmp_pair = chr_szs[chr_idx];
        const index_t sofar_len = tmp_pair.first;
        const index_t szs_idx = tmp_pair.second;
        bool inside_Ns = false;
        index_t add_pos = 0;
        assert(szs[szs_idx].first);
        for(index_t i = szs_idx; i < szs.size(); i++) {
            if(i != szs_idx && szs[i].first) break;
            if(pos < szs[i].off) {
                inside_Ns = true;
                break;
            } else {
                pos -= szs[i].off;
                if(pos < szs[i].len) {
                    break;
                } else {
                    pos -= szs[i].len;
                    add_pos += szs[i].len;
                }
            }
        }
        
        if(inside_Ns) continue;
        pos = sofar_len + add_pos + pos;
        if(chr_idx + 1 < chr_szs.size()) {
            if(pos >= chr_szs[chr_idx + 1].first) continue;
        } else {
            if(pos >= jlen) continue;
        }

        snps.expand();
        SNP<index_t>& snp = snps.back();
        snp.pos = pos;
        snp.diff.clear();
        if(type == "single") {
            snp.type = SNP_SGL;
            if(diff.size() != 1) {
                cerr << "Error: single type takes only one base" << endl;
                throw 1;
            }
            char bp = diff[0];
            if(bp == "ACGTN"[(int)s[pos]]) {
                cerr << "Error: single type (" << bp
                     << ") should have a different base than " << "ACGTN"[(int)s[pos]]
                     << " (" << snp_id << ")" << endl;
                throw 1;
            }
            snp.diff.push_back(bp);
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
    lastNode = nodes.size() - 1;
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
            nodes.back().value = nodes.size();
            
            edges.expand();
            edges.back().from = snp.pos;
            edges.back().to = nodes.size() - 1;
            
            edges.expand();
            edges.back().from = nodes.size() - 1;
            edges.back().to = snp.pos + 2;
        } else if(snp.type == SNP_DEL) {
            index_t deletionLen = snp.diff.size();
            assert_gt(deletionLen, 0);
            edges.expand();
            edges.back().from = snp.pos;
            edges.back().to = snp.pos + deletionLen + 1;
        } else if(snp.type == SNP_INS) {
            assert_gt(snp.diff.size(), 0);
            for(size_t j = 0; j < snp.diff.size(); j++) {
                char inserted_bp = snp.diff[j];
                nodes.expand();
                nodes.back().label = inserted_bp;
                nodes.back().value = nodes.size();
                
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
    
    if(debug) {
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
        reverseDeterminize(debug);
        assert(isReverseDeterministic());
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
    EList<bool> seen; seen.resize(4); seen.fillZero();
    for(index_t i = 1; i < edges.size(); i++) {
        index_t from = edges[i].from;
        assert_lt(from, nodes.size());
        char nt = nodes[from].label;
        nt = asc2dna[(int)nt];
        assert_lt(nt, seen.size());
        if(curr_to != edges[i].to) {
            curr_to = edges[i].to;
            seen.fillZero();
            seen[nt] = true;
        } else {
            if(seen[nt]) return false;
            seen[nt] = true;
        }
    }
    
    return true;
}


template <typename index_t>
void RefGraph<index_t>::reverseDeterminize(bool debug)
{
    EList<CompositeNode> cnodes; cnodes.ensure(nodes.size());
    map<basic_string<index_t>, index_t> cnode_map;
    deque<index_t> active_cnodes;
    EList<CompositeEdge> cedges; cedges.ensure(edges.size());
    
    // Start from the final node ('$')
    assert_lt(lastNode, nodes.size());
    const Node& last_node = nodes[lastNode];
    cnodes.expand();
    cnodes.back().reset();
    cnodes.back().label = last_node.label;
    cnodes.back().value = last_node.value;
    cnodes.back().backbone = true;
    cnodes.back().nodes.push_back(lastNode);
    active_cnodes.push_back(0);
    cnode_map[cnodes.back().nodes] = 0;
    
    sortEdgesTo();
  
    index_t firstNode = 0; // Z -> ... -> $
    EList<index_t> predecessors;
    while(!active_cnodes.empty()) {
        index_t cnode_id = active_cnodes.front(); active_cnodes.pop_front();
        assert_lt(cnode_id, cnodes.size());
        CompositeNode& cnode = cnodes[cnode_id];
        
        // Find predecessors of this composite node
        predecessors.clear();
        for(size_t i = 0; i < cnode.nodes.size(); i++) {
            index_t node_id = cnode.nodes[i];
            pair<index_t, index_t> edge_range = findEdgesTo(node_id);
            assert_leq(edge_range.first, edge_range.second);
            assert_leq(edge_range.second, edges.size());
            for(index_t j = edge_range.first; j < edge_range.second; j++) {
                assert_eq(edges[j].to, node_id);
                predecessors.push_back(edges[j].from);
            }
        }
        
        if(predecessors.size() >= 2) {
            // Remove redundant nodes
            predecessors.sort();
            index_t new_size = unique(predecessors.begin(), predecessors.end()) - predecessors.begin();
            predecessors.resize(new_size);
            
            // Create composite nodes by labels
            stable_sort(predecessors.begin(), predecessors.end(), TempNodeLabelCmp(nodes));
        }
        
        for(size_t i = 0; i < predecessors.size();) {
            index_t node_id = predecessors[i];
            assert_lt(node_id, nodes.size());
            const Node& node = nodes[node_id]; i++;
            cnodes.expand();
            cnodes.back().reset();
            cnodes.back().label = node.label;
            cnodes.back().value = node.value;
            cnodes.back().nodes.push_back(node_id);
            cnodes.back().potential = node.value < lastNode;
            
            if(node.label == 'Z' && firstNode == 0) {
                firstNode = cnodes.size() - 1;
                cnodes.back().backbone = true;
            }
            
            while(i < predecessors.size()) {
                index_t next_node_id = predecessors[i];
                assert_lt(next_node_id, nodes.size());
                const Node& next_node = nodes[next_node_id];
                if(next_node.label != node.label) break;
                cnodes.back().nodes.push_back(next_node_id);

                if((cnode.value < lastNode) != (next_node.value < lastNode)) {
                    cnode.value = min(cnode.value, next_node.value);
                } else {
                    cnode.value = max(cnode.value, next_node.value);
                }
                cnode.potential = cnode.value < lastNode;
                i++;
            }
            
            // Create edges from this new composite node to current composite node
            typename map<basic_string<index_t>, index_t>::iterator existing = cnode_map.find(cnodes.back().nodes);
            if(existing == cnode_map.end()) {
                cnode_map[cnodes.back().nodes] = cnodes.size() - 1;
                active_cnodes.push_back(cnodes.size() - 1);
                cedges.push_back(CompositeEdge(cnodes.size() - 1, cnode_id));
            } else {
                cnodes.pop_back();
                cedges.push_back(CompositeEdge((*existing).second, cnode_id));
            }
            
            // Increment indegree
            cnode.id++;
        }
    }
    
    // Specify backbone nodes
    // Interchange from and to
    for(index_t i = 0; i < cedges.size(); i++) {
        index_t tmp = cedges[i].from;
        cedges[i].from = cedges[i].to;
        cedges[i].to = tmp;
    }
    sort(cedges.begin(), cedges.end());
    active_cnodes.push_back(0);
    while(!active_cnodes.empty()) {
        index_t cnode_id = active_cnodes.front(); active_cnodes.pop_front();
        assert_lt(cnode_id, cnodes.size());
        const CompositeNode& cnode = cnodes[cnode_id];
        index_t i = cedges.bsearchLoBound(CompositeEdge(cnode_id, 0));
        while(i < cedges.size()) {
            assert_geq(cedges[i].from, cnode_id);
            if(cedges[i].from != cnode_id) break;
            index_t predecessor_cnode_id = cedges[i].to;
            assert_lt(predecessor_cnode_id, cnodes.size());
            CompositeNode& predecessor_cnode = cnodes[predecessor_cnode_id];
            if(predecessor_cnode.potential && cnode.value == predecessor_cnode.value + 1) {
                predecessor_cnode.backbone = true;
                predecessor_cnode.potential = false;
                active_cnodes.push_back(predecessor_cnode_id);
                break;
            }
            i++;
        }
    }
    // Restore from and to by interchanging them
    for(index_t i = 0; i < cedges.size(); i++) {
        index_t tmp = cedges[i].from;
        cedges[i].from = cedges[i].to;
        cedges[i].to = tmp;
    }
    
    // Create new nodes
    lastNode = 0; // Invalidate lastNode
    nodes.resizeExact(cnodes.size()); nodes.clear();
    assert_neq(firstNode, 0);
    assert_lt(firstNode, cnodes.size());
    CompositeNode& first_node = cnodes[firstNode];
    first_node.id = 0;
    nodes.expand();
    nodes.back() = first_node.getNode();
    active_cnodes.push_back(firstNode);    
    sort(cedges.begin(), cedges.end());
    while(!active_cnodes.empty()) {
        index_t cnode_id = active_cnodes.front(); active_cnodes.pop_front();
        assert_lt(cnode_id, cnodes.size());
        index_t i = cedges.bsearchLoBound(CompositeEdge(cnode_id, 0));
        while(i < cedges.size()) {
            assert_geq(cedges[i].from, cnode_id);
            if(cedges[i].from != cnode_id) break;
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
                if(nodes.back().label == '$') {
                    assert_eq(lastNode, 0);
                    assert_gt(nodes.size(), 1);
                    lastNode = nodes.size() - 1;
                }
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
    
    if(debug) {
        cerr << "Nodes:" << endl;
        for(size_t i = 0; i < nodes.size(); i++) {
            const Node& node = nodes[i];
            cerr << "\t" << i << "\t" << node.label << "\t" << node.value << (node.backbone ? "\tbackbone" : "") << endl;
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

static const uint8_t WORD_BITS = sizeof(TIndexOffU) * 8;

//--------------------------------------------------------------------------
template <typename index_t>
class PathGraph {
public:
    struct PathNode {
        index_t                from;
        index_t                to;
        pair<index_t, index_t> key;
        
        void setSorted()      { to = from; }
        bool isSorted() const { return to == from; }
        
        index_t value() const { return to; }
        index_t outdegree() const { return key.first; }
        
        bool  isBackbone() const {
            return (key.first & (((index_t)1) << (WORD_BITS - 1)));
        }
        void  setBackbone() {
            key.first |= ((index_t)1) << (WORD_BITS - 1);
        }
    };
    
    struct PathEdge {
        index_t from;
        index_t rank;
        char    label;
        
        PathEdge(index_t from_, index_t rank_, char label_) : from(from_), rank(rank_), label(label_) {}
    };
    
public:
    // Create a new graph in which paths are represented using nodes
    PathGraph(RefGraph<index_t>& parent);
    
    // Construct a next 2^(j+1) path graph from a 2^j path graph
    PathGraph(PathGraph<index_t>& previous);
    
    ~PathGraph() {}
    
    void printInfo();
    
    // status must be sorted. Calling this invalidates parent and
    // sets status to ready.
    // Writes outdegree to PathNode.key.second, value to PathNode.to, and
    // predecessor labels to PathNode.key.first.
    // Restores the labels of parent.
    bool generateEdges(RefGraph<index_t>& parent);
    
    // status must be ready.
    EList<pair<index_t, index_t> >* getSamples(index_t sample_rate, index_t& max_sample, const RefGraph<index_t>& base);
    
    EList<PathNode> nodes;
    EList<PathEdge> edges;
    index_t         ranks;
    index_t         max_label;  // Node label of initial nodes.
    index_t         temp_nodes; // Total number of nodes created before sorting.
    
    uint8_t         generation; // Sorted by paths of length 2^generation.
    
    enum status_t { error, ok, sorted, ready, edges_sorted } status;
    bool            has_stabilized;  // The number of nodes will probably not explode in subsequent doublings.
    
    // Can create an index by using key.second in PathNodes.
    // If the graph is not ready, its status becomes error.
    // Sorting edges by from actually sorts them by (from, to).
    void      sortEdges(bool by_from, bool create_index);
    pair<index_t, index_t> getEdges(index_t node, bool by_from); // Create index first.
    
private:
    void      createPathNode(const PathNode& left, const PathNode& right);
    void      sort();
    pair<index_t, index_t> nextMaximalSet(pair<index_t, index_t> node_range);
    
    void      sortByKey();
    void      sortEdges();  // By (from.label, to.rank)
    
    
    // Can create an index by using key.second in PathNodes.
    void      sortByFrom(bool create_index);
    pair<index_t, index_t> getNodesFrom(index_t node);        // Use sortByFrom(true) first.
    pair<index_t, index_t> getNextRange(pair<index_t, index_t> range);  // Use sortByFrom() first.
    
    // As in Graph.
    void      restoreLabels();
};


template <typename index_t>
PathGraph<index_t>::PathGraph(RefGraph<index_t>& base) :
ranks(0), max_label('Z'), temp_nodes(0), generation(0),
status(error), has_stabilized(false)
{
    if(!base.repOk()) return;

#if 0
    base.changeLabels(this->max_label);
#endif
    
    // Create a path node per edge with a key set to from node's label
    temp_nodes = base.edges.size() + 1;
    nodes.reserveExact(temp_nodes);
    for(index_t i = 0; i < base.edges.size(); i++) {
        const typename RefGraph<index_t>::Edge& e = base.edges[i];
        nodes.expand();
        nodes.back().from = e.from;
        nodes.back().to = e.to;
        nodes.back().key = pair<index_t, index_t>(base.nodes[e.from].label, 0);
    }
    // Final node.
    assert_lt(base.lastNode, base.nodes.size());
    assert_eq(base.nodes[base.lastNode].label, '$');
    nodes.expand();
    nodes.back().from = nodes.back().to = base.lastNode;
    nodes.back().key = pair<index_t, index_t>(0, 0);

#if 0
    base.restoreLabels(this->max_label);
#endif
   
    status = ok;
#if 0
    sort();
#endif
}

#if 0
template <typename index_t>
PathGraph<index_t>::PathGraph(PathGraph<index_t>& previous) :
ranks(0), max_label(previous.max_label), temp_nodes(0), generation(previous.generation + 1),
status(error), has_stabilized(false)
{
    if(previous.status != ok) { return; }
    
    previous.sortByFrom(true);
    
    // A heuristic to determine, whether the number of new nodes should be counted.
    usint new_nodes = previous.node_count + previous.node_count / 8;
    if(previous.ranks >= previous.node_count / 2 && !(previous.has_stabilized)) {
        new_nodes = 0;
        for(nvector::iterator left = previous.nodes.begin(); left != previous.nodes.end(); ++left) {
            if(left->isSorted()) { new_nodes++; continue; }
            pair_type pn_range = previous.getNodesFrom(left->to);
            new_nodes += length(pn_range);
        }
        std::cout << "Trying to allocate space for " << new_nodes << " nodes." << std::endl;
    }
    this->nodes.reserve(new_nodes);
    
    for(nvector::iterator left = previous.nodes.begin(); left != previous.nodes.end(); ++left) {
        if(left->isSorted()) { this->nodes.push_back(*left); continue; }
        
        pair_type pn_range = previous.getNodesFrom(left->to);
        for(usint pn = pn_range.first; pn <= pn_range.second; pn++) {
            this->createPathNode(*left, previous.nodes[pn]);
        }
    }
    this->temp_nodes = this->node_count = this->nodes.size();
    
    this->status = ok;
    this->sort();
    
    if(previous.has_stabilized || (this->generation >= 11 && this->ranks >= 0.8 * new_nodes)) {
        this->has_stabilized = true;
    }
}

#endif

#if 0

void
PathGraph::printInfo()
{
    std::cout << "Generation " << this->generation << " (" << this->temp_nodes << " -> " << this->node_count << " nodes, " << this->ranks << " ranks)" << std::endl;
}

//--------------------------------------------------------------------------

bool PathGraph::generateEdges(Graph& parent) {
    if(this->status != sorted) { return false; }
    
    this->sortByFrom(false);  // Sort nodes by from, do not create index.
    parent.sortEdges(false);  // Sort parent edges by to.
    parent.changeLabels(this->automata, this->max_label);
    
    // Create edges (from, to) as (from.from, to = rank(to)).
    pair_type pn_range = this->getNextRange(EMPTY_PAIR);
    pair_type ge_range = parent.getNextEdgeRange(EMPTY_PAIR, false);
    this->edges.reserve(this->node_count + this->node_count / 4);
    while(!isEmpty(pn_range) && !isEmpty(ge_range)) {
        if(this->nodes[pn_range.first].from == parent.edges[ge_range.first].to) {
            for(usint node = pn_range.first; node <= pn_range.second; node++) {
                for(usint edge = ge_range.first; edge <= ge_range.second; edge++) {
                    uint from = parent.edges[edge].from;
                    this->edges.push_back(PathEdge(from, this->nodes[node].key.first, parent.nodes[from].label));
                }
            }
            pn_range = this->getNextRange(pn_range);
            ge_range = parent.getNextEdgeRange(ge_range, false);
        }
        else if(this->nodes[pn_range.first].from < parent.edges[ge_range.first].to) {
            pn_range = this->getNextRange(pn_range);
        }
        else {
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
    while(node != this->nodes.end() && edge != this->edges.end()) {
        if(edge->from == node->from) {
            edge->from = node - this->nodes.begin(); ++edge;
            node->key.first++;
        } else {
            node->to = parent.nodes[node->from].value;
            ++node; node->key.first = 0;
        }
    }
    if(node != this->nodes.end()) {
        node->to = parent.nodes[node->from].value;
    }
    
    this->sortEdges(false, true);
    
    // daehwan - now we get BWT sequence from the labels in this->edges
    parent.restoreLabels(this->automata, this->max_label);
    this->restoreLabels();
    this->status = ready;
    return true;
}

//--------------------------------------------------------------------------

std::vector<pair_type>* PathGraph::getSamples(usint sample_rate, usint& max_sample, const Graph& parent) {
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

void PathGraph::createPathNode(const PathNode& left, const PathNode& right) {
    PathNode pn;
    pn.from = left.from; pn.to = right.to;
    pn.key = pair_type(left.key.first, right.key.first);
    this->nodes.push_back(pn);
}

void PathGraph::sort() {
    this->sortByKey();
    
    // Update ranks.
    usint rank = 0;
    pair_type key = this->nodes[0].key;
    for(nvector::iterator iter = this->nodes.begin(); iter != this->nodes.end(); ++iter) {
        if(iter->key != key) { rank++; key = iter->key; }
        iter->key = pair_type(rank, 0);
    }
    
    // Merge equivalent nodes.
    usint top = 0;
    pair_type node_range = EMPTY_PAIR;
    while(true) {
        node_range = this->nextMaximalSet(node_range);
        if(isEmpty(node_range)) { break; }
        this->nodes[top] = this->nodes[node_range.first]; top++;
    }
    this->node_count = top; this->nodes.resize(this->node_count);
    
    // Check if sorted.
    bool is_sorted = true;
    PathNode* candidate = &(this->nodes[0]);
    key = candidate->key; this->ranks = 1;
    for(usint i = 1; i < this->node_count; i++) {
        if(this->nodes[i].key != key) {
            if(candidate != 0) { candidate->setSorted(); }
            candidate = &(this->nodes[i]);
            key = candidate->key; this->ranks++;
        }
        else { candidate = 0; is_sorted = false; }
    }
    if(candidate != 0) { candidate->setSorted(); }
    
    // Replace the ranks of a sorted graph so that rank(i) = i.
    // Merges may otherwise leave gaps in the ranks.
    if(is_sorted) {
        for(usint i = 0; i < this->node_count; i++) {
            this->nodes[i].key.first = i;
        }
        this->status = sorted;
    }
}

// Returns the next maximal mergeable set of PathNodes.
// A set of PathNodes sharing adjacent keys is mergeable, if each of the
// PathNodes begins in the same GraphNode, and no other PathNode shares
// the key. If the maximal set is empty, returns the next PathNode.
pair_type PathGraph::nextMaximalSet(pair_type range) {
    if(range.second >= this->node_count - 1) { return EMPTY_PAIR; }
    
    if(isEmpty(range)) { range = pair_type(0, 0); }
    else {
        range.first = range.second + 1; range.second = range.first;
        if(this->nodes[range.first - 1].key == this->nodes[range.first].key) { return range; }
    }
    
    for(usint i = range.second + 1; i < this->node_count; i++) {
        if(this->nodes[i - 1].key != this->nodes[i].key) { range.second = i - 1; }
        if(this->nodes[i].from != this->nodes[range.first].from) { return range; }
    }
    
    range.second = this->node_count - 1;
    return range;
}

//--------------------------------------------------------------------------

struct PathNodeComparator {
    bool operator() (const PathNode& a, const PathNode& b) const
    {
        return (a.key < b.key);
    }
} pn_comparator;

void PathGraph::sortByKey() {
    parallelSort(this->nodes.begin(), this->nodes.end(), pn_comparator);
}

struct PathEdgeComparator {
    bool operator() (const PathEdge& a, const PathEdge& b) const
    {
        return ((a.label < b.label) ||
                ((a.label == b.label) && (a.rank < b.rank)));
    }
} pe_comparator;

void PathGraph::sortEdges() {
    parallelSort(this->edges.begin(), this->edges.end(), pe_comparator);
}

//--------------------------------------------------------------------------

struct PathEdgeFromComparator {
    bool operator() (const PathEdge& a, const PathEdge& b) const
    {
        return ((a.from < b.from) || (a.from == b.from && a.rank < b.rank));
    }
} pe_from_comparator;

struct PathEdgeToComparator {
    bool operator() (const PathEdge& a, const PathEdge& b) const
    {
        return (a.rank < b.rank);
    }
} pe_to_comparator;

void PathGraph::sortEdges(bool by_from, bool create_index) {
    if(by_from) {
        if(this->status != edges_sorted)
        {
            parallelSort(this->edges.begin(), this->edges.end(), pe_from_comparator);
            if(this->status == ready) { this->status = edges_sorted; }
            else                      { this->status = error; }
        }
    } else {
        parallelSort(this->edges.begin(), this->edges.end(), pe_to_comparator);
        this->status = error;
    }
    
    if(create_index) {
        for(nvector::iterator iter = this->nodes.begin(); iter != this->nodes.end(); ++iter) {
            iter->key.second = 0;
        }
        if(by_from) {
            for(evector::iterator iter = this->edges.begin(); iter != this->edges.end(); ++iter) {
                this->nodes[iter->from].key.second++;
            }
        } else {
            for(evector::iterator iter = this->edges.begin(); iter != this->edges.end(); ++iter) {
                this->nodes[iter->rank].key.second++;
            }
        }
        for(usint i = 1; i < this->node_count; i++) {
            this->nodes[i].key.second += this->nodes[i - 1].key.second;
        }
    }
}

pair_type PathGraph::getEdges(usint node, bool by_from) {
    if(node >= this->node_count) {
        std::cout << "Error: Trying to get edges " << (by_from ? "from " : "to ") << node << std::endl;
    }
    if(this->nodes[node].key.second == 0) { return EMPTY_PAIR; }
    if(node == 0) { return pair_type(0, this->nodes[node].key.second - 1); }
    return pair_type(this->nodes[node - 1].key.second, this->nodes[node].key.second - 1);
}

//--------------------------------------------------------------------------

struct PathNodeFromComparator {
    bool operator() (const PathNode& a, const PathNode& b) const {
        return (a.from < b.from);
    }
} pn_from_comparator;

void PathGraph::sortByFrom(bool create_index) {
    parallelSort(this->nodes.begin(), this->nodes.end(), pn_from_comparator);
    this->status = error;
    
    if(create_index) {
        usint current = 0;
        this->nodes[0].key.second = 0;
        for(nvector::iterator iter = this->nodes.begin(); iter != this->nodes.end(); ++iter) {
            while(current < iter->from) {
                current++;
                this->nodes[current].key.second = iter - this->nodes.begin();
            }
        }
        for(current++; current < this->node_count; current++) {
            this->nodes[current].key.second = this->nodes.size();
        }
    }
}

pair_type PathGraph::getNodesFrom(uint node) {
    if(node + 1 >= this->node_count) { return pair_type(this->nodes[node].key.second, this->node_count - 1); }
    if(this->nodes[node].key.second == this->nodes[node + 1].key.second) { return EMPTY_PAIR; }
    return pair_type(this->nodes[node].key.second, this->nodes[node + 1].key.second - 1);
}

pair_type PathGraph::getNextRange(pair_type range) {
    if(range.second >= this->node_count - 1) { return EMPTY_PAIR; }
    
    if(isEmpty(range)) { range = pair_type(0, 0); } // Create the first range.
    else { range.second++; range.first = range.second; }
    
    while(range.second < this->node_count - 1 && this->nodes[range.second + 1].from == this->nodes[range.first].from) {
        range.second++;
    }
    
    return range;
}

//--------------------------------------------------------------------------

void PathGraph::restoreLabels() {
    for(usint i = 0; i < this->edge_count; i++) {
        if(this->edges[i].label < this->automata) { this->edges[i].label = 0; }
        else {
            this->edges[i].label -= this->automata;
            if(this->edges[i].label >= this->max_label) { this->edges[i].label = this->max_label; }
        }
    }
}

//--------------------------------------------------------------------------

#endif

#endif /*GBWT_GRAPH_H_*/
