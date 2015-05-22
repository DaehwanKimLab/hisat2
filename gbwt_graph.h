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
    pair<index_t, index_t> getNextEdgeRange(pair<index_t, index_t> range, bool from) {
        if(range.second >= edges.size()) {
            return pair<index_t, index_t>(0, 0);
        }
        range.first = range.second; range.second++;
        
        if(from) {
            while(range.second < edges.size() && edges[range.second].from == edges[range.first].from) {
                range.second++;
            }
        } else {
            while(range.second < edges.size() && edges[range.second].to == edges[range.first].to) {
                range.second++;
            }
        }
        return range;
    }
    
private:
    EList<Node> nodes;
    EList<Edge> edges;
    index_t     lastNode; // $
    
#ifndef NDEBUG
    bool        debug;
#endif
    
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
    
#ifndef NDEBUG
    debug = (jlen <= 10);
#endif
    
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
            ASSERT_ONLY(RefRecord rec =) fastaRefReadAppend(*fb, first, s, distoff, rpcp, &refnames.back());
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
    
#ifndef NDEBUG
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
#endif
    
    if(!isReverseDeterministic()) {
        cerr << "\tis not reverse-deterministic, so reverse-determinize..." << endl;
        reverseDeterminize();
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
            if(mid == 0) {
                return pair<index_t, index_t>(0, 0);
            }
            
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
void RefGraph<index_t>::reverseDeterminize()
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
    
#ifndef NDEBUG
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
#endif
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
        
        bool operator< (const PathNode& o) const {
            return key < o.key;
        };
    };
    
    struct PathNodeFromCmp {
        bool operator() (const PathNode& a, const PathNode& b) const {
            return a.from < b.from;
        }
    };
    
    struct PathEdge {
        index_t from;
        index_t ranking;
        char    label;
        
        PathEdge() { reset(); }
        
        PathEdge(index_t from_, index_t ranking_, char label_) : from(from_), ranking(ranking_), label(label_) {}
        
        void reset() {
            from = 0;
            ranking = 0;
            label = 0;
        }
        
        bool operator< (const PathEdge& o) const {
            return label < o.label || (label == o.label && ranking < o.ranking);
        };
    };
    
    struct PathEdgeFromCmp {
        bool operator() (const PathEdge& a, const PathEdge& b) const {
            return a.from < b.from || (a.from == b.from && a.ranking < b.ranking);
        }
    };
    
    struct PathEdgeToCmp {
        bool operator() (const PathEdge& a, const PathEdge& b) const {
            return a.ranking < b.ranking;
        }
    };
    
public:
    // Create a new graph in which paths are represented using nodes
    PathGraph(RefGraph<index_t>& parent);
    
    // Construct a next 2^(j+1) path graph from a 2^j path graph
    PathGraph(PathGraph<index_t>& previous);
    
    ~PathGraph() {}
    
    void printInfo();
    
    bool IsSorted() const { return status != sorted; }
    bool repOk() const { return status != ok; }
    
    
    // status must be sorted. Calling this invalidates parent and
    // sets status to ready.
    // Writes outdegree to PathNode.key.second, value to PathNode.to, and
    // predecessor labels to PathNode.key.first.
    // Restores the labels of parent.
    bool generateEdges(RefGraph<index_t>& parent);

private:
    void      createPathNode(const PathNode& left, const PathNode& right);
    void      updateRank_and_merge();
    pair<index_t, index_t> nextMaximalSet(pair<index_t, index_t> node_range);
    
    void sortByKey() { sort(nodes.begin(), nodes.end()); } // by key

    // Can create an index by using key.second in PathNodes.
    void sortByFrom(bool create_index = true);
    pair<index_t, index_t> getNodesFrom(index_t node);        // Use sortByFrom(true) first.
    pair<index_t, index_t> getNextRange(pair<index_t, index_t> range);  // Use sortByFrom() first.

    void sortEdges() { sort(edges.begin(), edges.end()); } // by (from.label, to.rank)
    void sortEdgesFrom() {
        if(this->status == edges_sorted)
            return;
        sort(edges.begin(), edges.end(), PathEdgeFromCmp());
        status = edges_sorted;
    }
    void sortEdgesTo(bool create_index = false) {
        sort(edges.begin(), edges.end(), PathEdgeToCmp());
        status = error;
        
        if(create_index) {
            for(PathNode* node = nodes.begin(); node != nodes.end(); node++) {
                node->key.second = 0;
            }
            for(PathEdge* edge = edges.begin(); edge != edges.end(); edge++) {
                nodes[edge->ranking].key.second++;
            }
            for(index_t i = 1; i < nodes.size(); i++) {
                nodes[i].key.second += nodes[i - 1].key.second;
            }
        }
    }
    
private:
    // status must be ready.
    EList<pair<index_t, index_t> >* getSamples(index_t sample_rate, index_t& max_sample, const RefGraph<index_t>& base);
    
    EList<PathNode> nodes;
    EList<PathEdge> edges;
    index_t         ranks;
    index_t         max_label;  // Node label of initial nodes.
    index_t         temp_nodes; // Total number of nodes created before sorting.
    
    index_t         generation; // Sorted by paths of length 2^generation.
    
    enum status_t { error, ok, sorted, ready, edges_sorted } status;
    bool            has_stabilized;  // The number of nodes will probably not explode in subsequent doublings.
    
    // Can create an index by using key.second in PathNodes.
    // If the graph is not ready, its status becomes error.
    // Sorting edges by from actually sorts them by (from, to).
    void      sortEdges(bool by_from, bool create_index);
    pair<index_t, index_t> getEdges(index_t node, bool by_from); // Create index first.
    
#ifndef NDEBUG
    bool               debug;
    
    EList<char>        bwt_string;
    EList<char>        F_array;
    EList<char>        M_array;
    EList<index_t>     bwt_counts;
    
    // brute-force implementations
    index_t select(const EList<char>& array, index_t p, char c) {
        if(p <= 0) return 0;
        for(index_t i = 0; i < array.size(); i++) {
            if(array[i] == c) {
                assert_gt(p, 0);
                p--;
                if(p == 0)
                    return i;
            }
        }
        return array.size();
    }
    
    index_t select1(const EList<char>& array, index_t p) {
        return select(array, p, 1);
    }
    
    index_t rank(const EList<char>& array, index_t p, char c) {
        index_t count = 0;
        assert_lt(p, array.size());
        for(index_t i = 0; i <= p; i++) {
            if(array[i] == c)
                count++;
        }
        return count;
    }
    
    index_t rank1(const EList<char>& array, index_t p) {
        return rank(array, p, 1);
    }
#endif
};


template <typename index_t>
PathGraph<index_t>::PathGraph(RefGraph<index_t>& base) :
ranks(0), max_label('Z'), temp_nodes(0), generation(0),
status(error), has_stabilized(false)
{
    if(!base.repOk()) return;

#ifndef NDEBUG
    debug = base.nodes.size() <= 20;
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
   
    status = ok;
    updateRank_and_merge();
}

template <typename index_t>
PathGraph<index_t>::PathGraph(PathGraph<index_t>& previous) :
ranks(0), max_label(previous.max_label), temp_nodes(0), generation(previous.generation + 1),
status(error), has_stabilized(false)
{
    if(previous.status != ok)
        return;
    
#ifndef NDEBUG
    debug = previous.debug;
#endif
    
    previous.sortByFrom(true);
    
#ifndef NDEBUG
    if(debug) {
        cerr << "Path nodes (" << generation << "-generation) - soryByFrom" << endl;
        for(size_t i = 0; i < previous.nodes.size(); i++) {
            const PathNode& node = previous.nodes[i];
            cerr << "\t" << i << "\t(" << node.key.first << ", " << node.key.second << ")\t"
            << node.from << " --> " << node.to << (node.isSorted() ? "\tsorted" : "") << endl;
        }
    }
#endif
    
    // A heuristic to determine, whether the number of new nodes should be counted
    index_t new_nodes = previous.nodes.size() + previous.nodes.size() / 8;
    if(previous.ranks >= previous.nodes.size() / 2 && !(previous.has_stabilized)) {
        new_nodes = 0;
        for(index_t i = 0; i < previous.nodes.size(); i++) {
            if(previous.nodes[i].isSorted()) { new_nodes++; continue; }
            pair<index_t, index_t> pn_range = previous.getNodesFrom(previous.nodes[i].to);
            assert_lt(pn_range.first, pn_range.second);
            new_nodes += (pn_range.second - pn_range.first);
        }
        cerr << "Trying to allocate space for " << new_nodes << " nodes." << endl;
    }
    nodes.resizeExact(new_nodes);
    nodes.clear();
    
    for(index_t i = 0; i < previous.nodes.size(); i++) {
        if(previous.nodes[i].isSorted()) {
            nodes.push_back(previous.nodes[i]);
            continue;
        }
        pair<index_t, index_t> pn_range = previous.getNodesFrom(previous.nodes[i].to);
        for(index_t pn = pn_range.first; pn < pn_range.second; pn++) {
            createPathNode(previous.nodes[i], previous.nodes[pn]);
        }
    }
    temp_nodes = nodes.size();
    
#ifndef NDEBUG
    if(debug) {
        cerr << "Path nodes (" << generation << "-generation) - combine" << endl;
        for(size_t i = 0; i < nodes.size(); i++) {
            const PathNode& node = nodes[i];
            cerr << "\t" << i << "\t(" << node.key.first << ", " << node.key.second << ")\t"
                 << node.from << " --> " << node.to << (node.isSorted() ? "\tsorted" : "") << endl;
        }
    }
#endif
    
    status = ok;
    updateRank_and_merge();
    
    if(previous.has_stabilized || (generation >= 11 && ranks >= 0.8 * new_nodes)) {
        has_stabilized = true;
    }
}

template <typename index_t>
void PathGraph<index_t>::printInfo()
{
    cerr << "Generation " << generation
         << " (" << temp_nodes << " -> " << nodes.size() << " nodes, "
         << ranks << " ranks)" << endl;
}

//--------------------------------------------------------------------------
template <typename index_t>
bool PathGraph<index_t>::generateEdges(RefGraph<index_t>& base)
{
    if(status != sorted)
        return false;
    
    sortByFrom(false);
    base.sortEdgesTo();
    
    // Create edges (from, to) as (from.from, to = rank(to))
    pair<index_t, index_t> pn_range = getNextRange(pair<index_t, index_t>(0, 0));
    pair<index_t, index_t> ge_range = base.getNextEdgeRange(pair<index_t, index_t>(0, 0), false);
    edges.resizeExact(nodes.size() + nodes.size() / 4); edges.clear();
    while(pn_range.first < pn_range.second && ge_range.first < ge_range.second) {
        if(nodes[pn_range.first].from == base.edges[ge_range.first].to) {
            for(index_t node = pn_range.first; node < pn_range.second; node++) {
                for(index_t edge = ge_range.first; edge < ge_range.second; edge++) {
                    index_t from = base.edges[edge].from;
                    edges.push_back(PathEdge(from, nodes[node].key.first, base.nodes[from].label));
                }
            }
            pn_range = getNextRange(pn_range);
            ge_range = base.getNextEdgeRange(ge_range, false);
        } else if(nodes[pn_range.first].from < base.edges[ge_range.first].to) {
            pn_range = getNextRange(pn_range);
        } else {
            ge_range = base.getNextEdgeRange(ge_range, false);
        }
    }
    
#ifndef NDEBUG
    if(debug) {
        cerr << "just after creating path edges" << endl;
        cerr << "Ref edges" << endl;
        for(size_t i = 0; i < base.edges.size(); i++) {
            const typename RefGraph<index_t>::Edge& edge = base.edges[i];
            cerr << "\t" << i << "\t" << edge.from << " --> " << edge.to << endl;
        }
        
        cerr << "Path nodes" << endl;
        for(size_t i = 0; i < nodes.size(); i++) {
            const PathNode& node = nodes[i];
            cerr << "\t" << i << "\t(" << node.key.first << ", " << node.key.second << ")\t"
            << node.from << " --> " << node.to << endl;
        }
        
        cerr << "Path edges" << endl;
        for(size_t i = 0; i < edges.size(); i++) {
            const PathEdge& edge = edges[i];
            cerr << "\t" << i << "\tfrom: " << edge.from << "\tranking: " << edge.ranking << "\t" << edge.label << endl;
        }
    }
#endif
    
    sortByKey(); // Restore correct node order
    sortEdges(); // Sort edges by (from.label, to.rank)
    
#ifndef NDEBUG
    if(debug) {
        cerr << "after sorting nodes by ranking and edges by label and ranking" << endl;
        cerr << "Path nodes" << endl;
        for(size_t i = 0; i < nodes.size(); i++) {
            const PathNode& node = nodes[i];
            cerr << "\t" << i << "\t(" << node.key.first << ", " << node.key.second << ")\t"
            << node.from << " --> " << node.to << endl;
        }
        
        cerr << "Path edges" << endl;
        for(size_t i = 0; i < edges.size(); i++) {
            const PathEdge& edge = edges[i];
            cerr << "\t" << i << "\tfrom: " << edge.from << "\tranking: " << edge.ranking << "\t" << edge.label << endl;
        }
    }
#endif
    
    // Sets PathNode.to = GraphNode.value and PathNode.key.first to outdegree
    // Replaces (from.from, to) with (from, to)
    PathNode* node = nodes.begin(); node->key.first = 1;
    PathEdge* edge = edges.begin();
    while(node != nodes.end() && edge != edges.end()) {
        if(edge->from == node->from) {
            edge->from = node - nodes.begin(); edge++;
            node->key.first++;
        } else {
            node->to = base.nodes[node->from].value;
            node++; node->key.first = 0;
        }
    }
    if(node != nodes.end()) {
        node->to = base.nodes[node->from].value;
    }
    
    sortEdgesTo(true);
    status = ready;
    
#ifndef NDEBUG
    if(debug) {
        cerr << "Path nodes (final)" << endl;
        for(size_t i = 0; i < nodes.size(); i++) {
            const PathNode& node = nodes[i];
            cerr << "\t" << i << "\t(" << node.key.first << ", " << node.key.second << ")\t"
            << node.from << " --> " << node.to << endl;
        }
        
        cerr << "Path edges (final)" << endl;
        for(size_t i = 0; i < edges.size(); i++) {
            const PathEdge& edge = edges[i];
            cerr << "\t" << i << "\tfrom: " << edge.from << "\tranking: " << edge.ranking << "\t" << edge.label << endl;
        }
        
        nodes.pop_back(); // Remove 'Z' node
        bwt_string.clear();
        F_array.clear();
        M_array.clear();
        bwt_counts.resizeExact(5); bwt_counts.fillZero(); bwt_counts.front() = 1;
        for(index_t node = 0; node < nodes.size(); node++) {
            pair<index_t, index_t> edge_range = getEdges(node, false);
            for(index_t i = edge_range.first; i < edge_range.second; i++) {
                assert_lt(i, edges.size());
                char label = edges[i].label;
                if(label == 'Z') {
                    label = '$';
                }
                bwt_string.push_back(label);
                F_array.push_back(i == edge_range.first ? 1 : 0);
                
                if(label != '$') {
                    char nt = asc2dna[(int)label];
                    assert_lt(nt + 1, bwt_counts.size());
                    bwt_counts[nt + 1]++;
                }
            }
            for(index_t i = 0; i < nodes[node].key.first; i++) {
                M_array.push_back(i == 0 ? 1 : 0);
            }
        }
        
        assert_gt(bwt_string.size(), 0);
        assert_eq(bwt_string.size(), F_array.size());
        assert_eq(bwt_string.size(), M_array.size());
        cerr << "i\tBWT\tF\tM" << endl;
        for(index_t i = 0; i < bwt_string.size(); i++) {
            cerr << i << "\t" << bwt_string[i] << "\t"  // BWT char
                 << (int)F_array[i] << "\t"             // F bit value
                 << (int)M_array[i] << endl;            // M bit value
        }
        
        for(size_t i = 0; i < bwt_counts.size(); i++) {
            if(i > 0) bwt_counts[i] += bwt_counts[i - 1];
            cerr << i << "\t" << bwt_counts[i] << endl;
        }

        // Test searches, based on paper_example
        EList<string> queries;  EList<index_t> answers;
        queries.push_back("GACGT"); answers.push_back(9);
        queries.push_back("GATGT"); answers.push_back(9);
        queries.push_back("GACT");  answers.push_back(9);
        queries.push_back("ATGT");  answers.push_back(4);
        queries.push_back("GTAC");  answers.push_back(10);
        queries.push_back("ACTG");  answers.push_back(3);
        
        for(size_t q = 0; q < queries.size(); q++) {
            const string& query = queries[q];
            assert_gt(query.length(), 0);
            index_t top = 0, bot = nodes.size();
            cerr << "Aligning " << query <<  endl;
            
            for(size_t i = 0; i < query.length(); i++) {
                if(top >= bot) break;
                
                int nt = query[query.length() - i - 1];
                nt = asc2dna[nt];
                assert_lt(nt, 4);
                cerr << "\t" << i << ": " << "ACGT"[nt] << endl;
                cerr << "\t\tnode range: [" << top << ", " << bot << ")" << endl;
                
                top = select1(F_array, top + 1);
                bot = select1(F_array, bot + 1);
                cerr << "\t\tBWT range: [" << top << ", " << bot << ")" << endl;
                
                top = bwt_counts[(int)nt] + (top <= 0 ? 0 : rank(bwt_string, top - 1, "ACGT"[nt]));
                bot = bwt_counts[(int)nt] + rank(bwt_string, bot - 1, "ACGT"[nt]);
                cerr << "\t\tLF BWT range: [" << top << ", " << bot << ")" << endl;
                
                top = rank1(M_array, top) - 1;
                bot = rank1(M_array, bot - 1);
                cerr << "\t\tnode range: [" << top << ", " << bot << ")" << endl;
            }
            assert_eq(top, answers[q]);
            cerr << "finished... ";
            if(top < nodes.size()) {
                cerr << "being aligned at " << nodes[top].to;
            }
            cerr << endl << endl;
        }
    }
#endif
    
    return true;
}


template <typename index_t>
EList<pair<index_t, index_t> >* PathGraph<index_t>::getSamples(index_t sample_rate, index_t& max_sample, const RefGraph<index_t>& base)
{
    if(status != ready && status != edges_sorted)
        return NULL;

    EList<pair<index_t, index_t> >* sample_pairs = new EList<pair<index_t, index_t> >();
#if 0
    sortEdges(true, true);
    WriteBuffer* found_nodes = new WriteBuffer(nodes.size(), 1);
    ReadBuffer* found = found_nodes->getReadBuffer();
    WriteBuffer* sampled_nodes = new WriteBuffer(nodes.size(), 1);
    ReadBuffer* sampled = sampled_nodes->getReadBuffer();
    
    EList<index_t> active; // Push the initial nodes into stack.
    active.push_back(nodes.size() - 1);
    index_t unsampled = 0, current = 0;
    max_sample = 0;
    index_t total_found = 0;
    while(!active.empty()) {
        current = active.top(); active.pop();
        unsampled = 1;
        while(true) {
            // Nodes with multiple predecessors must be sampled.
            if(found->readItem(current)) {
                if(!sampled->readItem(current)) {
                    sample_pairs->push_back(pair_type(current, this->nodes[current].value()));
                    sampled_nodes->goToItem(current); sampled_nodes->writeItem(1);
                    max_sample = std::max(max_sample, (usint)(this->nodes[current].value()));
                }
                break;
            }
            found_nodes->goToItem(current); found_nodes->writeItem(1); total_found++;
            pair_type successors = getEdges(current, true);
            // No nearby samples.
            // Multiple or no outgoing edges.
            // Discontinuity in values.
            if(unsampled >= sample_rate || length(successors) != 1 ||
               this->nodes[this->edges[successors.first].rank].value() != this->nodes[current].value() + 1) {
                sample_pairs->push_back(pair_type(current, this->nodes[current].value()));
                sampled_nodes->goToItem(current); sampled_nodes->writeItem(1);
                max_sample = std::max(max_sample, (usint)(this->nodes[current].value()));
                for(usint suc = successors.first + 1; suc <= successors.second; suc++) {
                    active.push(this->edges[suc].rank);
                }
                unsampled = 0;
            }
            if(isEmpty(successors)) break;
            
            current = this->edges[successors.first].rank;
            unsampled++;
        }
    }
    delete found_nodes; delete found;
    delete sampled_nodes; delete sampled;
#endif
    
    //  cerr << "Found " << total_found << " nodes" << endl;
    return sample_pairs;
}

//--------------------------------------------------------------------------
template <typename index_t>
void PathGraph<index_t>::createPathNode(const PathNode& left, const PathNode& right)
{
    nodes.expand();
    nodes.back().from = left.from;
    nodes.back().to = right.to;
    nodes.back().key = pair<index_t, index_t>(left.key.first, right.key.first);
}

template <typename index_t>
void PathGraph<index_t>::updateRank_and_merge()
{
    sortByKey();
    
#ifndef NDEBUG
    if(debug) {
        cerr << "Path nodes (" << generation << "-generation) before merge" << endl;
        for(size_t i = 0; i < nodes.size(); i++) {
            const PathNode& node = nodes[i];
            cerr << "\t" << i << "\t(" << node.key.first << ", " << node.key.second << ")\t"
                 << node.from << " --> " << node.to << endl;
        }
    }
#endif
    
    // Update ranks
    index_t rank = 0;
    pair<index_t, index_t> key = nodes.front().key;
    for(index_t i = 0; i < nodes.size(); i++) {
        PathNode& node = nodes[i];
        if(node.key != key) {
            key = node.key;
            rank++;
        }
        node.key = pair<index_t, index_t>(rank, 0);
    }

    // Merge equivalent nodes
    index_t curr = 0;
    pair<index_t, index_t> range(0, 0); // Empty range
    while(true) {
        range = nextMaximalSet(range);
        if(range.first >= range.second)
            break;
        nodes[curr] = nodes[range.first]; curr++;
    }
    nodes.resize(curr);
    
    PathNode* candidate = &nodes.front();
    key = candidate->key; ranks = 1;
    for(index_t i = 1; i < nodes.size(); i++) {
        if(nodes[i].key != key) {
            if(candidate != NULL) {
                candidate->setSorted();
            }
            candidate = &nodes[i];
            key = candidate->key; ranks++;
        } else {
            candidate = NULL;
        }
    }
    if(candidate != NULL) {
        candidate->setSorted();
    }
    
    // Replace the ranks of a sorted graph so that rank(i) = i
    // Merges may otherwise leave gaps in the ranks
    if(ranks == nodes.size()) {
        for(index_t i = 0; i < nodes.size(); i++) {
            nodes[i].key.first = i;
        }
        status = sorted;
    }
    
#ifndef NDEBUG
    if(debug) {
        cerr << "Path nodes (" << generation << "-generation) after merge" << endl;
        for(size_t i = 0; i < nodes.size(); i++) {
            const PathNode& node = nodes[i];
            cerr << "\t" << i << "\t(" << node.key.first << ", " << node.key.second << ")\t"
                 << node.from << " --> " << node.to << (node.isSorted() ? "\tsorted" : "") << endl;
        }
    }
#endif
}

// Returns the next maximal mergeable set of PathNodes.
// A set of PathNodes sharing adjacent keys is mergeable, if each of the
// PathNodes begins in the same GraphNode, and no other PathNode shares
// the key. If the maximal set is empty, returns the next PathNode.
template <typename index_t>
pair<index_t, index_t> PathGraph<index_t>::nextMaximalSet(pair<index_t, index_t> range) {
    if(range.second >= nodes.size()) {
        return pair<index_t, index_t>(0, 0);
    }
    
    range.first = range.second;
    range.second = range.first + 1;
    if(range.first > 0 && nodes[range.first - 1].key == nodes[range.first].key) {
        return range;
    }
    
    for(index_t i = range.second; i < nodes.size(); i++) {
        if(nodes[i - 1].key != nodes[i].key) {
            range.second = i;
        }
        if(nodes[i].from != nodes[range.first].from) {
            return range;
        }
    }
    
    range.second = nodes.size();
    return range;
}

template <typename index_t>
pair<index_t, index_t> PathGraph<index_t>::getEdges(index_t node, bool by_from) {
    if(node >= nodes.size()) {
        cerr << "Error: Trying to get edges " << (by_from ? "from " : "to ") << node << endl;
    }
    if(nodes[node].key.second == 0) {
        return pair<index_t, index_t>(0, 0);
    }
    if(node == 0) {
        return pair<index_t, index_t>(0, nodes[node].key.second);
    } else {
        return pair<index_t, index_t>(nodes[node - 1].key.second, nodes[node].key.second);
    }
}

template <typename index_t>
void PathGraph<index_t>::sortByFrom(bool create_index) {
    sort(nodes.begin(), nodes.end(), PathNodeFromCmp());
    status = error;
    
    if(create_index) {
        index_t current = 0;
        nodes.front().key.second = 0;
        for(index_t i = 0; i < nodes.size(); i++) {
            while(current < nodes[i].from) {
                current++;
                assert_lt(current, nodes.size());
                nodes[current].key.second = i;
            }
        }
        for(current++; current < nodes.size(); current++) {
            nodes[current].key.second = nodes.size();
        }
    }
}

template <typename index_t>
pair<index_t, index_t> PathGraph<index_t>::getNodesFrom(index_t node) {
    if(node + 1 >= nodes.size()) {
        assert_eq(node + 1, nodes.size());
        return pair<index_t, index_t>(nodes.back().key.second, nodes.size());
    }
    if(nodes[node].key.second == nodes[node + 1].key.second) {
        return pair<index_t, index_t>(0, 0);
    }
    return pair<index_t, index_t>(nodes[node].key.second, nodes[node + 1].key.second);
}

template <typename index_t>
pair<index_t, index_t> PathGraph<index_t>::getNextRange(pair<index_t, index_t> range) {
    if(range.second >= nodes.size())
        return pair<index_t, index_t>(0, 0);
    
    if(range.first < range.second) {
        range.first = range.second; range.second++;
    }
    
    while(range.second < nodes.size() && nodes[range.second].from == nodes[range.first].from) {
        range.second++;
    }
    
    return range;
}

#endif /*GBWT_GRAPH_H_*/
