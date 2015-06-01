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
    
    bool operator< (const SNP& o) const {
        if(pos != o.pos) return pos < o.pos;
        if(type != o.type) return type < o.type;
        if(diff.size() != o.diff.size()) return diff.size() < o.diff.size();
        for(index_t i = 0; i < diff.size(); i++) {
            if(diff[i] != o.diff[i]) return diff[i] < o.diff[i];
        }
        return false;
    }
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
        
        bool write(ofstream& f_out, bool bigEndian) const {
            writeIndex<index_t>(f_out, value, bigEndian);
            writeU16(f_out, label, bigEndian);
            writeU16(f_out, backbone, bigEndian);
            return true;
        }
        
        bool read(ifstream& f_in, bool bigEndian) {
            value = readIndex<index_t>(f_in, bigEndian);
            label = (char)readU16(f_in, bigEndian);
            backbone = (bool)readU16(f_in, bigEndian);
            return true;
        }
        
        bool operator== (const Node& o) const {
            if(value != o.value) return false;
            if(label != o.label) return false;
            return true;
        }
        
        bool operator< (const Node& o) const {
            if(value != o.value) return value < o.value;
            return label < o.label;
        }
    };
    
    struct Edge {
        index_t from; // from Node
        index_t to;   // to Node
        
        Edge() {}
        Edge(index_t from_, index_t to_) : from(from_), to(to_) {}
        
        bool write(ofstream& f_out, bool bigEndian) const {
            writeIndex<index_t>(f_out, from, bigEndian);
            writeIndex<index_t>(f_out, to, bigEndian);
            return true;
        }
        
        bool read(ifstream& f_in, bool bigEndian) {
            from = readIndex<index_t>(f_in, bigEndian);
            to = readIndex<index_t>(f_in, bigEndian);
            return true;
        }
        
        bool operator< (const Edge& o) const {
            if(from != o.from) return from < o.from;
            return to < o.to;
        }
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
    
    // for debugging purposes
    struct TempEdgeNodeCmp {
        TempEdgeNodeCmp(const EList<Node>& nodes_) : nodes(nodes_) {}
        bool operator() (const Edge& a, const Edge& b) const {
            // Compare "from" nodes
            {
                assert_lt(a.from, nodes.size());
                const Node& a_node = nodes[a.from];
                assert_lt(b.from, nodes.size());
                const Node& b_node = nodes[b.from];
                if(!(a_node == b_node))
                    return a_node < b_node;
            }

            // Compare "to" nodes
            assert_lt(a.to, nodes.size());
            const Node& a_node = nodes[a.to];
            assert_lt(b.to, nodes.size());
            const Node& b_node = nodes[b.to];
            return a_node < b_node;
        }
        
        const EList<Node>& nodes;
    };
    
public:
    RefGraph(const SString<char>& s,
             const EList<RefRecord>& szs,
             const EList<SNP<index_t> >& snps,
             const string& out_fname,
             bool verbose);
    
    bool repOk() { return true; }
    
    void write(const std::string& base_name) {}
    void printInfo() {}
    
private:
    bool isReverseDeterministic();
    void reverseDeterminize(index_t lastNode_add = 0);
    
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
    EList<RefRecord> szs;
    
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
RefGraph<index_t>::RefGraph(const SString<char>& s,
                            const EList<RefRecord>& szs,
                            const EList<SNP<index_t> >& snps,
                            const string& out_fname,
                            bool verbose)
: lastNode(0)
{
    const bool bigEndian = false;
    assert_gt(szs.size(), 0);
    index_t jlen = s.length();
    
#ifndef NDEBUG
    debug = (jlen <= 10);
#endif
    
    // a memory-efficient way to create a population graph with known SNPs
    bool frag_automaton = true;
    if(frag_automaton) {
        const string rg_fname = out_fname + ".rf";
        ofstream rg_out_file(rg_fname.c_str(), ios::binary);
        if(!rg_out_file.good()) {
            cerr << "Could not open file for writing a reference graph: \"" << rg_fname << "\"" << endl;
            throw 1;
        }
        
        EList<RefRecord> tmp_szs;
        {
            index_t relax = 10;
            EList<pair<index_t, index_t> > snp_ranges; // each range inclusive
            for(index_t i = 0; i < snps.size(); i++) {
                const SNP<index_t>& snp = snps[i];
                pair<index_t, index_t> range;
                range.first = snp.pos > relax ? snp.pos - relax - 1 : 0;
                if(snp.type == SNP_SGL) {
                    range.second = snp.pos + 1;
                } else if(snp.type == SNP_DEL) {
                    range.second = snp.pos + snp.diff.size();
                } else if(snp.type == SNP_INS) {
                    range.second = snp.pos;
                } else assert(false);
                range.second += relax;
                
                if(snp_ranges.empty() || snp_ranges.back().second + 1 < range.first) {
                    snp_ranges.push_back(range);
                } else {
                    assert_leq(snp_ranges.back().first, range.first);
                    if(snp_ranges.back().second < range.second) {
                        snp_ranges.back().second = range.second;
                    }
                }
            }
            
            const index_t chunk_size = 1 << 20;
            index_t pos = 0, range_idx = 0;
            for(index_t i = 0; i < szs.size(); i++) {
                if(szs[i].len == 0) continue;
                if(szs[i].len <= chunk_size) {
                    tmp_szs.push_back(szs[i]);
                    pos += szs[i].len;
                } else {
                    index_t num_chunks = (szs[i].len + chunk_size - 1) / chunk_size;
                    assert_gt(num_chunks, 1);
                    index_t modified_chunk_size = szs[i].len / num_chunks;
                    index_t after_pos = pos + szs[i].len;
                    ASSERT_ONLY(index_t sum_len = 0);
                    while(pos < after_pos) {
                        index_t target_pos = pos + modified_chunk_size;
                        if(target_pos < after_pos) {
                            for(; range_idx < snp_ranges.size(); range_idx++) {
                                if(target_pos < snp_ranges[range_idx].first) break;
                            }
                            pair<index_t, index_t> snp_free_range;
                            if(range_idx == 0) {
                                snp_free_range.first = 0;
                            } else {
                                snp_free_range.first = snp_ranges[range_idx - 1].second + 1;
                            }
                        
                            if(range_idx == snp_ranges.size()) {
                                snp_free_range.second = jlen - 1;
                            } else {
                                snp_free_range.second = snp_ranges[range_idx].first - 1;
                            }
                            
                            assert_leq(snp_free_range.first, snp_free_range.second);
                            if(target_pos < snp_free_range.first) target_pos = snp_free_range.first;
                            if(target_pos > after_pos) target_pos = after_pos;
                        } else {
                            target_pos = after_pos;
                        }
                        
                        tmp_szs.expand();
                        tmp_szs.back().len = target_pos - pos;
                        tmp_szs.back().off = 0;
                        pos = target_pos;
                        ASSERT_ONLY(sum_len += tmp_szs.back().len);
                    }
                    assert_eq(pos, after_pos);
                    assert_eq(sum_len, szs[i].len);
                }
            }
#ifndef NDEBUG
            index_t modified_jlen = 0;
            for(index_t i = 0; i < tmp_szs.size(); i++) {
                modified_jlen += tmp_szs[i].len;
            }
            assert_eq(modified_jlen, jlen);
#endif
        }

        index_t num_nodes = 0, num_edges = 0, curr_pos = 0;
        index_t szs_idx = 0, snp_idx = 0, num_snp_nodes = 0;
        assert(tmp_szs[szs_idx].first);
        EList<index_t> prev_tail_nodes;
        for(; szs_idx < tmp_szs.size(); szs_idx++) {
            index_t curr_len = tmp_szs[szs_idx].len;
            if(curr_len <= 0) continue;
            
            index_t num_predicted_nodes = (index_t)(curr_len * 1.2);
            nodes.clear(); nodes.reserveExact(num_predicted_nodes);
            edges.clear(); edges.reserveExact(num_predicted_nodes);
            
            // Created head node
            nodes.expand();
            nodes.back().label = 'Z';
            nodes.back().value = 0;
            
            // Create nodes and edges corresponding to a reference genome
            assert_leq(curr_pos + curr_len, s.length());
            for(size_t i = curr_pos; i < curr_pos + curr_len; i++) {
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
            for(; snp_idx < snps.size(); snp_idx++) {
                const SNP<index_t>& snp = snps[snp_idx];
                assert_geq(snp.pos, curr_pos);
                if(snp.pos >= curr_pos + curr_len) break;
                if(snp.type == SNP_SGL) {
                    assert_eq(snp.diff.size(), 1);
                    nodes.expand();
                    num_snp_nodes++;
                    nodes.back().label = snp.diff[0];
                    nodes.back().value = jlen + 2 + num_snp_nodes;
                    
                    edges.expand();
                    edges.back().from = snp.pos - curr_pos;
                    edges.back().to = nodes.size() - 1;
                    
                    edges.expand();
                    edges.back().from = nodes.size() - 1;
                    edges.back().to = snp.pos - curr_pos + 2;
                } else if(snp.type == SNP_DEL) {
                    index_t deletionLen = snp.diff.size();
                    assert_gt(deletionLen, 0);
                    edges.expand();
                    edges.back().from = snp.pos - curr_pos;
                    edges.back().to = snp.pos + - curr_pos + deletionLen + 1;
                } else if(snp.type == SNP_INS) {
                    assert_gt(snp.diff.size(), 0);
                    for(size_t j = 0; j < snp.diff.size(); j++) {
                        char inserted_bp = snp.diff[j];
                        nodes.expand();
                        num_snp_nodes++;
                        nodes.back().label = inserted_bp;
                        nodes.back().value = jlen + 2 + num_snp_nodes;
                        
                        edges.expand();
                        edges.back().from = (j == 0 ? snp.pos - curr_pos : nodes.size() - 2);
                        edges.back().to = nodes.size() - 1;
                    }
                    edges.expand();
                    edges.back().from = nodes.size() - 1;
                    edges.back().to = snp.pos - curr_pos + 1;
                } else {
                    assert(false);
                }
            }
            
            if(!isReverseDeterministic()) {
                reverseDeterminize(curr_pos > 0 ? curr_pos + 1 : 0);
                assert(isReverseDeterministic());
            }

            // Identify head
            index_t head_node = nodes.size();
            for(index_t i = 0; i < nodes.size(); i++) {
                if(nodes[i].label == 'Z') {
                    head_node = i;
                    break;
                }
            }
            assert_lt(head_node, nodes.size());
            index_t tail_node = lastNode; assert_lt(tail_node, nodes.size());
            
            // Update edges
            const index_t invalid = std::numeric_limits<index_t>::max();
            bool head_off = curr_pos > 0, tail_off = curr_pos + curr_len < jlen;
            for(index_t i = 0; i < edges.size(); i++) {
                index_t from = edges[i].from;
                from = from + num_nodes;
                if(head_off && edges[i].from > head_node) from -= 1;
                if(tail_off && edges[i].from > tail_node) from -= 1;
                if(head_off && edges[i].from == head_node) {
                    edges[i].from = invalid;
                } else {
                    edges[i].from = from;
                }
                
                index_t to = edges[i].to;
                to = to + num_nodes;
                if(head_off && edges[i].to > head_node) to -= 1;
                if(tail_off && edges[i].to > tail_node) to -= 1;
                if(tail_off && edges[i].to == tail_node) {
                    edges[i].to = invalid;
                } else {
                    edges[i].to = to;
                }
            }
            head_node = tail_node = invalid;
            // Also update lastNode
            if(!tail_off) {
                lastNode += num_nodes;
                if(head_off) lastNode -= 1;
            }

            // Connect head nodes with tail nodes in the previous automaton
            index_t num_head_nodes = 0;
            index_t tmp_num_edges = edges.size();
            if(head_off) {
                assert_gt(prev_tail_nodes.size(), 0);
                for(index_t i = 0; i < tmp_num_edges; i++) {
                    if(edges[i].from == head_node) {
                        num_head_nodes++;
                        for(index_t j = 0; j < prev_tail_nodes.size(); j++) {
                            edges.expand();
                            edges.back().from = prev_tail_nodes[j];
                            edges.back().to = edges[i].to;
                            assert_lt(edges.back().from, edges.back().to);
                        }
                    }
                }
            }
            
            // List tail nodes
            prev_tail_nodes.clear();
            if(tail_off) {
                for(index_t i = 0; i < tmp_num_edges; i++) {
                    if(edges[i].to == tail_node) {
                        prev_tail_nodes.push_back(edges[i].from);
                    }
                }
            }
            
            // Write nodes and edges
            index_t tmp_num_nodes = nodes.size();
            assert_gt(tmp_num_nodes, 2);
            if(head_off) tmp_num_nodes--;
            if(tail_off) tmp_num_nodes--;
            writeIndex<index_t>(rg_out_file, tmp_num_nodes, bigEndian);
            ASSERT_ONLY(index_t num_nodes_written = 0);
            for(index_t i = 0; i < nodes.size(); i++) {
                if(head_off && nodes[i].label == 'Z') continue;
                if(tail_off && nodes[i].label == '$') continue;
                nodes[i].write(rg_out_file, bigEndian);
                ASSERT_ONLY(num_nodes_written++);
            }
            assert_eq(tmp_num_nodes, num_nodes_written);
            tmp_num_edges = edges.size();
            assert_gt(tmp_num_edges, num_head_nodes + prev_tail_nodes.size());
            if(head_off) tmp_num_edges -= num_head_nodes;
            if(tail_off) tmp_num_edges -= prev_tail_nodes.size();
            writeIndex<index_t>(rg_out_file, tmp_num_edges, bigEndian);
            ASSERT_ONLY(index_t num_edges_written = 0);
            for(index_t i = 0; i < edges.size(); i++) {
                if(head_off && edges[i].from == head_node) continue;
                if(tail_off && edges[i].to == tail_node) continue;
                edges[i].write(rg_out_file, bigEndian);
                ASSERT_ONLY(num_edges_written++);
            }
            assert_eq(tmp_num_edges, num_edges_written);
            
            // Clear nodes and edges
            nodes.clear(); edges.clear();
            
            curr_pos += curr_len;
            num_nodes += tmp_num_nodes;
            num_edges += tmp_num_edges;
        }
        
        // Close out file handle
        rg_out_file.close();
        
        // Read all the nodes and edges
        ifstream rg_in_file(rg_fname.c_str(), ios::binary);
        if(!rg_in_file.good()) {
            cerr << "Could not open file for reading a reference graph: \"" << rg_fname << "\"" << endl;
            throw 1;
        }
        nodes.resizeExact(num_nodes); nodes.clear();
        edges.resizeExact(num_edges); edges.clear();
        while(!rg_in_file.eof()) {
            index_t tmp_num_nodes = readIndex<index_t>(rg_in_file, bigEndian);
            for(index_t i = 0; i < tmp_num_nodes; i++) {
                nodes.expand();
                nodes.back().read(rg_in_file, bigEndian);
            }
            index_t tmp_num_edges = readIndex<index_t>(rg_in_file, bigEndian);
            for(index_t i = 0; i < tmp_num_edges; i++) {
                edges.expand();
                edges.back().read(rg_in_file, bigEndian);
            }
            
            if(nodes.size() >= num_nodes) {
                assert_eq(nodes.size(), num_nodes);
                assert_eq(edges.size(), num_edges);
                break;
            }
        }
        rg_in_file.close();
        std::remove(rg_fname.c_str());
    } else { // this is memory-consuming, but simple to implement
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
    
    // debugging purposes
#if 0
    if(frag_automaton) {
        cout << "frag automaton" << endl;
    } else {
        cout << "one automaton" << endl;
    }
    EList<Node> tmp_nodes = nodes;
    sort(tmp_nodes.begin(), tmp_nodes.end());
    cout << "num of nodes: " << tmp_nodes.size() << endl;
    for(index_t i = 0; i < tmp_nodes.size(); i++) {
        cout << tmp_nodes[i].label << "\t" << tmp_nodes[i].value << endl;
    }
    
    cout << "num of edges: " << edges.size() << endl;
    sort(edges.begin(), edges.end(), TempEdgeNodeCmp(nodes));
    for(index_t i = 0; i < edges.size(); i++) {
        const Node& from = nodes[edges[i].from];
        const Node& to = nodes[edges[i].to];
        cout << from.label << " " << from.value << " --> " << to.label << " " << to.value << endl;
    }
    
    exit(1);
#endif
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
void RefGraph<index_t>::reverseDeterminize(index_t lastNode_add)
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
            cnodes.back().potential = node.value < (lastNode + lastNode_add);
            
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

                if((cnode.value < (lastNode + lastNode_add)) != (next_node.value < (lastNode + lastNode_add))) {
                    cnode.value = min(cnode.value, next_node.value);
                } else {
                    cnode.value = max(cnode.value, next_node.value);
                }
                cnode.potential = cnode.value < (lastNode + lastNode_add);
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
        index_t                outdegree_;
        
        void setSorted()      { to = from; }
        bool isSorted() const { return to == from; }
        
        index_t value() const { return to; }
        index_t outdegree() const { return outdegree_; }
        
        /*
        bool  isBackbone() const {
            return (key.first & (((index_t)1) << (WORD_BITS - 1)));
        }
        void  setBackbone() {
            key.first |= ((index_t)1) << (WORD_BITS - 1);
        }
         */
        
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
    bool generateEdges(RefGraph<index_t>& parent, index_t ftabChars = 10);
    
    index_t getNumEdges() const { return edges.size(); }
    
    //
    bool nextRow(int& gbwtChar, int& F, int& M, index_t& pos, index_t& firstSeq, index_t& lastSeq, index_t ftabChars = 10) {
        if(report_node_idx >= nodes.size()) return false;
        bool firstOutEdge = false;
        if(report_edge_range.first >= report_edge_range.second) {
            report_edge_range = getEdges(report_node_idx, false /* from? */);
            firstOutEdge = true;
            if(report_node_idx == 0) {
                report_M = pair<index_t, index_t>(0, 0);
            }
        }
        assert_lt(report_edge_range.first, report_edge_range.second);
        assert_lt(report_edge_range.first, edges.size());
        const PathEdge& edge = edges[report_edge_range.first];
        gbwtChar = edge.label;
        if(gbwtChar == 'Z') gbwtChar = '$';
        assert_lt(report_node_idx, nodes.size());
        const PathNode& node = nodes[report_node_idx];
        pos = node.to;
        F = (firstOutEdge ? 1 : 0);
        
        // daehwan - for debugging purposes
        if(report_node_idx == 2) {
            int kk = 0;
            kk += 20;
        }
        
        assert_lt(edge.ranking, nodes.size());
        const PathNode& next_node = nodes[edge.ranking];
        if(gbwtChar != '$' && next_node.key.first != numeric_limits<index_t>::max()) {
            assert_neq(next_node.key.second, numeric_limits<index_t>::max());
            int nt = asc2dna[(int)gbwtChar];
            firstSeq = (nt << (ftabChars - 1)) | next_node.key.first;
            lastSeq = (nt << (ftabChars - 1)) | next_node.key.second;
        } else {
            firstSeq = numeric_limits<index_t>::max();
            lastSeq = numeric_limits<index_t>::max();
        }
        
        report_edge_range.first++;
        if(report_edge_range.first >= report_edge_range.second) {
            report_node_idx++;
        }
        assert_lt(report_M.first, nodes.size());
        M = (report_M.second == 0 ? 1 : 0);
        report_M.second++;
        if(report_M.second >= nodes[report_M.first].outdegree_) {
            report_M.first++;
            report_M.second = 0;
        }
        return true;
    }
    
    //
    index_t nextFLocation() {
        if(report_F_node_idx >= nodes.size()) return std::numeric_limits<index_t>::max();
        pair<index_t, index_t> edge_range = getEdges(report_F_node_idx, false /* from? */);
        report_F_node_idx++;
        assert_lt(edge_range.first, edge_range.second);
        report_F_location += (edge_range.second - edge_range.first);
        return report_F_location;
    }

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
        sort(edges.begin(), edges.end(), PathEdgeFromCmp());
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
    
    // For reporting GBWT char, F, and M values
    index_t                report_node_idx;
    pair<index_t, index_t> report_edge_range;
    pair<index_t, index_t> report_M;
    // For reporting location in F corresponding to 1 bit in M
    index_t                report_F_node_idx;
    index_t                report_F_location;
    
    // Can create an index by using key.second in PathNodes.
    // If the graph is not ready, its status becomes error.
    // Sorting edges by from actually sorts them by (from, to).
    void      sortEdges(bool by_from, bool create_index);
    pair<index_t, index_t> getEdges(index_t node, bool by_from); // Create index first.
    
#ifndef NDEBUG
    bool               debug;
#endif
    
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
    
private:
    // Used to construct FCharTable
    struct FNode {
        index_t id;
        index_t key;
        index_t depth;
    };
    
    pair<index_t, index_t> findEdges(index_t node, bool from) {
        pair<index_t, index_t> range(0, 0);
        assert_gt(edges.size(), 0);
        
        // Find lower bound
        index_t low = 0, high = edges.size() - 1;
        index_t temp;
        while(low < high) {
            index_t mid = low + (high - low) / 2;
            temp = (from ? edges[mid].from : edges[mid].ranking);
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
        temp = (from ? edges[low].from : edges[low].ranking);
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
            temp = (from ? edges[mid].from : edges[mid].ranking);
            if(node == temp) {
                low = mid;
            } else {
                assert_lt(node, temp);
                high = mid - 1;
            }
        }
#ifndef NDEBUG
        temp = (from ? edges[high].from : edges[high].ranking);
        assert_eq(node, temp);
#endif
        range.second = high + 1;
        return range;
    }
};


template <typename index_t>
PathGraph<index_t>::PathGraph(RefGraph<index_t>& base) :
ranks(0), max_label('Z'), temp_nodes(0), generation(0),
status(error), has_stabilized(false),
report_node_idx(0), report_edge_range(pair<index_t, index_t>(0, 0)), report_M(pair<index_t, index_t>(0, 0)),
report_F_node_idx(0), report_F_location(0)
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
status(error), has_stabilized(false),
report_node_idx(0), report_edge_range(pair<index_t, index_t>(0, 0)), report_M(pair<index_t, index_t>(0, 0)),
report_F_node_idx(0), report_F_location(0)
{
    if(previous.status != ok)
        return;
    
#ifndef NDEBUG
    debug = previous.debug;
#endif
    
//    previous.sortByFrom(true);

	//create a list of all unsorted nodes
	index_t hash;

	EList<PathNode> unsorted_nodes;
	for(index_t i = 0; i < previous.nodes.size(); i++) {
		if(!previous.nodes[i].isSorted()) {
			unsorted_nodes.push_back(previous.nodes[i]);
		}
	}
	//create an array with these pathnodes binned by hash
	EList<PathNode> nodes_table [unsorted_nodes.size() + unsorted_nodes.size() / 2];


	for(index_t i = 0; i < unsorted_nodes.size(); i++) {
		hash = (149 * unsorted_nodes[i].to % 10009729) % (unsorted_nodes.size() + unsorted_nodes.size() / 2);
		nodes_table[hash].push_back(unsorted_nodes[i]);
	}

#if 0
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
			if(previous.nodes[i].isSorted()) {
			new_nodes++;
			}
			index_t hash = (149 * previous.nodes[i].from % 10009729) % (unsorted_nodes.size() + unsorted_nodes.size() / 2);
       		for(index_t j = 0; j < nodes_table[hash].size(); j++) {
				if(previous.nodes[i].from == nodes_table[hash].get(j).to) {
            	new_nodes++;
     	   		}
        	}
		}
    }
	cerr << "Trying to allocate space for " << new_nodes << " nodes." << endl;
    nodes.resizeExact(new_nodes);
    nodes.clear();
    
    for(index_t i = 0; i < previous.nodes.size(); i++) {
		if(previous.nodes[i].isSorted()) {
			nodes.push_back(previous.nodes[i]);
		}
		index_t hash = (149 * previous.nodes[i].from % 10009729) % (unsorted_nodes.size() + unsorted_nodes.size() / 2);
        for(index_t j = 0; j < nodes_table[hash].size(); j++) {
			if(previous.nodes[i].from == nodes_table[hash].get(j).to) {
            createPathNode(nodes_table[hash].get(j), previous.nodes[i]);
     	   }
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
bool PathGraph<index_t>::generateEdges(RefGraph<index_t>& base, index_t ftabChars)
{
	return 1;
#if 0
    if(status != sorted)
        return false;
    
    sortByFrom(false);
    base.sortEdgesTo();
    
    // Create edges (from, to) as (from.from, to = rank(to))
    pair<index_t, index_t> pn_range = getNextRange(pair<index_t, index_t>(0, 0));
    pair<index_t, index_t> ge_range = base.getNextEdgeRange(pair<index_t, index_t>(0, 0), false /* from? */);
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
    PathNode* node = nodes.begin(); node->outdegree_ = 1;
    PathEdge* edge = edges.begin();
    while(node != nodes.end() && edge != edges.end()) {
        if(edge->from == node->from) {
            edge->from = node - nodes.begin(); edge++;
            node->outdegree_++;
        } else {
            node->to = base.nodes[node->from].value;
            node++; node->outdegree_ = 0;
        }
    }
    if(node != nodes.end()) {
        node->to = base.nodes[node->from].value;
    }
    
    // Pre-calculate prefix for each node
#if 0
    EList<FNode> tmp_nodes; // a list of (node_id, depth)
    for(index_t i = 0; i < nodes.size(); i++) {
        PathNode& node = nodes[i];
        tmp_nodes.clear();
        tmp_nodes.expand();
        tmp_nodes.back().id = i;
        tmp_nodes.back().key = 0;
        tmp_nodes.back().depth = 0;
        while(tmp_nodes.size() > 0) {
            FNode tmp_node = tmp_nodes.back(); tmp_nodes.pop_back();
            pair<index_t, index_t> edge_range = findEdges(tmp_node.id, true /* from? */);
            for(; edge_range.first < edge_range.second; edge_range.first++) {
                assert_lt(edge_range.first, edges.size());
                const PathEdge& edge = edges[edge_range.first];
                if(edge.label == '$') continue;
                index_t key = (tmp_node.key << 2) | asc2dna[(int)edge.label];
                if(tmp_node.depth + 2 == ftabChars) {
                    if(node.key.first == numeric_limits<index_t>::max()) {
                        assert_eq(node.key.second, numeric_limits<index_t>::max());
                        node.key.first = key;
                        node.key.second = key;
                    } else if(key < node.key.first) {
                        assert_neq(node.key.second, numeric_limits<index_t>::max());
                        node.key.first = key;
                    } else if(key > node.key.second) {
                        assert_neq(node.key.first, numeric_limits<index_t>::max());
                        node.key.second = key;
                    }
                } else {
                    if(node.key.first != numeric_limits<index_t>::max()) {
                        assert_neq(node.key.second, numeric_limits<index_t>::max());
                        index_t shift = ftabChars - 2 - tmp_node.depth;
                        index_t tmp_first = node.key.first >> shift;
                        index_t tmp_second = node.key.second >> shift;
                        if(tmp_first < key && tmp_second > key) continue;
                    }
                    
                    tmp_nodes.expand();
                    tmp_nodes.back().id = edge.ranking;
                    tmp_nodes.back().key = key;
                    tmp_nodes.back().depth = tmp_node.depth + 1;
                }
            }
        }
    }
    
#endif
    
    
    sortEdgesTo(true);
    status = ready;
    nodes.pop_back(); // Remove 'Z' node
    
    return true;
    
    bwt_string.clear();
    F_array.clear();
    M_array.clear();
    bwt_counts.resizeExact(5); bwt_counts.fillZero(); bwt_counts.front() = 1;
    for(index_t node = 0; node < nodes.size(); node++) {
        pair<index_t, index_t> edge_range = getEdges(node, false /* from? */);
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
        for(index_t i = 0; i < nodes[node].outdegree_; i++) {
            M_array.push_back(i == 0 ? 1 : 0);
        }
    }
    assert_gt(bwt_string.size(), 0);
    assert_eq(bwt_string.size(), F_array.size());
    assert_eq(bwt_string.size(), M_array.size());
    
    for(size_t i = 0; i < bwt_counts.size(); i++) {
        if(i > 0) bwt_counts[i] += bwt_counts[i - 1];
    }
 
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
    
        cerr << "i\tBWT\tF\tM" << endl;
        for(index_t i = 0; i < bwt_string.size(); i++) {
            cerr << i << "\t" << bwt_string[i] << "\t"  // BWT char
            << (int)F_array[i] << "\t"             // F bit value
            << (int)M_array[i] << endl;            // M bit value
        }
        
        for(size_t i = 0; i < bwt_counts.size(); i++) {
            cerr << i << "\t" << bwt_counts[i] << endl;
        }
    }
#endif
#endif
    
    // Test searches, based on paper_example
#if 0
    EList<string> queries;  EList<index_t> answers;
#   if 1
#      if 0
    queries.push_back("GACGT"); answers.push_back(9);
    queries.push_back("GATGT"); answers.push_back(9);
    queries.push_back("GACT");  answers.push_back(9);
    queries.push_back("ATGT");  answers.push_back(4);
    queries.push_back("GTAC");  answers.push_back(10);
    queries.push_back("ACTG");  answers.push_back(3);
#      else
    // rs55902548, at 402, ref, alt, unknown alt
    // queries.push_back("GGCAGCTCCCATGGGTACACACTGGGCCCAGAACTGGGATGGAGGATGCA");
    queries.push_back("GGCAGCTCCCATGGGTACACACTGGTCCCAGAACTGGGATGGAGGATGCA");
    // queries.push_back("GGCAGCTCCCATGGGTACACACTGGACCCAGAACTGGGATGGAGGATGCA");
    
    // rs5759268, at 926787, ref, alt, unknown alt
    // queries.push_back("AAATTGCTCAGCCTTGTGCTGTGCACACCTGGTTCTCTTTCCAGTGTTAT");
    // queries.push_back("AAATTGCTCAGCCTTGTGCTGTGCATACCTGGTTCTCTTTCCAGTGTTAT");
    // queries.push_back("AAATTGCTCAGCCTTGTGCTGTGCAGACCTGGTTCTCTTTCCAGTGTTAT");
#      endif
    
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
        // assert_eq(top, answers[q]);
        cerr << "finished... ";
        if(top < bot && top < nodes.size()) {
            index_t pos = nodes[top].to;
            index_t gpos = pos;
            const EList<RefRecord>& szs = base.szs;
            for(index_t i = 0; i < szs.size(); i++) {
                gpos += szs[i].off;
                if(pos < szs[i].len) break;
                pos -= szs[i].len;
            }
            
            cerr << "being aligned at " << gpos;
        }
        cerr << endl << endl;
    }
#   endif
    
    // See inconsistencies between F and M arrays
#   if 0
    cerr << endl << endl;
    EList<index_t> tmp_F;
    for(index_t i = 0; i < F_array.size(); i++) {
        if(F_array[i] == 1) tmp_F.push_back(i);
    }
    
    EList<index_t> tmp_M;
    for(index_t i = 0; i < M_array.size(); i++) {
        if(M_array[i] == 1) tmp_M.push_back(i);
    }
    
    index_t max_diff = 0;
    assert_eq(tmp_F.size(), tmp_M.size());
    for(index_t i = 0; i < tmp_F.size(); i++) {
        index_t diff = (tmp_F[i] >= tmp_M[i] ? tmp_F[i] - tmp_M[i] : tmp_M[i] - tmp_F[i]);
        if(diff > max_diff) {
            max_diff = diff;
            cerr << i << "\tdiff: " << max_diff << "\t" << (tmp_F[i] >= tmp_M[i] ? "+" : "-") << endl;
        }
    }
    cerr << "Final: " << tmp_F.back() << " vs. " << tmp_M.back() << endl;
#   endif
    
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
