/*
 * Copyright 2015, Daehwan Kim <infphilo@gmail.com>, Joe Paggi <jpaggi@mit.edu>
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
#include <time.h>
#include "alt.h"
#include "radix_sort.h"

// Reference:
// Jouni Sirén, Niko Välimäki, and Veli Mäkinen: Indexing Graphs for Path Queries with Applications in Genome Research.
// IEEE/ACM Transactions on Computational Biology and Bioinformatics 11(2):375-388, 2014.
// http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6698337

//--------------------------------------------------------------------------

struct NongraphException : public exception
{
    const char* what () const throw ()
    {
        return "Nongraph exception";
    }
};

struct ExplosionException : public exception
{
    const char* what () const throw ()
    {
        return "Explosion exception";
    }
};

template <typename index_t> class PathGraph;

// Note: I wrote the following codes based on Siren's work, gcsa (see the reference above).
template <typename index_t>
class RefGraph {
    friend class PathGraph<index_t>;
public:
    struct Node {
        char    label; // ACGTN + Y(head) + Z(tail)
        index_t value; // location in a whole genome
        
        Node() { reset(); }
        Node(char label_, index_t value_) : label(label_), value(value_) {}
        void reset() { label = 0; value = 0; }

        bool write(ofstream& f_out, bool bigEndian) const {
            writeIndex<index_t>(f_out, value, bigEndian);
            writeU16(f_out, label, bigEndian);
            return true;
        }

        bool read(ifstream& f_in, bool bigEndian) {
            value = readIndex<index_t>(f_in, bigEndian);
            label = (char)readU16(f_in, bigEndian);
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

    static index_t EdgeTo (Edge& a) {
        return a.to;
    }

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
    RefGraph(const SString<char>& s,
             const EList<RefRecord>& szs,
             const EList<ALT<index_t> >& alts,
             const EList<Haplotype<index_t> >& haplotypes,
             const string& out_fname,
             int nthreads_,
             bool verbose);

    bool repOk() { return true; }

    void write(const string& fname, bool bigEndian) {
        ofstream rg_file(fname.c_str(), ios::binary);
        if(!rg_file.good()) {
            cerr << "Could not open file for writing a reference graph: \"" << fname << "\"" << endl;
            throw 1;
        }
        writeIndex<index_t>(rg_file, (index_t)nodes.size(), bigEndian);
        for(index_t i = 0; i < nodes.size(); i++) {
            nodes[i].write(rg_file, bigEndian);
        }
        writeIndex<index_t>(rg_file, (index_t)edges.size(), bigEndian);
        for(index_t i = 0; i < edges.size(); i++) {
            edges[i].write(rg_file, bigEndian);
        }
        rg_file.close();
    }
    
    void nullify() {
        nodes.nullify();
        edges.nullify();
    }
    
    void read(const string& fname, bool bigEndian) {
        ifstream rg_file(fname.c_str(), ios::binary);
        if(!rg_file.good()) {
            cerr << "Could not open file for reading a reference graph: \"" << fname << "\"" << endl;
            throw 1;
        }
        index_t num_nodes = readIndex<index_t>(rg_file, bigEndian);
        nodes.resizeNoCopyExact(num_nodes);
        for(index_t i = 0; i < num_nodes; i++) {
            nodes[i].read(rg_file, bigEndian);
        }
        index_t num_edges = readIndex<index_t>(rg_file, bigEndian);
        edges.resizeNoCopyExact(num_edges);
        for(index_t i = 0; i < num_edges; i++) {
            edges[i].read(rg_file, bigEndian);
        }
        rg_file.close();
    }

private:
    static bool isReverseDeterministic(EList<Node>& nodes, EList<Edge>& edges);
    static void reverseDeterminize(EList<Node>& nodes, EList<Edge>& edges, index_t& lastNode, index_t lastNode_add = 0);

    static void sortEdgesFrom(EList<Edge>& edges) {
        std::sort(edges.begin(), edges.end(), EdgeFromCmp());
    }
    static void sortEdgesTo(EList<Edge>& edges) {
        std::sort(edges.begin(), edges.end(), EdgeToCmp());
    }

    // Return edge ranges [begin, end)
    static pair<index_t, index_t> findEdges(const EList<Edge>& edges, index_t node, bool from);
    static pair<index_t, index_t> findEdgesFrom(const EList<Edge>& edges, index_t node) {
        return findEdges(edges, node, true);
    }
    static pair<index_t, index_t> findEdgesTo(const EList<Edge>& edges, index_t node) {
        return findEdges(edges, node, false);
    }
    static pair<index_t, index_t> getNextEdgeRange(const EList<Edge>& sep_edges, pair<index_t, index_t> range, bool from) {
        if(range.second >= sep_edges.size()) {
            return pair<index_t, index_t>(0, 0);
        }
        range.first = range.second; range.second++;

        if(from) {
            while(range.second < sep_edges.size() && sep_edges[range.second].from == sep_edges[range.first].from) {
                range.second++;
            }
        } else {
            while(range.second < sep_edges.size() && sep_edges[range.second].to == sep_edges[range.first].to) {
                range.second++;
            }
        }
        return range;
    }

private:
    struct ThreadParam {
        // input
        index_t                           thread_id;
        RefGraph<index_t>*                refGraph;
        const SString<char>*              s;
        const EList<ALT<index_t> >*       alts;
        const EList<Haplotype<index_t> >* haplotypes;
        string                            out_fname;
        bool                              bigEndian;

        // output
        index_t                           num_nodes;
        index_t                           num_edges;
        index_t                           lastNode;
        bool                              multipleHeadNodes;
    };
    static void buildGraph_worker(void* vp);

private:
    EList<RefRecord> szs;
    EList<RefRecord> tmp_szs;

    EList<Node> nodes;
    EList<Edge> edges;
    index_t     lastNode; // Z

    int         nthreads;

#ifndef NDEBUG
    bool        debug;
#endif

private:
    // Following composite nodes and edges are used to reverse-determinize an automaton.
    struct CompositeNodeIDs {
        index_t id;
        EList<index_t> add_ids;
        
        CompositeNodeIDs() {
            id = (index_t)INDEX_MAX;
        }      
        
        bool operator<(const CompositeNodeIDs& o) const {
            if(id != o.id) return id < o.id;
            if(add_ids.size() != o.add_ids.size()) return add_ids.size() < o.add_ids.size();
            for(index_t i = 0; i < add_ids.size(); i++) {
                assert_lt(i, o.add_ids.size());
                if(add_ids[i] != o.add_ids[i]) return add_ids[i] < o.add_ids[i];
            }
            return false;
        }
        
        index_t size() const {
            if(id == (index_t)INDEX_MAX) return 0;
            return (index_t)add_ids.size() + 1;
        }
        index_t getID(index_t i) const {
            if(i == 0) return id;
            else {
                i -= 1;
                assert_lt(i, add_ids.size());
                return add_ids[i];
            }
        }
        void push_back(index_t node_id) {
            if(id == (index_t)INDEX_MAX) id = node_id;
            else add_ids.push_back(node_id);
        }
    };
    
    struct CompositeNode {
        CompositeNodeIDs           nodes;
        index_t                    id;
        char                       label;
        index_t                    value;

        CompositeNode() { reset(); }

        CompositeNode(char label_, index_t node_id) :
        id(0), label(label_)
        {
            nodes.push_back(node_id);
        }

        Node getNode() const {
            return Node(label, value);
        }

        void reset() {
            nodes.id = (index_t)INDEX_MAX;
            nodes.add_ids.clear();
            id = 0;
            label = 0;
            value = 0;
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
 * Load reference sequence file and alt information.
 * Construct a reference graph
 */
template <typename index_t>
RefGraph<index_t>::RefGraph(const SString<char>& s,
                            const EList<RefRecord>& szs,
                            const EList<ALT<index_t> >& alts,
                            const EList<Haplotype<index_t> >& haplotypes,
                            const string& out_fname,
                            int nthreads_,
                            bool verbose)
: lastNode(0), nthreads(nthreads_)
{
    const bool bigEndian = false;

    assert_gt(nthreads, 0);
    assert_gt(szs.size(), 0);
    index_t jlen = (index_t)s.length();

#ifndef NDEBUG
    debug = (jlen <= 20);
#endif

    // a memory-efficient way to create a population graph with known ALTs
    bool frag_automaton = jlen >= (1 << 16);
    if(frag_automaton) {
        {
            EList<pair<index_t, index_t> > alt_ranges; // each range inclusive
            for(index_t i = 0; i < alts.size(); i++) {
                const ALT<index_t>& alt = alts[i];
                index_t left_relax = 128, right_relax = 128;
                pair<index_t, index_t> range;
                range.first = alt.pos > left_relax ? alt.pos - left_relax - 1 : 0;
                if(alt.type == ALT_SNP_SGL) {
                    range.second = alt.pos + 1;
                } else if(alt.type == ALT_SNP_DEL) {
                    assert_gt(alt.len, 0);
                    range.second = alt.pos + alt.len;
                } else if(alt.type == ALT_SNP_INS) {
                    assert_gt(alt.len, 0);
                    range.second = alt.pos;
                } else if (alt.type == ALT_SPLICESITE) {
                    assert_lt(alt.left, alt.right);
                    range.second = alt.right + 1;
                } else {
                    assert(alt.exon());
                    continue;
                }
                range.second += right_relax;

                if(alt_ranges.empty() || alt_ranges.back().second + 1 < range.first) {
                    alt_ranges.push_back(range);
                } else {
                    assert_leq(alt_ranges.back().first, range.first);
                    if(alt_ranges.back().second < range.second) {
                        alt_ranges.back().second = range.second;
                    }
                }
            }

            index_t chunk_size = 1 << 20;
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
                            for(; range_idx < alt_ranges.size(); range_idx++) {
                                if(target_pos < alt_ranges[range_idx].first) break;
                            }
                            pair<index_t, index_t> alt_free_range;
                            if(range_idx == 0) {
                                alt_free_range.first = 0;
                            } else {
                                alt_free_range.first = alt_ranges[range_idx - 1].second + 1;
                                if(alt_free_range.first >= jlen) {
                                    alt_free_range.first = jlen - 1;
                                }
                            }

                            if(range_idx == alt_ranges.size()) {
                                alt_free_range.second = jlen - 1;
                            } else {
                                alt_free_range.second = alt_ranges[range_idx].first - 1;
                            }

                            assert_leq(alt_free_range.first, alt_free_range.second);
                            if(target_pos < alt_free_range.first) target_pos = alt_free_range.first;
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

        if(nthreads > (int)tmp_szs.size()) {
            nthreads = (int)tmp_szs.size();
        }
        assert_gt(nthreads, 0);
        AutoArray<tthread::thread*> threads(nthreads);
        EList<ThreadParam> threadParams;
        for(index_t i = 0; i < (index_t)nthreads; i++) {
            threadParams.expand();
            threadParams.back().thread_id = i;
            threadParams.back().refGraph = this;
            threadParams.back().s = &s;
            threadParams.back().alts = &alts;
            threadParams.back().haplotypes = &haplotypes;
            threadParams.back().out_fname = out_fname;
            threadParams.back().bigEndian = bigEndian;
            threadParams.back().num_nodes = 0;
            threadParams.back().num_edges = 0;
            threadParams.back().lastNode = 0;
            threadParams.back().multipleHeadNodes = false;
            if(nthreads == 1) {
                buildGraph_worker((void*)&threadParams.back());
            } else {
                threads[i] = new tthread::thread(buildGraph_worker, (void*)&threadParams.back());
            }
        }

        if(nthreads > 1) {
            for(index_t i = 0; i < (index_t)nthreads; i++)
                threads[i]->join();
        }

        index_t num_nodes = 0, num_edges = 0;
        for(index_t i = 0; i < threadParams.size(); i++) {
            num_nodes += threadParams[i].num_nodes;
            num_edges += threadParams[i].num_edges;
            // Make room for edges spanning graphs by different threads
            if(i > 0) {
                num_edges += 16;
            }
        }
        nodes.resizeExact(num_nodes); nodes.clear();
        edges.resizeExact(num_edges); edges.clear();

        // Read all the nodes and edges
        EList<index_t> tail_nodes;
        bool multipleHeadNodes = false;
        for(index_t i = 0; i < threadParams.size(); i++) {
            if(threadParams[i].multipleHeadNodes) multipleHeadNodes = true;
            std::ostringstream number; number << i;
            const string rg_fname = out_fname + "." + number.str() + ".rf";
            ifstream rg_in_file(rg_fname.c_str(), ios::binary);
            if(!rg_in_file.good()) {
                cerr << "Could not open file for reading a reference graph: \"" << rg_fname << "\"" << endl;
                throw 1;
            }
            index_t curr_num_nodes = (index_t)nodes.size();
            ASSERT_ONLY(index_t curr_num_edges = (index_t)edges.size());
            ASSERT_ONLY(index_t num_spanning_edges = 0);
            // Read nodes to be connected to last nodes in a previous thread
            if(i > 0) {
                assert_gt(tail_nodes.size(), 0)
                index_t num_head_nodes = readIndex<index_t>(rg_in_file, bigEndian);
                for(index_t j = 0; j < num_head_nodes; j++) {
                    index_t head_node = readIndex<index_t>(rg_in_file, bigEndian);
                    for(index_t k = 0; k < tail_nodes.size(); k++) {
                        edges.expand();
                        edges.back().from = tail_nodes[k];
                        edges.back().to = head_node + curr_num_nodes;
                        ASSERT_ONLY(num_spanning_edges++);
                    }
                }
            }
            while(!rg_in_file.eof()) {
                index_t tmp_num_nodes = readIndex<index_t>(rg_in_file, bigEndian);
                for(index_t j = 0; j < tmp_num_nodes; j++) {
                    nodes.expand();
                    nodes.back().read(rg_in_file, bigEndian);
                }
                index_t tmp_num_edges = readIndex<index_t>(rg_in_file, bigEndian);
                for(index_t j = 0; j < tmp_num_edges; j++) {
                    edges.expand();
                    edges.back().read(rg_in_file, bigEndian);
                    edges.back().from += curr_num_nodes;
                    edges.back().to += curr_num_nodes;
                }

                if(nodes.size() >= curr_num_nodes + threadParams[i].num_nodes) {
                    assert_eq(nodes.size(), curr_num_nodes + threadParams[i].num_nodes);
                    assert_eq(edges.size(), curr_num_edges + num_spanning_edges + threadParams[i].num_edges);
                    // Read last nodes in this thread
                    tail_nodes.clear();
                    if(i + 1 < (index_t)nthreads) {
                        index_t num_tail_nodes = readIndex<index_t>(rg_in_file, bigEndian);
                        for(index_t j = 0; j < num_tail_nodes; j++) {
                            index_t tail_node = readIndex<index_t>(rg_in_file, bigEndian);
                            tail_nodes.push_back(tail_node + curr_num_nodes);
                        }
                    }
                    break;
                }
            }
            rg_in_file.close();
            std::remove(rg_fname.c_str());
            if(i + 1 == (index_t)nthreads) {
                lastNode = threadParams.back().lastNode + curr_num_nodes;
                assert_lt(lastNode, nodes.size());
                assert_eq(nodes[lastNode].label, 'Z');
            }
        }
        
        if(s.length() + 2 == nodes.size() && nodes.size() == edges.size() + 1) {
            cerr << "Warning: no variants or splice sites in this graph" << endl;
            throw NongraphException();
        }

        if(multipleHeadNodes) {
            if(!isReverseDeterministic(nodes, edges)) {
                if(verbose) cerr << "\tis not reverse-deterministic, so reverse-determinize..." << endl;
                reverseDeterminize(nodes, edges, lastNode);
            }
        }
        assert(isReverseDeterministic(nodes, edges));
    } else { // this is memory-consuming, but simple to implement
        index_t num_predicted_nodes = (index_t)(jlen * 1.2);
        nodes.reserveExact(num_predicted_nodes);
        edges.reserveExact(num_predicted_nodes);

        // Created head node
        nodes.expand();
        nodes.back().label = 'Y';
        nodes.back().value = 0;
        // Create nodes and edges corresponding to a reference genome
        for(size_t i = 0; i < s.length(); i++) {
            nodes.expand();
            nodes.back().label = "ACGT"[(int)s[i]];
            nodes.back().value = (index_t)i;

            assert_geq(nodes.size(), 2);
            edges.expand();
            edges.back().from = (index_t)nodes.size() - 2;
            edges.back().to = (index_t)nodes.size() - 1;
        }

        // Create tail node
        nodes.expand();
        nodes.back().label = 'Z';
        nodes.back().value = (index_t)s.length();
        lastNode = (index_t)nodes.size() - 1;
        edges.expand();
        edges.back().from = (index_t)nodes.size() - 2;
        edges.back().to = (index_t)nodes.size() - 1;
        
        // Create nodes and edges for haplotypes
        for(index_t i = 0; i < haplotypes.size(); i++) {
            const Haplotype<index_t>& haplotype = haplotypes[i];
            const EList<index_t, 1>& snpIDs = haplotype.alts;
            assert_gt(snpIDs.size(), 0);
            assert_lt(haplotype.right, s.length());
            bool pass = true;
            for(index_t s = 0; s < snpIDs.size(); s++) {
                index_t snpID = snpIDs[s];
                assert_lt(snpID, alts.size());
                const ALT<index_t>& snp = alts[snpID];
                assert(snp.snp());
                if(s + 1 >= snpIDs.size()) break;
                index_t snpID2 = snpIDs[s+1];
                assert_lt(snpID2, alts.size());
                const ALT<index_t>& snp2 = alts[snpID2];
                assert(snp2.snp());
                if(snp.type == ALT_SNP_INS) {
                    if(snp.pos > snp2.pos) {
                        pass = false;
                        break;
                    }
                } else if(snp.type == ALT_SNP_DEL) {
                    if(snp2.type == ALT_SNP_DEL) {
                        if(snp.pos + snp.len >= snp2.pos) {
                            pass = false;
                            break;
                        }
                    } else {
                        if(snp.pos + snp.len - 1 >= snp2.pos) {
                            pass = false;
                            break;
                        }
                    }
                } else {
                    if(snp.pos >= snp2.pos) {
                        pass = false;
                        break;
                    }
                }
            }
            
            if(!pass) continue;
            
            index_t prev_ALT_type = ALT_NONE;
            index_t ID_i = 0;
            for(index_t j = haplotype.left; j <= haplotype.right; j++) {
                if(prev_ALT_type == ALT_SNP_INS) j--;
                const ALT<index_t>* altp = (ID_i < snpIDs.size() ? &(alts[snpIDs[ID_i]]) : NULL);                
                assert(altp == NULL || altp->pos >= j);
                if(altp != NULL && altp->pos == j) {
                    const ALT<index_t>& alt = *altp;
                    assert_lt(alt.pos, s.length());
                    assert(alt.snp());
                    if(alt.type == ALT_SNP_SGL) {
                        assert_eq(alt.len, 1);
                        nodes.expand();
                        assert_lt(alt.seq, 4);
                        assert_neq(alt.seq & 0x3, s[alt.pos]);
                        nodes.back().label = "ACGT"[alt.seq];
                        nodes.back().value = alt.pos;
                        if(prev_ALT_type != ALT_SNP_DEL) {
                            edges.expand();
                            if(j == haplotype.left) {
                                edges.back().from = alt.pos;
                            } else {
                                assert_gt(nodes.size(), 2);
                                edges.back().from = (index_t)nodes.size() - 2;
                            }
                            edges.back().to = (index_t)nodes.size() - 1;
                        }
                        if(j == haplotype.right) {
                            edges.expand();
                            edges.back().from = (index_t)nodes.size() - 1;
                            edges.back().to = alt.pos + 2;
                        }
                    }
                    else if(alt.type == ALT_SNP_DEL) {
                        assert_gt(alt.len, 0);
                        assert_leq(alt.pos + alt.len, s.length());
                        edges.expand();
                        if(j == haplotype.left) {
                            edges.back().from = alt.pos;
                        } else {
                            edges.back().from = (index_t)nodes.size() - 1;
                        }
                        j += (alt.len - 1);
                        assert_leq(j, haplotype.right);
                        if(j == haplotype.right) {
                            edges.back().to = alt.pos + alt.len + 1;
                        } else {
                            edges.back().to = (index_t)nodes.size();
                        }
                    } else {
                        assert_eq(alt.type, ALT_SNP_INS)
                        assert_gt(alt.len, 0);
                        for(size_t k = 0; k < alt.len; k++) {
                            uint64_t bp = alt.seq >> ((alt.len - k - 1) << 1);
                            bp &= 0x3;
                            char ch = "ACGT"[bp];
                            nodes.expand();
                            nodes.back().label = ch;
                            nodes.back().value = (index_t)INDEX_MAX;
                            if(prev_ALT_type == ALT_SNP_DEL && k == 0) continue;
                            edges.expand();
                            edges.back().from = ((k == 0 && j == haplotype.left) ? alt.pos : (index_t)nodes.size() - 2);
                            edges.back().to = (index_t)nodes.size() - 1;
                        }
                        if(j == haplotype.right) {
                            edges.expand();
                            edges.back().from = (index_t)nodes.size() - 1;
                            edges.back().to = alt.pos + 1;
                        }
                    }
                    ID_i++;
                    prev_ALT_type = alt.type;
                } else {
                    int nt = s[j];
                    assert_lt(nt, 4);
                    nodes.expand();
                    nodes.back().label = "ACGT"[nt];
                    nodes.back().value = j;
                    if(prev_ALT_type != ALT_SNP_DEL) {
                        edges.expand();
                        if(j == haplotype.left && prev_ALT_type == ALT_NONE) {
                            edges.back().from = j;
                        } else {
                            edges.back().from = (index_t)nodes.size() - 2;
                        }
                        edges.back().to = (index_t)nodes.size() - 1;
                    }
                    if(j == haplotype.right) {
                        edges.expand();
                        edges.back().from = (index_t)nodes.size() - 1;
                        edges.back().to = j + 2;
                    }
                    prev_ALT_type = ALT_SNP_SGL;
                }
            }
        }
        
        // Create nodes and edges for splice sites
        for(size_t i = 0; i < alts.size(); i++) {
            const ALT<index_t>& alt = alts[i];
            if(alt.pos >= s.length()) break;
            if(alt.type != ALT_SPLICESITE) continue;
            if(alt.excluded) continue;
            assert_lt(alt.left, alt.right);
            edges.expand();
            edges.back().from = alt.left;
            edges.back().to = alt.right + 2;
        }
        
        if(s.length() + 2 == nodes.size() && nodes.size() == edges.size() + 1) {
            throw NongraphException();
        }

        if(!isReverseDeterministic(nodes, edges)) {
            if(verbose) cerr << "\tis not reverse-deterministic, so reverse-determinize..." << endl;
            reverseDeterminize(nodes, edges, lastNode);
            assert(isReverseDeterministic(nodes, edges));
        }
    }
    
#ifndef NDEBUG
    if(debug) {
        cout << "num nodes: " << nodes.size() << endl;
        for(index_t i = 0; i < nodes.size(); i++) {
            const Node& n = nodes[i];
            cout << i << "\t" << n.label << "\t" << n.value << endl;
        }
        
        sort(edges.begin(), edges.end());
        cout << "num edges: " << edges.size() << endl;
        for(index_t i = 0; i < edges.size(); i++) {
            const Edge& e = edges[i];
            cout << i << "\t" << e.from << " --> " << e.to << endl;
        }
    }
#endif
}

template <typename index_t>
pair<index_t, index_t> RefGraph<index_t>::findEdges(const EList<Edge>& edges, index_t node, bool from) {
    pair<index_t, index_t> range(0, 0);
    assert_gt(edges.size(), 0);

    // Find lower bound
    index_t low = 0, high = (index_t)edges.size() - 1;
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
    high = (index_t)edges.size() - 1;
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
void RefGraph<index_t>::buildGraph_worker(void* vp) {
    ThreadParam* threadParam = (ThreadParam*)vp;
    RefGraph<index_t>& refGraph = *(threadParam->refGraph);

    const SString<char>& s = *(threadParam->s);
    index_t jlen = (index_t)s.length();

    const EList<ALT<index_t> >& alts = *(threadParam->alts);
    const EList<Haplotype<index_t> >& haplotypes = *(threadParam->haplotypes);

    EList<Node> nodes; EList<Edge> edges;
    const EList<RefRecord>& tmp_szs = refGraph.tmp_szs;

    index_t thread_id = threadParam->thread_id;
    index_t nthreads = refGraph.nthreads;
    std::ostringstream number; number << thread_id;
    const string rg_fname = threadParam->out_fname + "." + number.str() + ".rf";
    ofstream rg_out_file(rg_fname.c_str(), ios::binary);
    if(!rg_out_file.good()) {
        cerr << "Could not open file for writing a reference graph: \"" << rg_fname << "\"" << endl;
        throw 1;
    }
    
#ifndef NDEBUG
    set<index_t> snp_set;
#endif

    const bool bigEndian = threadParam->bigEndian;

    index_t& lastNode = threadParam->lastNode;

    index_t& num_nodes = threadParam->num_nodes;
    index_t& num_edges = threadParam->num_edges;
    index_t szs_idx = 0, szs_idx_end = (index_t)tmp_szs.size();
    if(threadParam->thread_id != 0) {
        szs_idx = (index_t)((tmp_szs.size() / nthreads) * thread_id);
    }
    if(thread_id + 1 < nthreads) {
        szs_idx_end = (index_t)((tmp_szs.size() / nthreads) * (thread_id + 1));
    }

    index_t curr_pos = 0;
    for(index_t i = 0; i < szs_idx; i++) {
        curr_pos += tmp_szs[i].len;
    }
    EList<index_t> prev_tail_nodes;
    index_t alt_idx = 0, haplotype_idx = 0;
    for(; szs_idx < szs_idx_end; szs_idx++) {
        index_t curr_len = tmp_szs[szs_idx].len;
        if(curr_len <= 0) continue;
        index_t num_predicted_nodes = (index_t)(curr_len * 1.2);
        nodes.resizeExact(num_predicted_nodes); nodes.clear();
        edges.resizeExact(num_predicted_nodes); edges.clear();

        // Created head node
        nodes.expand();
        nodes.back().label = 'Y';
        nodes.back().value = 0;

        // Create nodes and edges corresponding to a reference genome
        assert_leq(curr_pos + curr_len, s.length());
        for(size_t i = curr_pos; i < curr_pos + curr_len; i++) {
            nodes.expand();
            nodes.back().label = "ACGT"[(int)s[i]];
            nodes.back().value = (index_t)i;
            assert_geq(nodes.size(), 2);
            edges.expand();
            edges.back().from = (index_t)nodes.size() - 2;
            edges.back().to = (index_t)nodes.size() - 1;
        }

        // Create tail node
        nodes.expand();
        nodes.back().label = 'Z';
        nodes.back().value = (index_t)s.length();
        lastNode = (index_t)nodes.size() - 1;
        edges.expand();
        edges.back().from = (index_t)nodes.size() - 2;
        edges.back().to = (index_t)nodes.size() - 1;
        ASSERT_ONLY(index_t backbone_nodes = (index_t)nodes.size());
     
        // Create nodes and edges for haplotypes
        for(; haplotype_idx < haplotypes.size(); haplotype_idx++) {
            const Haplotype<index_t>& haplotype = haplotypes[haplotype_idx];
            if(haplotype.left < curr_pos) continue;
            if(haplotype.right >= curr_pos + curr_len) break;
            const EList<index_t, 1>& snpIDs = haplotype.alts;
            assert_gt(snpIDs.size(), 0);
            bool pass = true;
            for(index_t s = 0; s < snpIDs.size(); s++) {
                index_t snpID = snpIDs[s];
                assert_lt(snpID, alts.size());
                const ALT<index_t>& snp = alts[snpID];
                assert(snp.snp());
                if(s + 1 >= snpIDs.size()) break;
                index_t snpID2 = snpIDs[s+1];
                assert_lt(snpID2, alts.size());
                const ALT<index_t>& snp2 = alts[snpID2];
                assert(snp2.snp());
                if(snp.type == ALT_SNP_INS) {
                    if(snp.pos > snp2.pos) {
                        pass = false;
                        break;
                    }
                } else if(snp.type == ALT_SNP_DEL) {
                    if(snp2.type == ALT_SNP_DEL) {
                        if(snp.pos + snp.len >= snp2.pos) {
                            pass = false;
                            break;
                        }
                    } else {
                        if(snp.pos + snp.len - 1 >= snp2.pos) {
                            pass = false;
                            break;
                        }
                    }
                } else {
                    if(snp.pos >= snp2.pos) {
                        pass = false;
                        break;
                    }
                }
            }

            if(!pass) continue;
            
            index_t prev_ALT_type = ALT_NONE;
            index_t ID_i = 0;
            for(index_t j = haplotype.left; j <= haplotype.right; j++) {
                if(prev_ALT_type == ALT_SNP_INS) j--;
                const ALT<index_t>* altp = (ID_i < snpIDs.size() ? &(alts[snpIDs[ID_i]]) : NULL);
                assert(altp == NULL || altp->pos >= j);
                if(altp != NULL && altp->pos == j) {
                    const ALT<index_t>& alt = *altp;
                    assert_lt(alt.pos, s.length());
                    assert(alt.snp());
                    if(alt.type == ALT_SNP_SGL) {
                        assert_eq(alt.len, 1);
                        nodes.expand();
                        assert_lt(alt.seq, 4);
                        assert_neq(alt.seq & 0x3, s[alt.pos]);
                        nodes.back().label = "ACGT"[alt.seq];
                        nodes.back().value = alt.pos;
                        if(prev_ALT_type != ALT_SNP_DEL) {
                            edges.expand();
                            if(j == haplotype.left) {
                                edges.back().from = alt.pos - curr_pos;
                                assert_lt(edges.back().from, backbone_nodes);
                            } else {
                                assert_gt(nodes.size(), 2);
                                edges.back().from = (index_t)nodes.size() - 2;
                            }                            
                            edges.back().to = (index_t)nodes.size() - 1;
                        }
                        if(j == haplotype.right) {
                            edges.expand();
                            edges.back().from = (index_t)nodes.size() - 1;
                            edges.back().to = alt.pos - curr_pos + 2;
                            assert_lt(edges.back().to, backbone_nodes);
                        }
                    }
                    else if(alt.type == ALT_SNP_DEL) {
                        assert_gt(alt.len, 0);
                        assert_leq(alt.pos - curr_pos + alt.len, s.length());
                        edges.expand();
                        if(j == haplotype.left) {
                            edges.back().from = alt.pos - curr_pos;
                            assert_lt(edges.back().from, backbone_nodes);
                        } else {
                            edges.back().from = (index_t)nodes.size() - 1;
                        }
                        j += (alt.len - 1);
                        assert_leq(j, haplotype.right);
                        if(j == haplotype.right) {
                            edges.back().to = alt.pos - curr_pos + alt.len + 1;
                            assert_lt(edges.back().to, backbone_nodes);
                        } else {
                            edges.back().to = (index_t)nodes.size();
                        }
                    } else {
                        assert_eq(alt.type, ALT_SNP_INS)
                        assert_gt(alt.len, 0);
                        for(size_t k = 0; k < alt.len; k++) {
                            uint64_t bp = alt.seq >> ((alt.len - k - 1) << 1);
                            bp &= 0x3;
                            char ch = "ACGT"[bp];
                            nodes.expand();
                            nodes.back().label = ch;
                            nodes.back().value = (index_t)INDEX_MAX;
                            if(prev_ALT_type == ALT_SNP_DEL && k == 0) continue;
                            edges.expand();
                            edges.back().from = ((k == 0 && j == haplotype.left) ? alt.pos - curr_pos : (index_t)nodes.size() - 2);
                            edges.back().to = (index_t)nodes.size() - 1;
                        }
                        if(j == haplotype.right) {
                            edges.expand();
                            edges.back().from = (index_t)nodes.size() - 1;
                            edges.back().to = alt.pos - curr_pos + 1;
                        }
                    }
#ifndef NDEBUG
                    snp_set.insert(snpIDs[ID_i]);
#endif
                    ID_i++;
                    prev_ALT_type = alt.type;
                } else {
                    int nt = s[j];
                    assert_lt(nt, 4);
                    nodes.expand();
                    nodes.back().label = "ACGT"[nt];
                    nodes.back().value = j;
                    if(prev_ALT_type != ALT_SNP_DEL) {
                        edges.expand();
                        if(j == haplotype.left && prev_ALT_type == ALT_NONE) {
                            edges.back().from = j - curr_pos;
                            assert_lt(edges.back().from, backbone_nodes);
                        } else {
                            edges.back().from = (index_t)nodes.size() - 2;
                        }
                        edges.back().to = (index_t)nodes.size() - 1;
                    }
                    if(j == haplotype.right) {
                        edges.expand();
                        edges.back().from = (index_t)nodes.size() - 1;
                        edges.back().to = j - curr_pos + 2;
                        assert_lt(edges.back().to, backbone_nodes);
                    }
                    prev_ALT_type = ALT_SNP_SGL;
                }
            }
        }
        
        // Create nodes and edges for splice sites
        for(; alt_idx < alts.size(); alt_idx++) {
            const ALT<index_t>& alt = alts[alt_idx];
            if(alt.pos < curr_pos) continue;
            if(alt.pos >= curr_pos + curr_len) break;
            if(!alt.splicesite()) continue;
            if(alt.excluded) continue;
            assert_lt(alt.left, alt.right);
            edges.expand();
            edges.back().from = alt.left - curr_pos;
            edges.back().to = alt.right - curr_pos + 2;
            assert_lt(edges.back().from, backbone_nodes);
            assert_lt(edges.back().to, backbone_nodes);
        }

#ifndef NDEBUG
        if(refGraph.debug) {
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

        if(!isReverseDeterministic(nodes, edges)) {
            reverseDeterminize(nodes, edges, lastNode, curr_pos > 0 ? curr_pos + 1 : 0);
            assert(isReverseDeterministic(nodes, edges));
        }

        // Identify head
        index_t head_node = (index_t)nodes.size();
        for(index_t i = 0; i < nodes.size(); i++) {
            if(nodes[i].label == 'Y') {
                head_node = i;
                break;
            }
        }
        assert_lt(head_node, nodes.size());
        index_t tail_node = lastNode; assert_lt(tail_node, nodes.size());

        // Update edges
        const index_t invalid = (index_t)INDEX_MAX;
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
        index_t tmp_num_edges = (index_t)edges.size();
        if(head_off) {
            EList<index_t> nodes_to_head;
            for(index_t i = 0; i < tmp_num_edges; i++) {
                if(edges[i].from == head_node) {
                    num_head_nodes++;
                    if(prev_tail_nodes.size() > 0) {
                        for(index_t j = 0; j < prev_tail_nodes.size(); j++) {
                            edges.expand();
                            edges.back().from = prev_tail_nodes[j];
                            edges.back().to = edges[i].to;
                            assert_lt(edges.back().from, edges.back().to);
                        }
                    } else {
                        nodes_to_head.push_back(edges[i].to);
                    }
                }
            }

            if(nodes_to_head.size() > 0) {
                assert_gt(thread_id, 0);
                assert_eq(prev_tail_nodes.size(), 0);
                writeIndex<index_t>(rg_out_file, (index_t)nodes_to_head.size(), bigEndian);
                for(index_t i = 0; i < nodes_to_head.size(); i++) {
                    writeIndex<index_t>(rg_out_file, nodes_to_head[i], bigEndian);
                }
            }
        }

        // Need to check if it's reverse-deterministic
        if(num_head_nodes > 1) {
            threadParam->multipleHeadNodes = true;
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
        index_t tmp_num_nodes = (index_t)nodes.size();
        assert_gt(tmp_num_nodes, 2);
        if(head_off) tmp_num_nodes--;
        if(tail_off) tmp_num_nodes--;
        writeIndex<index_t>(rg_out_file, tmp_num_nodes, bigEndian);
        ASSERT_ONLY(index_t num_nodes_written = 0);
        for(index_t i = 0; i < nodes.size(); i++) {
            if(head_off && nodes[i].label == 'Y') continue;
            if(tail_off && nodes[i].label == 'Z') continue;
            nodes[i].write(rg_out_file, bigEndian);
            ASSERT_ONLY(num_nodes_written++);
        }
        assert_eq(tmp_num_nodes, num_nodes_written);
        tmp_num_edges = (index_t)edges.size();
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

    if(nthreads > 1 && thread_id + 1 < (index_t)nthreads && prev_tail_nodes.size() > 0) {
        writeIndex<index_t>(rg_out_file, (index_t)prev_tail_nodes.size(), bigEndian);
        for(index_t i = 0; i < prev_tail_nodes.size(); i++) {
            writeIndex<index_t>(rg_out_file, prev_tail_nodes[i], bigEndian);
        }
    }

    // Close out file handle
    rg_out_file.close();
}

template <typename index_t>
bool RefGraph<index_t>::isReverseDeterministic(EList<Node>& nodes, EList<Edge>& edges)
{
    if(edges.size() <= 0) return true;

    // Sort edges by "to" nodes
    sortEdgesTo(edges);

    index_t curr_to = (index_t)INDEX_MAX;
    EList<bool> seen; seen.resize(5); seen.fillZero();
    for(index_t i = 0; i < edges.size(); i++) {
        index_t from = edges[i].from;
        assert_lt(from, nodes.size());
        char nt = nodes[from].label;
        assert(nt == 'A' || nt == 'C' || nt == 'G' || nt == 'T' || nt == 'Y');
        if(nt == 'Y') nt = 4;
        else          nt = asc2dna[(int)nt];
        assert_lt(nt, seen.size());
        if(curr_to != edges[i].to) {
            curr_to = edges[i].to;
            seen.fillZero();
            seen[nt] = true;
        } else {
            if(seen[nt]) {
                return false;
            }
            seen[nt] = true;
        }
    }

    return true;
}


template <typename index_t>
void RefGraph<index_t>::reverseDeterminize(EList<Node>& nodes, EList<Edge>& edges, index_t& lastNode, index_t lastNode_add)
{
    EList<CompositeNode> cnodes; cnodes.ensure(nodes.size());
    map<CompositeNodeIDs, index_t> cnode_map;
    deque<index_t> active_cnodes;
    EList<CompositeEdge> cedges; cedges.ensure(edges.size());

    // Start from the final node ('Z')
    assert_lt(lastNode, nodes.size());
    const Node& last_node = nodes[lastNode];
    cnodes.expand();
    cnodes.back().reset();
    cnodes.back().label = last_node.label;
    cnodes.back().value = last_node.value;
    cnodes.back().nodes.push_back(lastNode);
    active_cnodes.push_back(0);
    cnode_map[cnodes.back().nodes] = 0;
    sortEdgesTo(edges);

    index_t firstNode = 0; // Y -> ... -> Z
    EList<index_t> predecessors;
    while(!active_cnodes.empty()) {
        index_t cnode_id = active_cnodes.front(); active_cnodes.pop_front();
        assert_lt(cnode_id, cnodes.size());

        // Find predecessors of this composite node
        predecessors.clear();
        for(size_t i = 0; i < cnodes[cnode_id].nodes.size(); i++) {
            index_t node_id = cnodes[cnode_id].nodes.getID((index_t)i);
            pair<index_t, index_t> edge_range = findEdgesTo(edges, node_id);
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
            index_t new_size = (index_t)(unique(predecessors.begin(), predecessors.end()) - predecessors.begin());
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
            
            if(node.label == 'Y' && firstNode == 0) {
                firstNode = (index_t)cnodes.size() - 1;
            }

            while(i < predecessors.size()) {
                index_t next_node_id = predecessors[i];
                assert_lt(next_node_id, nodes.size());
                const Node& next_node = nodes[next_node_id];
                if(next_node.label != node.label) break;
                cnodes.back().nodes.push_back(next_node_id);
                if(next_node.value != (index_t)INDEX_MAX) {
                    if(cnodes.back().value == (index_t)INDEX_MAX) {
                        cnodes.back().value = next_node.value;
                    } else {
                        cnodes.back().value = max(cnodes.back().value, next_node.value);
                    }
                }
                i++;
            }

            // Create edges from this new composite node to current composite node
            typename map<CompositeNodeIDs, index_t>::iterator existing = cnode_map.find(cnodes.back().nodes);
            if(existing == cnode_map.end()) {
                cnode_map[cnodes.back().nodes] = (index_t)cnodes.size() - 1;
                active_cnodes.push_back((index_t)cnodes.size() - 1);
                cedges.push_back(CompositeEdge((index_t)cnodes.size() - 1, cnode_id));
            } else {
                cnodes.pop_back();
                cedges.push_back(CompositeEdge((*existing).second, cnode_id));
            }

            // Increment indegree
            cnodes[cnode_id].id++;
        }
    }

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
        index_t i = (index_t)cedges.bsearchLoBound(CompositeEdge(cnode_id, 0));
        while(i < cedges.size()) {
            assert_geq(cedges[i].from, cnode_id);
            if(cedges[i].from != cnode_id) break;
            index_t predecessor_cnode_id = cedges[i].to;
            assert_lt(predecessor_cnode_id, cnodes.size());
            CompositeNode& predecessor_cnode = cnodes[predecessor_cnode_id];
            if(cnode.value == predecessor_cnode.value + 1) {
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
        index_t i = (index_t)cedges.bsearchLoBound(CompositeEdge(cnode_id, 0));
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
                successor_cnode.id = (index_t)nodes.size();
                nodes.expand();
                nodes.back() = successor_cnode.getNode();
                if(nodes.back().label == 'Z') {
                    assert_eq(lastNode, 0);
                    assert_gt(nodes.size(), 1);
                    lastNode = (index_t)nodes.size() - 1;
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
    sortEdgesFrom(edges);

#if 0
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
#endif
}

template <typename index_t>
class PathGraph {
public:
    struct PathNode {
        index_t                from;
        index_t                to;
        pair<index_t, index_t> key;

        void setSorted()      { to = (index_t)INDEX_MAX; }
        bool isSorted() const { return to == (index_t)INDEX_MAX; }

        index_t value() const { return to; }
        index_t outdegree() const { return key.first; }

        bool operator< (const PathNode& o) const {
            return key < o.key;
        };
    };

    static inline index_t PathNodeFrom (PathNode& a) {
        return a.from;
    }

    static inline index_t PathNodeKey (PathNode& a) {
        return a.key.first;
    }

    struct PathNodeKeySecondCmp {
        bool operator() (const PathNode& a, const PathNode& b) const {
            return a.key.second < b.key.second;
        }
    };

    struct PathNodeFromCmp {
        bool operator() (const PathNode& a, const PathNode& b) const {
            return a.from < b.from;
        }
    };

    struct PathNodeToCmp {
        bool operator() (const PathNode& a, const PathNode& b) const {
            return a.to < b.to;
        }
    };

    struct PathEdge {
        index_t from;
        union {
            index_t to;
            index_t ranking;
        };
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
    
    static inline index_t PathEdgeTo (PathEdge& a) {
        return a.to;
    }

    struct PathEdgeFromCmp {
        bool operator() (const PathEdge& a, const PathEdge& b) const {
            return a.from < b.from || (a.from == b.from && a.to < b.to);
        }
    };

    struct PathEdgeToCmp {
        bool operator() (const PathEdge& a, const PathEdge& b) const {
            return a.to < b.to || (a.to == b.to && a.from < b.from);
        }
    };

public:
    // Create a new graph in which paths are represented using nodes
    PathGraph(
              RefGraph<index_t>& parent,
              const string& base_fname,
              size_t max_num_nodes_ = std::numeric_limits<size_t>::max(),
              int nthreads_ = 1,
              bool verbose_ = false);

    ~PathGraph() {}

    void printInfo();

    bool generateEdges(RefGraph<index_t>& parent);

    index_t getNumNodes() const { return (index_t)nodes.size(); }
    index_t getNumEdges() const { return (index_t)edges.size(); }
    
    bool isSorted() const { return sorted; }

    bool nextRow(int& gbwtChar, int& F, int& M, index_t& pos) {
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
        assert_lt(report_node_idx, nodes.size());
        F = (firstOutEdge ? 1 : 0);

        report_edge_range.first++;
        if(report_edge_range.first >= report_edge_range.second) {
            report_node_idx++;
        }
        assert_lt(report_M.first, nodes.size());
        pos = nodes[report_M.first].to;
        M = (report_M.second == 0 ? 1 : 0);
        report_M.second++;
        if(report_M.second >= nodes[report_M.first].key.first) {
            report_M.first++;
            report_M.second = 0;
        }
        return true;
    }

    index_t nextFLocation() {
        if(report_F_node_idx >= nodes.size()) return (index_t)INDEX_MAX;
        index_t ret = report_F_location;
        pair<index_t, index_t> edge_range = getEdges(report_F_node_idx, false /* from? */);
        report_F_node_idx++;
        assert_lt(edge_range.first, edge_range.second);
        report_F_location += (edge_range.second - edge_range.first);
        return ret;
    }

private:
    void makeFromRef(RefGraph<index_t>& base);
    void generationOne();
    void earlyGeneration();
    void firstPruneGeneration();
    void lateGeneration();

    void mergeUpdateRank();
    pair<index_t, index_t> nextMaximalSet(pair<index_t, index_t> range);
    pair<index_t, index_t> getEdges(index_t node, bool by_from); // Create index first.

    struct CreateNewNodesParams {
        PathNode*           st;
        PathNode*           en;
        PathNode*           curr;
        index_t*            sub_temp_nodes;
        PathGraph<index_t>* graph;
    };
    static void createNewNodesCounter(void* vp);
    static void createNewNodesMaker(void* vp);
    void createNewNodes();

    struct GenEdgesParams {
        typename RefGraph<index_t>::Edge*        st;
        typename RefGraph<index_t>::Edge*        en;
        EList<index_t, 6>*                       label_index;
        EList<PathNode>*                         nodes;
        EList<PathEdge>*                         edges;
        EList<typename RefGraph<index_t>::Node>* ref_nodes;
    };
    static void generateEdgesCounter(void* vp);
    static void generateEdgesMaker(void* vp);

private:
    int             nthreads;
    bool            verbose;
    EList<PathNode> from_table;
    EList<PathNode> past_nodes;
    EList<PathNode> nodes;
    EList<PathEdge> edges;
    index_t         ranks;
    index_t         max_from; //number of nodes in RefGraph
    index_t         temp_nodes; // Total number of nodes created before sorting.
    index_t         generation; // Sorted by paths of length 2^generation.
    bool            sorted;
    
    // For reporting GBWT char, F, and M values
    index_t                report_node_idx;
    pair<index_t, index_t> report_edge_range;
    pair<index_t, index_t> report_M;
    // For reporting location in F corresponding to 1 bit in M
    index_t                report_F_node_idx;
    index_t                report_F_location;
    
    size_t          max_num_nodes;

    // following variables are for debugging purposes
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
        return (index_t)array.size();
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

    // for debugging purposes
#ifndef NDEBUG
public: EList<pair<index_t, index_t> > ftab;
#endif
};

//creates prefix-sorted PathGraph Nodes given a reverse determinized RefGraph
//outputs nodes sorted by their from attribute
template <typename index_t>
PathGraph<index_t>::PathGraph(
                              RefGraph<index_t>& base,
                              const string& base_fname,
                              size_t max_num_nodes_,
                              int nthreads_,
                              bool verbose_) :
nthreads(nthreads_), verbose(verbose_),
ranks(0), temp_nodes(0), generation(0), sorted(false),
report_node_idx(0), report_edge_range(pair<index_t, index_t>(0, 0)), report_M(pair<index_t, index_t>(0, 0)),
report_F_node_idx(0), report_F_location(0),
max_num_nodes(max_num_nodes_)
{
#ifndef NDEBUG
    debug = base.nodes.size() <= 20;
#endif
    // Fill nodes with a PathNode for each edge in base.edges.
    // Set max_from.
    makeFromRef(base);
    
    // Write RefGraph into a file
    const bool file_rf = base.nodes.size() > (1 << 22);
    const bool bigEndian = false;
    const string rf_fname = base_fname + ".rf";
    if(file_rf) {
        base.write(rf_fname, bigEndian);
        base.nullify();
    }
    
    // In the first generation the nodes enter, not quite sorted by from.
    // We use a counting sort to sort the nodes, otherwise same as early generation.
    generationOne();
    // In early generations no nodes become sorted.
    // Therefore, we skip the pruning step and leave the
    //   nodes sorted by from.
    while(generation < 3) {
        earlyGeneration();
    }
    // On the first generation we perform a pruning step,
    //   we are forced to sort the entire list of nodes by rank
    //   in order to perform pruning step.
    firstPruneGeneration();
    // In later generations, most nodes are already sorted, so we
    //   perform a more expensive random access join with nodes in rank order
    //   in return for avoiding having to sort by rank in order to prune nodes.
    while(!isSorted()) {
        lateGeneration();
    }
    // In the generateEdges method it is convenient to begin with nodes sorted by from.
    // We perform this action here, while we still have past_nodes allocated to avoid
    //   an in-place sort.
    nodes.resizeNoCopyExact(past_nodes.size());
    radix_sort_copy<PathNode, PathNodeFromCmp, index_t>(past_nodes.begin(), past_nodes.end(), nodes.ptr(),
            &PathNodeFrom, max_from, nthreads);
    past_nodes.nullify();
    from_table.nullify();
    
    if(file_rf) {
        base.read(rf_fname, bigEndian);
        std::remove(rf_fname.c_str());
    }
}

//make original unsorted PathNodes given a RefGraph
template <typename index_t>
void PathGraph<index_t>::makeFromRef(RefGraph<index_t>& base) {
    // Create a path node per edge with a key set to from node's label
    temp_nodes = (index_t)base.edges.size() + 1;
    max_from = 0;
    nodes.reserveExact(temp_nodes);
    for(index_t i = 0; i < base.edges.size(); i++) {
        const typename RefGraph<index_t>::Edge& e = base.edges[i];
        nodes.expand();
        nodes.back().from = e.from;
        if(e.from > max_from) max_from = e.from;
        nodes.back().to = e.to;

        switch(base.nodes[e.from].label) {
        case 'A':
            nodes.back().key = pair<index_t, index_t>(0, 0);
            break;
        case 'C':
            nodes.back().key = pair<index_t, index_t>(1, 0);
            break;
        case 'G':
            nodes.back().key = pair<index_t, index_t>(2, 0);
            break;
        case 'T':
            nodes.back().key = pair<index_t, index_t>(3, 0);
            break;
        case 'Y':
            nodes.back().key = pair<index_t, index_t>(4, 0);
            break;
        default:
            assert(false);
            throw 1;
        }
    }
    // Final node.
    assert_lt(base.lastNode, base.nodes.size());
    assert_eq(base.nodes[base.lastNode].label, 'Z');
    nodes.expand();
    nodes.back().from = nodes.back().to = base.lastNode;
    if(base.lastNode > max_from) max_from = base.lastNode;
    nodes.back().key = pair<index_t, index_t>(5, 0);
    printInfo();
}

template <typename index_t>
void PathGraph<index_t>::generationOne() {
    //nodes enter almost sorted by from
    //this is only generation method that whose
    // incoming nodes are in the nodes EList
    generation++;
    //Sort nodes by from using counting sort
    //Copy into past_nodes in the process
    //first count number with each from value
    for(PathNode* node = nodes.begin(); node != nodes.end(); node++) {
        nodes[node->from].key.second++;
    }
    //convert into an index
    index_t tot = nodes[0].key.second;
    nodes[0].key.second = 0;
    for(index_t i = 1; i < max_from + 2; i++) {
        tot += nodes[i].key.second;
        nodes[i].key.second = tot - nodes[i].key.second;
    }
    // use past_nodes as from_table
    past_nodes.resizeExact(nodes.size());
    for(PathNode* node = nodes.begin(); node != nodes.end(); node++) {
        past_nodes[nodes[node->from].key.second++] = *node;
    }
    //reset index
    for(index_t i = max_from + 1; i > 0; i--) {
        past_nodes[i].key.second = nodes[i - 1].key.second;
    }
    past_nodes[0].key.second = 0;
    //Now query direct-access table
    createNewNodes();
    printInfo();
    past_nodes.swap(nodes);
}

template <typename index_t>
void PathGraph<index_t>::earlyGeneration() {
    //past_nodes enter sorted by from
    //do not yet need to perform pruning step
    generation++;
    for(index_t i = 0; i < past_nodes.size(); i++) {
        past_nodes[past_nodes[i].from + 1].key.second = i + 1;
    }
    createNewNodes();
    printInfo();
    past_nodes.swap(nodes);
}

template <typename index_t>
void PathGraph<index_t>::firstPruneGeneration() {
    //past_nodes enter sorted by from
    //first generation where we need to perform pruning step
    // results in us needing to sort entirety of nodes after they are made
    generation++;
    //here past_nodes is already sorted by .from
    // first count where to start each from value
    time_t start = time(0);
    //Build from_index
    for(index_t i = 0; i < past_nodes.size(); i++) {
        past_nodes[past_nodes[i].from + 1].key.second = i + 1;
    }
    if(verbose) cerr << "BUILT FROM_INDEX: " << time(0) - start << endl;
    start = time(0);
    // Now query against direct-access table
    createNewNodes();
    past_nodes.resizeNoCopyExact(nodes.size());

    if(verbose) cerr << "RESIZE NODES: " << time(0) - start << endl;
    start = time(0);

    //max_rank always corresponds to repeated Z's
    // Z is mapped to 0x101
    // therefore max rank = 101101101101101101101101 = (101) 8 times
    index_t max_rank = 11983725;
    radix_sort_copy<PathNode, less<PathNode>, index_t>(nodes.begin(), nodes.end(), past_nodes.ptr(),
            &PathNodeKey, max_rank, nthreads);

    if(verbose) cerr << "SORT NODES: " << time(0) - start << endl;
    start = time(0);

    nodes.swap(past_nodes);
    mergeUpdateRank();

    if(verbose) cerr << "MERGE, UPDATE RANK: " << time(0) - start << endl;
    start = time(0);

    printInfo();
    past_nodes.swap(nodes);
}

template <typename index_t>
void PathGraph<index_t>::lateGeneration() {
    //past_nodes enter sorted by rank
    //build direct-access table sorted by from,
    //but query with original nodes sorted by rank
    //since nodes we query with are sorted by rank,
    // the nodes produced are automatically sorted by key.first
    // therefore we only need to sort clusters with same key.first
    generation++;
    time_t overall = time(0);
    time_t indiv = time(0);
    assert_gt(nthreads, 0);
    assert_neq(past_nodes.size(), ranks);
    from_table.resizeNoCopy(past_nodes.size());

    if(verbose) cerr << "ALLOCATE FROM_TABLE: " << time(0) - indiv << endl;
    indiv = time(0);

    radix_sort_copy<PathNode, PathNodeFromCmp, index_t>(past_nodes.begin(), past_nodes.end(), from_table.ptr(),
            &PathNodeFrom, max_from, nthreads);

    if(verbose) cerr << "BUILD TABLE: " << time(0) - indiv << endl;
    indiv = time(0);

    //Build from_index
    for(index_t i = 0; i < from_table.size(); i++) {
        from_table[from_table[i].from + 1].key.second = i + 1;
    }

    if(verbose) cerr << "BUILD INDEX: " << time(0) - indiv << endl;

    createNewNodes();

    indiv = time(0);

    mergeUpdateRank();

    if(verbose) cerr << "MERGEUPDATERANK: " << time(0) - indiv << endl;
    if(verbose) cerr << "TOTAL TIME: " << time(0) - overall << endl;
    
    if(ranks >= (index_t)max_num_nodes) {
        throw ExplosionException();
    }

    printInfo();
    past_nodes.swap(nodes);
}

//-----------------------------------------------------------------------------------------------

template <typename index_t>
void PathGraph<index_t>::createNewNodesCounter(void* vp) {
    CreateNewNodesParams* params = (CreateNewNodesParams*)vp;
    PathNode*           st       = params->st;
    PathNode*           en       = params->en;
    PathGraph<index_t>& graph    = *(params->graph);

    size_t count = 0;
    if(graph.generation > 4) {
        for(PathNode* node = st; node != en; node++) {
            if(node->isSorted()) {
                count++;
            } else {
                count += graph.from_table[node->to + 1].key.second - graph.from_table[node->to].key.second;
            }
        }
    } else {
        for(PathNode* node = st; node != en; node++) {
            count += graph.past_nodes[node->to + 1].key.second - graph.past_nodes[node->to].key.second;
        }
    }
    *(params->sub_temp_nodes) = (index_t)count;

    //check for overflow
    if(count > (index_t)-1) {
        cerr << "exceeded integer bounds, remove adjacent SNPs, use haplotypes, or switch to a large index (--large-index)" << endl;
        throw 1;
    }
}
template <typename index_t>
void PathGraph<index_t>::createNewNodesMaker(void* vp) {
    CreateNewNodesParams* params = (CreateNewNodesParams*)vp;
    PathNode*           st       = params->st;
    PathNode*           en       = params->en;
    PathNode*           curr     = params->curr;
    PathGraph<index_t>& graph    = *(params->graph);
    if(graph.generation > 4) {
        for(PathNode* node = st; node != en; node++) {
            if(node->isSorted()) {
                *curr++ = *node;
            } else {
                for(index_t j = graph.from_table[node->to].key.second; j < graph.from_table[node->to + 1].key.second; j++) {
                    curr->from = node->from;
                    curr->to = graph.from_table[j].to;
                    (curr++)->key  = pair<index_t, index_t>(node->key.first, graph.from_table[j].key.first);
                }
            }
        }
    } else if(graph.generation == 4) {
        for(PathNode* node = st; node != en; node++) {
            for(index_t j = graph.past_nodes[node->to].key.second; j < graph.past_nodes[node->to + 1].key.second; j++) {
                curr->from = node->from;
                curr->to = graph.past_nodes[j].to;
                (curr++)->key  = pair<index_t, index_t>(node->key.first, graph.past_nodes[j].key.first);
            }
        }
    } else {
        for(PathNode* node = st; node != en; node++) {
            for(index_t j = graph.past_nodes[node->to].key.second; j < graph.past_nodes[node->to + 1].key.second; j++) {
                curr->from = node->from;
                curr->to = graph.past_nodes[j].to;
                index_t bit_shift = 1 << (graph.generation - 1);
                bit_shift = (bit_shift << 1) + bit_shift;
                (curr++)->key  = pair<index_t, index_t>((node->key.first << bit_shift) + graph.past_nodes[j].key.first, 0);
            }
        }
    }
}

template <typename index_t>
void PathGraph<index_t>::createNewNodes() {
    time_t indiv = time(0);
    AutoArray<tthread::thread*> threads(nthreads);
    EList<CreateNewNodesParams> params; params.resizeExact(nthreads);
    EList<index_t> sub_temp_nodes; sub_temp_nodes.resizeExact(nthreads); sub_temp_nodes.fillZero();
    PathNode* st = past_nodes.begin();
    PathNode* en = st + past_nodes.size() / nthreads;
    for(int i = 0; i < nthreads; i++) {
        params[i].sub_temp_nodes = &sub_temp_nodes[i];
        params[i].st = st;
        params[i].en = en;
        params[i].graph = this;
        if(nthreads == 1) {
            createNewNodesCounter((void*)&params[0]);
        } else {
            threads[i] = new tthread::thread(&createNewNodesCounter, (void*)&params[i]);
        }
        st = en;
        if(i + 2 == nthreads) {
            en = past_nodes.end();
        } else {
            en = st + past_nodes.size() / nthreads;
        }
    }

    if(nthreads > 1) {
        for(int i = 0; i < nthreads; i++)
            threads[i]->join();
    }
    if(verbose) cerr << "COUNTED NEW NODES: " << time(0) - indiv << endl;
    indiv = time(0);
    //update all label indexes
    temp_nodes = 0;
    for(int i = 0; i < nthreads; i++) {
        // done to check if we exceed index_t range
        size_t val = (size_t)temp_nodes + (size_t)sub_temp_nodes[i];
        if(val > (index_t)-1) {
            cerr << "exceeded integer bounds, remove adjacent SNPs, use haplotypes, or switch to a large index (--large-index)" << endl;
            throw 1;
        }
        temp_nodes = (index_t)val;
    }
    if(verbose) cerr << "COUNTED TEMP NODES: " << time(0) - indiv << endl;
    indiv = time(0);
    nodes.resizeNoCopyExact(temp_nodes);
    if(verbose) cerr << "RESIZED NODES: " << time(0) - indiv << endl;
    indiv = time(0);
    temp_nodes = 0;
    for(int i = 0; i < nthreads; i++) {
        params[i].curr = nodes.begin() + temp_nodes;
        temp_nodes += sub_temp_nodes[i];
    }
    if(verbose) cerr << "RESIZED NODES: " << time(0) - indiv << endl;
    indiv = time(0);
    //make new nodes
    for(int i = 0; i < nthreads; i++) {
        if(nthreads == 1) {
            createNewNodesMaker((void*)&params[0]);
        } else {
            threads[i] = new tthread::thread(&createNewNodesMaker, (void*)&params[i]);
        }
    }

    if(nthreads > 1) {
        for(int i = 0; i < nthreads; i++)
            threads[i]->join();
    }
    if(verbose) cerr << "MADE NEW NODES: " << time(0) - indiv << endl;
    indiv = time(0);
}

//------------------------------------------------------------------------------------

template <typename index_t>
void PathGraph<index_t>::mergeUpdateRank()
{
    if(generation == 4) {
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

        // Set nodes that become sorted as sorted
        PathNode* candidate = &nodes.front();
        pair<index_t, index_t> key = candidate->key; ranks = 1;
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
        ranks = 0;
        key = nodes.front().key;
        for(index_t i = 0; i < nodes.size(); i++) {
            PathNode& node = nodes[i];
            if(node.key != key) {
                key = node.key;
                ranks++;
            }
            node.key = pair<index_t, index_t>(ranks, 0);
        }
        ranks++;
    } else {
        PathNode* block_start = nodes.begin();
        PathNode* curr = nodes.begin();
        PathNode* node = nodes.begin();
        ranks = 0;
        do {
            node++;
            if(node == nodes.end() || node->key.first != block_start->key.first) {
                if(node - block_start == 1) {
                    block_start->key.first = ranks++;
                    *curr++ = *block_start;
                } else {
                    sort(block_start, node, PathNodeKeySecondCmp());
                    while(block_start != node) {
                        //extend shift while share same key
                        index_t shift = 1;
                        while(block_start + shift != node && block_start->key == (block_start + shift)->key) {
                            shift++;
                        }
                        //check if all share .from
                        //if they share same from, then they are a mergable set
                        bool merge = true;
                        for(PathNode* n = block_start; n != (block_start + shift); n++) {
                            if(n->from != block_start->from) {
                                merge = false;
                                break;
                            }
                        }
                        //if not mergable, just write all to array
                        if(!merge) {
                            for(PathNode* n = block_start; n != (block_start + shift); n++) {
                                n->key.first = ranks;
                                *curr++ = *n;
                            }
                            ranks++;
                        } else if(curr == nodes.begin() || !(curr - 1)->isSorted() || (curr - 1)->from != block_start->from) {
                            block_start->setSorted();
                            block_start->key.first = ranks++;
                            *curr++ = *block_start;
                        }
                        block_start += shift;
                    }
                    // if we are at the last node or the last node is mergable into the previous node, we are done
                    if(node == nodes.end()) break;
                    if(node + 1 == nodes.end()) {
                        assert(curr >= nodes.begin() + 1);
                        if((curr - 1)->isSorted() && node->from == (curr - 1)->from)
                            break;
                    }
                    // check if we can safely merge the node immediately following the unsorted cluster into the previous node
                    // must be that:
                    // 1) node is not itself part of an unsorted cluster
                    // 2) the previous node is sorted
                    // 3) the nodes share the same from attribute
                    assert(node + 1 < nodes.end());
                    if(node->key.first != (node + 1)->key.first) {
                        assert(curr >= nodes.begin() + 1);
                        if((curr - 1)->isSorted() && node->from == (curr - 1)->from)
                            node++;
                    }
                }
                block_start = node;
            }
        } while(node != nodes.end());
        nodes.resizeExact((index_t)(curr - nodes.begin()));
    }
    // if all nodes have unique rank we are done!
    if(ranks == nodes.size()) sorted = true;
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
    range.second = (index_t)nodes.size();
    return range;
}

//-----------------------------------------------------------------------------------------

template <typename index_t>
void PathGraph<index_t>::printInfo()
{
    if(verbose) {
        cerr << "Generation " << generation
        << " (" << temp_nodes << " -> " << nodes.size() << " nodes, "
        << ranks << " ranks)" << endl;
    }
}

//------------------------------------------------------------------------------------------

template <typename index_t>
void PathGraph<index_t>::generateEdgesCounter(void* vp) {
    GenEdgesParams* params = (GenEdgesParams*)vp;
    typename RefGraph<index_t>::Edge*        st          = params->st;
    typename RefGraph<index_t>::Edge*        en          = params->en;
    EList<index_t, 6>&                       label_index = *(params->label_index);
    EList<typename RefGraph<index_t>::Node>& ref_nodes   = *(params->ref_nodes);
    EList<PathNode>&                         nodes       = *(params->nodes);
    //first count edges, fill out label_index
    for(typename RefGraph<index_t>::Edge* edge = st; edge != en; edge++) {
        char curr_label = ref_nodes[edge->from].label;
        int curr_label_index;
        switch(curr_label) {
            case 'A': curr_label_index = 0; break;
            case 'C': curr_label_index = 1; break;
            case 'G': curr_label_index = 2; break;
            case 'T': curr_label_index = 3; break;
            case 'Y': curr_label_index = 4; break;
            case 'Z': curr_label_index = 5; break;
            default: assert(false); throw 1;
        }
        assert_lt(edge->to + 1, nodes.size());
        assert_lt(nodes[edge->to].key.second, nodes[edge->to + 1].key.second);
        label_index[curr_label_index] += nodes[edge->to + 1].key.second - nodes[edge->to].key.second;
    }
}

template <typename index_t>
void PathGraph<index_t>::generateEdgesMaker(void* vp) {
    GenEdgesParams* params = (GenEdgesParams*)vp;
    typename RefGraph<index_t>::Edge*        st          = params->st;
    typename RefGraph<index_t>::Edge*        en          = params->en;
    EList<index_t, 6>&                       label_index = *(params->label_index);
    EList<typename RefGraph<index_t>::Node>& ref_nodes   = *(params->ref_nodes);
    EList<PathEdge>&                         edges       = *(params->edges);
    EList<PathNode>&                         nodes       = *(params->nodes);
    for(typename RefGraph<index_t>::Edge* edge = st; edge != en; edge++) {
        char curr_label = ref_nodes[edge->from].label;
        int curr_label_index;
        switch(curr_label) {
            case 'A': curr_label_index = 0; break;
            case 'C': curr_label_index = 1; break;
            case 'G': curr_label_index = 2; break;
            case 'T': curr_label_index = 3; break;
            case 'Y': curr_label_index = 4; break;
            case 'Z': curr_label_index = 5; break;
            default: assert(false); throw 1;
        }
        for(index_t j = nodes[edge->to].key.second; j < nodes[edge->to + 1].key.second; j++) {
            edges[label_index[curr_label_index]++] = PathEdge(edge->from, nodes[j].key.first, curr_label);
        }
    }
}

template <typename index_t>
bool PathGraph<index_t>::generateEdges(RefGraph<index_t>& base)
{
    //entering we have:
    // nodes - sorted by from
    // edges - empty
    // base.nodes - almost sorted by from/to
    // base.edges - almost sorted by from/to

    //need to join:
    // nodes.from -> base.nodes[]
    // nodes.from -> base.edges.to
    // nodes.from -> edges.from

    if(!sorted) return false;

    time_t indiv = time(0);
    time_t overall = time(0);

    //replace nodes.to with genomic position
    //fast because both roughly ordered by from
    for(PathNode* node = nodes.begin(); node != nodes.end(); node++) {
        node->to = base.nodes[node->from].value;
    }

    if(verbose) cerr << "NODE.TO -> GENOME POS: "  << time(0) - indiv << endl;
    indiv = time(0);

    // build an index for nodes
    for(index_t i = 0; i < nodes.size(); i++) {
        nodes[nodes[i].from + 1].key.second = i + 1;
    }

    if(verbose) cerr << "BUILD FROM_INDEX "  << time(0) - indiv << endl;
    indiv = time(0);

    // Now join nodes.from to edges.to
    // fast because base.edges roughly sorted by to

    //count number of edges
    AutoArray<tthread::thread*> threads(nthreads);
    EList<GenEdgesParams> params; params.resizeExact(nthreads);
    ELList<index_t, 6> label_index; label_index.resize(nthreads);
    typename RefGraph<index_t>::Edge* st = base.edges.begin();
    typename RefGraph<index_t>::Edge* en = st + base.edges.size() / nthreads;
    for(int i = 0; i < nthreads; i++) {
        label_index[i].resizeExact(6);
        label_index[i].fillZero();
        params[i].label_index = &label_index[i];
        params[i].st = st;
        params[i].en = en;
        params[i].nodes = &nodes;
        params[i].edges = &edges;
        params[i].ref_nodes = &base.nodes;
        if(nthreads == 1) {
            generateEdgesCounter((void*)&params[0]);
        } else {
            threads[i] = new tthread::thread(&generateEdgesCounter, (void*)&params[i]);
        }
        st = en;
        if(i + 2 == nthreads) {
            en = base.edges.end();
        } else {
            en = st + base.edges.size() / nthreads;
        }
    }

    if(nthreads > 1) {
        for(int i = 0; i < nthreads; i++)
            threads[i]->join();
    }
    
    if(verbose) cerr << "COUNTED NEW EDGES: " << time(0) - indiv << endl;
    indiv = time(0);
    //update all label indexes
    index_t tot = label_index[0][0];
    label_index[0][0] = 0;
    for(int i = 1; i < nthreads; i++) {
        tot += label_index[i][0];
        label_index[i][0] =  tot - label_index[i][0];
    }
    for(int j = 1; j < 6; j++) {
        for(int i = 0; i < nthreads; i++) {
            tot += label_index[i][j];
            label_index[i][j] =  tot - label_index[i][j];
        }
    }
    edges.resizeExact(tot);
    //make new edges
    for(int i = 0; i < nthreads; i++) {
        if(nthreads == 1) {
            generateEdgesMaker((void*)&params[0]);
        } else {
            threads[i] = new tthread::thread(&generateEdgesMaker, (void*)&params[i]);
        }
    }

    if(nthreads > 1) {
        for(int i = 0; i < nthreads; i++) {
            threads[i]->join();
        }
    }
    base.nullify();

    if(verbose) cerr << "MADE NEW EDGES: " << time(0) - indiv << endl;
    indiv = time(0);

    EList<index_t, 6>& index = label_index[nthreads - 1];

    EList<PathEdge> temp_edges; temp_edges.resizeExact(edges.size());

    radix_sort_copy<PathEdge, less<PathEdge>, index_t>(edges.begin()           , edges.begin() + index[0], temp_edges.ptr(),
            &PathEdgeTo, (index_t)nodes.size(), nthreads);
    radix_sort_copy<PathEdge, less<PathEdge>, index_t>(edges.begin() + index[0], edges.begin() + index[1], temp_edges.ptr() + index[0],
            &PathEdgeTo, (index_t)nodes.size(), nthreads);
    radix_sort_copy<PathEdge, less<PathEdge>, index_t>(edges.begin() + index[1], edges.begin() + index[2], temp_edges.ptr() + index[1],
            &PathEdgeTo, (index_t)nodes.size(), nthreads);
    radix_sort_copy<PathEdge, less<PathEdge>, index_t>(edges.begin() + index[2], edges.begin() + index[3], temp_edges.ptr() + index[2],
            &PathEdgeTo, (index_t)nodes.size(), nthreads);
    for(index_t i = index[3]; i < edges.size(); i++) {
        temp_edges[i] = edges[i];
    }
    sort(temp_edges.begin() + index[3], temp_edges.begin() + index[4]);
    sort(temp_edges.begin() + index[4], temp_edges.begin() + index[5]);
    edges.xfer(temp_edges);

    if(verbose) cerr << "SORTED NEW EDGES: " << time(0) - indiv << endl;
    indiv = time(0);

    EList<PathNode> past_nodes; past_nodes.resizeExact(nodes.size());
    radix_sort_copy<PathNode, less<PathNode>, index_t>(nodes.begin(), nodes.end(), past_nodes.ptr(),  &PathNodeKey, ranks, nthreads);
    nodes.xfer(past_nodes);

    if(verbose) cerr << "RE-SORTED NODES: " << time(0) - indiv << endl;
    indiv = time(0);

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

#ifndef NDEBUG

    // Switch char array[x][y]; to char** array;
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
    {
      PathNode* node = nodes.begin(); node->key.first = 0;
      PathEdge* edge = edges.begin();
      while(node != nodes.end() && edge != edges.end()) {
          if(edge->from == node->from) {
              edge->from = (index_t)(node - nodes.begin()); edge++;
              node->key.first++;
          } else {
              node++; node->key.first = 0;
          }
      }
    }

    if(verbose) cerr << "PROCESS EDGES: " << time(0) - indiv << endl;
    indiv = time(0);

    // Remove 'Y' node
    assert_gt(nodes.size(), 2);
    nodes.back().key.first = nodes[nodes.size() - 2].key.first;
    nodes[nodes.size() - 2] = nodes.back();
    nodes.pop_back();
    // Adjust edges accordingly
    for(size_t i = 0; i < edges.size(); i++) {
        PathEdge& edge = edges[i];
        if(edge.label == 'Y') {
            edge.label = 'Z';
        } else if(edge.ranking >= nodes.size()) {
            assert_eq(edge.ranking, nodes.size());
            edge.ranking -= 1;
        }
    }
    if(verbose) cerr << "REMOVE Y: " << time(0) - indiv << endl;
    indiv = time(0);

#ifndef NDEBUG
    if(debug) {
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
    temp_edges.resizeExact(edges.size());
    radix_sort_copy<PathEdge, PathEdgeToCmp, index_t>(edges.begin(), edges.end(), temp_edges.ptr(), &PathEdgeTo, (index_t)nodes.size(), nthreads);
    edges.xfer(temp_edges);
    for(index_t i = 0; i < edges.size(); i++) {
        nodes[edges[i].ranking].key.second = i + 1;
    }

    if(verbose) cerr << "SORT, Make index: " << time(0) - indiv << endl;
    if(verbose) cerr << "TOTAL: " << time(0) - overall << endl;
    return true;

//-----------------------------------------------------------------------------------------------------
    bwt_string.clear();
    F_array.clear();
    M_array.clear();
    bwt_counts.resizeExact(5); bwt_counts.fillZero();
    for(index_t node = 0; node < nodes.size(); node++) {
        pair<index_t, index_t> edge_range = getEdges(node, false /* from? */);
        for(index_t i = edge_range.first; i < edge_range.second; i++) {
            assert_lt(i, edges.size());
            char label = edges[i].label;
            if(label == 'Y') {
                label = 'Z';
            }
            bwt_string.push_back(label);
            F_array.push_back(i == edge_range.first ? 1 : 0);

            if(label != 'Z') {
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

    // Test searches, based on paper_example
#if 1
    EList<string> queries;  EList<index_t> answers;
#   if 1
#      if 1
    queries.push_back("GACGT"); answers.push_back(9);
    queries.push_back("GATGT"); answers.push_back(9);
    queries.push_back("GACT");  answers.push_back(9);
    queries.push_back("ATGT");  answers.push_back(4);
    queries.push_back("GTAC");  answers.push_back(10);
    queries.push_back("ACTG");  answers.push_back(3);
#      else
    // rs55902548, at 402, ref, alt, unknown alt
    queries.push_back("GGCAGCTCCCATGGGTACACACTGGGCCCAGAACTGGGATGGAGGATGCA");
    // queries.push_back("GGCAGCTCCCATGGGTACACACTGGTCCCAGAACTGGGATGGAGGATGCA");
    // queries.push_back("GGCAGCTCCCATGGGTACACACTGGACCCAGAACTGGGATGGAGGATGCA");

    // rs5759268, at 926787, ref, alt, unknown alt
    // queries.push_back("AAATTGCTCAGCCTTGTGCTGTGCACACCTGGTTCTCTTTCCAGTGTTAT");
    // queries.push_back("AAATTGCTCAGCCTTGTGCTGTGCATACCTGGTTCTCTTTCCAGTGTTAT");
    // queries.push_back("AAATTGCTCAGCCTTGTGCTGTGCAGACCTGGTTCTCTTTCCAGTGTTAT");
#      endif

    for(size_t q = 0; q < queries.size(); q++) {
        const string& query = queries[q];
        assert_gt(query.length(), 0);
        index_t top = 0, bot = edges.size();
        index_t node_top = 0, node_bot = 0;
        cerr << "Aligning " << query <<  endl;
        index_t i = 0;
        for(; i < query.length(); i++) {
            if(top >= bot) break;
            int nt = query[query.length() - i - 1];
            nt = asc2dna[nt];
            assert_lt(nt, 4);

            cerr << "\t" << i << "\tBWT range: [" << top << ", " << bot << ")" << endl;

            top = bwt_counts[(int)nt] + (top <= 0 ? 0 : rank(bwt_string, top - 1, "ACGT"[nt]));
            bot = bwt_counts[(int)nt] + rank(bwt_string, bot - 1, "ACGT"[nt]);
            cerr << "\t\tLF BWT range: [" << top << ", " << bot << ")" << endl;

            node_top = rank1(M_array, top) - 1;
            node_bot = rank1(M_array, bot - 1);
            cerr << "\t\tnode range: [" << node_top << ", " << node_bot << ")" << endl;

            top = select1(F_array, node_top + 1);
            bot = select1(F_array, node_bot + 1);
        }
        cerr << "\t" << i << "\tBWT range: [" << top << ", " << bot << ")" << endl;
        // assert_eq(top, answers[q]);
        cerr << "finished... ";
        if(node_top < node_bot && node_top < nodes.size()) {
            index_t pos = nodes[node_top].to;
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

    // See inconsistencies between F and M arraystimy thread
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

//--------------------------------------------------------------------------

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

#endif /*GBWT_GRAPH_H_*/
