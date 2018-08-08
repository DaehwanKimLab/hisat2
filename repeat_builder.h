 /*
 * Copyright 2018, Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>
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

#ifndef __REPEAT_BUILDER_H__
#define __REPEAT_BUILDER_H__

#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include "assert_helpers.h"
#include "word_io.h"
#include "mem_ids.h"
#include "ref_coord.h"
#include "ref_read.h"
#include "edit.h"
#include "ds.h"
#include "repeat.h"
#include "blockwise_sa.h"
#include "simple_func.h"
#include "scoring.h"
#include "aligner_sw.h"

//#define DEBUGLOG

using namespace std;

/**
 * Encapsulates repeat parameters.
 */
class RepeatParameter {
public:
    TIndexOffU seed_len;         // seed length
    TIndexOffU seed_count;       // seed count
    TIndexOffU seed_mm;          // maximum edit distance allowed during initial seed extension
    TIndexOffU repeat_count;     // repeat count
    TIndexOffU min_repeat_len;   // minimum repeat length
    TIndexOffU max_repeat_len;   // maximum repeat length
    TIndexOffU max_edit;         // maximum edit distance allowed
    bool       symmetric_extend; // extend symmetrically
    TIndexOffU extend_unit_len;  // extend seeds no longer than this length at a time
};

typedef pair<TIndexOffU, TIndexOffU> Range;

struct Fragments {
	bool contain(TIndexOffU pos) {
		if (pos >= joinedOff && pos < (joinedOff + length)) {
			return true;
		}
		return false;
	}

    TIndexOffU joinedOff;       // index within joined text
	TIndexOffU length;

	int frag_id;
	int seq_id;
	TIndexOffU seqOff;          // index within sequence 
	bool first;
};

class CoordHelper {
public:
    CoordHelper(TIndexOffU length,
                TIndexOffU forward_length,
                const EList<RefRecord>& szs,
                const EList<string>& ref_names);
    ~CoordHelper();
    int mapJoinedOffToSeq(TIndexOffU joined_pos);
    int getGenomeCoord(TIndexOffU joined_pos,
                       string& chr_name,
                       TIndexOffU& pos_in_chr);

    TIndexOffU getEnd(TIndexOffU e);
    TIndexOffU getStart(TIndexOffU e);
    TIndexOffU forward_length() const { return forward_length_; }
    TIndexOffU length() const { return length_; }
    
private:
    void buildNames();
    void buildJoinedFragment();
    
private:
    TIndexOffU       length_;
    TIndexOffU       forward_length_;
    EList<RefRecord> szs_;
    EList<string>    ref_names_;
    EList<string>    ref_namelines_;
    
    // mapping info from joined string to genome
    EList<Fragments> fraglist_;
    
    // Fragments Cache
#define CACHE_SIZE_JOINEDFRG    10
    Fragments cached_[CACHE_SIZE_JOINEDFRG];
    int num_cached_ = 0;
    int victim_ = 0;    /* round-robin */
};

struct SeedHP {
    bool operator==(const SeedHP &rhs) const {
        return range == rhs.range &&
               snpIDs == rhs.snpIDs;
    }

    bool operator!=(const SeedHP &rhs) const {
        return !(rhs == *this);
    }

    bool operator<(const SeedHP &rhs) const {
        return range < rhs.range;
    }

    bool operator>(const SeedHP &rhs) const {
        return rhs < *this;
    }

    bool operator<=(const SeedHP &rhs) const {
        return !(rhs < *this);
    }

    bool operator>=(const SeedHP &rhs) const {
        return !(*this < rhs);
    }

    Range range;
    EList<string> snpIDs;
};

struct SeedSNP {
    int type; // EDIT_TYPE_MM, EDIT_TYPE_READ_GAP, EDIT_TYPE_REF_GAP
    TIndexOffU pos;
    size_t len;
    string base; // valid when ty = MM, REF_GAP
    TIndexOffU id;

    SeedSNP() {}

    void init(int ty, TIndexOffU po, size_t ln, string bs) {
        type = ty;
        pos = po;
        len = ln;
        base = bs;
    }
    void init(int ty, TIndexOffU po, size_t ln, char bs) {
        type = ty;
        pos = po;
        len = ln;
        base.assign(1, bs);
    }

    friend ostream &operator<<(ostream &os, const SeedSNP &snp) {

        if(snp.type == EDIT_TYPE_MM) {
            os << "single";
        } else if(snp.type == EDIT_TYPE_REF_GAP) {
            os << "insertion";
        } else if(snp.type == EDIT_TYPE_READ_GAP) {
            os << "deletion";
        }
        os << " ";
        os << snp.pos << " ";

        if(snp.type == EDIT_TYPE_MM || snp.type == EDIT_TYPE_REF_GAP) {
            os << snp.base;
        } else if(snp.type == EDIT_TYPE_READ_GAP) {
            os << snp.len;
        }
        return os;
    }

    bool operator==(const SeedSNP &rhs) const {
        return type == rhs.type &&
               pos == rhs.pos &&
               len == rhs.len &&
               base == rhs.base;
    }

    bool operator!=(const SeedSNP &rhs) const {
        return !(rhs == *this);
    }

    static bool cmpSeedSNPByPos( SeedSNP * const& a,  SeedSNP * const& b)  {
        return a->pos < b->pos;
    }

};

class SeedExt {
public:
    SeedExt() {
        reset();
    };

    void reset() {
        done = false;
        curr_ext_len = 0;
        ed = total_ed = 0;
        orig_pos.first = 0;
        orig_pos.second = 0;
        pos.first = 0;
        pos.second = 0;
        bound.first = 0;
        bound.second = 0;
        consensus_pos.first = 0;
        consensus_pos.second = 0;
        left_gaps.clear();
        right_gaps.clear();
        aligned = true;
    };

    TIndexOffU getLeftExtLength() const {
        assert_leq(pos.first, orig_pos.first);
        TIndexOffU len = orig_pos.first - pos.first;
        return len;
    }

    TIndexOffU getRightExtLength() const {
        assert_geq(pos.second, orig_pos.second);
        return pos.second - orig_pos.second;
    }

    TIndexOffU getLength() const {
        TIndexOffU len = orig_pos.second - orig_pos.first;
        len += getLeftExtLength();
        len += getRightExtLength();
        return len;
    }

    template<typename TStr>
    void getLeftExtString(const TStr& s, string& seq)
    {
        seq.clear();
        seq = getString(s, pos.first, orig_pos.first - pos.first);
    }

    template<typename TStr>
    void getRightExtString(const TStr& s, string& seq)
    {
        seq.clear();
        seq = getString(s, orig_pos.second, pos.second - orig_pos.second);
    }

    template<typename TStr>
    void getExtendedSeedSequence(const TStr& s,
                                 string& seq) const;

    template<typename TStr>
    void generateSNPs(const TStr& s, const string& consensus, EList<SeedSNP *>& repeat_snps);

    static bool isSameConsensus(const SeedExt& a, const SeedExt& b) {
        return (a.consensus_pos == b.consensus_pos);
            //&& (a.getLength() == b.getLength());
    }

    static bool isSameSNPs(const SeedExt& a, const SeedExt& b) {
        const EList<SeedSNP *>& a_edit = a.snps;
        const EList<SeedSNP *>& b_edit = b.snps;

        if(a_edit.size() != b_edit.size()) {
            return false;
        }

        for(size_t i = 0; i < a_edit.size(); i++) {
            if(!(*a_edit[i] == *b_edit[i])) {
                return false;
            }
        }
        return true;
    }

    static bool isSameAllele(const SeedExt& a, const SeedExt& b) {
        const EList<Edit>& a_edit = a.edits;
        const EList<Edit>& b_edit = b.edits;

        if(a_edit.size() != b_edit.size()) {
            return false;
        }

        for(size_t i = 0; i < a_edit.size(); i++) {
            if(!(a_edit[i] == b_edit[i])) {
                return false;
            }
        }
        return true;
    }

    static bool sort_by_edits(const SeedExt& a, const SeedExt& b) {
        const EList<Edit>& a_edit = a.edits;
        const EList<Edit>& b_edit = b.edits;

        size_t i = 0;

        while((i < a_edit.size())
                && (i < b_edit.size())) {

            if(!(a_edit[i] == b_edit[i])) {
                if(a_edit[i] < b_edit[i]) {
                    return true;
                } else {
                    return false;
                }
            }

            i++;
        }

        if(a_edit.size() < b_edit.size()) {
            return true;
        }

        if((a_edit.size() == b_edit.size())
                && (a.pos.first < b.pos.first)) {
            return true;
        }
        return false;
    }

    void generateEdits(const string& consensus_merged, const string& seed_ext)
    {
        size_t ed_done = 0;
        size_t ext_len = getLength();
        assert_eq(ext_len, seed_ext.length());

        edits.clear();

        for(size_t i = 0; i < ext_len && ed_done < total_ed; i++) {
            char con_base = consensus_merged[consensus_pos.first + i];
            char seed_base = seed_ext[i];

            if (con_base != seed_base) {
                edits.expand();
                edits.back().init(consensus_pos.first + i,
                                  con_base, seed_base, EDIT_TYPE_MM);
                ed_done++;
            }
        }
    }
    
    Range getExtendedRange(size_t consensus_len) const
    {
        assert_leq(consensus_pos.second, consensus_len);
        Range range;
        range.first = (pos.first < consensus_pos.first ? 0 : pos.first - consensus_pos.first);
        range.second = pos.second + (consensus_len - consensus_pos.second);
        return range;
    }
    
#ifndef NDEBUG
    bool valid() const {
        assert_leq(consensus_pos.first, consensus_pos.second);
        TIndexOffU constr_len = consensus_pos.second - consensus_pos.first;
        TIndexOffU allele_len = getLength();
        
        TIndexOffU cur_off = 0;
        for(size_t i = 0; i < left_gaps.size(); i++) {
            assert_geq(left_gaps[i].first, cur_off);
            cur_off = left_gaps[i].first;
            int gap_len = left_gaps[i].second;
            assert_neq(gap_len, 0);
            if(gap_len > 0) { // deletion
                allele_len += gap_len;
            } else {
                allele_len += gap_len;
                cur_off += (-gap_len);
            }
        }
        cur_off = 0;
        for(size_t i = 0; i < right_gaps.size(); i++) {
            assert_geq(right_gaps[i].first, cur_off);
            cur_off = right_gaps[i].first;
            int gap_len = right_gaps[i].second;
            assert_neq(gap_len, 0);
            if(gap_len > 0) { // deletion
                allele_len += gap_len;
            } else {
                allele_len += gap_len;
                cur_off += (-gap_len);
            }
        }
        assert_eq(constr_len, allele_len);
        return true;
    }
#endif
    
public:
    // seed extended position [first, second)
    pair<TIndexOffU, TIndexOffU> orig_pos;
    pair<TIndexOffU, TIndexOffU> pos;
    
    // extension bound. the seed must be placed on same fragment
    // [first, second)
    pair<TIndexOffU, TIndexOffU> bound;
    
    // positions relative to consensus sequence
    pair<TIndexOffU, TIndexOffU> consensus_pos;
    
    // a list of gaps (deletions and insertions) in both directions
    // offsets from seed's left and right ("pos" above)
    // positive and negative values indicate deletions and insertions, resp.
    EList<pair<TIndexOffU, int> > left_gaps;
    EList<pair<TIndexOffU, int> > right_gaps;
    
    uint32_t ed;            // edit distance
    uint32_t total_ed;      // total edit distance
    bool done;              // done flag
    uint32_t curr_ext_len;  //
    
    bool aligned;
    
    EList<Edit> edits;      // edits w.r.t. consensus_merged

    EList<SeedSNP *> snps;
};

class RB_AlleleCoord {
public:
    RB_AlleleCoord() :
    left(0),
    right(0),
    idx(0)
    {}
    
    RB_AlleleCoord(TIndexOffU l, TIndexOffU r, size_t i) :
    left(l),
    right(r),
    idx(i)
    {}
    
    TIndexOffU len() const { return right - left; }
    
    bool operator<(const RB_AlleleCoord& o) const
    {
        if(left != o.left)
            return left < o.left;
        if(right != o.right)
            return right > o.right;
        return false;
    }
    
    bool contain(const RB_AlleleCoord& o, size_t relax = 0) const
    {
        if(o.left + relax >= left && o.right <= right + relax)
            return true;
        else
            return false;
    }
    
    bool contained(const RB_AlleleCoord& o) const
    {
        return o.contain(*this);
    }
    
    TIndexOffU overlap_len(const RB_AlleleCoord& o) const
    {
        if(contain(o)) return o.len();
        else if(o.contain(*this)) return len();
        else if(left < o.right && right > o.left) {
            if(left <= o.left) {
                return right - o.left;
            } else {
                return o.right - left;
            }
        }
        return 0;
    }
    
public:
    TIndexOffU left;
    TIndexOffU right;
    size_t     idx;
};

class RB_RepeatManager;
class RB_SWAligner;
class RB_SubSA;

class RB_RepeatBase {
public:
    string seq;
    EList<TIndexOffU> nodes;
};

class RB_Repeat {
public:
    RB_Repeat() {}
    ~RB_Repeat() {
        if(snps_.size() > 0) {
            for(size_t i = 0; i < snps_.size(); i++) {
                delete snps_[i];
            }
        }
    }
    
    void repeat_id(size_t repeat_id) { repeat_id_ = repeat_id; }
    size_t repeat_id() const { return repeat_id_; }
    
    void parent_id(size_t parent_id) { parent_id_ = parent_id; }
    size_t parent_id() const { return parent_id_; }
    
    string& consensus() { return consensus_; }
    const string& consensus() const { return consensus_; }
    
    EList<SeedExt>& seeds() { return seeds_; }
    const EList<SeedExt>& seeds() const { return seeds_; }
    
    EList<RB_AlleleCoord>& seed_ranges() { return seed_ranges_; }
    const EList<RB_AlleleCoord>& seed_ranges() const { return seed_ranges_; }
    
    EList<SeedSNP *>& snps() { return snps_; }
    const EList<SeedSNP *>& snps() const { return snps_; }
    
    template<typename TStr>
    void init(const RepeatParameter& rp,
              const TStr& s,
              CoordHelper& coordHelper,
              const RB_SubSA& subSA,
              const RB_RepeatBase& repeatBase);
    
    template<typename TStr>
    void extendConsensus(const RepeatParameter& rp,
                         const TStr& s);
    
    template<typename TStr>
    void getNextRepeat(const RepeatParameter& rp,
                       const TStr& s,
                       RB_Repeat& o);
    
    
    template<typename TStr>
    void saveSeedExtension(const RepeatParameter& rp,
                           const TStr& s,
                           CoordHelper& coordHelper,
                           TIndexOffU seed_grp_id,
                           ostream& fp,
                           size_t& total_repeat_seq_len,
                           size_t& total_allele_seq_len) const;
    
    void saveSNPs(ofstream& fp,
                  TIndexOffU grp_id,
                  TIndexOffU& snp_id_base);
    void saveConsensus(ofstream& fp,
                       TIndexOffU grp_id);
    
    bool overlap(const RB_Repeat& o,
                 bool& contain,
                 bool& left,
                 size_t& seed_i,
                 size_t& seed_j,
                 bool debug = false) const;
    
    bool satisfy(const RepeatParameter& rp) const
    {
        if(consensus_.length() < rp.min_repeat_len)
            return false;
        if(seeds_.size() < rp.repeat_count)
            return false;
        if(seeds_[rp.repeat_count - 1].getLength() < rp.min_repeat_len)
            return false;
        return true;
    }
    
    void reset() {
        consensus_.clear();
        for(size_t i = 0; i < seeds_.size(); i++)
            seeds_[i].reset();
        seeds_.clear();
        seed_ranges_.clear();
    }
    
    void showInfo(const RepeatParameter& rp,
                  CoordHelper& coordHelper) const;
    
    template<typename TStr>
    void generateSNPs(const RepeatParameter&, const TStr& s, TIndexOffU grp_id);
    
    bool self_repeat() const { return self_repeat_; }
    
    void update() { internal_update(); }
    void addSeed(const SeedExt& seed)
    {
        seeds_.expand();
        seeds_.back() = seed;
    }
    
    bool contain(TIndexOffU left, TIndexOffU right) const;
    
protected:
    template<typename TStr>
    void get_consensus_seq(const TStr& s,
                           EList<SeedExt>& seeds,
                           size_t sb,
                           size_t se,
                           size_t min_left_ext,
                           size_t min_right_ext,
                           size_t max_ed,
                           const RepeatParameter& rp,
                           EList<size_t>& ed_seed_nums,
                           EList<string>* left_consensuses,
                           EList<string>* right_consensuses) const;
    
    void internal_update();
    
    
protected:
    size_t                repeat_id_;
    size_t                parent_id_;
    string                consensus_;
    EList<SeedExt>        seeds_;
    EList<RB_AlleleCoord> seed_ranges_;
    EList<SeedSNP*>       snps_;
    bool                  self_repeat_;
    
    static EList<size_t>  ca_ed_;
    static EList<size_t>  ca_ed2_;
    static string         ca_s_;
    static string         ca_s2_;
    
public:
    static size_t         seed_merge_tried;
    static size_t         seed_merged;
};

class RB_RepeatExt : public RB_Repeat {
public:
    RB_RepeatExt() {}
    ~RB_RepeatExt() {}

    template<typename TStr>
    void extendConsensus(const RepeatParameter& rp,
                         const TStr& s);
    
    float mergeable(const RB_Repeat& o) const;
    
    template<typename TStr>
    bool merge(const RepeatParameter& rp,
               const TStr& s,
               RB_SWAligner& swalginer,
               const RB_RepeatExt& o,
               bool contain,
               size_t seed_i,
               size_t seed_j,
               bool debug = false);
    
protected:
    template<typename TStr>
    void get_consensus_seq(const TStr& s,
                           EList<SeedExt>& seeds,
                           size_t sb,
                           size_t se,
                           size_t min_left_ext,
                           size_t min_right_ext,
                           size_t max_ed,
                           const RepeatParameter& rp,
                           EList<size_t>& ed_seed_nums,
                           EList<string>* left_consensuses,
                           EList<string>* right_consensuses) const;
    
    template<typename TStr>
    bool align(const RepeatParameter& rp,
               const TStr& ref,
               const string& s,
               const EList<pair<size_t, size_t> >& s_kmer_table,
               const string& s2,
               EList<int>& offsets,
               size_t k,
               SeedExt& seed,
               int consensus_approx_left,
               int consensus_approx_right,
               size_t left,
               size_t right,
               bool debug);
    
    bool isSelfRepeat(const RepeatParameter& rp,
                      const string& s,
                      const EList<pair<size_t, size_t> >& s_kmer_table,
                      EList<int>& offsets,
                      size_t k,
                      bool debug);
    
};

// check if a set of seeds are already processed
class RB_RepeatManager {
public:
    size_t numCoords() const { return range_to_repeats_.size(); }
    
    bool checkRedundant(const RepeatParameter& rp,
                        const map<size_t, RB_Repeat*>& repeat_map,
                        const EList<TIndexOffU>& positions,
                        EList<size_t>& to_remove) const;
    
    void addRepeat(const RB_Repeat* repeat);
    void addRepeat(Range range, size_t repeat_id);
    void removeRepeat(const RB_Repeat* repeat);
    void removeRepeat(Range range, size_t repeat_id);
    
public:
    void showInfo(const RepeatParameter& rp,
                  CoordHelper& coordHelper,
                  const map<size_t, RB_Repeat*>& repeat_map,
                  size_t num_case = 5) const;
    
private:
    map<Range, EList<size_t> > range_to_repeats_;
};

class RB_SWAligner {
public:
    RB_SWAligner();
    ~RB_SWAligner();
    
    void init_dyn(const RepeatParameter& rp);
    
    int alignStrings(const string &ref,
                     const string &read,
                     EList<Edit>& edits,
                     Coord& coord);
    
    void makePadString(const string& ref,
                       const string& read,
                       string& pad,
                       size_t len);
    
    void doTest(const RepeatParameter& rp,
                const string& refstr,
                const string& readstr);
    
    void doTestCase1(const string& refstr,
                     const string& readstr,
                     TIndexOffU rpt_edit);
    
private:
    //
    SimpleFunc scoreMin_;
    SimpleFunc nCeil_;
    SimpleFunc penCanIntronLen_;
    SimpleFunc penNoncanIntronLen_;
    
    Scoring *sc_;
    SwAligner swa;
    LinkedEList<EList<Edit> > rawEdits_;
    RandomSource rnd_;
};

// SA Subset
class RB_SubSA {
public:
    RB_SubSA() {}
    ~RB_SubSA();

    void writeFile(ofstream& fp);
    void readFile(ifstream &fp);
    
    void init(TIndexOffU sa_size,
              TIndexOffU seed_len,
              TIndexOffU seed_count);
    
    template<typename TStr>
    void push_back(const TStr& s,
                   CoordHelper& coordHelper,
                   TIndexOffU saElt,
                   bool lastInput = false);
    
    TIndexOffU get(size_t i) const;
    inline TIndexOffU operator[](size_t i) const { return get(i); }
    
    /**
     * Return true iff there are no elements
     * @return
     */
    inline bool empty() const { return cur_ == 0; }
    inline size_t size() const { return cur_; }
    inline TIndexOffU seed_len() const { return seed_len_; }
    inline TIndexOffU seed_count() const { return seed_count_; }
    
    const EList<TIndexOffU>& getRepeatIndex() const { return repeat_index_; }
    
    template<typename TStr>
    Range find(const TStr& s,
               const string& seq) const;
    
    template<typename TStr>
    TIndexOffU find_repeat_idx(const TStr& s,
                               const string& seq) const;
    
    void setDone(TIndexOffU off, TIndexOffU len = 1);
    bool isDone(TIndexOffU off, TIndexOffU len = 1) const;
    
    void dump() const
    {
        cerr << "seed length: " << seed_len_ << endl;
        cerr << "minimum seed count: " << seed_count_ << endl;
        cerr << "number of seed groups: " << repeat_index_.size() << endl;
        cerr << "item_bit_size_: " << item_bit_size_ << endl;
        cerr << "block_size_: " << block_size_ << endl;
        cerr << "items_per_block_: " << items_per_block_ << endl;
        cerr << "cur_: " << cur_ << endl;
        cerr << "sz_: " << sz_ << endl;
        cerr << "number of blocks: " << blocks.size() << endl;
    }
    
    size_t getMemUsage() const {
        size_t tot = blocks.size() * block_size_;
        tot += blocks.totalCapacityBytes();
        return tot;
    }
    
    template<typename TStr>
    void buildRepeatBase(const TStr& s,
                         CoordHelper& coordHelper,
                         const size_t max_len,
                         EList<RB_RepeatBase>& repeatBases);
    
private:
    void put(size_t i, TIndexOffU);
    void push_back(TIndexOffU t);
    
    inline uint32_t bit_to_mask(size_t bit) const
    {
        return (uint32_t)((1L << bit) - 1);
    }
    
    void expand(size_t sz = 1);
    // increase size
    void allocSize(size_t sz);
    void allocItems(size_t count);
    
private:
    TIndexOffU getItem(uint32_t *block, size_t idx, size_t offset) const;
    void setItem(uint32_t *block, size_t idx, size_t offset, TIndexOffU val);
    
    pair<size_t, size_t> index_to_addr(size_t index) const;
    pair<size_t, size_t> col_to_pos(size_t col) const;
    
private:
    TIndexOffU sa_size_;
    TIndexOffU seed_len_;
    TIndexOffU seed_count_;
    EList<TIndexOffU> temp_suffixes_;
    
    size_t item_bit_size_; // item bit size(e.g. 33bit)
    
    size_t block_bit_size_; // 8bit
    size_t items_per_block_bit_;
    size_t items_per_block_bit_mask_;
    size_t items_per_block_; // number of items in block
    
    size_t cur_;        // current size (in count)
    size_t sz_;         // maximum (in count)
    
    size_t block_size_;       // block size in Byte
    
    // list of packed array
    EList<uint32_t *> blocks;
    
    // repeat index
    EList<TIndexOffU> repeat_index_;
    EList<uint32_t>   done_;
};


// find and write repeats
template<typename TStr>
class RepeatBuilder {

public:
	RepeatBuilder(TStr& s,
                  const EList<RefRecord>& szs,
                  const EList<string>& ref_names,
                  bool forward_only,
                  const string& filename);
    ~RepeatBuilder();


public:
    void readSA(const RepeatParameter& rp,
                BlockwiseSA<TStr>& sa);

    void readSA(const RepeatParameter& rp,
                const string& filename);

    void writeSA(const RepeatParameter& rp,
                 const string& filename);
    
	void build(const RepeatParameter& rp);

	void sortRepeatGroup();

    void saveRepeats(const RepeatParameter& rp);
    void saveConsensus(const RepeatParameter& rp);
    void saveFile(const RepeatParameter& rp);
    
    void addRepeatGroup(const RepeatParameter& rp,
                        size_t& repeat_id,
                        RB_RepeatManager& repeat_manger,
                        const RB_RepeatBase& repeatBase,
                        ostream& fp);
    
    void reassignSeeds(const RepeatParameter& rp,
                       size_t repeat_bid,
                       size_t repeat_eid,
                       EList<SeedExt>& seeds);

    bool checkSequenceMergeable(const string& ref,
                                const string& read,
                                EList<Edit>& edits,
                                Coord& coord,
                                TIndexOffU rpt_len,
                                TIndexOffU max_edit = 10);

    void doTest(const RepeatParameter& rp,
                const string& refstr,
                const string& readstr)
    {
        swaligner_.doTest(rp, refstr, readstr);
    }
    
private:
    void saveAlleles(const RepeatParameter& rp,
                     RB_Repeat& repeat,
                     ofstream&,
                     ofstream&,
                     TIndexOffU grp_id,
                     TIndexOffU&);

    void writeAllele(TIndexOffU grp_id, 
                     TIndexOffU allele_id,
                     Range range,
                     const string& seq_name,
                     const EList<SeedExt>& seeds,
                     ostream &fp);


    void writeSNPs(ostream& fp, 
                   const string& rep_chr_name,
                   const ESet<Edit>& editset);

    void writeHaploType(const EList<SeedHP>& haplo_list,
                        const EList<SeedExt>& seeds,
                        TIndexOffU& hapl_id_base,
                        const string& seq_name,
                        ostream& fp);

    void generateHaploType(Range range,
                           const EList<SeedExt>& seeds,
                           EList<SeedHP>& haplo_list);

private:
    const int output_width = 60;
    
    TStr& s_;
    bool forward_only_;
    string filename_;
    
    RB_SubSA subSA_;
    RB_SubSA test_subSA_;
    
    TIndexOffU forward_length_;
    
    CoordHelper coordHelper_;
    
    //
    RB_SWAligner swaligner_;

    // Seeds
    EList<string> consensus_all_;
    ELList<SeedExt> seeds_;
    map<size_t, RB_Repeat*> repeat_map_;
};


int strcmpPos(const string&, const string&, TIndexOffU&);
template<typename TStr> void dump_tstr(TStr& s);

#endif /* __REPEAT_BUILDER_H__ */
