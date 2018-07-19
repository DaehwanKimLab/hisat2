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
#include "ds.h"
#include "repeat.h"
#include "blockwise_sa.h"

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

// Dump
//
// to_string
static string to_string(int val)
{
	stringstream ss;
	ss << val;
	return ss.str();
}

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
    int getGenomeCoord(TIndexOffU joined_pos, string& chr_name, TIndexOffU& pos_in_chr);
    
    TIndexOffU getEnd(TIndexOffU e);
    TIndexOffU getStart(TIndexOffU e);
    TIndexOffU forward_length() const { return forward_length_; }
    
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

struct RepeatGroup {
	string seq;
    
    EList<RepeatCoord<TIndexOffU> > positions;

    Coord coord; 
	EList<Edit> edits; 
    EList<string> snpIDs;
    
    EList<RepeatGroup> alt_seq;
    size_t base_offset;

	void merge(const RepeatGroup& rg)
	{
        alt_seq.push_back(rg);
#if 0
		alt_seq.push_back(rg.seq);
		for (int i = 0; i < rg.alt_seq.size(); i++) {
			alt_seq.push_back(rg.alt_seq[i]);
		}

		for (int i = 0; i < rg.positions.size(); i++) {
			positions.push_back(rg.positions[i]);
		}
#endif
	}

    void merge(const RepeatGroup& rg, const EList<Edit>& ed, const Coord& coord)
    {
        merge(rg);
        alt_seq.back().edits = ed;
        alt_seq.back().coord = coord;
    }

	bool empty(void) 
	{ 
		return positions.size() == 0;
	}

	void set_empty(void) 
	{ 
		positions.clear();
        assert(positions.size() == 0);
	}

    void writeSNPs(ostream& fp, const string& rep_chr_name)
    {
        size_t ref_base = base_offset + coord.off();
        int rd_gaps = 0;    // Read Gaps
        int rf_gaps = 0;    // Ref Gaps

        for(size_t i = 0; i < edits.size(); i++) {
            Edit& ed = edits[i];

            assert_geq(edits[i].pos + rd_gaps - rf_gaps, 0);
            if(ed.isMismatch()) {

                fp << snpIDs[i];
                fp << "\t" << "single";
                fp << "\t" << rep_chr_name;
                fp << "\t" << ref_base + edits[i].pos + rd_gaps - rf_gaps;
                fp << "\t" << edits[i].qchr;
                fp << endl;
            } else if(ed.isReadGap()) {

                fp << snpIDs[i];
                fp << "\t" << "deletion";
                fp << "\t" << rep_chr_name;
                fp << "\t" << ref_base + edits[i].pos + rd_gaps - rf_gaps;
                fp << "\t" << 1;    // TODO
                fp << endl;

                rd_gaps++;
            } else if(ed.isRefGap()) {
                fp << snpIDs[i];
                fp << "\t" << "insertion";
                fp << "\t" << rep_chr_name;
                fp << "\t" << ref_base + edits[i].pos + rd_gaps - rf_gaps;
                fp << "\t" << (char)edits[i].qchr;  // TODO
                fp << endl;

                rf_gaps++;
            } else {
                assert(false);
            }
        }
    }

    void buildSNPs(size_t& base_idx)
    {
        snpIDs.clear();
        for(size_t i = 0; i < edits.size(); i++) {
            snpIDs.push_back("rps" + to_string(base_idx++));
        }
    }

    void writeHaploType(ofstream& fp, const string& rep_chr_name, size_t& base_idx)
    {
        assert_gt(edits.size(), 0);
        assert_eq(edits.size(), snpIDs.size());

        size_t ref_base = base_offset + coord.off();

        int rd_gaps = 0;
        int rf_gaps = 0;

        for(size_t i = 0; i < edits.size(); i++) {
            if(edits[i].isReadGap()) {
                rd_gaps++;
            } else if(edits[i].isRefGap()) {
                rf_gaps++;
            }
        }

        size_t left = edits[0].pos;
        size_t right = edits[edits.size() - 1].pos + rd_gaps - rf_gaps;

        fp << "rpht" << base_idx++;
        fp << "\t" << rep_chr_name;
        fp << "\t" << ref_base + left;
        fp << "\t" << ref_base + right;
        fp << "\t";
        for(size_t i = 0; i < edits.size(); i++) {
            if (i != 0) {
                fp << ",";
            }
            fp << snpIDs[i];
        }
        fp << endl;
    }

};

struct SeedExt {
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

    EList<Edit> edits;      // edits w.r.t. consensus_merged
    
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

    static bool isSameConsensus(const SeedExt& a, const SeedExt& b) {
        return (a.consensus_pos == b.consensus_pos)
            && (a.getLength() == b.getLength());
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
};

class RB_Repeat {
public:
    RB_Repeat() {}
    ~RB_Repeat() {}
    
    string& consensus() { return consensus_; }
    EList<SeedExt>& seeds() { return seeds_; }
    EList<Range>& seed_ranges() { return seed_ranges_; }
    
    template<typename TStr>
    void getExtendedSeedSequence(const TStr& s,
                                 const SeedExt& seed,
                                 string& seq);
    
    template<typename TStr>
    void extendConsensus(const RepeatParameter& rp,
                         const TStr& s);
    
    template<typename TStr>
    void saveSeedExtension(const RepeatParameter& rp,
                           const TStr& s,
                           CoordHelper& coordHelper,
                           TIndexOffU seed_grp_id,
                           ostream& fp,
                           size_t& total_repeat_seq_len,
                           size_t& total_allele_seq_len);
    
    bool satisfy(const RepeatParameter& rp)
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
    
private:
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
                           EList<string>* right_consensuses);

private:
    string         consensus_;
    EList<SeedExt> seeds_;
    EList<Range>   seed_ranges_;
};

// find and write repeats
template<typename TStr>
class RepeatBuilder {

public:
	RepeatBuilder();
	RepeatBuilder(TStr& s,
                  const EList<RefRecord>& szs,
                  const EList<string>& ref_names,
                  bool forward_only,
                  BlockwiseSA<TStr>& sa,
                  const string& filename);
    ~RepeatBuilder();


public:
	void build(const RepeatParameter& rp);
    
	static bool compareRepeatGroupByJoinedOff(const RepeatGroup& a, const RepeatGroup& b)
	{
		return a.positions[0].joinedOff < b.positions[0].joinedOff;
	}
	void sortRepeatGroup();

    void saveRepeatPositions(ofstream& fp, RepeatGroup& rg);
	void saveFile(const RepeatParameter& rp);
	void saveRepeatSequence();
	void saveRepeatGroup();

    void addRepeatGroup(const RepeatParameter& rp,
                        map<Range, EList<TIndexOffU> >& range_to_repeatgroup,
                        const string& seed_str,
                        const EList<RepeatCoord<TIndexOffU> >& positions,
                        ostream& fp);
    
    void mergeRepeatGroup();
    void groupRepeatGroup(TIndexOffU rpt_edit);
	void adjustRepeatGroup(bool flagGrouping = false);
    RepeatGroup* findRepeatGroup(const string&);

	TIndexOffU getLCP(TIndexOffU a, TIndexOffU b);

	void repeat_masking();
    void init_dyn(const RepeatParameter& rp);

    bool checkSequenceMergeable(const string& ref,
                                const string& read,
                                EList<Edit>& edits,
                                Coord& coord,
                                TIndexOffU rpt_len,
                                TIndexOffU max_edit = 10);
    int alignStrings(const string&, const string&, EList<Edit>&, Coord&);
    void makePadString(const string&, const string&, string&, size_t);

    void saveSeedEdit(const string& consensus_merged,
                      EList<SeedExt>& seeds,
                      TIndexOffU min_rpt_len,
                      ostream& fp);

    void refineConsensus(const string& seed_string,
                         EList<SeedExt>& seeds,
                         const RepeatParameter& rp,
                         const string& old_consensus,
                         string& refined_consensus,
                         ostream& fp);

    void seedGrouping(const RepeatParameter& rp);

    void doTest(const RepeatParameter& rp,
                const string& refstr,
                const string& readstr);
    void doTestCase1(const string&, const string&, TIndexOffU);
    
private:
    void get_alleles(TIndexOffU grp_id,
                     const string& seq_name, 
                     TIndexOffU baseoff,
                     TIndexOffU& allele_id,
                     TIndexOffU& hapl_id_base,
                     const Range range, 
                     EList<SeedExt>& seeds, 
                     ostream& info_fp,
                     ostream& hapl_fp);

    void updateSeedBaseoff(EList<SeedExt>&, Range, size_t);

    void writeAllele(TIndexOffU grp_id, 
                     TIndexOffU allele_id,
                     TIndexOffU baseoff,
                     Range range,
                     const string& seq_name,
                     EList<SeedExt>& seeds,
                     ostream &fp);

    void writeSNPs(ostream& fp, 
                   TIndexOffU baseoff,
                   const string& rep_chr_name, 
                   const ESet<Edit>& editset);

    void writeHaploType(TIndexOffU& hapl_id_base,
                        TIndexOffU baseoff,
                        Range range,
                        const string& seq_name,
                        EList<SeedExt>& seeds,
                        ostream &fp);

    void saveSeedSNPs(TIndexOffU seed_group_id,
                      TIndexOffU& snp_id_base, 
                      TIndexOffU& hapl_id_base,
                      TIndexOffU baseoff,
                      EList<SeedExt>& seeds,
                      ostream& info_fp,
                      ostream& snp_fp,
                      ostream& hapl_fp);

    void saveSeeds(const RepeatParameter& rp);
    void saveConsensusSequence();
    
private:
    const int output_width = 60;
    TIndexOffU min_repeat_len_;
    
    TStr& s_;
    bool forward_only_;
    string filename_;
    
    BlockwiseSA<TStr>& bsa_;
    
    //
    EList<RepeatGroup> rpt_grp_;
    
    TIndexOffU forward_length_;
    
    CoordHelper coordHelper_;
    
    //
    SimpleFunc scoreMin_;
    SimpleFunc nCeil_;
    SimpleFunc penCanIntronLen_;
    SimpleFunc penNoncanIntronLen_;
    
    Scoring *sc_;
    SwAligner swa;
    LinkedEList<EList<Edit> > rawEdits_;
    RandomSource rnd_;


    // Seeds
    EList<string> consensus_all_;
    ELList<SeedExt> seeds_;
    EList<RB_Repeat> repeats_;

};

int strcmpPos(const string&, const string&, TIndexOffU&);
template<typename TStr> void dump_tstr(TStr& s);

#endif /* __REPEAT_BUILDER_H__ */
