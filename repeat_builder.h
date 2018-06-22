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

using namespace std;

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

struct RepeatGroup {
	string seq;

    EList<RepeatCoord<TIndexOffU> > positions;

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

    void merge(const RepeatGroup& rg, const EList<Edit>& ed)
    {
        merge(rg);
        alt_seq.back().edits = ed;
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

    void writeSNPs(ofstream& fp, const string& rep_chr_name)
    {
        for(size_t i = 0; i < edits.size(); i++) {
            Edit& ed = edits[i];

            assert(ed.isMismatch());

            fp << snpIDs[i];
            fp << "\t" << "single";
            fp << "\t" << rep_chr_name;
            fp << "\t" << base_offset + edits[i].pos;
            fp << "\t" << edits[i].qchr;
            fp << endl;
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

        size_t left = edits[0].pos;
        size_t right = edits[edits.size() - 1].pos;

        fp << "rpht" << base_idx++;
        fp << "\t" << rep_chr_name;
        fp << "\t" << base_offset + left;
        fp << "\t" << base_offset + right;
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

// build Non-repetitive Genome
template<typename TStr>
class NRG {

public:
	NRG();
	NRG(
		const EList<RefRecord>& szs,
		EList<string>& ref_names,
		TStr& s,
		const string& filename,
		BlockwiseSA<TStr>& sa);
    ~NRG();

public:
	const int output_width = 60;

	const EList<RefRecord>& szs_;
	EList<string>& ref_namelines_;
	EList<string> ref_names_;
	TStr& s_;
	string filename_;

	BlockwiseSA<TStr>& bsa_;
    
	// mapping info from joined string to genome
	EList<Fragments> fraglist_;

	//
	EList<RepeatGroup> rpt_grp_;
    
    TIndexOffU half_length_;

	// Fragments Cache
#define CACHE_SIZE_JOINEDFRG	10
	Fragments cached_[CACHE_SIZE_JOINEDFRG];
	int num_cached_ = 0;
	int victim_ = 0;	/* round-robin */

    //
    SimpleFunc scoreMin_;
    SimpleFunc nCeil_;
    SimpleFunc penCanIntronLen_;
    SimpleFunc penNoncanIntronLen_;

    Scoring *sc_;
    SwAligner swa;
    LinkedEList<EList<Edit> > rawEdits_;
    RandomSource rnd_;

    TIndexOffU rpt_matchlen_;
    TIndexOffU rpt_edit_;


public:

	void build(TIndexOffU rpt_len,
               TIndexOffU rpt_cnt,
               bool flagGrouping,
               TIndexOffU rpt_edit,
               TIndexOffU rpt_matchlen);
	void buildNames();
	int mapJoinedOffToSeq(TIndexOffU joined_pos);
	int getGenomeCoord(TIndexOffU joined_pos, string& chr_name, TIndexOffU& pos_in_chr);

	void buildJoinedFragment();

	static bool compareRepeatGroupByJoinedOff(const RepeatGroup& a, const RepeatGroup& b)
	{
		return a.positions[0].joinedOff < b.positions[0].joinedOff;
	}
	void sortRepeatGroup();

    void saveRepeatPositions(ofstream& fp, RepeatGroup& rg);
	void saveFile();
	void saveRepeatSequence();
	void saveRepeatGroup();

	void addRepeatGroup(const string&, const EList<RepeatCoord<TIndexOffU> >&);
    void mergeRepeatGroup();
    void groupRepeatGroup(TIndexOffU rpt_edit);
	void adjustRepeatGroup(bool flagGrouping = false);
    RepeatGroup* findRepeatGroup(const string&);

    TIndexOffU getEnd(TIndexOffU e);
	TIndexOffU getLCP(TIndexOffU a, TIndexOffU b);

	void repeat_masking();
    void init_dyn();

    bool checkSequenceMergeable(const string& s1, const string& s2, 
            EList<Edit>&, TIndexOffU max_edit = 10);
    int alignStrings(const string&, const string&, EList<Edit>&);
};

int strcmpPos(const string&, const string&, TIndexOffU&);
template<typename TStr> void dump_tstr(TStr& s);

#endif /* __REPEAT_BUILDER_H__ */
