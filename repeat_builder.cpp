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

#include <iostream>
#include <vector>
#include <algorithm>
#include "timer.h"
#include "aligner_sw.h"
#include "aligner_result.h"
#include "scoring.h"
#include "sstring.h"

#include "repeat_builder.h"

unsigned int levenshtein_distance(const std::string& s1, const std::string& s2)
{
    const std::size_t len1 = s1.size(), len2 = s2.size();
    std::vector<unsigned int> col(len2+1), prevCol(len2+1);

    for (unsigned int i = 0; i < prevCol.size(); i++)
        prevCol[i] = i;
    for (unsigned int i = 0; i < len1; i++) {
        col[0] = i+1; 
        for (unsigned int j = 0; j < len2; j++)
            // note that std::min({arg1, arg2, arg3}) works only in C++11,
            // for C++98 use std::min(std::min(arg1, arg2), arg3)
            col[j+1] = std::min( std::min(prevCol[1 + j] + 1, col[j] + 1), prevCol[j] + (s1[i]==s2[j] ? 0 : 1) );
        col.swap(prevCol);
    }
    return prevCol[len2];
}


unsigned int hamming_distance(const string& s1, const string& s2)
{
    assert_eq(s1.length(), s2.length());

    uint32_t cnt = 0;

    for(size_t i = 0; i < s1.length(); i++) {
        if(s1[i] != s2[i]) {
            cnt++;
        }
    }

    return cnt;
}

static const Range EMPTY_RANGE = Range(1, 0);

struct RepeatRange {
	RepeatRange() {
        forward = true;
    };
	RepeatRange(Range r, int id) : 
		range(r), rg_id(id), forward(true) {};
    RepeatRange(Range r, int id, bool fw) :
        range(r), rg_id(id), forward(fw) {};

	Range range;
	int rg_id;
    bool forward;
};

Range reverseRange(const Range& r, TIndexOffU size)
{
	size_t len = r.second - r.first;
	Range rc;

	rc.first = size - r.second;
	rc.second = rc.first + len;

	return rc;
}

string reverseComplement(const string& str)
{
	string rc;
	for(TIndexOffU si = str.length(); si > 0; si--) {
		rc.push_back(asc2dnacomp[str[si - 1]]);
	}
	return rc;
}


string toMDZ(const EList<Edit>& edits, const string& read)
{
    StackedAln stk;
    BTDnaString btread;
    BTString buf;

    btread.install(read.c_str(), true);

    stk.init(btread, edits, 0, 0, 0, 0);
    stk.buildMdz();
    stk.writeMdz(&buf, NULL);

    return string(buf.toZBuf());
}

string applyEdits(const string& ref, size_t read_len, const EList<Edit>& edits, const Coord& coord)
{
    string read;
    size_t ref_pos = coord.off();
    size_t read_pos = 0;

    for(size_t i = 0; i < edits.size(); i++) {
        for(; read_pos < edits[i].pos; read_pos++, ref_pos++) {
            read.push_back(ref[ref_pos]);
        }

        if(edits[i].type == EDIT_TYPE_READ_GAP) {
            // delete on read
            ref_pos++;
        } else if(edits[i].type == EDIT_TYPE_REF_GAP) {
            // insert on read
            read.push_back(edits[i].qchr);
            read_pos++;
        } else if(edits[i].type == EDIT_TYPE_MM) {
            assert_eq(ref[ref_pos], edits[i].chr);
            read.push_back(edits[i].qchr);

            read_pos++;
            ref_pos++;
        } else {
            assert(false);
        }
    }

    for(; read_pos < read_len; read_pos++, ref_pos++) {
        read.push_back(ref[ref_pos]);
    }

    return read;
}

/**
 * @brief return true iff a U b = a 
 *
 * @param a
 * @param b
 *
 * @return 
 */
static bool checkRangeMergeable(const Range& a, const Range& b)
{
	if(a.first <= b.first && a.second >= b.second) {
		return true;
	}

	return false;
}


static bool compareRepeatRangeByRange(const RepeatRange& a, const RepeatRange& b)
{
	if((a.range.second > b.range.second) ||
			(a.range.second == b.range.second && a.range.first < b.range.first)) {
		return true;
	}

	return false;
}

static bool compareRepeatRangeByRgID(const RepeatRange& a, const RepeatRange& b)
{
	return a.rg_id < b.rg_id;
}


template<typename index_t>
static bool compareRepeatCoordByJoinedOff(const RepeatCoord<index_t>& a, const RepeatCoord<index_t>& b)
{
	return a.joinedOff < b.joinedOff;
}


template<typename TStr>
string getString(const TStr& ref, TIndexOffU start, size_t len)
{
    const size_t ref_len = ref.length();
    assert_leq(start + len, ref_len);

    string s;
	for(size_t i = 0; i < len; i++) {
		char nt = "ACGT"[ref[start + i]];
		s.push_back(nt);
	}

	return s;
}


template<typename TStr>
inline uint8_t getSequenceBase(const TStr& ref, TIndexOffU pos)
{
    assert_lt(pos, ref.length());
    return (uint8_t)ref[pos];
}


template<typename TStr>
void dump_tstr(const TStr& s)
{
	static int print_width = 60;

	size_t s_len = s.length();

	for(size_t i = 0; i < s_len; i += print_width) {
		string buf;
		for(size_t j = 0; (j < print_width) && (i + j < s_len); j++) {
			buf.push_back("ACGTN"[s[i + j]]);
		}
		cerr << buf << endl;
	}
	cerr << endl;
}


size_t getMaxMatchLen(const EList<Edit>& edits, const size_t read_len)
{
    size_t max_length = 0;
    uint32_t last_edit_pos = 0;
    uint32_t len = 0;

    if (edits.size() == 0) {
        // no edits --> exact match
        return read_len;
    }

    for(size_t i = 0; i < edits.size(); i++) {
        if (last_edit_pos > edits[i].pos) {
            continue;
        }

        len = edits[i].pos - last_edit_pos;
        if(len > max_length) {
            max_length = len;
        }

        last_edit_pos = edits[i].pos + 1;
    }

    if (last_edit_pos < read_len) {
        len = read_len - last_edit_pos;
        if(len > max_length) {
            max_length = len;
        }
    }

    return max_length;
}

template<typename TStr>
RepeatGenerator<TStr>::RepeatGenerator(TStr& s,
                                       const EList<RefRecord>& szs,
                                       EList<string>& ref_names,
                                       bool forward_only,
                                       BlockwiseSA<TStr>& sa,
                                       const string& filename) :
	s_(s), szs_(szs), forward_only_(forward_only),
    ref_namelines_(ref_names),
    bsa_(sa), filename_(filename),
    forward_length_(forward_only ? s.length() : s.length() / 2)
{
	cerr << "RepeatGenerator: " << filename_ << endl;

    rnd_.init(0);

	// build ref_names_ from ref_namelines_
    buildNames();
    buildJoinedFragment();
}

template<typename TStr>
RepeatGenerator<TStr>::~RepeatGenerator()
{
    if (sc_) {
        delete sc_;
    }
}

#define DMAX std::numeric_limits<double>::max()

template<typename TStr>
void RepeatGenerator<TStr>::init_dyn(const RepeatParameter& rp)
{
    const int MM_PEN = 3;
    // const int MM_PEN = 6;
    const int GAP_PEN_LIN = 2;
    // const int GAP_PEN_LIN = (((MM_PEN) * rpt_edit_ + 1) * 1.0);
    const int GAP_PEN_CON = 4;
    // const int GAP_PEN_CON = (((MM_PEN) * rpt_edit_ + 1) * 1.0);
    const int MAX_PEN = MAX_I16;

    scoreMin_.init(SIMPLE_FUNC_LINEAR, rp.max_edit * MM_PEN * -1.0, 0.0);
    nCeil_.init(SIMPLE_FUNC_LINEAR, 0.0, 0.0);

    penCanIntronLen_.init(SIMPLE_FUNC_LOG, -8, 1);
    penNoncanIntronLen_.init(SIMPLE_FUNC_LOG, -8, 1);

    sc_ = new Scoring(
            DEFAULT_MATCH_BONUS,     // constant reward for match
            DEFAULT_MM_PENALTY_TYPE,     // how to penalize mismatches
            MM_PEN,      // max mm penalty
            MM_PEN,      // min mm penalty
            MAX_PEN,      // max sc penalty
            MAX_PEN,      // min sc penalty
            scoreMin_,       // min score as function of read len
            nCeil_,          // max # Ns as function of read len
            DEFAULT_N_PENALTY_TYPE,      // how to penalize Ns in the read
            DEFAULT_N_PENALTY,           // constant if N pelanty is a constant
            DEFAULT_N_CAT_PAIR,    // whether to concat mates before N filtering


            GAP_PEN_CON, // constant coeff for read gap cost
            GAP_PEN_CON, // constant coeff for ref gap cost
            GAP_PEN_LIN, // linear coeff for read gap cost
            GAP_PEN_LIN, // linear coeff for ref gap cost
            1 /* gGapBarrier */    // # rows at top/bot only entered diagonally
            );
}


template<typename TStr>
void RepeatGenerator<TStr>::doTest(const RepeatParameter& rp,
                                   const string& refstr,
                                   const string& readstr)
{
	TIndexOffU count = 0;

    init_dyn(rp);

    doTestCase1(refstr, readstr, rp.max_edit);
}

template<typename TStr>
void RepeatGenerator<TStr>::build(const RepeatParameter& rp)
{
	TIndexOffU count = 0;
    min_repeat_len_ = rp.min_repeat_len;
    
    init_dyn(rp);

    map<TIndexOffU, TIndexOffU>* seedpos_to_repeatgroup = new map<TIndexOffU, TIndexOffU>;
    
	EList<RepeatCoord<TIndexOffU> > rpt_positions;
	TIndexOffU min_lcp_len = s_.length();

	while(count < s_.length() + 1) {
		TIndexOffU saElt = bsa_.nextSuffix();
		count++;
		
		if(count && (count % 10000000 == 0)) {
			cerr << "SA count " << count << endl;
		}
        
        if(saElt == s_.length()) {
            assert_eq(count, s_.length() + 1);
            break;
        }

		if(rpt_positions.size() == 0) {
            rpt_positions.expand();
			rpt_positions.back().joinedOff = saElt;
		} else {
			TIndexOffU prev_saElt = rpt_positions.back().joinedOff;

			// calculate common prefix length between two text.
			//   text1 is started from prev_saElt and text2 is started from saElt
			int lcp_len = getLCP(prev_saElt, saElt);

			if(lcp_len >= rp.seed_len) {
				rpt_positions.expand();
				rpt_positions.back().joinedOff = saElt;

				if(min_lcp_len > lcp_len) {
					min_lcp_len = lcp_len;
				}
			} else {
				if (rpt_positions.size() >= rp.seed_count) {
                    sort(rpt_positions.begin(), 
                            rpt_positions.begin() + rpt_positions.size(),
                            compareRepeatCoordByJoinedOff<TIndexOffU>);

                    string ss = getString(s_, prev_saElt, min_lcp_len);
					addRepeatGroup(*seedpos_to_repeatgroup, ss, rpt_positions);
				}

				// flush previous positions 
                // clear and expand
				rpt_positions.resize(1);
				rpt_positions.back().joinedOff = saElt;
				min_lcp_len = s_.length();
			}
		}
	}
    
    cerr << "number of seed positions is " << seedpos_to_repeatgroup->size() << endl;
    delete seedpos_to_repeatgroup;
    seedpos_to_repeatgroup = NULL;

    //mergeRepeatGroup();
    
#if 0
    if(rp.max_edit > 0) {
        groupRepeatGroup(rp.max_edit);
    }
#endif

	// we found repeat_group
	cerr << rpt_grp_.size() << " groups found" << endl;


    // Seed extension
    seedGrouping(rp);

}

template<typename TStr>
void RepeatGenerator<TStr>::buildNames()
{
	ref_names_.resize(ref_namelines_.size());
	for(size_t i = 0; i < ref_namelines_.size(); i++) {
		string& nameline = ref_namelines_[i];

		for(size_t j = 0; j < nameline.length(); j++) {
			char n = nameline[j];
			if(n == ' ') {
				break;
			}
			ref_names_[i].push_back(n);
		}
	}
}


string reverse(const string& str)
{
    string rev = str;
    size_t str_len = str.length();

    for(size_t i = 0; i < str_len; i++) {
        rev[i] = str[str_len - i - 1];
    }

    return rev;
}

template<typename TStr>
int RepeatGenerator<TStr>::mapJoinedOffToSeq(TIndexOffU joinedOff)
{

	/* search from cached_list */
	if(num_cached_ > 0) {
		for(size_t i = 0; i < num_cached_; i++) {
			Fragments *frag = &cached_[i];
			if(frag->contain(joinedOff)) {
				return frag->frag_id;
			}
		}
		/* fall through */
	}

	/* search list */
	int top = 0;
	int bot = fraglist_.size() - 1; 
	int pos = 0;

	Fragments *frag = &fraglist_[pos];
	while((bot - top) > 1) {
		pos = top + ((bot - top) >> 1);
		frag = &fraglist_[pos];

		if (joinedOff < frag->joinedOff) {
			bot = pos;
		} else {
			top = pos;
		}
	}

	frag = &fraglist_[top];
	if (frag->contain(joinedOff)) {
		// update cache
		if (num_cached_ < CACHE_SIZE_JOINEDFRG) {
			cached_[num_cached_] = *frag;
			num_cached_++;
		} else {
			cached_[victim_] = *frag;
			victim_++; // round-robin
			victim_ %= CACHE_SIZE_JOINEDFRG;
		}

		return top;
	}

	return -1;
}

template<typename TStr>
int RepeatGenerator<TStr>::getGenomeCoord(TIndexOffU joinedOff,
		string& chr_name, TIndexOffU& pos_in_chr)
{
	int seq_id = mapJoinedOffToSeq(joinedOff);
	if (seq_id < 0) {
		return -1;
	}

	Fragments *frag = &fraglist_[seq_id];
	TIndexOffU offset = joinedOff - frag->joinedOff;

	pos_in_chr = frag->seqOff + offset;
	chr_name = ref_names_[frag->seq_id];

	return 0;
}

template<typename TStr>
void RepeatGenerator<TStr>::buildJoinedFragment()
{
    TIndexOffU acc_joinedOff = 0;
    TIndexOffU acc_seqOff = 0;
    TIndexOffU seq_id = 0;
    TIndexOffU frag_id = 0;
    for(size_t i = 0; i < szs_.size(); i++) {
        const RefRecord& rec = szs_[i];
        if(rec.len == 0) {
            continue;
        }
        if (rec.first) {
            acc_seqOff = 0;
            seq_id++;
        }

        fraglist_.expand();

        fraglist_.back().joinedOff = acc_joinedOff;
        fraglist_.back().length = rec.len;

        fraglist_.back().frag_id = frag_id++;
        fraglist_.back().seq_id = seq_id - 1;
        fraglist_.back().first = rec.first;

        acc_seqOff += rec.off;
        fraglist_.back().seqOff = acc_seqOff;

        acc_joinedOff += rec.len;
        acc_seqOff += rec.len;
    }

	// Add Last Fragment(empty)
    fraglist_.expand();
    fraglist_.back().joinedOff = acc_joinedOff;
    fraglist_.back().length = 0;
	fraglist_.back().seqOff = acc_seqOff;
    fraglist_.back().first = false;
    fraglist_.back().frag_id = frag_id; 
    fraglist_.back().seq_id = seq_id;
}

template<typename TStr>
void RepeatGenerator<TStr>::sortRepeatGroup()
{
	if(rpt_grp_.size() > 0) {
		sort(rpt_grp_.begin(), rpt_grp_.begin() + rpt_grp_.size(), 
                compareRepeatGroupByJoinedOff);
	}
}


template<typename TStr>
void RepeatGenerator<TStr>::saveRepeatPositions(ofstream& fp, RepeatGroup& rg)
{
    EList<RepeatCoord<TIndexOffU> >& positions = rg.positions;

#ifndef NDEBUG
    for(size_t j = 0; j < positions.size(); j += 10) {
        string seq_cmp = getString(s_, positions[j].joinedOff, rg.seq.length());
        assert_eq(rg.seq, seq_cmp);
    }
#endif

    for(size_t j = 0; j < positions.size(); j++) {
        if(positions[j].joinedOff < forward_length_) {
            positions[j].fw = true;
        } else {
            positions[j].joinedOff = s_.length() - positions[j].joinedOff - rg.seq.length();
            positions[j].fw = false;
        }
    }
    positions.sort();

    // Positions
    for(size_t j = 0; j < positions.size(); j++) {
        if(j > 0 && (j % 10 == 0)) {
            fp << endl;
        }

        if(j % 10 != 0) {
            fp << " ";
        }            
        string chr_name;
        TIndexOffU pos_in_chr;
        getGenomeCoord(positions[j].joinedOff, chr_name, pos_in_chr);

        char direction = (positions[j].fw ? '+' : '-');
        fp << chr_name << ":" << pos_in_chr << ":" << direction;
    }
    fp << endl;
}


int max_index(size_t array[4])
{
    int max_idx = 0;

    for(size_t i = 1; i < 4; i++) {
        if(array[max_idx] < array[i]) {
            max_idx = i;
        }
    }

    return max_idx;
}

#if 0
string makeAverageString(EList<SeedExt>& seeds, bool flag_left)
{
    string ext = "";

    size_t ext_len = flag_left ? seeds[0].left_ext.length() : seeds[0].right_ext.length();

    for(size_t i = 0; i < ext_len; i++) {

        char cand = 'A';

        size_t count[4] = {0,};

        for(size_t si = 0; si < seeds.size(); si++) {
            string& seed_ext = flag_left ? seeds[si].left_ext : seeds[si].right_ext;

            uint8_t code = asc2dna[seed_ext[i]];
            count[code]++;
        }

        int max_code = 0;
        cerr << count[0];
        cerr << ", " << count[1];
        cerr << ", " << count[2];
        cerr << ", " << count[3];
        cerr << endl;
        for(size_t si = 1; si < 4; si++) {
            if (count[max_code] < count[si]) {
                max_code = si;
            }
        }

        cand = "ACGT"[max_code];

        ext = ext + cand;
    }

    return ext;
}
#endif

template<typename TStr>
void RepeatGenerator<TStr>::seedGrouping(const RepeatParameter& rp)
{
    string seed_filename = filename_ + ".rep.seed";
    ofstream fp(seed_filename.c_str());
    
    size_t total_rep_seq_len = 0;

    seeds_.clear();
    
    for(size_t i = 0; i < rpt_grp_.size(); i++)
    {
        EList<RepeatCoord<TIndexOffU> >& positions = rpt_grp_[i].positions;
        string seed_str = rpt_grp_[i].seq.substr(0, rp.seed_len);
        
        // DK - debugging purposes
#if 1
        if(positions.size() < 50)
            continue;
        
        if(i == 1172) {
            int kk = 20;
            kk += 20;
        } else {
            // continue;
        }
#endif

        seeds_.expand();
        EList<SeedExt>& seeds = seeds_.back();

        seeds.reserveExact(positions.size());

        for(size_t pi = 0; pi < positions.size(); pi++) {
            TIndexOffU left = positions[pi].joinedOff;
            TIndexOffU right = positions[pi].joinedOff + rp.seed_len;

            seeds.expand();
            seeds.back().reset();
            seeds.back().orig_pos = pair<TIndexOffU, TIndexOffU>(left, right);
            seeds.back().pos = seeds.back().orig_pos;
            seeds.back().bound = pair<TIndexOffU, TIndexOffU>(getStart(left), getEnd(left));
            

        }
        assert(seeds.size() > 0);
        string consensus;
        consensus_all_.expand();
        string& refined_consensus = consensus_all_.back();
        seedExtension(seed_str,
                      seeds,
                      consensus,
                      rp);

        refineConsensus(seed_str,
                        seeds,
                        rp,
                        consensus,
                        refined_consensus,
                        fp);

        saveSeedExtension(seed_str,
                          seeds,
                          rp,
                          i,
                          fp,
                          refined_consensus,
                          total_rep_seq_len);
    }

    assert_eq(consensus_all_.size(), seeds_.size());
#ifndef NDEBUG
    {
        size_t tot_len = 0;
        for(size_t i = 0; i < consensus_all_.size(); i++) {
            tot_len += consensus_all_[i].length();
        }
        assert_eq(total_rep_seq_len, tot_len); 
    }
#endif

    fp << "total repeat sequence length: " << total_rep_seq_len << endl;

    fp.close();
}

template<typename TStr>
size_t calc_edit_dist(const TStr& s,
                      const SeedExt& seed,
                      const SeedExt& seed2,
                      size_t left_ext,
                      size_t right_ext,
                      size_t max_ed)
{
    if(seed.bound.first + left_ext > seed.pos.first ||
       seed.pos.second + right_ext > seed.bound.second ||
       seed2.bound.first + left_ext > seed2.pos.first ||
       seed2.pos.second + right_ext > seed2.bound.second)
        return max_ed + 1;
    
    size_t ed = 0;
    for(size_t i = 0; i < left_ext; i++) {
        int ch = getSequenceBase(s, seed.pos.first - i - 1);
        assert_range(0, 3, ch);
        int ch2 = getSequenceBase(s, seed2.pos.first - i - 1);
        assert_range(0, 3, ch2);
        if(ch != ch2) ed++;
        if(ed > max_ed)
            return ed;
        
    }
    for(size_t i = 0; i < right_ext; i++) {
        int ch = getSequenceBase(s, seed.pos.second + i);
        assert_range(0, 3, ch);
        int ch2 = getSequenceBase(s, seed2.pos.second + i);
        assert_range(0, 3, ch2);
        if(ch != ch2) ed++;
        if(ed > max_ed)
            return ed;
    }
    
    return ed;
}

template<typename TStr>
void calc_edit_dist(const TStr& s,
                    EList<SeedExt>& seeds,
                    size_t sb,
                    size_t se,
                    const string& left_consensus,
                    const string& right_consensus,
                    uint32_t max_ed)
{
    size_t left_ext = left_consensus.length();
    size_t right_ext = right_consensus.length();
    for(size_t i = sb; i < se; i++) {
        SeedExt& seed = seeds[i];
        assert(!seed.done);
        uint32_t ed = 0;
        if(seed.bound.first + left_ext > seed.pos.first ||
           seed.pos.second + right_ext > seed.bound.second) {
            ed = max_ed + 1;
            seed.ed = ed;
            continue;
        }
            
        for(size_t j = 0; j < left_ext; j++) {
            int ch = getSequenceBase(s, seed.pos.first - j - 1);
            assert_range(0, 3, ch);
            if("ACGT"[ch] != left_consensus[j]) ed++;
            if(ed > max_ed) break;
        }
        
        for(size_t j = 0; j < right_ext; j++) {
            int ch = getSequenceBase(s, seed.pos.second + j);
            assert_range(0, 3, ch);
            if("ACGT"[ch] != right_consensus[j]) ed++;
            if(ed > max_ed) break;
        }
        
        seed.ed = ed;
    }
}

template<typename TStr>
void RepeatGenerator<TStr>::get_consensus_seq(EList<SeedExt>& seeds,
                                              size_t sb,
                                              size_t se,
                                              size_t min_left_ext,
                                              size_t min_right_ext,
                                              size_t max_ed,
                                              const RepeatParameter& rp,
                                              EList<size_t>& ed_seed_nums,
                                              EList<string>* left_consensuses,
                                              EList<string>* right_consensuses)
{
    assert_lt(sb, se);
    assert_leq(se, seeds.size());
    if(left_consensuses != NULL) {
        left_consensuses->clear(); left_consensuses->resize(max_ed + 1);
        for(size_t i = 0; i < max_ed + 1; i++) {
            (*left_consensuses)[i].clear();
        }
    }
    if(right_consensuses != NULL) {
        right_consensuses->clear(); right_consensuses->resize(max_ed + 1);
        for(size_t i = 0; i < max_ed + 1; i++) {
            (*right_consensuses)[i].clear();
        }
    }
    
    // cluster sequences
    EList<size_t> belongto;
    belongto.resizeExact(se - sb);
    for(size_t i = 0; i < belongto.size(); i++) belongto[i] = i;
    
    for(size_t i = 0; i + 1 < belongto.size(); i++) {
        for(size_t j = i + 1; j < belongto.size(); j++) {
            if(belongto[j] != j) continue;
            size_t ed = calc_edit_dist(s_,
                                       seeds[sb + i],
                                       seeds[sb + j],
                                       min_left_ext,
                                       min_right_ext,
                                       max_ed + 1);
            if(ed <= max_ed + 1) {
                belongto[j] = belongto[i];
            }
            
        }
    }
    
    // find the maximum group that has the most sequences
    EList<size_t> vote;
    vote.resizeExact(seeds.size());
    vote.fillZero();
    size_t max_group = 0;
    for(size_t i = 0; i < belongto.size(); i++) {
        size_t cur_group = belongto[i];
        assert_lt(cur_group, vote.size());
        vote[cur_group]++;
        if(cur_group != max_group && vote[cur_group] > vote[max_group]) {
            max_group = cur_group;
        }
    }
    
    // reuse vote to store seeds
    EList<size_t>& consensus_group = vote;
    consensus_group.clear();
    for(size_t i = 0; i < belongto.size(); i++) {
        if(belongto[i] == max_group)
            consensus_group.push_back(i);
    }
    
    // update ed
    // extend consensus string
    ed_seed_nums.resize(max_ed + 1); ed_seed_nums.fillZero();
    EList<size_t> next_ed_seed_nums; next_ed_seed_nums.resize(max_ed + 1);
    for(size_t i = sb; i < se; i++) seeds[i].ed = 0;
    size_t seed_ext_len = 0;
    while(seed_ext_len < max(min_left_ext, min_right_ext)) {
        // count base
        size_t l_count[4] = {0, };
        size_t r_count[4] = {0, };

        for(size_t i = 0; i < consensus_group.size(); i++) {
            size_t si = sb + consensus_group[i];
            assert(!seeds[si].done);
            int ch;
            if(seed_ext_len < min_left_ext) {
                ch = getSequenceBase(s_, seeds[si].pos.first - seed_ext_len - 1);
                assert_range(0, 3, ch);
                l_count[ch]++;
            }
            if(seed_ext_len < min_right_ext) {
                ch = getSequenceBase(s_, seeds[si].pos.second + seed_ext_len);
                assert_range(0, 3, ch);
                r_count[ch]++;
            }
        }
        
        // select base to be used for extension
        uint8_t left_ext_base = 0;
        if(seed_ext_len < min_left_ext) left_ext_base = max_index(l_count);
        uint8_t right_ext_base = 0;
        if(seed_ext_len < min_right_ext) right_ext_base= max_index(r_count);
        
        // estimate extended ed
        next_ed_seed_nums.fillZero();
        for(size_t i = sb; i < se; i++) {
            assert(!seeds[i].done);
            uint32_t est_ed = seeds[i].ed;
            if(seed_ext_len < min_left_ext) {
                TIndexOffU left_pos = seeds[i].pos.first - seed_ext_len - 1;
                if(left_pos < seeds[i].bound.first) {
                    est_ed = max_ed + 1;
                } else {
                    int ch = getSequenceBase(s_, left_pos);
                    assert_range(0, 3, ch);
                    if (ch != left_ext_base) {
                        est_ed++;
                    }
                }
            }
            
            if(seed_ext_len < min_right_ext) {
                TIndexOffU right_pos = seeds[i].pos.second + seed_ext_len;
                if(right_pos >= seeds[i].bound.second) {
                    est_ed = max_ed + 1;
                } else {
                    int ch = getSequenceBase(s_, right_pos);
                    assert_range(0, 3, ch);
                    if (ch != right_ext_base) {
                        est_ed++;
                    }
                }
            }
            
            if(est_ed <= max_ed) {
                next_ed_seed_nums[est_ed]++;
            }
        }
        
        for(size_t i = 1; i < next_ed_seed_nums.size(); i++) next_ed_seed_nums[i] += next_ed_seed_nums[i-1];
        if(next_ed_seed_nums[max_ed] < rp.repeat_count) {
            break;
        }
        
        for(size_t i = sb; i < se; i++) {
            assert(!seeds[i].done);
            if(seed_ext_len < min_left_ext) {
                TIndexOffU left_pos = seeds[i].pos.first - seed_ext_len - 1;
                if(left_pos < seeds[i].bound.first) {
                    seeds[i].ed = max_ed + 1;
                } else {
                    int ch = getSequenceBase(s_, left_pos);
                    assert_range(0, 3, ch);
                    if (ch != left_ext_base) {
                        seeds[i].ed++;
                    }
                }
            }
            
            if(seed_ext_len < min_right_ext) {
                TIndexOffU right_pos = seeds[i].pos.second + seed_ext_len;
                if(right_pos >= seeds[i].bound.second) {
                    seeds[i].ed = max_ed + 1;
                } else {
                    int ch = getSequenceBase(s_, right_pos);
                    assert_range(0, 3, ch);
                    if (ch != right_ext_base) {
                        seeds[i].ed++;
                    }
                }
            }
        }
        
        for(size_t i = 0; i < next_ed_seed_nums.size(); i++) {
            if(next_ed_seed_nums[i] < rp.repeat_count)
                continue;
            ed_seed_nums[i] = next_ed_seed_nums[i];
            if(seed_ext_len < min_left_ext) {
                assert(left_consensuses != NULL);
                (*left_consensuses)[i] += (char)("ACGT"[left_ext_base]);
            }
            if(seed_ext_len < min_right_ext) {
                assert(right_consensuses != NULL);
                (*right_consensuses)[i] += (char)("ACGT"[right_ext_base]);
            }
        }
        
        seed_ext_len++;
    }
}

template<typename TStr>
void refineSeed(const TStr& s,
                const string& consensus,
                EList<SeedExt>& seeds,
                TIndexOffU backbone_seed_idx,
                TIndexOffU seed_idx,
                const RepeatParameter& rp)
{
    const SeedExt& backbone_seed = seeds[backbone_seed_idx];
    SeedExt& seed = seeds[seed_idx];
    
    TIndexOffU ext_len = seed.getExtLength();
    TIndexOffU backbone_ext_len = backbone_seed.getExtLength();
    if(ext_len >= backbone_ext_len ||
       seed.leftExtLength() > backbone_seed.leftExtLength() ||
       seed.rightExtLength() > backbone_seed.rightExtLength())
        return;

    string constr = consensus.substr(backbone_seed.baseoff, backbone_ext_len);
    TIndexOffU l_diff = backbone_seed.leftExtLength() - seed.leftExtLength();
    TIndexOffU r_diff = backbone_seed.rightExtLength() - seed.rightExtLength();
    
    if(seed.bound.first + l_diff > seed.pos.first ||
       seed.pos.second + r_diff > seed.bound.second)
        return;
    
    string seqstr = getString(s, seed.pos.first - l_diff, backbone_ext_len);
    assert_eq(constr.length(), seqstr.length());
    
    TIndexOffU l = 0;
    for(; l < l_diff; l++) {
        TIndexOffU str_i = l_diff - l - 1;
        if(constr[str_i] != seqstr[str_i])
            break;
    }
    
    TIndexOffU r = 0;
    for(; r < r_diff; r++) {
        TIndexOffU str_i = l_diff + ext_len + r;
        if(constr[str_i] != seqstr[str_i])
            break;
    }
    
    if(l <= 0 && r <= 0)
        return;
    
    seed.backbone = backbone_seed_idx;
    seed.baseoff = backbone_seed.baseoff + l_diff - l;
    seed.pos.first -= l;
    seed.pos.second += r;
}

struct SeedExtInfo {
    TIndexOffU sb;
    TIndexOffU se;
    TIndexOffU left_ext_len;
    TIndexOffU right_ext_len;
};

bool cmpSeed(const SeedExt& s, const SeedExt& s2) {
    if(s.backbone != s2.backbone)
        return s.backbone < s2.backbone;
    if(s.getExtLength() != s2.getExtLength())
        return s.getExtLength() > s2.getExtLength();
    assert_neq(s.orig_pos.first, s2.orig_pos.second);
    return s.orig_pos.first > s2.orig_pos.second;
}

template<typename TStr>
void RepeatGenerator<TStr>::seedExtension(string& seed_string,
                                          EList<SeedExt>& seeds,
                                          string& consensus_merged,
                                          const RepeatParameter& rp)
{
    const size_t seed_mm = rp.seed_mm;
    TIndexOffU baseoff = 0;
    
    // initialize the first group with the size of seeds
    EList<SeedExtInfo> seed_groups;
    seed_groups.expand();
    seed_groups.back().sb = 0;
    seed_groups.back().se = seeds.size();
    
    EList<string> left_consensuses, right_consensuses;
    EList<size_t> ed_seed_nums, ed_seed_nums2;
    bool first_ext = true;
    while(!seed_groups.empty()) {
        TIndexOffU orig_sb = seed_groups.back().sb;  // seed begin
        TIndexOffU sb = orig_sb;                     // seed begin
        TIndexOffU se = seed_groups.back().se;       // seed end
        bool left_extended = seed_groups.back().left_ext_len > 0;
        assert_lt(sb, se);
        seed_groups.pop_back();
        
        size_t max_group_rep = 0, max_group_num = 0;
        TIndexOffU prev_left_ext_len = seeds[sb].leftExtLength();
        TIndexOffU prev_right_ext_len = seeds[sb].rightExtLength();
    
        for(size_t i = sb; i < se; i++) {
#ifndef NDEBUG
            if(!first_ext) assert(seeds[i].done);
#endif
            seeds[i].done = false;
            seeds[i].ed = 0;
        }
        
        while(se - sb >= rp.repeat_count) {
#ifndef NDEBUG
            for(size_t i = sb; i < se; i++) {
                assert(!seeds[i].done);
            }
#endif
            static const string empty_string = "";
            const string* pleft_consensus = NULL;
            const string* pright_consensus = NULL;
            size_t allowed_seed_mm = seed_mm;
            if(first_ext) {
                get_consensus_seq(seeds,
                                  sb,
                                  se,
                                  rp.extend_unit_len,
                                  rp.extend_unit_len,
                                  rp.seed_mm,
                                  rp,
                                  ed_seed_nums,
                                  &left_consensuses,
                                  &right_consensuses);
                
                pleft_consensus = &left_consensuses[seed_mm];
                pright_consensus = &right_consensuses[seed_mm];
            } else {
                assert_lt(seeds[sb].pos.second - seeds[sb].pos.first, rp.max_repeat_len);
                TIndexOffU max_ext_len = rp.max_repeat_len - (seeds[sb].pos.second - seeds[sb].pos.first);
                if(rp.symmetric_extend) max_ext_len /= 2;
                if(max_ext_len > rp.extend_unit_len)
                    max_ext_len = rp.extend_unit_len;
                
                if(rp.symmetric_extend) {
                    get_consensus_seq(seeds,
                                      sb,
                                      se,
                                      max_ext_len,
                                      max_ext_len,
                                      seed_mm,
                                      rp,
                                      ed_seed_nums,
                                      &left_consensuses,
                                      &right_consensuses);
                  
                    pleft_consensus = &empty_string;
                    pright_consensus = &empty_string;
                    for(int i = (int)seed_mm; i >= 0; i--) {
                        size_t left_extlen = (ed_seed_nums[i] < rp.repeat_count ? 0 : left_consensuses[i].length());
                        ASSERT_ONLY(size_t right_extlen = (ed_seed_nums[i] < rp.repeat_count ? 0 : right_consensuses[i].length()));
                        assert_eq(left_extlen, right_extlen);
                        if(i > 0) {
                            if(left_extlen < max_ext_len)
                                continue;
                        } else {
                            if(left_extlen <= 0)
                                continue;
                        }
                        
                        pleft_consensus = &left_consensuses[i];
                        pright_consensus = &right_consensuses[i];
                        allowed_seed_mm = (size_t)i;
                        assert(pleft_consensus->length() > 0 || pright_consensus->length() > 0);
                        break;
                    }
                } else {
                    get_consensus_seq(seeds,
                                      sb,
                                      se,
                                      max_ext_len,
                                      0,
                                      seed_mm,
                                      rp,
                                      ed_seed_nums,
                                      &left_consensuses,
                                      NULL);
                    
                    get_consensus_seq(seeds,
                                      sb,
                                      se,
                                      0,
                                      max_ext_len,
                                      seed_mm,
                                      rp,
                                      ed_seed_nums2,
                                      NULL,
                                      &right_consensuses);
                    
                    pleft_consensus = &empty_string;
                    pright_consensus = &empty_string;
                    for(int i = (int)seed_mm; i >= 0; i--) {
                        size_t left_extlen = (ed_seed_nums[i] < rp.repeat_count ? 0 : left_consensuses[i].length());
                        size_t right_extlen = (ed_seed_nums2[i] < rp.repeat_count ? 0 : right_consensuses[i].length());
                        if(i > 0) {
                            if(max(left_extlen, right_extlen) < rp.extend_unit_len)
                                continue;
                        } else {
                            if(max(left_extlen, right_extlen) <= 0)
                                continue;
                        }
                        
                        if(left_extlen > right_extlen || (left_extlen == right_extlen && !left_extended)) {
                            pleft_consensus = &left_consensuses[i];
                        } else {
                            pright_consensus = &right_consensuses[i];
                        }
                        allowed_seed_mm = (size_t)i;
                        assert(pleft_consensus->length() > 0 || pright_consensus->length() > 0);
                        break;
                    }
                }
                
                calc_edit_dist(s_,
                               seeds,
                               sb,
                               se,
                               *pleft_consensus,
                               *pright_consensus,
                               allowed_seed_mm);
            }
            
            const string& left_consensus = *pleft_consensus;
            const string& right_consensus = *pright_consensus;
            
            if(left_consensus.empty() && right_consensus.empty())
                break;
            
            string consensus;
            if(!left_consensus.empty()) consensus = reverse(left_consensus);
            if(first_ext) {
                consensus += seed_string;
            } else {
                consensus += consensus_merged.substr(seeds[sb].baseoff,
                                                     seeds[sb].pos.second - seeds[sb].pos.first);
            }
            if(!right_consensus.empty()) consensus += right_consensus;
            consensus_merged += consensus;
            
            // Move up "done" seeds
            size_t num_passed_seeds = 0;
            size_t j = sb;
            for(size_t i = sb; i < se; i++) {
                assert(!seeds[i].done);
                if(seeds[i].ed > allowed_seed_mm) {
                    // reset
                    seeds[i].ed = 0;
                    continue;
                }
                
                // select
                seeds[i].done = true;
                seeds[i].baseoff = baseoff;
                seeds[i].pos.first = seeds[i].pos.first - left_consensus.length();
                seeds[i].pos.second = seeds[i].pos.second + right_consensus.length();
                seeds[i].total_ed += seeds[i].ed;
                assert_eq(seeds[i].pos.second - seeds[i].pos.first, consensus.length());
                seeds[i].backbone = sb;
                assert_geq(i, j);
                if(i > j) {
                    assert(!seeds[j].done);
                    SeedExt temp = seeds[j];
                    seeds[j] = seeds[i];
                    seeds[i] = temp;
                    // Find next "undone" seed
                    j++;
                    while(j < i && seeds[j].done) {
                        j++;
                    }
                    assert(j < seeds.size() && !seeds[j].done);
                } else {
                    j = i + 1;
                }
                num_passed_seeds++;
            }
            
            if(num_passed_seeds >= (first_ext ? rp.seed_count : rp.repeat_count) &&
               num_passed_seeds > max_group_num) {
                max_group_rep = sb;
                max_group_num = num_passed_seeds;
            }
            
            baseoff += consensus.length();
            assert_eq(baseoff, consensus_merged.length());
            sb += num_passed_seeds;
            assert_leq(sb, se);
        } // while(se - sb >= rp.repeat_count)
        
        if(se - sb > 0) {
            for(size_t i = sb; i < se; i++) {
                assert(!seeds[i].done);
                if(first_ext) {
                    seeds[i].baseoff = baseoff;
                } else {
                    if(max_group_num > 0) {
                        seeds[i].backbone = max_group_rep;
                        refineSeed(s_,
                                   consensus_merged,
                                   seeds,
                                   max_group_rep,
                                   i,
                                   rp);
                    }
                }
            }
            if(first_ext) {
                consensus_merged += seed_string;
                baseoff += seed_string.length();
            }
        }
        
        sort(seeds.begin() + orig_sb, seeds.begin() + se, cmpSeed);
        // update backbone
        for(size_t i = orig_sb; i < se;) {
            size_t j = i + 1;
            for(; j < se; j++) {
                if(seeds[i].backbone != seeds[j].backbone) break;
                seeds[j].backbone = i;
            }
            seeds[i].backbone = i;
            i = j;
        }
        
        const TIndexOffU repeat_count = (first_ext ? rp.seed_count : rp.repeat_count);
        for(size_t i = orig_sb; i + repeat_count < se;) {
            if(!seeds[i].done) {
                i++;
                continue;
            }
            size_t j = i + 1;
            for(; j < se; j++) {
                if(!seeds[j].done) break;
                assert_leq(seeds[i].backbone, seeds[j].backbone);
                if(seeds[i].backbone != seeds[j].backbone) break;
                assert_geq(seeds[i].getExtLength(), seeds[j].getExtLength());
                if(seeds[i].getExtLength() > seeds[j].getExtLength()) break;
            }
            if(j - i >= repeat_count) {
                bool further_extend = true;
                if(first_ext) {
                    if(seeds[i].getExtLength() < rp.min_repeat_len) {
                        further_extend = false;
                    }
                }
                if(seeds[i].getExtLength() >= rp.max_repeat_len) {
                    further_extend = false;
                }
                if(further_extend) {
                    seed_groups.expand();
                    seed_groups.back().sb = i;
                    seed_groups.back().se = j;
                    assert_leq(prev_left_ext_len, seeds[i].leftExtLength());
                    seed_groups.back().left_ext_len = seeds[i].leftExtLength() - prev_left_ext_len;
                    assert_leq(prev_right_ext_len, seeds[i].rightExtLength());
                    seed_groups.back().right_ext_len = seeds[i].rightExtLength() - prev_right_ext_len;
                    
                }
            }
            i = j;
        }
        
        first_ext = false;
    } // while(!seed_groups.empty())
    
#ifndef NDEBUG
    // make sure seed positions are unique
    EList<pair<TIndexOffU, TIndexOffU> > seed_poses;
    for(size_t i = 0; i < seeds.size(); i++) {
        seed_poses.expand();
        seed_poses.back() = seeds[i].orig_pos;
    }
    seed_poses.sort();
    for(size_t i = 0; i + 1 < seed_poses.size(); i++) {
        if(seed_poses[i].first == seed_poses[i+1].first) {
            assert_lt(seed_poses[i].second, seed_poses[i+1].second);
        } else {
            assert_lt(seed_poses[i].first, seed_poses[i+1].first);
        }
    }
#endif
    
    // update seed dependency
    for(size_t i = 0; i < seeds.size(); i++) {
        size_t idep = i;
        while(idep != seeds[idep].backbone) {
            idep = seeds[idep].backbone;
        }
        seeds[i].backbone = idep;
    }
}


#if 0
template<typename TStr>
void RepeatGenerator<TStr>::saveSeedEdit(
        const string& consensus_merged,
        EList<SeedExt>& seeds,
        TIndexOffU min_rpt_len,
        ostream& fp)
{
    // TODO
    // check ESet performance 
    ESet<Edit> editset;

    size_t org_edit_count = 0;

    // get unique edit list (wrt consensus_merged)
    for(size_t i = 0; i < seeds.size(); i++) {
        size_t ext_len = seeds[i].getExtLength();
        if (ext_len < min_rpt_len) continue;

        for(size_t j = 0; j < seeds[i].edits.size(); j++) {
            editset.insert(seeds[i].edits[j]);
            org_edit_count++;
        }
    }

    fp << "org_edit_count: " << org_edit_count << endl;
    fp << "uniq edit_count: " << editset.size() << endl;

    // assign id
    for(size_t i = 0; i < editset.size(); i++) {
        editset[i].snpID = i;
    }

    for(size_t i = 0; i < seeds.size(); i++) {
        size_t ext_len = seeds[i].getExtLength();
        if (ext_len < min_rpt_len) continue;

        for(size_t j = 0; j < seeds[i].edits.size(); j++) {
            const Edit& e = editset.get(seeds[i].edits[j]);
            seeds[i].edits[j].snpID = e.snpID;
        }
    }

    // write SNPs
    writeSNPs(fp, "rep", editset);

    // group by baseoff
    for(size_t sb = 0; sb < seeds.size(); ) {
        size_t se = sb + 1;
        size_t count = 1;
        for(; se < seeds.size(); se++) {
            if(!SeedExt::isSameConsensus(seeds[sb], seeds[se])) {
                break;
            }
            count++;
        }

        // [sb, se) has same baseoff
        size_t ext_len = seeds[sb].getExtLength();
        cerr << sb << ", " << se << ", " << count << endl;
        // DK
        if(ext_len < min_rpt_len) {
            // go baseoff
            sb = se;
            continue;
        }

        // sort by edits
        sort(seeds.begin() + sb, seeds.begin() + se, SeedExt::sort_by_edits);


        // write SNP, Haplotype
        // group by allele
        get_alleles(sb, se, seeds, fp);

        sb = se;
    }
}
#endif

template<typename TStr>
void RepeatGenerator<TStr>::get_alleles(TIndexOffU grp_id,
                                        const string& seq_name, 
                                        TIndexOffU baseoff,
                                        TIndexOffU& allele_id,
                                        TIndexOffU& hapl_id_base,
                                        const Range range, 
                                        EList<SeedExt>& seeds, 
                                        ostream& info_fp,
                                        ostream& hapl_fp)
{
    // [ragne.first, range.last)

    for(size_t sb = range.first; sb < range.second;) {
        size_t se = sb + 1;
        for(; se < range.second; se++) {
            if(!SeedExt::isSameAllele(seeds[sb], seeds[se])) {
                break;
            }
        }

        // [sb, se) are same allele

        writeAllele(grp_id, allele_id, baseoff, Range(sb, se), seq_name, seeds, info_fp);
        writeHaploType(hapl_id_base, baseoff, Range(sb, se), seq_name, seeds, hapl_fp);

        allele_id++;
        sb = se;
    }
    
}


template<typename TStr>
void RepeatGenerator<TStr>::updateSeedBaseoff(EList<SeedExt>& seeds, Range range, size_t diff)
{
    for(size_t i = range.first; i < range.second; i++) {
        seeds[i].baseoff += diff;
    }
}


#if 1
template<typename TStr>
void RepeatGenerator<TStr>::saveConsensusSequence(const EList<bool>& filter_out)
{
    string fname = filename_ + ".rep.fa";
    ofstream fp(fname.c_str());

	/* TODO */
	fp << ">" << "rep" << endl;

	int oskip = 0;

    size_t acc_len = 0;

    assert_eq(consensus_all_.size(), seeds_.size());

	for(size_t idx = 0; idx < consensus_all_.size(); idx++) {
        if(filter_out[idx]) continue;
        string& constr = consensus_all_[idx];
        size_t constr_len = constr.length();
        acc_len += constr_len;

		size_t si = 0;
		while(si < constr_len) {
			size_t out_len = std::min((size_t)(output_width - oskip), (size_t)(constr_len - si));

			fp << constr.substr(si, out_len);

			if((oskip + out_len) == output_width) {
				fp << endl;
				oskip = 0;
			} else {
				// last line
				oskip = oskip + out_len;
			}

			si += out_len;
		}
	}
	if(oskip) {
		fp << endl;
	}

	fp.close();
}


template<typename TStr>
void RepeatGenerator<TStr>::saveSeedSNPs(TIndexOffU seed_group_id,
                                         TIndexOffU& snp_id_base, 
                                         TIndexOffU& hapl_id_base,
                                         TIndexOffU baseoff,
                                         EList<SeedExt>& seeds,
                                         ostream& info_fp,
                                         ostream& snp_fp,
                                         ostream& hapl_fp)
{
    const string rep_basename = "rep";

    // TODO
    // check ESet performance 
    ESet<Edit> editset;

    size_t org_edit_count = 0;

    // get unique edit list
    for(size_t i = 0; i < seeds.size(); i++) {
        size_t ext_len = seeds[i].getExtLength();
        if (ext_len < min_repeat_len_) continue;

        for(size_t j = 0; j < seeds[i].edits.size(); j++) {
            editset.insert(seeds[i].edits[j]);
            org_edit_count++;
        }
    }

    //cerr << "org_edit_count: " << org_edit_count << endl;
    //cerr << "uniq edit_count: " << editset.size() << endl;

    // assign id
    for(size_t i = 0; i < editset.size(); i++) {
        editset[i].snpID = snp_id_base + i;
    }
    for(size_t i = 0; i < seeds.size(); i++) {
        size_t ext_len = seeds[i].getExtLength();
        if (ext_len < min_repeat_len_) continue;

        for(size_t j = 0; j < seeds[i].edits.size(); j++) {
            const Edit& e = editset.get(seeds[i].edits[j]);
            seeds[i].edits[j].snpID = e.snpID;
        }
    }

    // write SNPs
    writeSNPs(snp_fp, baseoff, rep_basename, editset);

    TIndexOffU allele_id = 0;

    // group by baseoff
    for(size_t sb = 0; sb < seeds.size(); ) {
        size_t se = sb + 1;
        size_t count = 1;
        for(; se < seeds.size(); se++) {
            if(!SeedExt::isSameConsensus(seeds[sb], seeds[se])) {
                break;
            }
            count++;
        }

        // [sb, se) has same baseoff
        size_t ext_len = seeds[sb].getExtLength();

        if(ext_len < min_repeat_len_) {
            // go baseoff
            sb = se;
            continue;
        }

        // sort by edits
        sort(seeds.begin() + sb, seeds.begin() + se, SeedExt::sort_by_edits);

        // write SNP, Haplotype
        // group by allele
        get_alleles(seed_group_id, 
                    rep_basename, 
                    baseoff, 
                    allele_id,
                    hapl_id_base, 
                    Range(sb, se), 
                    seeds, 
                    info_fp, 
                    hapl_fp);

        sb = se;
    }

    // Update snpIDBase
    snp_id_base += editset.size();
}

void filter_seeds(const ELList<SeedExt>& seedss,
                  EList<bool>& filter_out) {
    if(seedss.size() <= 0)
        return;
    filter_out.resizeExact(seedss.size());
    filter_out.fill(false);
    
    EList<TIndexOffU> positions, positions2;
    map<TIndexOffU, TIndexOffU> seedpos_to_group;
    
    // skip if a given set of seeds corresponds to an existing set of seeds
    const TIndexOffU pos_diff = 50;
    const size_t sampling = 10;
    for(size_t s = 0; s < seedss.size(); s++) {
        const EList<SeedExt>& seeds = seedss[s];
        bool add = true;
        positions.clear();
        for(size_t i = 0; i < seeds.size(); i++) {
            positions.push_back(seeds[i].pos.first);
        }
        positions.sort();
        
        for(size_t i = 0; i < seeds.size(); i += sampling) {
            TIndexOffU joinedOff = seeds[i].pos.first;
            map<TIndexOffU, TIndexOffU>::iterator it = seedpos_to_group.lower_bound(joinedOff >= pos_diff ? joinedOff - pos_diff : 0);
            for(; it != seedpos_to_group.end(); it++) {
                if(it->first > joinedOff + pos_diff) {
                    break;
                }
                assert_geq(it->first + pos_diff, joinedOff);
                assert_lt(it->second, seedss.size());
                size_t s2 = it->second;
                const EList<SeedExt>& seeds2 = seedss[s2];
                
                positions2.clear();
                for(size_t i2 = 0; i2 < seeds2.size(); i2++) {
                    positions2.push_back(seeds2[i2].pos.first);
                }
                positions2.sort();
                
                size_t num_match = 0;
                size_t p = 0, p2 = 0;
                while(p < positions.size() && p2 < positions2.size()) {
                    TIndexOffU pos = positions[p], pos2 = positions2[p2];
                    if(pos + pos_diff >= pos2 && pos2 + pos_diff >= pos) {
                        num_match++;
                    }
                    if(pos <= pos2) p++;
                    else            p2++;
                }
                
                // if the number of matches is >= 10% of positions in the smaller group
                if(num_match * 10 >= min(positions.size(), positions2.size()) * 1) {
                    if(seeds.size() <= seeds.size()) {
                        add = false;
                        filter_out[s] = true;
                        break;
                    } else {
                        filter_out[s2] = true;
                        for(size_t p2 = 0; p2 < seeds2.size(); p2++) {
                            if(seedpos_to_group[positions[p2]] == s2) {
                                seedpos_to_group.erase(positions[p2]);
                            }
                        }
                    }
                }
            }
        }
        if(add) {
            for(size_t i = 0; i < positions.size(); i++) {
                seedpos_to_group[positions[i]] = s;
            }
        }
    }
}


template<typename TStr>
void RepeatGenerator<TStr>::saveSeeds()
{
	string rptinfo_filename = filename_ + ".rep.info";
    string snp_filename = filename_ + ".rep.snp";
    string hapl_filename = filename_ + ".rep.haplotype";
    
	ofstream info_fp(rptinfo_filename.c_str());
    ofstream snp_fp(snp_filename.c_str());
    ofstream hapl_fp(hapl_filename.c_str());

    TIndexOffU snp_id_base = 0;
    TIndexOffU hapl_id_base = 0;
    TIndexOffU baseoff = 0;
    
    EList<bool> filter_out;
    filter_seeds(seeds_, filter_out);

    assert_eq(seeds_.size(), consensus_all_.size());
    for(size_t i = 0; i < seeds_.size(); i++) {
        if(filter_out[i]) continue;
        saveSeedSNPs(i, snp_id_base, hapl_id_base, baseoff, seeds_[i], info_fp, snp_fp, hapl_fp);
        baseoff += consensus_all_[i].length();
    }

    info_fp.close();
    snp_fp.close();
    hapl_fp.close();

    saveConsensusSequence(filter_out);
    
    // DK - debugging purposes
#if 1
    string test_filename = filename_ + ".rep.test";
    ofstream test_fp(test_filename.c_str());

    for(size_t i = 0; i < seeds_.size(); i++) {
        test_fp << i << '\t' << seeds_[i].size() << '\t' << seeds_[i][0].pos.first << '\t' << (filter_out[i] ? "filtered_out" : "filtered_in") << endl;
    }
    
    test_fp.close();
#endif
}

#endif

template<typename TStr>
void RepeatGenerator<TStr>::writeHaploType(TIndexOffU& hapl_id_base,
                                           TIndexOffU baseoff,
                                           Range range,
                                           const string& seq_name,
                                           EList<SeedExt>& seeds,
                                           ostream &fp)
{
    EList<Edit>& edits = seeds[range.first].edits;
    if(edits.size() == 0)
        return;
    
    // right-most position
    TIndexOffU max_right_pos = seeds[range.first].getExtLength() - 1;

    // create haplotypes of at least 16 bp long to prevent combinations of SNPs
    // break a list of SNPs into several haplotypes if a SNP is far from the next SNP in the list
    const TIndexOffU min_ht_len = 16;
    size_t eb = 0, ee = 1;
    while(ee < edits.size() + 1) {
        if(ee == edits.size() ||
           edits[eb].pos + (min_ht_len << 1) < edits[ee].pos) {
            TIndexOffU left_pos = edits[eb].pos + baseoff;
            TIndexOffU right_pos = (ee == edits.size() ? edits[ee-1].pos : edits[ee].pos);
            right_pos = min<TIndexOffU>(max_right_pos, right_pos + min_ht_len) + baseoff;
            fp << "rpht" << hapl_id_base++;
            fp << "\t" << seq_name;
            fp << "\t" << left_pos;
            fp << "\t" << right_pos;
            fp << "\t";
            for(size_t i = eb; i < ee; i++) {
                if(i > eb) {
                    fp << ",";
                }
                fp << "rps" << edits[i].snpID;
            }
            fp << endl;
            eb = ee;
        }
        ee++;
    }
}

template<typename TStr>
void RepeatGenerator<TStr>::writeAllele(TIndexOffU grp_id, 
                                        TIndexOffU allele_id,
                                        TIndexOffU baseoff,
                                        Range range,
                                        const string& seq_name,
                                        EList<SeedExt>& seeds,
                                        ostream &fp)
{
    // >rpt_name*0\trep\trep_pos\trep_len\tpos_count\t0
    // chr_name:pos:direction chr_name:pos:direction
    //
    // >rep1*0	rep	0	100	470	0
    // 22:112123123:+ 22:1232131113:+
    //
    size_t snp_size = seeds[range.first].edits.size();
    size_t pos_size = range.second - range.first;
    fp << ">";
    fp << "rpt_" << grp_id << "*" << allele_id;
    fp << "\t" << seq_name;
    fp << "\t" << seeds[range.first].baseoff + baseoff;
    fp << "\t" << seeds[range.first].getExtLength();
    fp << "\t" << range.second - range.first;
    fp << "\t" << snp_size;

    fp << "\t";
    for(size_t i = 0; i < snp_size; i++) {
        if(i > 0) {fp << ",";}
        fp << "rps" << seeds[range.first].edits[i].snpID;
    }
    fp << endl;

    // print positions
    for(size_t i = 0; i < pos_size; i++) {
        if(i > 0 && (i % 10 == 0)) {
            fp << endl;
        }

        if(i % 10 != 0) {
            fp << " ";
        }            
        string chr_name;
        TIndexOffU pos_in_chr;
        TIndexOffU joinedOff = seeds[range.first + i].pos.first;

        bool fw = true;
        if(joinedOff < forward_length_) {
            fw = true;
        } else {
            fw = false;
            joinedOff = s_.length() - joinedOff - seeds[range.first].getExtLength();
        }

        getGenomeCoord(joinedOff, chr_name, pos_in_chr);

        char direction = fw ? '+' : '-';
        fp << chr_name << ":" << pos_in_chr << ":" << direction;
    }

    fp << endl;
}

template<typename TStr>
void RepeatGenerator<TStr>::writeSNPs(ostream& fp, 
                                      TIndexOffU baseoff,
                                      const string& rep_chr_name, 
                                      const ESet<Edit>& editset)
{
    for(size_t i = 0; i < editset.size(); i++) {
        const Edit& ed = editset[i];

        if(ed.isMismatch()) {
            fp << "rps" << ed.snpID;
            fp << "\t" << "single";
            fp << "\t" << rep_chr_name;
            fp << "\t" << ed.pos + baseoff;
            fp << "\t" << ed.qchr;
            fp << endl;
        } else {
            assert(false);
        }
    }
}

template<typename TStr>
void RepeatGenerator<TStr>::refineConsensus(const string& seed_string,
                                            EList<SeedExt>& seeds,
                                            const RepeatParameter& rp,
                                            const string& old_consensus,
                                            string& refined_consensus,
                                            ostream& fp)
{
#ifndef NDEBUG
    size_t tot_exp_len = 0;
    for(size_t i = 0; i < seeds.size(); i++) {
        if(seeds[i].backbone == i) {
            size_t ext_len = seeds[i].getExtLength();
            if (ext_len < rp.min_repeat_len) continue;
            tot_exp_len += ext_len;
        }
    }
#endif
    
    refined_consensus.clear();
    TIndexOffU refined_baseoff = 0;

    // group by baseoff
    for(size_t sb = 0; sb < seeds.size(); ) {
        size_t se = sb + 1;
        size_t count = 1;
        for(; se < seeds.size(); se++) {
            if(!SeedExt::isSameConsensus(seeds[sb], seeds[se])) {
                break;
            }
            count++;
        }

        // [sb, se) has same baseoff
        size_t ext_len = seeds[sb].getExtLength();
        string constr = old_consensus.substr(seeds[sb].baseoff, ext_len);

        if(ext_len < rp.min_repeat_len) {
            // go baseoff
            sb = se;
            continue;
        }

        // check backbone
        if(seeds[sb].backbone == sb) {
            // this group has own consensus string

            // update baseoff
            for(size_t i = sb; i < se; i++) {
                seeds[i].baseoff = refined_baseoff;
            }

            refined_consensus += constr;
            refined_baseoff += ext_len;

        } else {
            // This group's consensus is located on backbone's one
            // Reuse backbone's consensus
            TIndexOffU backbone = seeds[sb].backbone;
            assert_leq(backbone, sb);

            TIndexOffU backbone_baseoff = seeds[backbone].baseoff;
            size_t backbone_leftext = seeds[backbone].leftExtLength();
            size_t leftext = seeds[sb].leftExtLength();
            assert_geq(backbone_leftext, leftext);
            size_t diff = backbone_leftext - leftext;

#ifndef NDEBUG
            {
                string bb_con = refined_consensus.substr(backbone_baseoff + diff, ext_len);
                string cons = old_consensus.substr(seeds[sb].baseoff, ext_len);
                assert_eq(bb_con, cons);
            }
#endif

            // update baseoff
            for(size_t i = sb; i < se; i++) {
                seeds[i].baseoff = backbone_baseoff + diff; 
            }
        }

        for(size_t i = sb; i < se; i++) {
            string deststr = getString(s_, seeds[i].pos.first, ext_len);
            seeds[i].generateEdits(refined_consensus, deststr);
        }

        sb = se;
    }

    assert_eq(refined_consensus.length(), tot_exp_len);
}

template<typename TStr>
void RepeatGenerator<TStr>::saveSeedExtension(const string& seed_string,
                                              const EList<SeedExt>& seeds,
                                              const RepeatParameter& rp,
                                              TIndexOffU rpt_grp_id,
                                              ostream& fp,
                                              const string& consensus_merged,
                                              size_t& total_repeat_seq_len)
{
    // apply color, which is compatible with linux commands such as cat and less -r
#if 1
    const string red = "\033[31m", reset = "\033[0m";
#else
    const string red = "", reset = "";
#endif
    
    TIndexOffU max_left_ext_len = 0, max_right_ext_len = 0;
    for(size_t i = 0; i < seeds.size(); i++) {
        TIndexOffU left_ext_len = seeds[i].leftExtLength();
        if(max_left_ext_len < left_ext_len) {
            max_left_ext_len = left_ext_len;
        }
        assert_leq(seeds[i].orig_pos.second, seeds[i].pos.second);
        TIndexOffU right_ext_len = seeds[i].rightExtLength();
        if(max_right_ext_len < right_ext_len) {
            max_right_ext_len = right_ext_len;
        }
    }
    
    TIndexOffU prev_consensus_baseoff = (TIndexOffU)-1;
    size_t count = 0, total_count = 0;
    for(size_t i = 0; i < seeds.size(); i++) {
        // assert(seeds[i].done);
        size_t ext_len = seeds[i].getExtLength();
        if(ext_len < rp.min_repeat_len) continue;
        total_count++;
        
        string deststr = getString(s_, seeds[i].pos.first, ext_len);
        if(prev_consensus_baseoff == seeds[i].baseoff) {
            count++;
        } else {
            if(count > 0) fp << setw(5) << count << "  " << setw(5) << total_count << endl << endl;
            count = 1;
            prev_consensus_baseoff = seeds[i].baseoff;
        }

        fp << rpt_grp_id << "  " << rpt_grp_[rpt_grp_id].positions.size();
        fp << "  " << setw(4) << i << " -> " << setw(4) << std::left << seeds[i].backbone << std::right;
        fp << "  " << setw(4) << ext_len;
        fp << "  " << setw(3) << seeds[i].total_ed;
        fp << "  " << setw(5) << seeds[i].baseoff;
        fp << "  " << setw(10) << seeds[i].pos.first << "  " << setw(10) << seeds[i].pos.second;
        fp << "  " << setw(10) << seeds[i].bound.first << "  " << setw(10) << seeds[i].bound.second;

        string chr_name;
        TIndexOffU pos_in_chr;
        if(seeds[i].pos.first < forward_length_) {
            getGenomeCoord(seeds[i].pos.first, chr_name, pos_in_chr);
        } else {
            getGenomeCoord(s_.length() - seeds[i].pos.first - (seeds[i].pos.second - seeds[i].pos.first), chr_name, pos_in_chr);
        }
        fp << "  " << chr_name << ":" << setw(10) << std::left << pos_in_chr << std::right;

        if(seeds[i].pos.second < forward_length_) {
            getGenomeCoord(seeds[i].pos.second, chr_name, pos_in_chr);
        } else {
            getGenomeCoord(s_.length() - seeds[i].pos.second - (seeds[i].pos.second - seeds[i].pos.first), chr_name, pos_in_chr);
        }
        fp << "  " << setw(5) << chr_name << ":" << pos_in_chr;
        string constr = consensus_merged.substr(seeds[i].baseoff, ext_len);
        if(seeds[i].backbone == i) total_repeat_seq_len += constr.length();

        // fp << "\t" << constr;
        assert_eq(constr.length(), deststr.length());
        fp << "  ";
        
        // print sequence w.r.t. the current group
        TIndexOffU left_ext_len = seeds[i].leftExtLength();
        assert_leq(left_ext_len, max_left_ext_len);
        TIndexOffU left_indent = max_left_ext_len - left_ext_len;
        for(size_t j = 0; j < left_indent; j++) fp << ' ';
        for(size_t j = 0; j < constr.length(); j++) {
            bool different = (constr[j] != deststr[j]);
            if(different) fp << red;
            fp << deststr[j];
            if(different) fp << reset;
        }

#if 0
        fp << "\t";
        for(size_t ei = 0; ei < seeds[i].edits.size(); ei++) {
            const Edit& edit = seeds[i].edits[ei];
            if(ei > 0) fp << ",";
            fp << edit;
            if (edit.snpID != std::numeric_limits<uint32_t>::max()) {
                fp << "@" << edit.snpID;
            }
        }
#endif
        
        // print sequence w.r.t. the backbone sequence
        if(rp.symmetric_extend) {
            TIndexOffU right_ext_len = seeds[i].rightExtLength();
            assert_leq(right_ext_len, max_right_ext_len);
            TIndexOffU right_indent = max_right_ext_len - right_ext_len;
            for(size_t j = 0; j < right_indent; j++) fp << ' ';
            fp << "  ";

            TIndexOffU max_backbone = 0, max_ext_len = 0;
            for(int j = (int)seeds[i].backbone; j >= 0; j--) {
                TIndexOffU backbone_ext_len = seeds[j].getExtLength();
                if(backbone_ext_len > max_ext_len) {
                    max_backbone = (TIndexOffU)j;
                    max_ext_len = backbone_ext_len;
                }
            }
            const SeedExt& backbone_seed = seeds[max_backbone];
            size_t backbone_ext_len = backbone_seed.getExtLength();
            string backbone_constr = consensus_merged.substr(backbone_seed.baseoff, backbone_ext_len);
            TIndexOffU backbone_left_ext_len = backbone_seed.leftExtLength();
            assert_leq(backbone_left_ext_len, max_left_ext_len);
            TIndexOffU backbone_left_indent = max_left_ext_len - backbone_left_ext_len;
            assert_leq(left_ext_len, backbone_left_ext_len);
            if(seeds[i].pos.first >= backbone_left_ext_len - left_ext_len) {
                string backbone_deststr = getString(s_, seeds[i].pos.first - (backbone_left_ext_len - left_ext_len), backbone_ext_len);
                assert_eq(backbone_constr.length(), backbone_deststr.length());
                for(size_t j = 0; j < backbone_left_indent; j++) fp << ' ';
                for(size_t j = 0; j < backbone_constr.length(); j++) {
                    bool different = (backbone_constr[j] != backbone_deststr[j]);
                    if(different) fp << red;
                    fp << backbone_deststr[j];
                    if(different) fp << reset;
                }
            }
        }
        
        fp << endl;
    }
    
    if(count > 0) fp << setw(5) << count << "  " << setw(5) << total_count << endl << endl;
}


template<typename TStr>
void RepeatGenerator<TStr>::saveRepeatGroup()
{
    const string rep_basename = "rep";
	string rptinfo_filename = filename_ + ".rep.info";
    string snp_filename = filename_ + ".rep.snp";
    string hapl_filename = filename_ + ".rep.haplotype";

	int rpt_count = rpt_grp_.size();
	TIndexOffU acc_pos = 0;
    size_t snp_base_idx = 0;
    size_t hapl_base_idx = 0;

	ofstream fp(rptinfo_filename.c_str());
    ofstream snp_fp(snp_filename.c_str());
    ofstream hapl_fp(hapl_filename.c_str());

	for(size_t i = 0; i < rpt_count; i++) {
		RepeatGroup& rg = rpt_grp_[i];
		EList<RepeatCoord<TIndexOffU> >& positions = rg.positions;

        size_t allele_count = rg.alt_seq.size();
        

		// >rpt_name*0\trep\trep_pos\trep_len\tpos_count\t0
		// chr_name:pos:direction chr_name:pos:direction
		//
		// >rep1*0	rep	0	100	470	0
		// 22:112123123:+ 22:1232131113:+
		//

		// Header line
		fp << ">" << "rpt_" << i << "*0";
		fp << "\t" << rep_basename; // TODO
		fp << "\t" << rg.base_offset;
		fp << "\t" << rg.seq.length();
		fp << "\t" << positions.size();
		fp << "\t" << "0"; 
        // debugging
        fp << "\t" << rg.seq.substr(0, 50);
		fp << endl;

        saveRepeatPositions(fp, rg);


        //save SNPs 
        for(size_t j = 0; j < allele_count; j++) {
            RepeatGroup& alt_rg = rg.alt_seq[j];

            alt_rg.base_offset = rg.base_offset;

            assert_gt(alt_rg.edits.size(), 0);
#ifndef NDEBUG
            {
                string seq_cmp = applyEdits(rg.seq, alt_rg.seq.length(), alt_rg.edits, alt_rg.coord); 
                assert_eq(alt_rg.seq, seq_cmp);
            }
#endif

            // save snps
            alt_rg.buildSNPs(snp_base_idx);
            alt_rg.writeSNPs(snp_fp, rep_basename);
            alt_rg.writeHaploType(hapl_fp, rep_basename, hapl_base_idx);

            // Header line
            fp << ">" << "rpt_" << i << "*" << (j + 1);
            fp << "\t" << "rep"; // TODO
            fp << "\t" << rg.base_offset;
            fp << "\t" << rg.seq.length();
            fp << "\t" << alt_rg.positions.size();
            fp << "\t" << alt_rg.edits.size(); 
            fp << "\t";
            for(size_t snp_idx = 0; snp_idx < alt_rg.edits.size(); snp_idx++) {
                if (snp_idx != 0) { 
                    fp << ",";
                }
                fp << alt_rg.snpIDs[snp_idx];
            }
            // debugging
            // fp << "\t" << rg.seq;
            fp << endl;

            saveRepeatPositions(fp, alt_rg);
        }
	}		
	fp.close();
    snp_fp.close();
    hapl_fp.close();
}
	
	
template<typename TStr>
void RepeatGenerator<TStr>::saveRepeatSequence()
{
	string fname = filename_ + ".rep.fa";

	ofstream fp(fname.c_str());

	/* TODO */
	fp << ">" << "rep" << endl;

	int oskip = 0;

    size_t acc_len = 0;

	for(TIndexOffU grp_idx = 0; grp_idx < rpt_grp_.size(); grp_idx++) {
		RepeatGroup& rg = rpt_grp_[grp_idx];

        rg.base_offset = acc_len;
        size_t seq_len = rg.seq.length();
        acc_len += seq_len;

		TIndexOffU si = 0;
		while(si < seq_len) {
			size_t out_len = std::min((size_t)(output_width - oskip), (size_t)(seq_len - si));

			fp << rg.seq.substr(si, out_len);

			if((oskip + out_len) == output_width) {
				fp << endl;
				oskip = 0;
			} else {
				// last line
				oskip = oskip + out_len;
			}

			si += out_len;
		}
	}
	if(oskip) {
		fp << endl;
	}

	fp.close();
}

template<typename TStr>
void RepeatGenerator<TStr>::saveFile()
{
    saveSeeds();
	//saveRepeatSequence();
	//saveRepeatGroup();
}


/**
 * TODO
 * @brief 
 *
 * @param rpt_seq
 * @param rpt_range
 */
template<typename TStr>
void RepeatGenerator<TStr>::addRepeatGroup(map<TIndexOffU, TIndexOffU>& seedpos_to_repeatgroup,
                                           const string& rpt_seq,
                                           const EList<RepeatCoord<TIndexOffU> >& positions)
{
#if 0
	// rpt_seq is always > 0
	//
	const int rpt_len = rpt_seq.length();

	for (int i = 0; i < rpt_grp_.size(); i++) {
		RepeatGroup& rg = rpt_grp_[i];
		string& rseq = rg.rpt_seq;
		const int rlen = rseq.length();
		if (rlen == 0) {
			// skip
			continue;
		}

		if (rlen > rpt_len) {
			// check if rpt_seq is substring of rpt_groups sequeuce
			if (rseq.find(rpt_seq) != string::npos) {
				// substring. exit
				return;
			}
		} else if (rlen <= rpt_len) {
			// check if rpt_groups sequeuce is substring of rpt_seq
			if (rpt_seq.find(rseq) != string::npos) {
				// remove rseq
				rg.rpt_seq = "";
			}
		}
	}
#endif
    
    // count the number of seeds on the sense strand
    size_t sense_mer_count = 0;
    for(size_t i = 0; i < positions.size(); i++) {
        if(positions[i].joinedOff < forward_length_)
            sense_mer_count++;
    }
    
    // skip if there is no sense seeds
    if(sense_mer_count <= 0)
        return;
    
    // skip if a given set of seeds corresponds to an existing set of seeds
    const TIndexOffU pos_diff = 5;
    const size_t sampling = 10;
    size_t add_idx = rpt_grp_.size();
    for(size_t i = 0; i < positions.size(); i += sampling) {
        TIndexOffU joinedOff = positions[i].joinedOff;
        map<TIndexOffU, TIndexOffU>::iterator it = seedpos_to_repeatgroup.lower_bound(joinedOff >= pos_diff ? joinedOff - pos_diff : 0);
        for(; it != seedpos_to_repeatgroup.end(); it++) {
            if(it->first > joinedOff + pos_diff) {
                break;
            }
            assert_geq(it->first + pos_diff, joinedOff);
            assert_lt(it->second, rpt_grp_.size());
            const EList<RepeatCoord<TIndexOffU> >& positions2 = rpt_grp_[it->second].positions;
            
            size_t num_match = 0;
            size_t p = 0, p2 = 0;
            while(p < positions.size() && p2 < positions2.size()) {
                TIndexOffU pos = positions[p].joinedOff, pos2 = positions2[p2].joinedOff;
                if(pos + pos_diff >= pos2 && pos2 + pos_diff >= pos) {
                    num_match++;
                }
                if(pos <= pos2) p++;
                else            p2++;
            }
            
            // if the number of matches is >= 90% of positions in the smaller group
            if(num_match * 10 >= min(positions.size(), positions2.size()) * 9) {
                if(positions.size() <= positions2.size()) {
                    return;
                } else {
                    add_idx = it->second;
                    for(size_t p2 = 0; p2 < positions2.size(); p2++) {
                        if(seedpos_to_repeatgroup[positions2[p2].joinedOff] == add_idx) {
                            seedpos_to_repeatgroup.erase(positions2[p2].joinedOff);
                        }
                    }
                    break;
                }
            }
        }
        if(add_idx < rpt_grp_.size())
            break;
    }

    assert_leq(add_idx, rpt_grp_.size());
    if(add_idx == rpt_grp_.size()) {
        rpt_grp_.expand();
    }
    rpt_grp_[add_idx].seq = rpt_seq;
    rpt_grp_[add_idx].positions = positions;
    
    for(size_t i = 0; i < positions.size(); i++) {
        seedpos_to_repeatgroup[positions[i].joinedOff] = add_idx;
    }
}


/**
 * @brief 
 *
 * @tparam TStr
 */
template<typename TStr> 
void RepeatGenerator<TStr>::mergeRepeatGroup()
{
	int range_count = 0;
	for(size_t i = 0; i < rpt_grp_.size(); i++) {
		range_count += rpt_grp_[i].positions.size();
	}

	cerr << "range_count " << range_count << endl;

	if(range_count == 0) {
		cerr << "no repeat sequeuce" << endl; 
		return;
	}

	EList<RepeatRange> rpt_ranges;
	rpt_ranges.reserveExact(range_count);

    for(size_t i = 0; i < rpt_grp_.size(); i++) {
        RepeatGroup& rg = rpt_grp_[i];
        size_t s_len = rg.seq.length();

        for(size_t j = 0; j < rg.positions.size(); j++) {
            rpt_ranges.push_back(RepeatRange(
                        make_pair(rg.positions[j].joinedOff, rg.positions[j].joinedOff + s_len),
                        i, rg.positions[j].fw));
        }
    }

    assert_gt(rpt_ranges.size(), 0);

	sort(rpt_ranges.begin(), rpt_ranges.begin() + rpt_ranges.size(), 
            compareRepeatRangeByRange);

	// Merge
	int merged_count = 0;
	for(size_t i = 0; i < rpt_ranges.size() - 1;) {
		size_t j = i + 1;
		for(; j < rpt_ranges.size(); j++) {
			// check i, j can be merged 
			//

			if(!checkRangeMergeable(rpt_ranges[i].range, rpt_ranges[j].range)) {
				break;
			}

			rpt_ranges[j].range = EMPTY_RANGE;
			rpt_ranges[j].rg_id = std::numeric_limits<int>::max();
			merged_count++;
		}
		i = j;
	}

	cerr << "merged_count: " << merged_count;
	cerr << endl;

#if 1
	{
		string fname = filename_ + ".rptinfo";
		ofstream fp(fname.c_str());

		for(size_t i = 0; i < rpt_ranges.size(); i++) {
			if(rpt_ranges[i].range == EMPTY_RANGE) {
				continue;
			}
			Range rc = reverseRange(rpt_ranges[i].range, s_.length());
			fp << i;
			fp << "\t" << rpt_ranges[i].range.first;
			fp << "\t" << rpt_ranges[i].range.second;
			fp << "\t" << rpt_grp_[rpt_ranges[i].rg_id].seq;
			fp << "\t" << rc.first;
			fp << "\t" << rc.second;
			fp << "\t" << reverseComplement(rpt_grp_[rpt_ranges[i].rg_id].seq);
			fp << "\t" << rpt_ranges[i].rg_id;
			fp << endl;
		}

		fp.close();
	}
#endif


    /* remake RepeatGroup from rpt_ranges */

	// sort by rg_id
	sort(rpt_ranges.begin(), rpt_ranges.begin() + rpt_ranges.size(), 
            compareRepeatRangeByRgID);

	EList<RepeatGroup> mgroup;

	mgroup.reserveExact(rpt_grp_.size());
	mgroup.swap(rpt_grp_);

	for(size_t i = 0; i < rpt_ranges.size() - 1;) {
		if(rpt_ranges[i].rg_id == std::numeric_limits<int>::max()) {
			break;
		}

		size_t j = i + 1;
		for(; j < rpt_ranges.size(); j++) {
			if(rpt_ranges[i].rg_id != rpt_ranges[j].rg_id) {
				break;
			}
		}

		/* [i, j) has a same rg_id */

		int rg_id = rpt_ranges[i].rg_id;
		rpt_grp_.expand();
		rpt_grp_.back().seq = mgroup[rg_id].seq;
		for (int k = i; k < j; k++) {
			rpt_grp_.back().positions.push_back(
                    RepeatCoord<TIndexOffU>(0, 0,
                        rpt_ranges[k].range.first,
                        rpt_ranges[k].forward));
		}

		// sort positions
        assert_gt(rpt_grp_.back().positions.size(), 0);
        rpt_grp_.back().positions.sort();
#if 0
        sort(rpt_grp_.back().positions.begin(), 
                rpt_grp_.back().positions.begin()
                + rpt_grp_.back().positions.size(),
                compareRepeatCoordByJoinedOff<TIndexOffU>);
#endif

		i = j;
	}

}

template<typename TStr>
void RepeatGenerator<TStr>::groupRepeatGroup(TIndexOffU rpt_edit)
{
    if (rpt_grp_.size() == 0) {
        cerr << "no repeat group" << endl;
        return;
    }

    cerr << "before grouping " << rpt_grp_.size() << endl;

    int step = rpt_grp_.size() >> 8;

    if(step == 0) {step = 1;}

    Timer timer(cerr, "Total time for grouping sequences: ", true);

    for(size_t i = 0; i < rpt_grp_.size() - 1; i++) {
        if(i % step == 0) {
            cerr << i << "/" << rpt_grp_.size() << endl;
        }

        if(rpt_grp_[i].empty()) {
            // empty -> skip
            continue;
        }

        string& str1 = rpt_grp_[i].seq;

        for(size_t j = i + 1; j < rpt_grp_.size(); j++) {
            if(rpt_grp_[j].empty()) {
                // empty -> skip
                continue;
            }
            string& str2 = rpt_grp_[j].seq;
            EList<Edit> edits;
            Coord coord;

            if(checkSequenceMergeable(str1, str2, edits, coord, rpt_edit)) {
                /* i, j merge into i */
                rpt_grp_[i].merge(rpt_grp_[j], edits, coord);

                rpt_grp_[j].set_empty();
            }
        }
    }

	EList<RepeatGroup> mgroup;
    mgroup.reserveExact(rpt_grp_.size());
    mgroup.swap(rpt_grp_);

    for(size_t i = 0; i < mgroup.size(); i++) {
        if (!mgroup[i].empty()) {
            rpt_grp_.expand();
            rpt_grp_.back() = mgroup[i];
        }
    }

    cerr << "after merge " << rpt_grp_.size() << endl;

#if 1
    {
        string fname = filename_ + ".altseq";
        ofstream fp(fname.c_str());

        for(size_t i = 0; i < rpt_grp_.size(); i++) {
            RepeatGroup& rg = rpt_grp_[i];
            if(rg.empty()) {
                continue;
            }
            fp << i;
            fp << "\t" << rg.alt_seq.size();
            fp << "\t" << rg.seq;
            for(size_t j = 0; j < rg.alt_seq.size(); j++) {
                fp << "\t" << toMDZ(rg.alt_seq[j].edits, rg.alt_seq[j].seq /* read */);
                fp << "\t" << rg.alt_seq[j].seq;
            }
            fp << endl;
        }

        fp.close();
    }
#endif
}

template<typename TStr>
TIndexOffU RepeatGenerator<TStr>::getEnd(TIndexOffU e) {
    assert_lt(e, s_.length())
    
    TIndexOffU end = 0;
    if(e < forward_length_) {
        int frag_id = mapJoinedOffToSeq(e);
        assert_geq(frag_id, 0);
        end = fraglist_[frag_id].joinedOff + fraglist_[frag_id].length;
    } else {
        // ReverseComplement
        // a, b are positions w.r.t reverse complement string.
        // fragment map is based on forward string
        int frag_id = mapJoinedOffToSeq(s_.length() - e - 1);
        assert_geq(frag_id, 0);
        end = s_.length() - fraglist_[frag_id].joinedOff;
    }
    
    assert_leq(end, s_.length());
    return end;
}

template<typename TStr>
TIndexOffU RepeatGenerator<TStr>::getStart(TIndexOffU e) {
    assert_lt(e, s_.length())
    
    TIndexOffU start = 0;
    if(e < forward_length_) {
        int frag_id = mapJoinedOffToSeq(e);
        assert_geq(frag_id, 0);
        start = fraglist_[frag_id].joinedOff;
    } else {
        // ReverseComplement
        // a, b are positions w.r.t reverse complement string.
        // fragment map is based on forward string
        int frag_id = mapJoinedOffToSeq(s_.length() - e - 1);
        assert_geq(frag_id, 0);
        start = s_.length() - (fraglist_[frag_id].joinedOff + fraglist_[frag_id].length);
    }
    
    assert_leq(start, s_.length());
    return start;
}

template<typename TStr>
TIndexOffU RepeatGenerator<TStr>::getLCP(TIndexOffU a, TIndexOffU b)
{
    size_t a_end = getEnd(a);
    size_t b_end = getEnd(b);

    assert_leq(a_end, s_.length());
    assert_leq(b_end, s_.length());
    
    TIndexOffU k = 0;
    while((a + k) < a_end && (b + k) < b_end) {
        if(s_[a + k] != s_[b + k]) {
            break;
        }
        k++;
    }
    
    return k;
}

template<typename TStr>
void RepeatGenerator<TStr>::makePadString(const string& ref, const string& read,
        string& pad, size_t len)
{
    pad.resize(len);

    for(size_t i = 0; i < len; i++) {
        // shift A->C, C->G, G->T, T->A
        pad[i] = "CGTA"[asc2dna[ref[i]]];

        if(read[i] == pad[i]) {
            // shift
            pad[i] = "CGTA"[asc2dna[pad[i]]];
        }
    }

    int head_len = len / 2;
    size_t pad_start = len - head_len;

    for(size_t i = 0; i < head_len; i++) {
        if(read[i] == pad[pad_start + i]) {
            // shift
            pad[pad_start + i] = "CGTA"[asc2dna[pad[pad_start + i]]];
        }
    }
}

template<typename TStr>
bool RepeatGenerator<TStr>::checkSequenceMergeable(const string& ref,
                                                   const string& read,
                                                   EList<Edit>& edits,
                                                   Coord& coord,
                                                   TIndexOffU rpt_len,
                                                   TIndexOffU max_edit)
{
    size_t max_matchlen = 0;
    EList<Edit> ed;

    string pad;
    makePadString(ref, read, pad, 5);

    string ref2 = pad + ref;
    string read2 = pad + read;

    alignStrings(ref2, read2, ed, coord);

    // match should start from pad string
    if(coord.off() != 0) {
        return false;
    }

    // no edits on pad string
    if(ed.size() > 0 && ed[0].pos < pad.length()) {
        return false;
    }

    size_t left = pad.length();
    size_t right = left + read.length();

    edits.clear();
    edits.reserveExact(ed.size());
    for(size_t i = 0; i < ed.size(); i++) {
        if(ed[i].pos >= left && ed[i].pos <= right) {
            edits.push_back(ed[i]);
            edits.back().pos -= left;
        }
    }

    max_matchlen = getMaxMatchLen(edits, read.length());

#ifdef DEBUGLOG
    {
        cerr << "After pad removed" << endl;
        BTDnaString btread;
        btread.install(read.c_str(), true);
        Edit::print(cerr, edits); cerr << endl;
        Edit::printQAlign(cerr, btread, edits);
    }
#endif

    return (max_matchlen >= rpt_len);
}

template<typename TStr>
int RepeatGenerator<TStr>::alignStrings(const string &ref, const string &read, EList<Edit>& edits, Coord& coord)
{
    // Prepare Strings

    // Read -> BTDnaString
    // Ref -> bit-encoded string

    //SwAligner swa;

    BTDnaString btread;
    BTString btqual;
    BTString btref;
    BTString btref2;

    BTDnaString btreadrc;
    BTString btqualrc;


    string qual = "";
    for(size_t i = 0; i < read.length(); i++) {
        qual.push_back('I');
    }

#if 0
    cerr << "REF : " << ref << endl;
    cerr << "READ: " << read << endl;
    cerr << "QUAL: " << qual << endl;
#endif

    btread.install(read.c_str(), true);
    btreadrc = btread;
    btreadrc.reverseComp();

    btqual.install(qual.c_str());
    btqualrc = btqual;

    btref.install(ref.c_str());

    TAlScore min_score = sc_->scoreMin.f<TAlScore >((double)btread.length());

    btref2 = btref;

    size_t nceil = 0;
    size_t nrow = btread.length();

    // Convert reference string to mask
    for(size_t i = 0; i < btref2.length(); i++) {
        if(toupper(btref2[i]) == 'N') {
            btref2.set(16, i);
        } else {
            int num = 0;
            int alts[] = {4, 4, 4, 4};
            decodeNuc(toupper(btref2[i]), num, alts);
            assert_leq(num, 4);
            assert_gt(num, 0);
            btref2.set(0, i);
            for(int j = 0; j < num; j++) {
                btref2.set(btref2[i] | (1 << alts[j]), i);
            }
        }
    }


    bool fw = true;
    uint32_t refidx = 0;
    
    swa.initRead(
            btread,     // read sequence
            btreadrc,
            btqual,     // read qualities
            btqualrc,
            0,          // offset of first character within 'read' to consider
            btread.length(), // offset of last char (exclusive) in 'read' to consider
            *sc_);      // local-alignment score floor

    DynProgFramer dpframe(false);
    size_t readgaps = 0;
    size_t refgaps = 0;
    size_t maxhalf = 0;

    DPRect rect;
    dpframe.frameSeedExtensionRect(
            0,              // ref offset implied by seed hit assuming no gaps
            btread.length(),    // length of read sequence used in DP table
            btref2.length(),    // length of reference
            readgaps,   // max # of read gaps permitted in opp mate alignment
            refgaps,    // max # of ref gaps permitted in opp mate alignment
            (size_t)nceil,  // # Ns permitted
            maxhalf, // max width in either direction
            rect);      // DP Rectangle

    assert(rect.repOk());

    size_t cminlen = 2000, cpow2 = 4;

    swa.initRef(
            fw,                 // whether to align forward or revcomp read
            refidx,             // reference ID
            rect,               // DP rectangle
            btref2.wbuf(),      // reference strings
            0,                  // offset of first reference char to align to
            btref2.length(),    // offset of last reference char to align to
            btref2.length(),    // length of reference sequence
            *sc_,               // scoring scheme
            min_score,          // minimum score
            true,               // use 8-bit SSE if positions
            cminlen,            // minimum length for using checkpointing scheme
            cpow2,              // interval b/t checkpointed diags; 1 << this
            false,              // triangular mini-fills?
            false               // is this a seed extension?
            );


    TAlScore best = std::numeric_limits<TAlScore>::min();
    bool found = swa.align(rnd_, best);
#ifdef DEBUGLOG 
    cerr << "found: " << found << "\t" << best << "\t" << "minsc: " << min_score << endl;
#endif

    if (found) {
#ifdef DEBUGLOG 
        //cerr << "CP " << "found: " << found << "\t" << best << "\t" << "minsc: " << min_score << endl;
        cerr << "REF : " << ref << endl;
        cerr << "READ: " << read << endl;
#endif

        SwResult res;
        int max_match_len = 0;
        res.reset();
        res.alres.init_raw_edits(&rawEdits_);

        found = swa.nextAlignment(res, best, rnd_);
        if (found) {
            edits = res.alres.ned();
            //const TRefOff ref_off = res.alres.refoff();
            //const Coord& coord = res.alres.refcoord();
            coord = res.alres.refcoord();
            //assert_geq(genomeHit._joinedOff + coord.off(), genomeHit.refoff());

#ifdef DEBUGLOG 
            cerr << "num edits: " << edits.size() << endl;
            cerr << "coord: " << coord.off();
            cerr << ", " << coord.ref();
            cerr << ", " << coord.orient();
            cerr << ", " << coord.joinedOff();
            cerr << endl;
            Edit::print(cerr, edits); cerr << endl;
            Edit::printQAlign(cerr, btread, edits);
#endif

            max_match_len = getMaxMatchLen(edits, btread.length());
#ifdef DEBUGLOG 
            cerr << "max match length: " << max_match_len << endl;
#endif
        }
#ifdef DEBUGLOG 
        cerr << "nextAlignment: " << found << endl;
        cerr << "-------------------------" << endl;
#endif
    }

    return 0;
}

template<typename TStr>
void RepeatGenerator<TStr>::doTestCase1(const string& refstr, const string& readstr, TIndexOffU rpt_edit)
{
    cerr << "doTestCase1----------------" << endl;
    EList<Edit> edits;
    Coord coord;

    if (refstr.length() == 0 ||
            readstr.length() == 0) {
        return;
    }

    EList<Edit> ed;

    string pad;
    makePadString(refstr, readstr, pad, 5);

    string ref2 = pad + refstr + pad;
    string read2 = pad + readstr + pad;
    alignStrings(refstr, readstr, edits, coord);

    size_t left = pad.length();
    size_t right = left + readstr.length();

    edits.reserveExact(ed.size());
    for(size_t i = 0; i < ed.size(); i++) {
        if(ed[i].pos >= left && ed[i].pos <= right) {
            edits.push_back(ed[i]);
            edits.back().pos -= left;
        }
    }


    RepeatGroup rg;

    rg.edits = edits;
    rg.coord = coord;
    rg.seq = readstr;
    rg.base_offset = 0;

    string chr_name = "rep";

    cerr << "REF : " << refstr << endl;
    cerr << "READ: " << readstr << endl;
    size_t snpids = 0;
    rg.buildSNPs(snpids);
    rg.writeSNPs(cerr, chr_name); cerr << endl;

}


/****************************/
template class RepeatGenerator<SString<char> >;
template void dump_tstr(const SString<char>& );
template bool compareRepeatCoordByJoinedOff(const RepeatCoord<TIndexOffU>& , const RepeatCoord<TIndexOffU>&);
