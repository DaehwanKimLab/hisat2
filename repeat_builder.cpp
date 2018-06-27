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

#include "opal_wrapper.h"

#include "repeat_builder.h"

typedef pair<TIndexOffU, TIndexOffU> Range;
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
        return 0;
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
NRG<TStr>::NRG(TStr& s,
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
	cerr << "NRG: " << filename_ << endl;

    rnd_.init(0);

	// build ref_names_ from ref_namelines_
    buildNames();
    buildJoinedFragment();
}

template<typename TStr>
NRG<TStr>::~NRG()
{
    if (sc_) {
        delete sc_;
    }
}

#define DMAX std::numeric_limits<double>::max()

template<typename TStr>
void NRG<TStr>::init_dyn()
{

//#define MM_PEN  3
#define MM_PEN  6
#define GAP_PEN_LIN 2
//#define GAP_PEN_LIN (((MM_PEN) * rpt_edit_ + 1) * 1.0)
#define GAP_PEN_CON 1
//#define GAP_PEN_CON (((MM_PEN) * rpt_edit_ + 1) * 1.0)
#define MAX_PEN (MAX_I16)

    scoreMin_.init(SIMPLE_FUNC_LINEAR, rpt_edit_ * MM_PEN * -1.0, 0.0);
    nCeil_.init(SIMPLE_FUNC_LINEAR, 0.0, 0.0);

    penCanIntronLen_.init(SIMPLE_FUNC_LOG, -8, 1);
    penNoncanIntronLen_.init(SIMPLE_FUNC_LOG, -8, 1);

    sc_ = new Scoring(
            DEFAULT_MATCH_BONUS,     // constant reward for match
            DEFAULT_MM_PENALTY_TYPE,     // how to penalize mismatches
            DEFAULT_MM_PENALTY_MAX, //MM_PEN,      // max mm penalty
            DEFAULT_MM_PENALTY_MAX, //MM_PEN,      // min mm penalty
            MAX_PEN,      // max sc penalty
            MAX_PEN,      // min sc penalty
            scoreMin_,       // min score as function of read len
            nCeil_,          // max # Ns as function of read len
            DEFAULT_N_PENALTY_TYPE,      // how to penalize Ns in the read
            DEFAULT_N_PENALTY,           // constant if N pelanty is a constant
            DEFAULT_N_CAT_PAIR,    // whether to concat mates before N filtering

            DEFAULT_READ_GAP_CONST, //GAP_PEN_CON, //MM_PEN,  // constant coeff for read gap cost
            DEFAULT_REF_GAP_CONST, //GAP_PEN_CON, //MM_PEN,  // constant coeff for ref gap cost
            DEFAULT_READ_GAP_LINEAR, //GAP_PEN_LIN, //MM_PEN, // linear coeff for read gap cost
            DEFAULT_REF_GAP_LINEAR, //GAP_PEN_LIN, //MM_PEN, // linear coeff for ref gap cost
            1 /* gGapBarrier */    // # rows at top/bot only entered diagonally
            );
}


template<typename TStr>
void NRG<TStr>::doTest(TIndexOffU rpt_len,
               TIndexOffU rpt_cnt,
               bool flagGrouping,
               TIndexOffU rpt_edit,
               TIndexOffU rpt_matchlen,
               const string& refstr,
               const string& readstr)
{
	TIndexOffU count = 0;

    rpt_matchlen_ = rpt_matchlen;
    rpt_edit_ = rpt_edit;

    init_dyn();

    doTestCase1(refstr, readstr, rpt_edit);
    doTestCase2(refstr, readstr, rpt_edit);
}



template<typename TStr>
void NRG<TStr>::build(TIndexOffU rpt_len,
                      TIndexOffU rpt_cnt,
                      bool flagGrouping,
                      TIndexOffU rpt_edit,
                      TIndexOffU rpt_matchlen)
{
	TIndexOffU count = 0;

    rpt_matchlen_ = rpt_matchlen;
    rpt_edit_ = rpt_edit;

    init_dyn();

	EList<RepeatCoord<TIndexOffU> > rpt_positions;
	TIndexOffU min_lcp_len = s_.length();

	while(count < s_.length() + 1) {
		TIndexOffU saElt = bsa_.nextSuffix();
		count++;
		
		if(count && (count % 1000000 == 0)) {
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

			if(lcp_len >= rpt_len) {
				rpt_positions.expand();
				rpt_positions.back().joinedOff = saElt;

				if(min_lcp_len > lcp_len) {
					min_lcp_len = lcp_len;
				}

			} else {
				if (rpt_positions.size() >= rpt_cnt) {
                    sort(rpt_positions.begin(), 
                            rpt_positions.begin() + rpt_positions.size(),
                            compareRepeatCoordByJoinedOff<TIndexOffU>);

                    string ss = getString(s_, prev_saElt, min_lcp_len);
					addRepeatGroup(ss, rpt_positions);
				}

				// flush previous positions 
				rpt_positions.resize(1);
				rpt_positions.back().joinedOff = saElt;
				min_lcp_len = s_.length();
			}
		}
	}

    mergeRepeatGroup();
    
    if(flagGrouping && rpt_edit > 0) {
        groupRepeatGroup(rpt_edit);
    }

	// we found repeat_group
	cerr << "CP " << rpt_grp_.size() << " groups found" << endl;

}

template<typename TStr>
void NRG<TStr>::buildNames()
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

template<typename TStr>
int NRG<TStr>::mapJoinedOffToSeq(TIndexOffU joinedOff)
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
int NRG<TStr>::getGenomeCoord(TIndexOffU joinedOff, 
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
void NRG<TStr>::buildJoinedFragment()
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
void NRG<TStr>::sortRepeatGroup()
{
	if(rpt_grp_.size() > 0) {
		sort(rpt_grp_.begin(), rpt_grp_.begin() + rpt_grp_.size(), 
                compareRepeatGroupByJoinedOff);
	}
}


template<typename TStr>
void NRG<TStr>::saveRepeatPositions(ofstream& fp, RepeatGroup& rg)
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


template<typename TStr>
void NRG<TStr>::saveRepeatGroup()
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
        // fp << "\t" << rg.seq;
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
void NRG<TStr>::saveRepeatSequence()
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
void NRG<TStr>::saveFile()
{
	saveRepeatSequence();
	saveRepeatGroup();
}


/**
 * TODO
 * @brief 
 *
 * @param rpt_seq
 * @param rpt_range
 */
template<typename TStr>
void NRG<TStr>::addRepeatGroup(const string& rpt_seq, const EList<RepeatCoord<TIndexOffU> >& positions)
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
    
    // DK - check this out
    size_t sense_mer_count = 0;
    for(size_t i = 0; i < positions.size(); i++) {
        if(positions[i].joinedOff < forward_length_)
            sense_mer_count++;
    }

	// add to last
    if(sense_mer_count > 0) {
        rpt_grp_.expand();
        rpt_grp_.back().seq = rpt_seq;
        rpt_grp_.back().positions = positions;
    }
}


/**
 * @brief 
 *
 * @tparam TStr
 */
template<typename TStr> 
void NRG<TStr>::mergeRepeatGroup()
{
	int range_count = 0;
	for(size_t i = 0; i < rpt_grp_.size(); i++) {
		range_count += rpt_grp_[i].positions.size();
	}

	cerr << "CP " << "range_count " << range_count << endl;

	if(range_count == 0) {
		cerr << "CP " << "no repeat sequeuce" << endl; 
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

	cerr << "CP ";
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
			fp << "CP " << i ;
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
void NRG<TStr>::groupRepeatGroup(TIndexOffU rpt_edit)
{
    if (rpt_grp_.size() == 0) {
        cerr << "CP " << "no repeat group" << endl;
        return;
    }

    cerr << "CP " << "before grouping " << rpt_grp_.size() << endl;

    int step = rpt_grp_.size() >> 8;

    if(step == 0) {step = 1;}
    cerr << "CP " << step << endl;

    Timer timer(cerr, "Total time for grouping sequences: ", true);

    for(size_t i = 0; i < rpt_grp_.size() - 1; i++) {
        if(i % step == 0) {
            cerr << "CP " << i << "/" << rpt_grp_.size() << endl;
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

    cerr << "CP " << "after merge " << rpt_grp_.size() << endl;

#if 1
    {
        string fname = filename_ + ".altseq";
        ofstream fp(fname.c_str());

        for(size_t i = 0; i < rpt_grp_.size(); i++) {
            RepeatGroup& rg = rpt_grp_[i];
            if(rg.empty()) {
                continue;
            }
            fp << "CP " << i ;
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
TIndexOffU NRG<TStr>::getEnd(TIndexOffU e) {
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
TIndexOffU NRG<TStr>::getLCP(TIndexOffU a, TIndexOffU b)
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
bool NRG<TStr>::checkSequenceMergeable(const string& ref, const string& read, 
        EList<Edit>& edits, Coord& coord, TIndexOffU max_edit)
{
    size_t max_matchlen = 0;

    // TODO:
#if 1
    // merge two strings if have same length
    if (ref.length() != read.length()) {
        return false;
    }

    alignStrings(ref, read, edits, coord);
#else

    coord.init(0, 0, true, 0);
    opal.alignStrings(ref, read, edits);
    //max_matchlen = getMaxMatchLen(edits, read.length());
    //cerr << "max_matchlen: " << max_matchlen << endl;
#endif

    max_matchlen = getMaxMatchLen(edits, read.length());
    return (max_matchlen >= rpt_matchlen_);

#if 0
	unsigned int ed = levenshtein_distance(s1, s2);
    return false;
	//return (ed <= max_edit);
    //return (alignStrings(s1, s2) >= 50);
#endif

}

template<typename TStr>
int NRG<TStr>::alignStrings(const string &ref, const string &read, EList<Edit>& edits, Coord& coord)
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
    TAlScore floor_score = sc_->scoreMin.f<TAlScore >((double)btread.length());

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
    cerr << "CP " << "found: " << found << "\t" << best << "\t" << "minsc: " << min_score << endl;
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
            cerr << "CP ";
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
        cerr << "CP " << "nextAlignment: " << found << endl;
        cerr << "CP -------------------------" << endl;
#endif
    }

    return 0;
}

template<typename TStr>
void NRG<TStr>::doTestCase1(const string& refstr, const string& readstr, TIndexOffU rpt_edit)
{
    cerr << "doTestCase1----------------" << endl;
    EList<Edit> edits;
    Coord coord;

    if (refstr.length() == 0 ||
            readstr.length() == 0) {
        return;
    }

    //checkSequenceMergeable(refstr, readstr, edits, coord, rpt_edit);
    alignStrings(refstr, readstr, edits, coord);


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

template<typename TStr>
void NRG<TStr>::doTestCase2(const string& refstr, const string& readstr, TIndexOffU rpt_edit)
{
    cerr << "doTestCase2----------------" << endl;
    EList<Edit> edits;
    Coord coord;

    if (refstr.length() == 0 ||
            readstr.length() == 0) {
        return;
    }

    OPAL opal;
    opal.setMinScore(rpt_edit);

    //checkSequenceMergeable(refstr, readstr, edits, coord, rpt_edit);
    opal.alignStrings(refstr, readstr, edits);
    coord.init(0, 0, true, 0);

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
template class NRG<SString<char> >;
template void dump_tstr(const SString<char>& );
template bool compareRepeatCoordByJoinedOff(const RepeatCoord<TIndexOffU>& , const RepeatCoord<TIndexOffU>&);
