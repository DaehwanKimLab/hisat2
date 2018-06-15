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

bool seq_mergeable(const string& s1, const string& s2, TIndexOffU max_edit = 10)
{
	unsigned int ed = levenshtein_distance(s1, s2);
	return (ed <= max_edit);
}

typedef pair<TIndexOffU, TIndexOffU> Range;
static const Range EMPTY_RANGE = Range(1, 0);

struct RepeatRange {
	RepeatRange() {};
	RepeatRange(Range r, int id) : 
		range(r), rg_id(id) {};

	Range range;
	int rg_id;
};

/**
 * @brief return true iff a U b = a 
 *
 * @param a
 * @param b
 *
 * @return 
 */
static bool range_mergeable(const Range& a, const Range& b)
{
	if (a.first <= b.first && a.second >= b.second) {
		return true;
	}

	return false;
}


static bool cmp_by_repeat_range(const RepeatRange& a, const RepeatRange& b)
{
	if ((a.range.second > b.range.second) ||
			(a.range.second == b.range.second && a.range.first < b.range.first)) {
		return true;
	}

	return false;
}

static bool cmp_by_rg_id(const RepeatRange& a, const RepeatRange& b)
{
	return a.rg_id < b.rg_id;
}



template<typename TStr>
string get_string(TStr& ref, TIndexOffU start, int len)
{
	string s;
	size_t ref_len = ref.length();

	for (int i = 0; (i < len) && (start + i < ref_len); i++) {
		char nt = "ACGT"[ref[start + i]];
		s.push_back(nt);
	}

	return s;
}

template<typename TStr>
void masking_with_N(TStr& s, TIndexOffU start, size_t length)
{
	size_t s_len = s.length();

	for (size_t pos = 0; (pos < length) && (start + pos < s_len); pos++) {
		s[start + pos] = 0x04;
	}
}

template<typename TStr>
void dump_tstr(TStr& s)
{
	static int print_width = 60;

	size_t s_len = s.length();

	for (size_t i = 0; i < s_len; i += print_width) {
		string buf;
		for (size_t j = 0; (j < print_width) && (i + j < s_len); j++) {
			buf.push_back("ACGTN"[s[i + j]]);
		}
		cerr << buf << endl;
	}
	cerr << endl;
}


template<typename TStr>
NRG<TStr>::NRG(
		EList<RefRecord>& szs,
		EList<string>& ref_names,
		TStr& s,
		string& filename,
		BlockwiseSA<TStr>& sa) :
	szs_(szs), ref_namelines_(ref_names), 
	s_(s), filename_(filename), bsa_(sa)
{
	cerr << "NRG: " << filename_ << endl;

	// build ref_names_ from ref_namelines_
	build_names();
	build_joined_fragment();
}

template<typename TStr>
void NRG<TStr>::build(TIndexOffU rpt_len,
                      TIndexOffU rpt_cnt,
                      bool flagGrouping,
                      TIndexOffU rpt_edit)
{
	TIndexOffU count = 0;

	EList<TIndexOffU> rpt_positions;
	TIndexOffU min_lcp_len = s_.length();

	while(count < s_.length() + 1) {
		TIndexOffU saElt = bsa_.nextSuffix();
		count++;

		if (count && (count % 1000000 == 0)) {
			cerr << "SA count " << count << endl;
		}

		if (rpt_positions.size() == 0) {
			rpt_positions.expand();
			rpt_positions.back() = saElt;
		} else {
			TIndexOffU prev_saElt = rpt_positions.back();

			// calculate common prefix length between two text.
			//   text1 is started from prev_saElt and text2 is started from saElt
			int lcp_len = get_lcp(prev_saElt, saElt);

			// check prev_saElt

			if (lcp_len >= rpt_len) {
				rpt_positions.expand();
				rpt_positions.back() = saElt;

				if (min_lcp_len > lcp_len) {
					min_lcp_len = lcp_len;
				}

			} else {
				if (rpt_positions.size() >= rpt_cnt) {
					// save ranges
					//cerr << "CP " << "we have " << rpt_positions.size() << " positions" << ", " << min_lcp_len << endl;
					//

					rpt_positions.sort();

					string ss = get_string(s_, prev_saElt, min_lcp_len);

					add_repeat_group(ss, rpt_positions);

#if 0
					rg.positions = repeat_range;
					rg.positions.sort();
					rg.rpt_seq = get_string(s, prev_saElt, min_lcp_len);
#endif
				}

				// flush previous positions 
				rpt_positions.resize(1);
				rpt_positions.back() = saElt;
				min_lcp_len = s_.length();
			}
		}

#if 0//{{{
		// DK - debugging purposes
		if(count < 100) {
			suffixes.expand();
			cerr << setw(12) << saElt << "\t";
			for(int k = 0; k < 100; k++) {
				char nt = "ACGT"[s[saElt+k]];
				cerr << nt;
				suffixes.back().push_back(nt);
			}
			cerr << endl;

			if(count > 1) {
				SwAligner al;
				SwResult res;

				SimpleFunc scoreMin;
				scoreMin.init(SIMPLE_FUNC_LINEAR, 0.0f, -0.2f);

				SimpleFunc nCeil;
				nCeil.init(SIMPLE_FUNC_LINEAR, 0.0f, 0, 2.0f, 0.1f);

				const string& str1 = suffixes[suffixes.size() - 2];
				const string& str2 = suffixes[suffixes.size() - 1];

				string qual = "";
				for(int i = 0; i < str1.length(); i++) {
					qual.push_back('I');
				}

				// Set up penalities
				Scoring sc(
						DEFAULT_MATCH_BONUS,     // constant reward for match
						DEFAULT_MM_PENALTY_TYPE,     // how to penalize mismatches
						30,        // constant if mm pelanty is a constant
						30,        // penalty for decoded SNP
						0,
						0,
						scoreMin,       // min score as function of read len
						nCeil,          // max # Ns as function of read len
						DEFAULT_N_PENALTY_TYPE,      // how to penalize Ns in the read
						DEFAULT_N_PENALTY,          // constant if N pelanty is a constant
						DEFAULT_N_CAT_PAIR,      // true -> concatenate mates before N filtering
						25,  // constant coeff for cost of gap in read
						25,  // constant coeff for cost of gap in ref
						15, // linear coeff for cost of gap in read
						15, // linear coeff for cost of gap in ref
						1,    // # rows at top/bot only entered diagonally
						0,   // canonical splicing penalty
						0,   // non-canonical splicing penalty
						0);  // conflicting splice site penalt

				doTestCase2(
						al,
						str1.c_str(),
						qual.c_str(),
						str2.c_str(),
						0,
						sc,
						DEFAULT_MIN_CONST,
						DEFAULT_MIN_LINEAR,
						res);
			}
		}
#endif//}}}

	}

	adjust_repeat_group(flagGrouping, rpt_edit);
	// we found repeat_group
	cerr << "CP " << rpt_grp_.size() << " groups found" << endl;

	//dump_tstr(s);

	// write to FA
	// sequence
	// repeat sequeuce
	savefile();
}

template<typename TStr>
void NRG<TStr>::build_names(void)
{
	ref_names_.resize(ref_namelines_.size());
	for (int i = 0; i < ref_namelines_.size(); i++) {
		string& nameline = ref_namelines_[i];

		for (int j = 0; j < nameline.length(); j++) {
			char n = nameline[j];
			if (n == ' ') {
				break;
			}
			ref_names_[i].push_back(n);
		}
	}
}

template<typename TStr>
int NRG<TStr>::map_joined_pos_to_seq(TIndexOffU joined_pos)
{

	/* search from cached_list */
	if (num_cached_ > 0) {
		for (int i = 0; i < num_cached_; i++) {
			Fragments *frag = &cached_[i];
			if (frag->contain(joined_pos)) {
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
	while ((bot - top) > 1) {
		pos = top + ((bot - top) >> 1);
		frag = &fraglist_[pos];

		if (joined_pos < frag->start) {
			bot = pos;
		} else {
			top = pos;
		}
	}

	frag = &fraglist_[top];
	if (frag->contain(joined_pos)) {
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
int NRG<TStr>::get_genome_coord(TIndexOffU joined_pos, 
		string& chr_name, TIndexOffU& pos_in_chr)
{
	int seq_id = map_joined_pos_to_seq(joined_pos);
	if (seq_id < 0) {
		return -1;
	}

	Fragments *frag = &fraglist_[seq_id];
	TIndexOffU offset = joined_pos - frag->start;

	pos_in_chr = frag->start_in_seq + offset;
	chr_name = ref_names_[frag->seq_id];

	return 0;
}

template<typename TStr>
void NRG<TStr>::build_joined_fragment(void)
{
	int n_seq = 0;
	int n_frag = 0;

	for (int i = 0; i < szs_.size(); i++) {
		if (szs_[i].len > 0) n_frag++;
		if (szs_[i].first && szs_[i].len > 0) n_seq++;
	}

	int npos = 0;
	int seq_id = -1;
	TIndexOffU acc_frag_length = 0;
	TIndexOffU acc_ref_length = 0;
	fraglist_.resize(n_frag + 1);

	for (int i = 0; i < szs_.size(); i++) {
		if (szs_[i].len == 0) {
			continue;
		}

		fraglist_[npos].start = acc_frag_length;
		fraglist_[npos].length = szs_[i].len;
		fraglist_[npos].start_in_seq = acc_ref_length + szs_[i].off;
		fraglist_[npos].frag_id = i;
		fraglist_[npos].frag_id = npos;
		if (szs_[i].first) {
			seq_id++;
			fraglist_[npos].first = true;
		}
		fraglist_[npos].seq_id = seq_id;

		acc_frag_length += szs_[i].len;
		acc_ref_length += szs_[i].off + szs_[i].len;

		npos++;
	}

	// Add Last Fragment(empty)
	fraglist_[npos].start = acc_frag_length;
	fraglist_[npos].length = 0;
	fraglist_[npos].start_in_seq = acc_ref_length + szs_.back().off;
}

template<typename TStr>
void NRG<TStr>::sort_rpt_grp(void)
{
	if (rpt_grp_.size() > 0) {
		sort(rpt_grp_.begin(), rpt_grp_.begin() + rpt_grp_.size(), repeat_group_cmp);
	}
}

template<typename TStr>
void NRG<TStr>::saveRepeatGroup(void)
{
	string rptinfo_filename = filename_ + ".rep.info";

	int rpt_count = rpt_grp_.size();
	TIndexOffU acc_pos = 0;

	ofstream fp(rptinfo_filename.c_str());

	for (int i = 0; i < rpt_count; i++) {
		RepeatGroup& rg = rpt_grp_[i];
		EList<TIndexOffU>& positions = rg.positions;

		// >rpt_name*0\trep\trep_pos\trep_len\tpos_count\t0
		// chr_name:pos:direction chr_name:pos:direction
		//
		// >rep1*0	rep	0	100	470	0
		// 22:112123123:+ 22:1232131113:+
		//

		// Header line
		fp << ">" << "rpt_" << i << "*0";
		fp << "\t" << "rep"; // TODO
		fp << "\t" << acc_pos;
		fp << "\t" << rg.seq.length();
		fp << "\t" << positions.size();
		fp << "\t" << "0"; 
		fp << endl;

		acc_pos += rg.seq.length();

		// Positions
		for (int j = 0; j < positions.size(); j++) {
			if (j && (j % 10 == 0)) {
				fp << endl;
			}

			if (j % 10) {
				fp << " ";
			}

			string chr_name;
			TIndexOffU pos_in_chr;

			get_genome_coord(positions[j], chr_name, pos_in_chr);
			char direction = '+';

			fp << chr_name << ":" << pos_in_chr << ":" << direction;
		}
		fp << endl;
	}		
	fp.close();
}
	
	
template<typename TStr>
void NRG<TStr>::saveRepeatSequence(void)
{
	string fname = filename_ + ".rep.fa";

	ofstream fp(fname.c_str());

	/* TODO */
	fp << ">" << "rep" << endl;

	int oskip = 0;

	for (TIndexOffU grp_idx = 0; grp_idx < rpt_grp_.size(); grp_idx++) {

		RepeatGroup& rg = rpt_grp_[grp_idx];
        size_t seq_len = rg.seq.length();

		TIndexOffU si = 0;
		while (si < seq_len) {
			size_t out_len = std::min((size_t)(output_width - oskip), (size_t)(seq_len - si));

			fp << rg.seq.substr(si, out_len);

			if ((oskip + out_len) == output_width) {
				fp << endl;
				oskip = 0;
			} else {
				// last line
				oskip = oskip + out_len;
			}

			si += out_len;
		}
	}
	if (oskip) {
		fp << endl;
	}

	fp.close();
}

template<typename TStr>
void NRG<TStr>::savefile(void)
{
	// save repeat sequences
	saveRepeatSequence();

	// save repeat infos
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
void NRG<TStr>::add_repeat_group(string& rpt_seq, EList<TIndexOffU>& rpt_range)
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

	// add to last
	rpt_grp_.expand();
	rpt_grp_.back().seq = rpt_seq;
	rpt_grp_.back().positions = rpt_range;
}

/**
 * @brief Remove empty repeat group
 */
template<typename TStr>
void NRG<TStr>::adjust_repeat_group(bool flagGrouping, TIndexOffU rpt_edit)
{
	cerr << "CP " << "repeat_group size " << rpt_grp_.size() << endl;

	int range_count = 0;
	for (int i = 0; i < rpt_grp_.size(); i++) {
		range_count += rpt_grp_[i].positions.size();
	}

	cerr << "CP " << "range_count " << range_count << endl;

	if (range_count == 0) {
		cerr << "CP " << "no repeat sequeuce" << endl; 
		return;
	}

	cerr << "Build RepeatRange" << endl;

	EList<RepeatRange> rpt_ranges;
	rpt_ranges.reserveExact(range_count);

	for (int i = 0; i < rpt_grp_.size(); i++) {
		RepeatGroup& rg = rpt_grp_[i];
		size_t s_len = rg.seq.length();

		for (int j = 0; j < rg.positions.size(); j++) {
			rpt_ranges.push_back(RepeatRange(make_pair(rg.positions[j], rg.positions[j] + s_len), i));
		}
	}
    assert_eq(rpt_ranges.size(), range_count);


	sort(rpt_ranges.begin(), rpt_ranges.begin() + rpt_ranges.size(), 
			cmp_by_repeat_range);

	// Dump
#if 0
	for (int i = 0; i < rpt_ranges.size(); i++) {
		cerr << "CP " << i ;
		cerr << "\t" << rpt_ranges[i].range.first;
		cerr << "\t" << rpt_ranges[i].range.second;
		cerr << "\t" << rpt_grp_[rpt_ranges[i].rg_id].rpt_seq;
		cerr << "\t" << rpt_ranges[i].rg_id;
		cerr << endl;
	}
#endif

#if 1	
	// Merge
	int merged_count = 0;
	for (int i = 0; i < rpt_ranges.size() - 1;) {
		int j = i + 1;
		for (; j < rpt_ranges.size(); j++) {
			// check i, j can be merged 
			//

			if (!range_mergeable(rpt_ranges[i].range, rpt_ranges[j].range)) {
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
#endif

	// Dump
#if 1
	{
		string fname = filename_ + ".rptinfo";
		ofstream fp(fname.c_str());

		for (int i = 0; i < rpt_ranges.size(); i++) {
			if (rpt_ranges[i].range == EMPTY_RANGE) {
				continue;
			}
			fp << "CP " << i ;
			fp << "\t" << rpt_ranges[i].range.first;
			fp << "\t" << rpt_ranges[i].range.second;
			fp << "\t" << rpt_grp_[rpt_ranges[i].rg_id].seq;
			fp << "\t" << rpt_ranges[i].rg_id;
			fp << endl;
		}

		fp.close();
	}
#endif

	// sort by rg_id
	sort(rpt_ranges.begin(), rpt_ranges.begin() + rpt_ranges.size(), 
			cmp_by_rg_id);


	/***********/

	/* rebuild rpt_grp_ */
	EList<RepeatGroup> mgroup;

	mgroup.reserveExact(rpt_grp_.size());
	mgroup.swap(rpt_grp_);

	for (int i = 0; i < rpt_ranges.size() - 1;) {
		if (rpt_ranges[i].rg_id == std::numeric_limits<int>::max()) {
			break;
		}

		int j = i + 1;
		for (; j < rpt_ranges.size(); j++) {
			if (rpt_ranges[i].rg_id != rpt_ranges[j].rg_id) {
				break;
			}
		}

		/* [i, j) has a same rg_id */


		int rg_id = rpt_ranges[i].rg_id;
		rpt_grp_.expand();
		rpt_grp_.back().seq = mgroup[rg_id].seq;
		for (int k = i; k < j; k++) {
			rpt_grp_.back().positions.push_back(rpt_ranges[k].range.first);
		}

		i = j;
	}


	if (flagGrouping) {
		cerr << "CP " << "before grouping " << rpt_grp_.size() << endl;
		int step = rpt_grp_.size() >> 8;

		if (step == 0) {step = 1;}
		cerr << "CP " << step << endl;
		Timer timer(cerr, "Total time for grouping sequences: ", true);

		for (int i = 0; i < rpt_grp_.size() - 1; i++) {

			if (i % step == 0) {
				cerr << "CP " << i << "/" << rpt_grp_.size() << endl;
			}

			if (rpt_grp_[i].empty()) {
				// empty -> skip
				continue;
			}

			string& str1 = rpt_grp_[i].seq;


			for (int j = i + 1; j < rpt_grp_.size(); j++) {
				string& str2 = rpt_grp_[j].seq;

				if (seq_mergeable(str1, str2, rpt_edit)) {
					/* i, j merge into i */
					rpt_grp_[i].merge(rpt_grp_[j]);

					rpt_grp_[j].set_empty();
				}
			}
		}

		mgroup.clear();
		mgroup.reserveExact(rpt_grp_.size());
		mgroup.swap(rpt_grp_);

		for (int i = 0; i < mgroup.size(); i++) {
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

			for (int i = 0; i < rpt_grp_.size(); i++) {
				RepeatGroup& rg = rpt_grp_[i];
				if (rg.empty()) {
					continue;
				}
				fp << "CP " << i ;
				fp << "\t" << rg.seq;
				for (int j = 0; j < rg.alt_seq.size(); j++) {
					fp << "\t" << rg.alt_seq[j];
				}
				fp << endl;
			}

			fp.close();
		}
#endif
	}


}

template<typename TStr>
void NRG<TStr>::repeat_masking(void)
{
	for (int i = 0; i < rpt_grp_.size(); i++) {
		RepeatGroup *rg = &rpt_grp_[i];

		size_t rpt_sqn_len = rg->seq.length();

		for (int j = 0; j < rg->positions.size(); j++) {
			TIndexOffU pos = rg->positions[j];

			// masking [pos, pos + rpt_sqn_len) to 'N'
			masking_with_N(s_, pos, rpt_sqn_len);
		}
	}
}

// Dump
//
// to_string
static string to_string(int val)
{
	stringstream ss;
	ss << val;
	return ss.str();
}


template<typename TStr>
int get_lcp_back(TStr& s, TIndexOffU a, TIndexOffU b)
{
	int k = 0;
	TIndexOffU s_len = s.length();

	if (a == s_len || b == s_len) {
		return 0;
	}

	while ((a - k) > 0 && (b - k) > 0) {
		if (s[a - k - 1] != s[b - k - 1]) {
			break;
		}
		k++;
	}

	return k;
}

#if 0
int get_lcp(string a, string b)
{
    int k = 0;
    int a_len = a.length();
    int b_len = b.length();
    
    while (k < a_len && k < b_len) {
        if (a[k] != b[k]) {
            break;
        }
        k++;
    }
    
    return k;
}
#endif

#if 0
template<typename TStr>
int NRG<TStr>::get_lcp(TIndexOffU a, TIndexOffU b)
{
    int k = 0;
    TIndexOffU s_len = s_.length();
    
    while ((a + k) < s_len && (b + k) < s_len) {
        if (s_[a + k] != s_[b + k]) {
            break;
        }
        k++;
    }
    
    return k;
}
#else
template<typename TStr>
int NRG<TStr>::get_lcp(TIndexOffU a, TIndexOffU b)
{
    int k = 0;
    TIndexOffU s_len = s_.length();
    
    if (a >= s_len || b >= s_len) {
        return 0;
    }
    
#if 1
    int a_frag_id = map_joined_pos_to_seq(a);
    int b_frag_id = map_joined_pos_to_seq(b);
    
    if (a_frag_id < 0 || b_frag_id < 0) {
        cerr << "CP " << a_frag_id << ", " << b_frag_id << endl;
        return 0;
    }
    
    size_t a_end = fraglist_[a_frag_id].start + fraglist_[a_frag_id].length;
    size_t b_end = fraglist_[b_frag_id].start + fraglist_[b_frag_id].length;
    
    assert_leq(a_end, s_len);
    assert_leq(b_end, s_len);
#else
    size_t a_end = s_len;
    size_t b_end = s_len;
#endif
    
    
    while ((a + k) < a_end && (b + k) < b_end) {
        if (s_[a + k] != s_[b + k]) {
            break;
        }
        k++;
    }
    
    return k;
}
#endif




/****************************/
template class NRG<SString<char> >;
