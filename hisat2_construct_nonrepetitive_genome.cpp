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
#include <fstream>
#include <string>
#include <cassert>
#include <getopt.h>
#include <algorithm>
#include "assert_helpers.h"
#include "endian_swap.h"
#include "formats.h"
#include "sequence_io.h"
#include "tokenize.h"
#include "timer.h"
#include "ref_read.h"
#include "filebuf.h"
#include "reference.h"
#include "ds.h"
#include "gfm.h"
#include "aligner_sw.h"
#include "aligner_result.h"
#include "search_globals.h"
#include "scoring.h"
#include "mask.h"

/**
 * \file Driver for the bowtie-build indexing tool.
 */

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

MemoryTally gMemTally;
// Build parameters
int verbose;
static int sanityCheck;
static int format;
static TIndexOffU bmax;
static TIndexOffU bmaxMultSqrt;
static uint32_t bmaxDivN;
static int dcv;
static int noDc;
static int entireSA;
static int seed;
static int showVersion;
//   GFM parameters
static int32_t lineRate;
static bool    lineRate_provided;
static int32_t linesPerSide;
static int32_t offRate;
static int32_t ftabChars;
static int32_t localOffRate;
static int32_t localFtabChars;
static int  bigEndian;
static bool autoMem;
static int nthreads;      // number of pthreads operating concurrently
static string wrapper;
static int repeat_count;
static int repeat_length;

static void resetOptions() {
	verbose        = true;  // be talkative (default)
	sanityCheck    = 0;     // do slow sanity checks
	format         = FASTA; // input sequence format
	bmax           = OFF_MASK; // max blockwise SA bucket size
	bmaxMultSqrt   = OFF_MASK; // same, as multplier of sqrt(n)
	bmaxDivN       = 4;          // same, as divisor of n
	dcv            = 1024;  // bwise SA difference-cover sample sz
	noDc           = 0;     // disable difference-cover sample
	entireSA       = 0;     // 1 = disable blockwise SA
	seed           = 0;     // srandom seed
	showVersion    = 0;     // just print version and quit?
	// GFM parameters
	lineRate       = GFM<TIndexOffU>::default_lineRate_gfm;
    lineRate_provided = false;
	linesPerSide   = 1;  // 1 64-byte line on a side
	offRate        = 4;  // sample 1 out of 16 SA elts
	ftabChars      = 10; // 10 chars in initial lookup table
    localOffRate   = 3;
    localFtabChars = 6;
	bigEndian      = 0;  // little endian
	autoMem        = true;  // automatically adjust memory usage parameters
	nthreads       = 1;
	repeat_length  = 50;
	repeat_count   = 5;
    wrapper.clear();
}

// Argument constants for getopts
enum {
	ARG_BMAX = 256,
	ARG_BMAX_MULT,
	ARG_BMAX_DIV,
	ARG_DCV,
	ARG_SEED,
	ARG_CUTOFF,
	ARG_PMAP,
	ARG_NTOA,
	ARG_USAGE,
	ARG_REVERSE_EACH,
    ARG_SA,
	ARG_REPEAT_LENGTH,
	ARG_REPEAT_CNT,
	ARG_WRAPPER
};

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "HISAT2 version " << string(HISAT2_VERSION).c_str() << " by Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>" << endl;
    
#ifdef BOWTIE_64BIT_INDEX
	string tool_name = "hisat2-construct-nonrepetitive-genome-l";
#else
	string tool_name = "hisat2-construct-nonrepetitive-genome-s";
#endif
	if(wrapper == "basic-0") {
		tool_name = "hisat2-build";
	}
    
	out << "Usage: " << tool_name << " [options]* <reference_in>" << endl
	    << "    reference_in            comma-separated list of files with ref sequences" << endl
        << "Options:" << endl
        << "    -c                      reference sequences given on cmd line (as" << endl
        << "                            <reference_in>)" << endl;
    if(wrapper == "basic-0") {
        out << "    --large-index           force generated index to be 'large', even if ref" << endl
		<< "                            has fewer than 4 billion nucleotides" << endl;
	}
    out << "    -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting" << endl
	    << "    -p <int>                number of threads" << endl
	    << "    --bmax <int>            max bucket sz for blockwise suffix-array builder" << endl
	    << "    --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)" << endl
	    << "    --dcv <int>             diff-cover period for blockwise (default: 1024)" << endl
	    << "    --nodc                  disable diff-cover (algorithm becomes quadratic)" << endl
	    << "    --seed <int>            seed for random number generator" << endl
	    << "    -q/--quiet              disable verbose output (for debugging)" << endl
	    << "    -h/--help               print detailed description of tool and its options" << endl
	    << "    --usage                 print this usage message" << endl
	    << "    --version               print version information and quit" << endl
	    ;
}

static const char *short_options = "qrap:h?nscfl:i:o:t:h:3C";

static struct option long_options[] = {
	{(char*)"quiet",          no_argument,       0,            'q'},
	{(char*)"sanity",         no_argument,       0,            's'},
	{(char*)"threads",        required_argument, 0,            'p'},
	{(char*)"little",         no_argument,       &bigEndian,   0},
	{(char*)"big",            no_argument,       &bigEndian,   1},
	{(char*)"bmax",           required_argument, 0,            ARG_BMAX},
	{(char*)"bmaxmultsqrt",   required_argument, 0,            ARG_BMAX_MULT},
	{(char*)"bmaxdivn",       required_argument, 0,            ARG_BMAX_DIV},
	{(char*)"dcv",            required_argument, 0,            ARG_DCV},
	{(char*)"nodc",           no_argument,       &noDc,        1},
	{(char*)"seed",           required_argument, 0,            ARG_SEED},
	{(char*)"entiresa",       no_argument,       &entireSA,    1},
	{(char*)"version",        no_argument,       &showVersion, 1},
	{(char*)"noauto",         no_argument,       0,            'a'},
	{(char*)"noblocks",       required_argument, 0,            'n'},
	{(char*)"linerate",       required_argument, 0,            'l'},
	{(char*)"linesperside",   required_argument, 0,            'i'},
	{(char*)"usage",          no_argument,       0,            ARG_USAGE},
	{(char*)"repeat-length",  required_argument, 0,            ARG_REPEAT_LENGTH},
	{(char*)"repeat-count",   required_argument, 0,            ARG_REPEAT_CNT},
    {(char*)"wrapper",        required_argument, 0,            ARG_WRAPPER},
	{(char*)0, 0, 0, 0} // terminator
};

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', then output the given error message and
 * exit with an error and a usage message.
 */
template<typename T>
static T parseNumber(T lower, const char *errmsg) {
	char *endPtr= NULL;
	T t = (T)strtoll(optarg, &endPtr, 10);
	if (endPtr != NULL) {
		if (t < lower) {
			cerr << errmsg << endl;
			printUsage(cerr);
			throw 1;
		}
		return t;
	}
	cerr << errmsg << endl;
	printUsage(cerr);
	throw 1;
	return -1;
}

/**
 * Read command-line arguments
 */
static void parseOptions(int argc, const char **argv) {
	int option_index = 0;
	int next_option;
	do {
		next_option = getopt_long(
			argc, const_cast<char**>(argv),
			short_options, long_options, &option_index);
		switch (next_option) {
            case ARG_WRAPPER:
				wrapper = optarg;
				break;
			case 'f': format = FASTA; break;
			case 'c': format = CMDLINE; break;
			//case 'p': packed = true; break;
			case 'C':
				cerr << "Error: -C specified but Bowtie 2 does not support colorspace input." << endl;
				throw 1;
				break;
			case 'l':
				lineRate = parseNumber<int>(3, "-l/--lineRate arg must be at least 3");
                lineRate_provided = true;
				break;
			case 'i':
				linesPerSide = parseNumber<int>(1, "-i/--linesPerSide arg must be at least 1");
				break;
			case 'o':
				offRate = parseNumber<int>(0, "-o/--offRate arg must be at least 0");
				break;
           case 'n':
				// all f-s is used to mean "not set", so put 'e' on end
				bmax = 0xfffffffe;
				break;
			case 'h':
			case ARG_USAGE:
				printUsage(cout);
				throw 0;
				break;
           case ARG_BMAX:
				bmax = parseNumber<TIndexOffU>(1, "--bmax arg must be at least 1");
				bmaxMultSqrt = OFF_MASK; // don't use multSqrt
				bmaxDivN = 0xffffffff;     // don't use multSqrt
				break;
			case ARG_BMAX_MULT:
				bmaxMultSqrt = parseNumber<TIndexOffU>(1, "--bmaxmultsqrt arg must be at least 1");
				bmax = OFF_MASK;     // don't use bmax
				bmaxDivN = 0xffffffff; // don't use multSqrt
				break;
			case ARG_BMAX_DIV:
				bmaxDivN = parseNumber<uint32_t>(1, "--bmaxdivn arg must be at least 1");
				bmax = OFF_MASK;         // don't use bmax
				bmaxMultSqrt = OFF_MASK; // don't use multSqrt
				break;
			case ARG_DCV:
				dcv = parseNumber<int>(3, "--dcv arg must be at least 3");
				break;
			case ARG_SEED:
				seed = parseNumber<int>(0, "--seed arg must be at least 0");
				break;
			case ARG_REPEAT_CNT:
				repeat_count = parseNumber<int>(2, "--seed arg must be at least 2");
				break;
			case ARG_REPEAT_LENGTH:
				repeat_length = parseNumber<int>(5, "--seed arg must be at least 5");
				break;
			case 'a': autoMem = false; break;
			case 'q': verbose = false; break;
			case 's': sanityCheck = true; break;
            case 'p':
                nthreads = parseNumber<int>(1, "-p arg must be at least 1");
                break;

			case -1: /* Done with options. */
				break;
			case 0:
				if (long_options[option_index].flag != 0)
					break;
			default:
				printUsage(cerr);
				throw 1;
		}
	} while(next_option != -1);
	if(bmax < 40) {
		cerr << "Warning: specified bmax is very small (" << bmax << ").  This can lead to" << endl
		     << "extremely slow performance and memory exhaustion.  Perhaps you meant to specify" << endl
		     << "a small --bmaxdivn?" << endl;
	}
}

static EList<bool> stbuf, enbuf;
static BTDnaString btread;
static BTString btqual;
static BTString btref;
static BTString btref2;

static BTDnaString readrc;
static BTString qualrc;

/**
 * Helper function for running a case consisting of a read (sequence
 * and quality), a reference string, and an offset that anchors the 0th
 * character of the read to a reference position.
 */
static void doTestCase(
                       SwAligner&         al,
                       const BTDnaString& read,
                       const BTString&    qual,
                       const BTString&    refin,
                       TRefOff            off,
                       EList<bool>        *en,
                       const Scoring&     sc,
                       TAlScore           minsc,
                       TAlScore           floorsc,
                       SwResult&          res,
                       bool               nsInclusive,
                       bool               filterns,
                       uint32_t           seed)
{
    RandomSource rnd(seed);
    btref2 = refin;
    assert_eq(read.length(), qual.length());
    size_t nrow = read.length();
    TRefOff rfi, rff;
    // Calculate the largest possible number of read and reference gaps given
    // 'minsc' and 'pens'
    size_t maxgaps;
    size_t padi, padf;
    {
        int readGaps = sc.maxReadGaps(minsc, read.length());
        int refGaps = sc.maxRefGaps(minsc, read.length());
        assert_geq(readGaps, 0);
        assert_geq(refGaps, 0);
        int maxGaps = max(readGaps, refGaps);
        padi = 2 * maxGaps;
        padf = maxGaps;
        maxgaps = (size_t)maxGaps;
    }
    size_t nceil = 1; // (size_t)sc.nCeil.f((double)read.length());
    size_t width = 1 + padi + padf;
    rfi = off;
    off = 0;
    // Pad the beginning of the reference with Ns if necessary
    if(rfi < padi) {
        size_t beginpad = (size_t)(padi - rfi);
        for(size_t i = 0; i < beginpad; i++) {
            btref2.insert('N', 0);
            off--;
        }
        rfi = 0;
    } else {
        rfi -= padi;
    }
    assert_geq(rfi, 0);
    // Pad the end of the reference with Ns if necessary
    while(rfi + nrow + padi + padf > btref2.length()) {
        btref2.append('N');
    }
    rff = rfi + nrow + padi + padf;
    // Convert reference string to masks
    for(size_t i = 0; i < btref2.length(); i++) {
        if(toupper(btref2[i]) == 'N' && !nsInclusive) {
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
    size_t solwidth = width;
    if(maxgaps >= solwidth) {
        solwidth = 0;
    } else {
        solwidth -= maxgaps;
    }
    if(en == NULL) {
        enbuf.resize(solwidth);
        enbuf.fill(true);
        en = &enbuf;
    }
    assert_geq(rfi, 0);
    assert_gt(rff, rfi);
    readrc = read;
    qualrc = qual;
    al.initRead(
                read,          // read sequence
                readrc,
                qual,          // read qualities
                qualrc,
                0,             // offset of first character within 'read' to consider
                read.length(), // offset of last char (exclusive) in 'read' to consider
                sc);        // local-alignment score floor
    
    DynProgFramer dpframe(false);  // trimToRef
    size_t readGaps = 10, refGaps = 10, maxhalf = 10;
    DPRect rect;
    dpframe.frameSeedExtensionRect(0,                // ref offset implied by seed hit assuming no gaps
                                   read.length(),    // length of read sequence used in DP table
                                   read.length(),    // length of reference
                                   readGaps,       // max # of read gaps permitted in opp mate alignment
                                   refGaps,        // max # of ref gaps permitted in opp mate alignment
                                   (size_t)nceil,  // # Ns permitted
                                   maxhalf,        // max width in either direction
                                   rect);          // DP rectangle
    assert(rect.repOk());
    
    size_t cminlen = 2000, cpow2 = 4;
    al.initRef(fw,                // whether to align forward or revcomp read
               refidx,            // reference aligned against
               rect,              // DP rectangle
               btref2.wbuf(),     // Reference strings
               rfi,
               rff,
               read.length(),
               sc,                // scoring scheme
               minsc,             // minimum score permitted
               true,              // use 8-bit SSE if possible?
               cminlen,           // minimum length for using checkpointing scheme
               cpow2,             // interval b/t checkpointed diags; 1 << this
               false,             // triangular mini-fills?
               true);              // this is a seed extension - not finding a mate
    
    TAlScore best = 0;
    al.align(rnd, best);
}

/**
 * Another interface for running a case.
 */
static void doTestCase2(
                        SwAligner&         al,
                        const char        *read,
                        const char        *qual,
                        const char        *refin,
                        TRefOff            off,
                        const Scoring&     sc,
                        float              costMinConst,
                        float              costMinLinear,
                        SwResult&          res,
                        bool               nsInclusive = false,
                        bool               filterns = false,
                        uint32_t           seed = 0)
{
    btread.install(read, true);
    TAlScore minsc = (TAlScore)(Scoring::linearFunc(
                                                    btread.length(),
                                                    costMinConst,
                                                    costMinLinear));
    TAlScore floorsc = (TAlScore)(Scoring::linearFunc(
                                                      btread.length(),
                                                      costMinConst,
                                                      costMinLinear));
    btqual.install(qual);
    btref.install(refin);
    doTestCase(
               al,
               btread,
               btqual,
               btref,
               off,
               NULL,
               sc,
               minsc,
               floorsc,
               res,
               nsInclusive,
               filterns,
               seed
               );
}

extern void initializeCntLut();
extern void initializeCntBit();

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


// Dump
//
// to_string
string to_string(int val)
{
	stringstream ss;
	ss << val;
	return ss.str();
}

template<typename TStr>
int get_lcp(TStr& s, TIndexOffU a, TIndexOffU b)
{
    int k = 0;
    TIndexOffU s_len = s.length();
    
    while ((a + k) < s_len && (b + k) < s_len) {
        if (s[a + k] != s[b + k]) {
            break;
        }
        k++;
    }
    
    return k;
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
template<typename TStr>
int get_lcp(TStr& s, TIndexOffU a, TIndexOffU b)
{
    int k = 0;
    TIndexOffU s_len = s.length();
    
    if (a >= s_len || b >= s_len) {
        return 0;
    }
    
#if 1
    int a_frag_id = map_joined_pos_to_seq(joined_fragments, a);
    int b_frag_id = map_joined_pos_to_seq(joined_fragments, b);
    
    if (a_frag_id < 0 || b_frag_id < 0) {
        cerr << "CP " << a_frag_id << ", " << b_frag_id << endl;
        return 0;
    }
    
    size_t a_remains = joined_fragments[a_frag_id].length + joined_fragments[a_frag_id].start;
    size_t b_remains = joined_fragments[b_frag_id].length + joined_fragments[b_frag_id].start;
    
    assert_leq(a_remains, s_len);
    assert_leq(b_remains, s_len);
#else
    size_t a_remains = s_len;
    size_t b_remains = s_len;
#endif
    
    
    while ((a + k) < a_remains && (b + k) < b_remains) {
        if (s[a + k] != s[b + k]) {
            break;
        }
        k++;
    }
    
    return k;
}
#endif

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
void masking_with_N(TStr& s, TIndexOffU start, size_t length)
{
	size_t s_len = s.length();

	for (size_t pos = 0; (pos < length) && (start + pos < s_len); pos++) {
		s[start + pos] = 0x04;
	}
}

// build Non-repetitive Genome
template<typename TStr>
class NRG {
	struct Fragments {
		bool contain(TIndexOffU pos) {
			if (pos >= start && pos < (start + length)) {
				return true;
			}
			return false;
		}

		TIndexOffU start;   // index within joined text
		TIndexOffU length;

		int frag_id;
		int seq_id;
		TIndexOffU start_in_seq;    // index within global 
		TIndexOffU start_in_block;  // index within Fasta Block
		bool first;

		string seq_name;
		string nameline;
	};

	struct RepeatGroup {
		string repeat_sequence;
		string repeat_sequence2;
		EList<TIndexOffU> positions;
	};

public:
	NRG(
		EList<RefRecord>& szs,
		EList<string>& ref_names,
		TStr& s,
		string& filename) :
		szs_(szs), ref_namelines_(ref_names), 
		s_(s), filename_(filename)
	{
		cerr << "NRG: " << filename_ << endl;

		// build ref_names_ from ref_namelines_
		build_names();
		build_joined_fragment();
	
	}


public:
	const int output_width = 60;

	EList<RefRecord>& szs_;
	EList<string>& ref_namelines_;
	EList<string> ref_names_;
	TStr& s_;
	string& filename_;

	// mapping info from joined string to genome
	EList<Fragments> fraglist_;

	//
	EList<RepeatGroup> rpt_grp_;

	// Fragments Cache
#define CACHE_SIZE_JOINEDFRG	10
	Fragments cached_[CACHE_SIZE_JOINEDFRG];
	int num_cached_ = 0;
	int victim_ = 0;	/* round-robin */

public:
	void build_names()
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

	int map_joined_pos_to_seq(TIndexOffU joined_pos)
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
				victim_++;
				victim_ %= CACHE_SIZE_JOINEDFRG;
			}

			return top;
		}

		return -1;
	}

	int get_genome_coord(TIndexOffU joined_pos, 
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

	void build_joined_fragment()
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

	static bool repeat_group_cmp(const RepeatGroup& a, const RepeatGroup& b)
	{
		return a.positions[0] < b.positions[0];
	}

	void sort_rpt_grp()
	{
		if (rpt_grp_.size() > 0) {
			sort(rpt_grp_.begin(), rpt_grp_.begin() + rpt_grp_.size(), repeat_group_cmp);
		}
	}

	void save_repeat_groups()
	{
		string rptinfo_filename = filename_ + ".rptinfo";
		int len = rpt_grp_.size();
		ofstream fp(rptinfo_filename.c_str());

		for (int i = 0; i < len; i++) {
			RepeatGroup& rg = rpt_grp_[i];
			EList<TIndexOffU>& positions = rg.positions;

			// @rpt_name length count
			// pos0 pos1 ... pos15
			// pos16 pos17 ... pos31
			//
			// where pos# = chr_name:pos

			// Header line
			fp << "@" << "rpt_" << i;
			fp << "\t" << rg.repeat_sequence.length();
			fp << "\t" << positions.size();
			fp << "\t" << rg.repeat_sequence;
			fp << endl;

			// Positions
			for (int j = 0; j < positions.size(); j++) {
				if (j && (j % 16 == 0)) {
					fp << endl;
				}

				if (j % 16) {
					fp << "\t";
				}

				string chr_name;
				TIndexOffU pos_in_chr;

				get_genome_coord(positions[j], chr_name, pos_in_chr);

				fp << chr_name << ":" << pos_in_chr;
			}
			fp << endl;
		}		
		fp.close();
	}

	void fill_with_N(ofstream& fp, TIndexOffU start, TIndexOffU end_in_seq)
	{
		// [pos_in_seq, end_in_seq) fill with N
		for (; start < end_in_seq; start++) {
			if (start && (start % output_width == 0)) {
				fp << endl;
			}
			fp << "N";
		}
	}


	// Save single repeatgroup to fp
	void savefile_rpt(ofstream& fp, RepeatGroup& rpt, int rpt_idx)
	{
		// >rpt_###
		// ACCAGGCATATA

		// Name
		stringstream rpt_name;
		rpt_name << "rpt_" << rpt_idx;
		fp << ">" << rpt_name.str() << endl;

		// Sequeunce
		size_t rpt_len = rpt.repeat_sequence.length();
		for (size_t i = 0; i < rpt_len; i += output_width) {
			size_t out_len = std::min((size_t)output_width, (size_t)(rpt_len - i));

			fp << rpt.repeat_sequence.substr(i, out_len) << endl;
		}

	}

	void savefile()
	{
		string nonrpt_name = filename_ + ".nonrpt";
		ofstream fp(nonrpt_name.c_str());

		TIndexOffU ref_name_idx = 0;
		size_t acc_sum_seq = 0;
		size_t acc_sum_joined = 0;

		for (int i = 0; i < szs_.size(); i++) {
			RefRecord *rec = &szs_[i];

			if (rec->first) {
				if (i) {
					fp << endl;
				}
				fp << ">" << ref_namelines_[ref_name_idx] << endl;
				ref_name_idx++;
				acc_sum_seq = 0;
			}

			if (rec->off) {
				fill_with_N(fp, acc_sum_seq, acc_sum_seq + rec->off);
				acc_sum_seq += rec->off;
			}

			if (rec->len) {
				for (TIndexOffU start = acc_sum_joined; start < (acc_sum_joined + rec->len); start++, acc_sum_seq++) {
					if (acc_sum_seq && (acc_sum_seq % output_width == 0)) {
						fp << endl;
					}
					fp << "ACGTN"[s_[start]];
				}
				acc_sum_joined += rec->len;
			}
		}

		fp << endl;

		// save repeat sequnce 
		for (int i = 0; i < rpt_grp_.size(); i++) {
			savefile_rpt(fp, rpt_grp_[i], i);
		}

		fp.close();

		// save rpt infos
		save_repeat_groups();
	}


	void add_repeat_group(string& rpt_seq, string& rpt_seq2, EList<TIndexOffU>& rpt_range)
	{
		// rpt_seq is always > 0
		//
		const int rpt_len = rpt_seq.length();

		for (int i = 0; i < rpt_grp_.size(); i++) {
			RepeatGroup& rg = rpt_grp_[i];
			string& rseq = rg.repeat_sequence;
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
					rg.repeat_sequence = "";
				}
			}
		}

		// add to last
		rpt_grp_.expand();
		rpt_grp_.back().repeat_sequence = rpt_seq;
		rpt_grp_.back().repeat_sequence2 = rpt_seq2;
		rpt_grp_.back().positions = rpt_range;
	}

	void adjust_repeat_group(void)
	{
		// remove empty repeat group
		EList<RepeatGroup> mgroup;


		mgroup.reserveExact(rpt_grp_.size());
		mgroup.swap(rpt_grp_);

		for (int i = 0; i < mgroup.size(); i++) {
			if (mgroup[i].repeat_sequence.length() > 0) {
				rpt_grp_.push_back(mgroup[i]);
			}
		}
	}

	void repeat_masking(void)
	{
		for (int i = 0; i < rpt_grp_.size(); i++) {
			RepeatGroup *rg = &rpt_grp_[i];

			size_t rpt_sqn_len = rg->repeat_sequence.length();

			for (int j = 0; j < rg->positions.size(); j++) {
				TIndexOffU pos = rg->positions[j];

				// masking [pos, pos + rpt_sqn_len) to 'N'
				masking_with_N(s_, pos, rpt_sqn_len);
			}
		}
	}
};

/**
 * Drive the index construction process and optionally sanity-check the
 * result.
 */
template<typename TStr>
static void driver(
                   const string& infile,
                   EList<string>& infiles,
                   const string& outfile,
                   bool packed,
                   int reverse)
{
    initializeCntLut();
    initializeCntBit();
	EList<FileBuf*> is(MISC_CAT);
	bool bisulfite = false;
    bool nsToAs = false;
	RefReadInParams refparams(false, reverse, nsToAs, bisulfite);
	assert_gt(infiles.size(), 0);
	if(format == CMDLINE) {
		// Adapt sequence strings to stringstreams open for input
		stringstream *ss = new stringstream();
		for(size_t i = 0; i < infiles.size(); i++) {
			(*ss) << ">" << i << endl << infiles[i].c_str() << endl;
		}
		FileBuf *fb = new FileBuf(ss);
		assert(fb != NULL);
		assert(!fb->eof());
		assert(fb->get() == '>');
		ASSERT_ONLY(fb->reset());
		assert(!fb->eof());
		is.push_back(fb);
	} else {
		// Adapt sequence files to ifstreams
		for(size_t i = 0; i < infiles.size(); i++) {
			FILE *f = fopen(infiles[i].c_str(), "r");
			if (f == NULL) {
				cerr << "Error: could not open "<< infiles[i].c_str() << endl;
				throw 1;
			}
			FileBuf *fb = new FileBuf(f);
			assert(fb != NULL);
			if(fb->peek() == -1 || fb->eof()) {
				cerr << "Warning: Empty fasta file: '" << infile.c_str() << "'" << endl;
				continue;
			}
			assert(!fb->eof());
			assert(fb->get() == '>');
			ASSERT_ONLY(fb->reset());
			assert(!fb->eof());
			is.push_back(fb);
		}
	}
	if(is.empty()) {
		cerr << "Warning: All fasta inputs were empty" << endl;
		throw 1;
	}
    // Vector for the ordered list of "records" comprising the input
	// sequences.  A record represents a stretch of unambiguous
	// characters in one of the input sequences.
	EList<RefRecord> szs(MISC_CAT);
	EList<string> ref_names;
	std::pair<size_t, size_t> sztot;
	{
		if(verbose) cerr << "Reading reference sizes" << endl;
		Timer _t(cerr, "  Time reading reference sizes: ", verbose);
        sztot = BitPairReference::szsFromFasta(is, outfile, bigEndian, refparams, szs, sanityCheck, &ref_names);
	}
	assert_gt(sztot.first, 0);
	assert_gt(sztot.second, 0);
	assert_gt(szs.size(), 0);

    // Compose text strings into single string
    cerr << "Calculating joined length" << endl;
    TIndexOffU jlen = 0;
    for(unsigned int i = 0; i < szs.size(); i++) {
        jlen += (TIndexOffU)szs[i].len;
    }
    // assert_geq(jlen, sztot);
    
    TStr s;
    {
        cerr << "Reserving space for joined string" << endl;
        s.resize(jlen);
        cerr << "Joining reference sequences" << endl;
        Timer timer(cerr, "  Time to join reference sequences: ", verbose);
        
        s = GFM<TIndexOffU>::join<TStr>(
                                                  is,
                                                  szs,
                                                  (TIndexOffU)sztot.first,
                                                  refparams,
                                                  seed);
    }
        
    // Succesfully obtained joined reference string
    assert_geq(s.length(), jlen);

	// NRG
	NRG<TStr> nrg(szs, ref_names, s, infiles[0]);

    if(bmax != (TIndexOffU)OFF_MASK) {
        // VMSG_NL("bmax according to bmax setting: " << bmax);
    }
    else if(bmaxDivN != (TIndexOffU)OFF_MASK) {
        bmax = max<uint32_t>(jlen / bmaxDivN, 1);
        // VMSG_NL("bmax according to bmaxDivN setting: " << bmax);
    }
    else {
        bmax = (uint32_t)sqrt(s.length());
        // VMSG_NL("bmax defaulted to: " << bmax);
    }
    int iter = 0;
    bool first = true;
    bool passMemExc = false, sanity = false;
    // Look for bmax/dcv parameters that work.
    while(true) {
        if(!first && bmax < 40 && passMemExc) {
            cerr << "Could not find approrpiate bmax/dcv settings for building this index." << endl;
            cerr << "Please try indexing this reference on a computer with more memory." << endl;
            if(sizeof(void*) == 4) {
                cerr << "If this computer has more than 4 GB of memory, try using a 64-bit executable;" << endl
                << "this executable is 32-bit." << endl;
            }
            throw 1;
        }
        if(dcv > 4096) dcv = 4096;
        if((iter % 6) == 5 && dcv < 4096 && dcv != 0) {
            dcv <<= 1; // double difference-cover period
        } else {
            bmax -= (bmax >> 2); // reduce by 25%
        }
        iter++;
        try {
            cerr << "Using parameters --bmax " << bmax << endl;
            if(dcv == 0) {
                cerr << " and *no difference cover*" << endl;
            } else {
                cerr << " --dcv " << dcv << endl;
            }
            {
                cerr << "  Doing ahead-of-time memory usage test" << endl;
                // Make a quick-and-dirty attempt to force a bad_alloc iff
                // we would have thrown one eventually as part of
                // constructing the DifferenceCoverSample
                dcv <<= 1;
                TIndexOffU sz = (TIndexOffU)DifferenceCoverSample<TStr>::simulateAllocs(s, dcv >> 1);
                if(nthreads > 1) sz *= (nthreads + 1);
                AutoArray<uint8_t> tmp(sz, EBWT_CAT);
                dcv >>= 1;
                // Likewise with the KarkkainenBlockwiseSA
                sz = (TIndexOffU)KarkkainenBlockwiseSA<TStr>::simulateAllocs(s, bmax);
                AutoArray<uint8_t> tmp2(sz, EBWT_CAT);
                // Grab another 20 MB out of caution
                AutoArray<uint32_t> extra(20*1024*1024, EBWT_CAT);
                // If we made it here without throwing bad_alloc, then we
                // passed the memory-usage stress test
                cerr << "  Passed!  Constructing with these parameters: --bmax " << bmax << " --dcv " << dcv << endl;
                cerr << "" << endl;
            }
            cerr << "Constructing suffix-array element generator" << endl;
            KarkkainenBlockwiseSA<TStr> bsa(s, bmax, nthreads, dcv, seed, sanity, passMemExc, verbose, outfile);
            assert(bsa.suffixItrIsReset());
            assert_eq(bsa.size(), s.length() + 1);
            cerr << "Converting suffix-array elements to index image" << endl;
            EList<string> suffixes;

            {
                // CP - this is a suffix array
                TIndexOffU count = 0;
                
				EList<TIndexOffU> repeat_range;
				TIndexOffU min_lcp_len = s.length();

                while(count < s.length() + 1) {
                    TIndexOffU saElt = bsa.nextSuffix();
                    count++;

					if (count && (count % 1000000 == 0)) {
						cerr << "SA count " << count << endl;
					}

					if (repeat_range.size() == 0) {
						repeat_range.expand();
						repeat_range.back() = saElt;
					} else {
						TIndexOffU prev_saElt = repeat_range.back();

						// calculate common prefix length between two text.
						//   text1 is started from prev_saElt and text2 is started from saElt
						int lcp_len = get_lcp<TStr>(s, prev_saElt, saElt);
                        
                        // check prev_saElt

						if (lcp_len >= repeat_length) {
							repeat_range.expand();
							repeat_range.back() = saElt;
							if (min_lcp_len > lcp_len) {
								min_lcp_len = lcp_len;
							}

						} else {
							if (repeat_range.size() >= repeat_count) {
								// save ranges
								//cerr << "CP " << "we have " << repeat_range.size() << " continuous range" << ", " << min_lcp_len << endl;
								//

#if 0
								repeat_groups.expand();
								RepeatGroup& rg = repeat_groups.back();
#endif

								repeat_range.sort();

								string ss = get_string(s, prev_saElt, min_lcp_len);

								string ss2 = get_string(s, prev_saElt + min_lcp_len, min_lcp_len );
								ss2 += " " + to_string(min_lcp_len);

								nrg.add_repeat_group(ss, ss2, repeat_range);

#if 0
								rg.positions = repeat_range;
								rg.positions.sort();
								rg.repeat_sequence = get_string(s, prev_saElt, min_lcp_len);
#endif
							}

							// flush previous range
							repeat_range.resize(1);
							repeat_range.back() = saElt;
							min_lcp_len = s.length();
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

				nrg.adjust_repeat_group();
				// we found repeat_group
				cerr << "CP " << nrg.rpt_grp_.size() << " groups found" << endl;

				//dump_tstr(s);
				// Masking repeat sequeuce to 'N'
				nrg.repeat_masking();

				//cerr << "After masking" << endl;
				//dump_tstr(nrg.s_);

				// write to FA
				// sequence
				// repeat sequeuce
				nrg.savefile();

				break;
            }

        } catch(bad_alloc& e) {
            if(passMemExc) {
                cerr << "  Ran out of memory; automatically trying more memory-economical parameters." << endl;
            } else {
                cerr << "Out of memory while constructing suffix array.  Please try using a smaller" << endl
                << "number of blocks by specifying a smaller --bmax or a larger --bmaxdivn" << endl;
                throw 1;
            }
        }
        first = false;
    }
    // assert(repOk());
}

static const char *argv0 = NULL;

extern "C" {
/**
 * main function.  Parses command-line arguments.
 */
int hisat2_construct_nonrepetitive_genome(int argc, const char **argv) {
    string outfile;
	try {
		// Reset all global state, including getopt state
		opterr = optind = 1;
		resetOptions();

		string infile;
		EList<string> infiles(MISC_CAT);

		parseOptions(argc, argv);
		argv0 = argv[0];
		if(showVersion) {
			cout << argv0 << " version " << string(HISAT2_VERSION).c_str() << endl;
			if(sizeof(void*) == 4) {
				cout << "32-bit" << endl;
			} else if(sizeof(void*) == 8) {
				cout << "64-bit" << endl;
			} else {
				cout << "Neither 32- nor 64-bit: sizeof(void*) = " << sizeof(void*) << endl;
			}
			cout << "Built on " << BUILD_HOST << endl;
			cout << BUILD_TIME << endl;
			cout << "Compiler: " << COMPILER_VERSION << endl;
			cout << "Options: " << COMPILER_OPTIONS << endl;
			cout << "Sizeof {int, long, long long, void*, size_t, off_t}: {"
				 << sizeof(int)
				 << ", " << sizeof(long) << ", " << sizeof(long long)
				 << ", " << sizeof(void *) << ", " << sizeof(size_t)
				 << ", " << sizeof(off_t) << "}" << endl;
			return 0;
		}

		// Get input filename
		if(optind >= argc) {
			cerr << "No input sequence or sequence file specified!" << endl;
			printUsage(cerr);
			return 1;
		}
		infile = argv[optind++];
        
		// Get output filename
#if 0
		if(optind >= argc) {
			cerr << "No output file specified!" << endl;
			printUsage(cerr);
			return 1;
		}
		outfile = argv[optind++];
#endif

		tokenize(infile, ",", infiles);
		if(infiles.size() < 1) {
			cerr << "Tokenized input file list was empty!" << endl;
			printUsage(cerr);
			return 1;
		}
        
   		// Optionally summarize
		if(verbose) {
			cerr << "Settings:" << endl
            << "  Output files: \"" << outfile.c_str() << ".*." << gfm_ext << "\"" << endl;
			cerr << "  Endianness: " << (bigEndian? "big":"little") << endl
				 << "  Actual local endianness: " << (currentlyBigEndian()? "big":"little") << endl
				 << "  Sanity checking: " << (sanityCheck? "enabled":"disabled") << endl;
	#ifdef NDEBUG
			cerr << "  Assertions: disabled" << endl;
	#else
			cerr << "  Assertions: enabled" << endl;
	#endif
			cerr << "  Random seed: " << seed << endl;
			cerr << "  Sizeofs: void*:" << sizeof(void*) << ", int:" << sizeof(int) << ", long:" << sizeof(long) << ", size_t:" << sizeof(size_t) << endl;
			cerr << "Input files DNA, " << file_format_names[format].c_str() << ":" << endl;
			for(size_t i = 0; i < infiles.size(); i++) {
				cerr << "  " << infiles[i].c_str() << endl;
			}
		}
		// Seed random number generator
        srand(seed);
        {
            Timer timer(cerr, "Total time for call to driver() for forward index: ", verbose);
            try {
                driver<SString<char> >(infile, infiles, outfile, false, REF_READ_FORWARD);
            } catch(bad_alloc& e) {
                if(autoMem) {
                    cerr << "Switching to a packed string representation." << endl;
                } else {
                    throw e;
                }
            }
        }
        return 0;
    } catch(std::exception& e) {
		cerr << "Error: Encountered exception: '" << e.what() << "'" << endl;
		cerr << "Command: ";
		for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
		cerr << endl;
		return 1;
	} catch(int e) {
		if(e != 0) {
			cerr << "Error: Encountered internal HISAT2 exception (#" << e << ")" << endl;
			cerr << "Command: ";
			for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
			cerr << endl;
		}
		return e;
	}
}
}
