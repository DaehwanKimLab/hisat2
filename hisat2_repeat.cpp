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
#include "repeat_builder.h"

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
static TIndexOffU repeat_count;
static TIndexOffU repeat_length;
static TIndexOffU max_repeat_edit;
static TIndexOffU max_repeat_matchlen;
static bool repeat_indel;
static bool forward_only;

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
    max_repeat_edit = 10;
    max_repeat_matchlen = repeat_length / 2; // half of repeat_length
    repeat_indel = false;
    forward_only = false;
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
    ARG_REPEAT_EDIT,
	ARG_REPEAT_MATCHLEN,
	ARG_REPEAT_INDEL,
	ARG_WRAPPER,
	ARG_FORWARD_ONLY
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
        << "    --repeat-length <int>   minimum repeat length (defaultL 50)" << endl
        << "    --repeat-count <int>    minimum repeat count (default: 5)" << endl
        << "    --repeat-edit <int>     maximum repeat edit distance (default: 10)" << endl
        << "    --repeat-matchlen <int>" << endl
        << "    --repeat-indel" << endl
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
    {(char*)"repeat-edit",    required_argument, 0,            ARG_REPEAT_EDIT},
    {(char*)"repeat-matchlen",required_argument, 0,            ARG_REPEAT_MATCHLEN},
	{(char*)"repeat-indel",   no_argument,       0,            ARG_REPEAT_INDEL},
    {(char*)"wrapper",        required_argument, 0,            ARG_WRAPPER},
	{(char*)"forward-only",   no_argument,       0,            ARG_FORWARD_ONLY},
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
    bool saw_max_repeat_matchlen = false;
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
            case ARG_REPEAT_LENGTH:
                repeat_length = parseNumber<TIndexOffU>(5, "--repeat-length arg must be at least 5");
                break;
			case ARG_REPEAT_CNT:
				repeat_count = parseNumber<TIndexOffU>(2, "--repeat-count arg must be at least 2");
				break;
            case ARG_REPEAT_EDIT:
                max_repeat_edit = parseNumber<TIndexOffU>(0, "--repeat-edit arg must be at least 0");
                break;
            case ARG_REPEAT_MATCHLEN:
                max_repeat_matchlen = parseNumber<TIndexOffU>(0, "--repeat-matchlen arg must be at least 0");
                saw_max_repeat_matchlen = true;
                break;
			case ARG_REPEAT_INDEL:
				repeat_indel = true;
				break;
			case ARG_FORWARD_ONLY:
				forward_only = true;
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

    if(!saw_max_repeat_matchlen) {
        max_repeat_matchlen = repeat_length / 2;
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
        sztot = BitPairReference::szsFromFasta(is, "", bigEndian, refparams, szs, sanityCheck, &ref_names);
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
        cerr << "Joining reference sequences" << endl;
        Timer timer(cerr, "  Time to join reference sequences: ", verbose);
        GFM<TIndexOffU>::join<TStr>(
                                    is,
                                    szs,
                                    (TIndexOffU)sztot.first,
                                    refparams,
                                    seed,
                                    s,
                                    true); // include reverse complemented sequence
    }

    // Successfully obtained joined reference string
    assert_geq(s.length(), jlen << 1);

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
            cerr << "Could not find appropriate bmax/dcv settings for building this index." << endl;
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
            KarkkainenBlockwiseSA<TStr> bsa(s, bmax, nthreads, dcv, seed, sanity, passMemExc, false /* verbose */, outfile);
            
            assert(bsa.suffixItrIsReset());
            assert_eq(bsa.size(), s.length() + 1);

			// NRG
			NRG<TStr> nrg(szs, ref_names, s, outfile, bsa);

			nrg.build(repeat_length,
                      repeat_count,
                      true,
                      max_repeat_edit,
                      max_repeat_matchlen);

            nrg.saveFile();

			break;

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
        
		if(optind >= argc) {
			cerr << "No output file specified!" << endl;
			printUsage(cerr);
			return 1;
		}
		outfile = argv[optind++];

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
