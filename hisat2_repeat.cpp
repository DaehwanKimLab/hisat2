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
static TIndexOffU seed_length;
static TIndexOffU repeat_count;
static TIndexOffU repeat_length;
static TIndexOffU max_repeat_edit;
static TIndexOffU max_repeat_matchlen;
static bool repeat_indel;
static bool forward_only;
static bool CGtoTG;
static string repeat_str1;
static string repeat_str2;
TIndexOffU max_seed_mm;
TIndexOffU max_seed_repeat;
TIndexOffU max_seed_extlen;

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
    seed_length    = 50;
	repeat_length  = 100;
	repeat_count   = 5;
    max_repeat_edit = 10;
    max_repeat_matchlen = repeat_length / 2; // half of repeat_length
    repeat_indel = false;
    forward_only = false;
    CGtoTG = false;
    max_seed_mm = 5;
    max_seed_repeat = 5;
    max_seed_extlen = 25;
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
    ARG_SEED_LENGTH,
	ARG_REPEAT_LENGTH,
	ARG_REPEAT_CNT,
    ARG_REPEAT_EDIT,
	ARG_REPEAT_MATCHLEN,
	ARG_REPEAT_INDEL,
	ARG_WRAPPER,
	ARG_FORWARD_ONLY,
    ARG_CGTOTG,
    ARG_REPEAT_STR1,
    ARG_REPEAT_STR2,
    ARG_MAX_SEED_MM,
    ARG_MAX_SEED_REPEAT,
    ARG_MAX_SEED_EXTLEN,
};

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "HISAT2 version " << string(HISAT2_VERSION).c_str() << " by Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>" << endl;
    
    string tool_name = "hisat2-repeat";
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
        << "    --seed-length <int>     seed length (default: 50)" << endl
        << "    --repeat-length <int>   minimum repeat length (default: 100)" << endl
        << "    --repeat-count <int>    minimum repeat count (default: 5)" << endl
        << "    --repeat-edit <int>     maximum repeat edit distance (default: 10)" << endl
        << "    --repeat-matchlen <int>" << endl
        << "    --repeat-indel" << endl
        << "    --repeat-str1" << endl
        << "    --repeat-str2" << endl
        << "    --forward-only          use forward strand only" << endl
        << "    --CGtoTG                change CG to TG" << endl
        << "    --max-seed-mm <int>" << endl
        << "    --max-seed-repeat <int>" << endl
        << "    --max-seed-extlen <int>" << endl
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
    {(char*)"seed-length",   required_argument, 0,             ARG_SEED_LENGTH},
	{(char*)"repeat-length",  required_argument, 0,            ARG_REPEAT_LENGTH},
	{(char*)"repeat-count",   required_argument, 0,            ARG_REPEAT_CNT},
    {(char*)"repeat-edit",    required_argument, 0,            ARG_REPEAT_EDIT},
    {(char*)"repeat-matchlen",required_argument, 0,            ARG_REPEAT_MATCHLEN},
	{(char*)"repeat-indel",   no_argument,       0,            ARG_REPEAT_INDEL},
    {(char*)"wrapper",        required_argument, 0,            ARG_WRAPPER},
	{(char*)"forward-only",   no_argument,       0,            ARG_FORWARD_ONLY},
	{(char*)"CGtoTG",   	  no_argument,       0,            ARG_CGTOTG},
	{(char*)"repeat-str1",    required_argument, 0,            ARG_REPEAT_STR1},
	{(char*)"repeat-str2",    required_argument, 0,            ARG_REPEAT_STR2},
	{(char*)"max-seed-mm",    required_argument, 0,            ARG_MAX_SEED_MM},
	{(char*)"max-seed-repeat",required_argument, 0,            ARG_MAX_SEED_REPEAT},
	{(char*)"max-seed-extlen",required_argument, 0,            ARG_MAX_SEED_EXTLEN},
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
            case ARG_SEED_LENGTH:
                seed_length = parseNumber<TIndexOffU>(1, "--seed-length arg must be at least 1");
                break;
            case ARG_REPEAT_LENGTH:
                repeat_length = parseNumber<TIndexOffU>(1, "--repeat-length arg must be at least 1");
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
			case ARG_CGTOTG:
				CGtoTG = true;
				break;
			case ARG_REPEAT_STR1:
				repeat_str1 = optarg;
				break;
			case ARG_REPEAT_STR2:
				repeat_str2 = optarg;
				break;
			case ARG_MAX_SEED_MM:
                max_seed_mm = parseNumber<TIndexOffU>(1, "--max_seed_mm arg must be at least 1");
				break;
			case ARG_MAX_SEED_REPEAT:
                max_seed_repeat = parseNumber<TIndexOffU>(5, "--max_seed_repeat arg must be at least 5");
				break;
			case ARG_MAX_SEED_EXTLEN:
                max_seed_extlen = parseNumber<TIndexOffU>(0, "--max_seed_extlen arg must be at least 0");
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
                   bool forward_only,
				   bool CGtoTG)
{
    initializeCntLut();
    initializeCntBit();
	EList<FileBuf*> is(MISC_CAT);
	bool bisulfite = false;
    bool nsToAs = false;
	RefReadInParams refparams(false, false /* reverse */, nsToAs, bisulfite);
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
        bool both_strand = forward_only ? false : true;
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
                                    both_strand, // include reverse complemented sequence
									CGtoTG); //Change CG to TG
    }

    // Successfully obtained joined reference string
#ifndef NDEBUG
    if(forward_only) {
        assert_geq(s.length(), jlen);
    } else {
        assert_geq(s.length(), jlen << 1);
    }
#endif

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
			NRG<TStr> nrg(s, szs, ref_names, forward_only, bsa, outfile);

            if (repeat_str1.length() && repeat_str2.length()) {
                nrg.doTest(repeat_length,
                           repeat_count,
                           true,
                           max_repeat_edit,
                           max_repeat_matchlen,
                           repeat_str1,
                           repeat_str2);
            } else {
                nrg.build(seed_length,
                          repeat_length,
                          repeat_count,
                          true,
                          max_repeat_edit,
                          max_repeat_matchlen);
                nrg.saveFile();
            }


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
                driver<SString<char> >(infile, infiles, outfile, false, forward_only, CGtoTG);
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
