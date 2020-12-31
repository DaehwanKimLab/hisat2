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
#include "filebuf.h"
#include "ds.h"
#include "quant.h"

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
static int showVersion;
static int nthreads;      // number of pthreads operating concurrently 
static string wrapper;        // type of wrapper script, so we can print correct usage
static int  bigEndian;
static int seed;

bool bTranscriptome;    // run Transcriptome alignment

static void resetOptions() {
  verbose        = true;  // be talkative (default)
  sanityCheck    = 0;     // do slow sanity checks
  wrapper.clear();
  bigEndian      = 0;  // little endian
  seed           = 0;     // srandom seed
}

// Argument constants for getopts
enum {
  ARG_USAGE = 256,
  ARG_SEED,
  ARG_WRAPPER,
};

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
  out << "HISAT2 version " << string(HISAT2_VERSION).c_str() << " by Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>" << endl;
  
  string tool_name = "hisat2-quant";
  out << "Usage: " << tool_name << " [options]* <sam_in>" << endl
      << "    reference_in            comma-separated list of SAM files" << endl
      << "Options:" << endl
      << "    -p <int>                number of threads" << endl
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
  {(char*)"version",        no_argument,       &showVersion, 1},
  {(char*)"seed",           required_argument, 0,            ARG_SEED},
  {(char*)"wrapper",        required_argument, 0,            ARG_WRAPPER},
  {(char*)0, 0, 0, 0} // terminator
};

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', then output the given error message and
 * exit with an error and a usage message.
 */
template<typename T>
static T parseNumber(const char *str, T lower, const char *errmsg)
{
    char *endPtr= NULL;
    T t = (T)strtoll(str, &endPtr, 10);
    if(endPtr != NULL) {
        if(t < lower) {
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

template<typename T>
static T parseNumber(T lower, const char *errmsg) {
    return parseNumber(optarg, lower, errmsg);
}

/**
 * Read command-line arguments
 */
static void parseOptions(int argc, const char **argv) {
  int option_index = 0;
  int next_option;
  do {
    next_option = getopt_long(argc, const_cast<char**>(argv),
			      short_options, long_options, &option_index);
    switch (next_option) {
    case ARG_WRAPPER:
      wrapper = optarg;
      break;
    case 'h':
    case ARG_USAGE:
      printUsage(cout);
      throw 0;
      break;
    case 'p':
      nthreads = parseNumber<int>(1, "-p arg must be at least 1");
      break;
    case ARG_SEED:
      seed = parseNumber<int>(0, "--seed arg must be at least 0");
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
}

/**
 * Drive the index construction process and optionally sanity-check the
 * result.
 */
template<typename TStr>
static void driver(
                   const string& infile,
                   vector<string>& infiles,
                   const string& outfile)
{
    assert_gt(infiles.size(), 0);

    Quant quant;
    quant.init(infiles);
}

static const char *argv0 = NULL;

extern "C" {
/**
 * main function.  Parses command-line arguments.
 */
int hisat2_quant(int argc, const char **argv) {
  string outfile;
  try {
    // Reset all global state, including getopt state
    opterr = optind = 1;
    resetOptions();
    
    string infile;
    vector<string> infiles;
    
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
      cerr << "Settings:" << endl;
      // << "  Output files: \"" << outfile.c_str() << ".*." << gfm_ext << "\"" << endl;
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
	driver<SString<char> >(infile, infiles, outfile);
      } catch(bad_alloc& e) {
	throw e;
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
