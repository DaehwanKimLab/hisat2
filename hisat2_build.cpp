/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <getopt.h>
#include "assert_helpers.h"
#include "endian_swap.h"
#include "bt2_idx.h"
#include "hier_idx.h"
#include "formats.h"
#include "sequence_io.h"
#include "tokenize.h"
#include "timer.h"
#include "ref_read.h"
#include "filebuf.h"
#include "reference.h"
#include "ds.h"

/**
 * \file Driver for the bowtie-build indexing tool.
 */

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include <misc/utils.h>

#include "graph.h"


using namespace CSA;


/*
 Builds an automaton from aligned DNA sequences. Accepted characters are the
 bases ('A', 'C', 'G', 'T'), 'N' denoting any base, and '-' denoting a gap.
 
 Input format: n lines of m characters each, separated by '\n'. Use
 clean_alignment first to reformat the alignment into supported format.
 
 The program considers a context of k bases after each position.
 Two bases are merged to a common state, if they
 
 1) are identical,
 2) occur in the same position in the alignment, and
 3) share a context
 
 The reason for using post-context is that the automaton must be
 deterministic in reverse direction.
 */

const std::string alphabet = "ACGNTZ";
const usint alphabet_size = 7;

struct Node
{
    char  label;  // Character value for the node.
    usint key;    // Label + context.
    usint context_filter;
    
    usint  seq, pos;     // Current coordinates.
    usint  context_end;  // Last pos of the context.
    char*  data;
    usint  sequences;
    
    usint* positions;
    
    void init(char* _data, usint _sequences, usint _seq, usint context_length, usint* _positions)
    {
        this->context_filter = alphabet_size;
        this->seq = _seq; this->pos = 0;
        this->context_end = 0;
        this->data = _data; this->sequences = _sequences;
        this->positions = _positions;
        
        usint ptr = this->pos * (this->sequences + 1) + this->seq;
        this->label = this->data[ptr];
        
        for(usint i = 0; i < context_length; i++) { this->context_filter *= alphabet_size; }
        this->key = this->rank(this->label);
        for(usint i = 0; i < context_length; i++) { this->updateKey(); }
    }
    
    void print()
    {
        std::cout << "(" << this->label << ", " << this->key << ", " << this->context_filter << ", " << this->seq << ", " << this->pos << ", " << this->context_end << ")" << std::endl;
    }
    
    usint getPtr()     { return this->pos * (this->sequences + 1) + this->seq; }
    usint getContext() { return this->context_end * (this->sequences + 1) + this->seq; }
    bool atEnd()       { return (this->label == '\0'); }
    bool isGap()       { return (this->label == '-'); }
    
    GraphNode getNode() { return GraphNode(this->label, this->positions[this->pos]); }
    
    usint rank(char c)
    {
        usint temp = alphabet.find(c);
        return ((temp != std::string::npos) ? temp + 1 : 0);
    }
    
    void advance()
    {
        if(this->atEnd()) { return; }
        this->pos++;
        this->label = this->data[this->getPtr()];
        if(this->isGap()) { return; }
        
        this->updateKey();
    }
    
    void updateKey()
    {
        usint ptr = this->getContext();
        this->key = (alphabet_size * this->key) % this->context_filter;
        while(this->data[ptr] != '\0')
        {
            ptr += this->sequences + 1; this->context_end++;
            if(this->data[ptr] != '-') { break; }
        }
        this->key += this->rank(this->data[ptr]);
    }
};


struct NodeComparator
{
    bool operator() (const Node* a, const Node* b)
    {
        return ((a->key < b->key) || (a->key == b->key && a->seq < b->seq));
    }
} node_context_comparator;


void
buildAutomaton(std::ifstream& input, const std::string& base_name, usint sequences, usint context_length, bool backbone)
{
    usint  lines = fileSize(input) / (sequences + 1) + 2;
    usint  size = lines * (sequences + 1);
    char* data = new char[size];
    char* real_data = data + sequences + 1;
    char* data_end = data + size;
    char* real_end = data_end - (sequences + 1);
    
    // Read the input file.
    input.clear();
    input.seekg(0, std::ios::beg);
    for(usint i = 0; i < sequences; i++) { data[i] = 'Z'; }
    input.read(real_data, real_end - real_data);
    for(char* ptr = real_end; ptr < data_end; ++ptr) { *ptr = '\0'; }
    input.close();
    
    // Map positions in the multiple alignment to positions in the first sequence (backbone).
    // Initial node 0, sequence 1..n, final node n+1.
    // Alignment position that correspond to no position in the backbone are mapped to
    // separate positions starting from n+2 (or more).
    usint* positions = new usint[lines];
    if(backbone)
    {
        for(usint i = 0, pos = 0, extra = 0; i < lines; i++)
        {
            if(data[i * (sequences + 1)] == '-') { positions[i] = lines + extra; extra++; }
            else { positions[i] = pos; pos++; }
        }
    }
    else { for(usint i = 0; i < lines; i++) { positions[i] = i; } }
    
    Node nodes[sequences];
    usint previous[sequences]; // Previous node id for every sequence.
    for(usint i = 0; i < sequences; i++)
    {
        nodes[i].init(data, sequences, i, context_length, positions);
        previous[i] = 0;
    }
    std::vector<GraphNode> node_vector;
    std::vector<pair_type> edge_buffer;
    std::vector<GraphEdge> edge_vector;
    node_vector.push_back(GraphNode('Z', positions[0])); // Initial node.
    
    for(usint line = 1; line < lines; line++)
    {
        // Advance all nodes.
        Node* active[sequences];
        usint active_nodes = 0;
        for(usint i = 0; i < sequences; i++)
        {
            nodes[i].advance();
            if(!(nodes[i].isGap())) { active[active_nodes] = nodes + i; active_nodes++; }
        }
        
        // Create nodes and temporary edges.
        sequentialSort(active, active + active_nodes, node_context_comparator);
        edge_buffer.clear();
        for(usint i = 0; i < active_nodes; i++)
        {
            if(i == 0 || active[i]->key != active[i - 1]->key)
            {
                GraphNode n = active[i]->getNode();
                if(backbone && active[i]->seq > 0) { n.label = tolower(n.label); }
                node_vector.push_back(n);
            }
            edge_buffer.push_back(pair_type(previous[active[i]->seq], node_vector.size() - 1));
            previous[active[i]->seq] = node_vector.size() - 1;
        }
        
        // Create edges.
        sequentialSort(edge_buffer.begin(), edge_buffer.end());
        for(usint i = 0; i < edge_buffer.size(); i++)
        {
            if(i == 0 || edge_buffer[i] != edge_buffer[i - 1])
            {
                edge_vector.push_back(GraphEdge(edge_buffer[i].first, edge_buffer[i].second));
            }
        }
    }
    
    std::cout << "Lines: " << (lines - 2) << std::endl;
    std::cout << "Nodes: " << node_vector.size() << std::endl;
    std::cout << "Edges: " << edge_vector.size() << std::endl;
    std::cout << std::endl;
    
    delete[] data;
    delete[] positions;
    Graph graph(node_vector, edge_vector);
    graph.write(base_name);
}

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
//   Ebwt parameters
static int32_t lineRate;
static int32_t linesPerSide;
static int32_t offRate;
static int32_t ftabChars;
static int32_t localOffRate;
static int32_t localFtabChars;
static int  bigEndian;
static bool nsToAs;
static bool autoMem;
static bool packed;
static bool writeRef;
static bool justRef;
static bool reverseEach;
static string wrapper;

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
	//   Ebwt parameters
	lineRate       = Ebwt<TIndexOffU>::default_lineRate;
	linesPerSide   = 1;  // 1 64-byte line on a side
	offRate        = 4;  // sample 1 out of 16 SA elts
	ftabChars      = 10; // 10 chars in initial lookup table
    localOffRate   = 3;
    localFtabChars = 6;
	bigEndian      = 0;  // little endian
	nsToAs         = false; // convert reference Ns to As prior to indexing
	autoMem        = true;  // automatically adjust memory usage parameters
	packed         = false; //
	writeRef       = true;  // write compact reference to .3.bt2/.4.bt2
	justRef        = false; // *just* write compact reference, don't index
	reverseEach    = false;
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
	ARG_WRAPPER,
    ARG_LOCAL_OFFRATE,
    ARG_LOCAL_FTABCHARS
};

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "HISAT2 version " << string(HISAT2_VERSION).c_str() << " by Daehwan Kim (infphilo@gmail.com, http://www.ccb.jhu.edu/people/infphilo)" << endl;
    
#ifdef BOWTIE_64BIT_INDEX
	string tool_name = "hisat-build-l";
#else
	string tool_name = "hisat-build-s";
#endif
	if(wrapper == "basic-0") {
		tool_name = "hisat-build";
	}
    
	out << "Usage: hisat2-build [options]* <reference_in> <bt2_index_base>" << endl
	    << "    reference_in            comma-separated list of files with ref sequences" << endl
	    << "    hisat_index_base          write " << gEbwt_ext << " data to files with this dir/basename" << endl
        << "Options:" << endl
        << "    -c                      reference sequences given on cmd line (as" << endl
        << "                            <reference_in>)" << endl;
    if(wrapper == "basic-0") {
        out << "    --large-index           force generated index to be 'large', even if ref" << endl
		<< "                            has fewer than 4 billion nucleotides" << endl;
	}
    out << "    -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting" << endl
	    << "    -p/--packed             use packed strings internally; slower, uses less mem" << endl
	    << "    --bmax <int>            max bucket sz for blockwise suffix-array builder" << endl
	    << "    --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)" << endl
	    << "    --dcv <int>             diff-cover period for blockwise (default: 1024)" << endl
	    << "    --nodc                  disable diff-cover (algorithm becomes quadratic)" << endl
	    << "    -r/--noref              don't build .3/.4.bt2 (packed reference) portion" << endl
	    << "    -3/--justref            just build .3/.4.bt2 (packed reference) portion" << endl
	    << "    -o/--offrate <int>      SA is sampled every 2^offRate BWT chars (default: 5)" << endl
	    << "    -t/--ftabchars <int>    # of chars consumed in initial lookup (default: 10)" << endl
        << "    --localoffrate <int>    SA (local) is sampled every 2^offRate BWT chars (default: 3)" << endl
        << "    --localftabchars <int>  # of chars consumed in initial lookup in a local index (default: 6)" << endl
	    << "    --seed <int>            seed for random number generator" << endl
	    << "    -q/--quiet              verbose output (for debugging)" << endl
	    << "    -h/--help               print detailed description of tool and its options" << endl
	    << "    --usage                 print this usage message" << endl
	    << "    --version               print version information and quit" << endl
	    ;
    
    if(wrapper.empty()) {
		cerr << endl
        << "*** Warning ***" << endl
        << "'" << tool_name << "' was run directly.  It is recommended "
        << "that you run the wrapper script 'bowtie2-build' instead."
        << endl << endl;
	}
}

static const char *short_options = "qraph?nscfl:i:o:t:h:3C";

static struct option long_options[] = {
	{(char*)"quiet",          no_argument,       0,            'q'},
	{(char*)"sanity",         no_argument,       0,            's'},
	{(char*)"packed",         no_argument,       0,            'p'},
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
	{(char*)"offrate",        required_argument, 0,            'o'},
	{(char*)"ftabchars",      required_argument, 0,            't'},
    {(char*)"localoffrate",   required_argument, 0,            ARG_LOCAL_OFFRATE},
	{(char*)"localftabchars", required_argument, 0,            ARG_LOCAL_FTABCHARS},
	{(char*)"help",           no_argument,       0,            'h'},
	{(char*)"ntoa",           no_argument,       0,            ARG_NTOA},
	{(char*)"justref",        no_argument,       0,            '3'},
	{(char*)"noref",          no_argument,       0,            'r'},
	{(char*)"color",          no_argument,       0,            'C'},
    {(char*)"sa",             no_argument,       0,            ARG_SA},
	{(char*)"reverse-each",   no_argument,       0,            ARG_REVERSE_EACH},
	{(char*)"usage",          no_argument,       0,            ARG_USAGE},
    {(char*)"wrapper",        required_argument, 0,            ARG_WRAPPER},
	{(char*)0, 0, 0, 0} // terminator
};

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', then output the given error message and
 * exit with an error and a usage message.
 */
template<typename T>
static int parseNumber(T lower, const char *errmsg) {
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
			case 'p': packed = true; break;
			case 'C':
				cerr << "Error: -C specified but Bowtie 2 does not support colorspace input." << endl;
				throw 1;
				break;
			case 'l':
				lineRate = parseNumber<int>(3, "-l/--lineRate arg must be at least 3");
				break;
			case 'i':
				linesPerSide = parseNumber<int>(1, "-i/--linesPerSide arg must be at least 1");
				break;
			case 'o':
				offRate = parseNumber<int>(0, "-o/--offRate arg must be at least 0");
				break;
            case ARG_LOCAL_OFFRATE:
                localOffRate = parseNumber<int>(0, "-o/--localoffrate arg must be at least 0");
                break;
			case '3':
				justRef = true;
				break;
			case 't':
				ftabChars = parseNumber<int>(1, "-t/--ftabChars arg must be at least 1");
				break;
            case ARG_LOCAL_FTABCHARS:
				localFtabChars = parseNumber<int>(1, "-t/--localftabchars arg must be at least 1");
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
			case ARG_REVERSE_EACH:
				reverseEach = true;
				break;
			case ARG_NTOA: nsToAs = true; break;
			case 'a': autoMem = false; break;
			case 'q': verbose = false; break;
			case 's': sanityCheck = true; break;
			case 'r': writeRef = false; break;

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

EList<string> filesWritten;

/**
 * Delete all the index files that we tried to create.  For when we had to
 * abort the index-building process due to an error.
 */
static void deleteIdxFiles(
	const string& outfile,
	bool doRef,
	bool justRef)
{
	
	for(size_t i = 0; i < filesWritten.size(); i++) {
		cerr << "Deleting \"" << filesWritten[i].c_str()
		     << "\" file written during aborted indexing attempt." << endl;
		remove(filesWritten[i].c_str());
	}
}

extern void initializeCntLut();

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
	EList<FileBuf*> is(MISC_CAT);
	bool bisulfite = false;
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
	std::pair<size_t, size_t> sztot;
	{
		if(verbose) cout << "Reading reference sizes" << endl;
		Timer _t(cout, "  Time reading reference sizes: ", verbose);
		if(!reverse && (writeRef || justRef)) {
			filesWritten.push_back(outfile + ".3." + gEbwt_ext);
			filesWritten.push_back(outfile + ".4." + gEbwt_ext);
			sztot = BitPairReference::szsFromFasta(is, outfile, bigEndian, refparams, szs, sanityCheck);
		} else {
			sztot = BitPairReference::szsFromFasta(is, string(), bigEndian, refparams, szs, sanityCheck);
		}
	}
	if(justRef) return;
	assert_gt(sztot.first, 0);
	assert_gt(sztot.second, 0);
	assert_gt(szs.size(), 0);
	// Construct index from input strings and parameters
	filesWritten.push_back(outfile + ".1." + gEbwt_ext);
	filesWritten.push_back(outfile + ".2." + gEbwt_ext);
	TStr s;
	HierEbwt<TIndexOffU> hierEbwt(
                                  s,
                                  packed,
                                  0,
                                  1,  // TODO: maybe not?
                                  lineRate,
                                  offRate,      // suffix-array sampling rate
                                  ftabChars,    // number of chars in initial arrow-pair calc
                                  localOffRate,
                                  localFtabChars,
                                  outfile,      // basename for .?.ebwt files
                                  reverse == 0, // fw
                                  !entireSA,    // useBlockwise
                                  bmax,         // block size for blockwise SA builder
                                  bmaxMultSqrt, // block size as multiplier of sqrt(len)
                                  bmaxDivN,     // block size as divisor of len
                                  noDc? 0 : dcv,// difference-cover period
                                  is,           // list of input streams
                                  szs,          // list of reference sizes
                                  (TIndexOffU)sztot.first,  // total size of all unambiguous ref chars
                                  refparams,    // reference read-in parameters
                                  seed,         // pseudo-random number generator seed
                                  -1,           // override offRate
                                  verbose,      // be talkative
                                  autoMem,      // pass exceptions up to the toplevel so that we can adjust memory settings automatically
                                  sanityCheck); // verify results and internal consistency
	// Note that the Ebwt is *not* resident in memory at this time.  To
	// load it into memory, call ebwt.loadIntoMemory()
	if(verbose) {
		// Print Ebwt's vital stats
		hierEbwt.eh().print(cout);
	}
	if(sanityCheck) {
		// Try restoring the original string (if there were
		// multiple texts, what we'll get back is the joined,
		// padded string, not a list)
		hierEbwt.loadIntoMemory(
								0,
								reverse ? (refparams.reverse == REF_READ_REVERSE) : 0,
								true,  // load SA sample?
								true,  // load ftab?
								true,  // load rstarts?
								false,
								false);
		SString<char> s2;
		hierEbwt.restore(s2);
		hierEbwt.evictFromMemory();
		{
			SString<char> joinedss = Ebwt<>::join<SString<char> >(
				is,          // list of input streams
				szs,         // list of reference sizes
				(TIndexOffU)sztot.first, // total size of all unambiguous ref chars
				refparams,   // reference read-in parameters
				seed);       // pseudo-random number generator seed
			if(refparams.reverse == REF_READ_REVERSE) {
				joinedss.reverse();
			}
			assert_eq(joinedss.length(), s2.length());
			assert(sstr_eq(joinedss, s2));
		}
		if(verbose) {
			if(s2.length() < 1000) {
				cout << "Passed restore check: " << s2.toZBuf() << endl;
			} else {
				cout << "Passed restore check: (" << s2.length() << " chars)" << endl;
			}
		}
	}
}

static const char *argv0 = NULL;

extern "C" {
/**
 * main function.  Parses command-line arguments.
 */
int hisat2_build(int argc, const char **argv) {
    
    std::cout << "Automaton builder" << std::endl;
    std::cout << std::endl;
    if(argc < 3)
    {
        std::cout << "Usage: build_automaton [-b] alignment_file context_length" << std::endl;
        std::cout << "  -b  create backbone information" << std::endl;
        return 1;
    }
    
    bool backbone = false;
    usint name_arg = 1, context_arg = 2;
    if(argc >= 4 && argv[1][0] == '-' && argv[1][1] == 'b')
    {
        backbone = true; name_arg = 2; context_arg = 3;
    }
    
    std::cout << "Alignment file: " << argv[name_arg] << std::endl;
    std::string base_name = argv[name_arg];
    std::ifstream alignment_file(argv[name_arg], std::ios_base::binary);
    if(!alignment_file)
    {
        std::cerr << "Error opening alignment file!" << std::endl;
        return 2;
    }
    
    usint context_length = atoi(argv[context_arg]);
    std::cout << "Context length: " << context_length << std::endl;
    
    std::string line;
    std::getline(alignment_file, line);
    usint sequences = line.size();
    std::cout << "Number of sequences: " << sequences << std::endl;
    std::cout << std::endl;
    
    double start = readTimer();
    buildAutomaton(alignment_file, base_name, sequences, context_length, backbone);
    double time = readTimer() - start;
    std::cout << "Used " << time << " seconds." << std::endl;
    std::cout << "Memory: " << memoryUsage() << " kB" << std::endl;
    std::cout << std::endl;
    
    return 0;
    
    
    ///////////////////////////
    
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
			cout << argv0 << " version " << string(HISAT_VERSION).c_str() << endl;
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
			cout << "Settings:" << endl
				 << "  Output files: \"" << outfile.c_str() << ".*." << gEbwt_ext << "\"" << endl
				 << "  Line rate: " << lineRate << " (line is " << (1<<lineRate) << " bytes)" << endl
				 << "  Lines per side: " << linesPerSide << " (side is " << ((1<<lineRate)*linesPerSide) << " bytes)" << endl
				 << "  Offset rate: " << offRate << " (one in " << (1<<offRate) << ")" << endl
				 << "  FTable chars: " << ftabChars << endl
				 << "  Strings: " << (packed? "packed" : "unpacked") << endl
                 << "  Local offset rate: " << localOffRate << " (one in " << (1<<localOffRate) << ")" << endl
                 << "  Local fTable chars: " << localFtabChars << endl
                 << "  Local sequence length: " << local_index_size << endl
                 << "  Local sequence overlap between the two consecutive indexes: " << local_index_overlap << endl
				 ;
			if(bmax == OFF_MASK) {
				cout << "  Max bucket size: default" << endl;
			} else {
				cout << "  Max bucket size: " << bmax << endl;
			}
			if(bmaxMultSqrt == OFF_MASK) {
				cout << "  Max bucket size, sqrt multiplier: default" << endl;
			} else {
				cout << "  Max bucket size, sqrt multiplier: " << bmaxMultSqrt << endl;
			}
			if(bmaxDivN == 0xffffffff) {
				cout << "  Max bucket size, len divisor: default" << endl;
			} else {
				cout << "  Max bucket size, len divisor: " << bmaxDivN << endl;
			}
			cout << "  Difference-cover sample period: " << dcv << endl;
			cout << "  Endianness: " << (bigEndian? "big":"little") << endl
				 << "  Actual local endianness: " << (currentlyBigEndian()? "big":"little") << endl
				 << "  Sanity checking: " << (sanityCheck? "enabled":"disabled") << endl;
	#ifdef NDEBUG
			cout << "  Assertions: disabled" << endl;
	#else
			cout << "  Assertions: enabled" << endl;
	#endif
			cout << "  Random seed: " << seed << endl;
			cout << "  Sizeofs: void*:" << sizeof(void*) << ", int:" << sizeof(int) << ", long:" << sizeof(long) << ", size_t:" << sizeof(size_t) << endl;
			cout << "Input files DNA, " << file_format_names[format].c_str() << ":" << endl;
			for(size_t i = 0; i < infiles.size(); i++) {
				cout << "  " << infiles[i].c_str() << endl;
			}
		}
		// Seed random number generator
		srand(seed);
		{
			Timer timer(cout, "Total time for call to driver() for forward index: ", verbose);
			if(!packed) {
				try {
					driver<SString<char> >(infile, infiles, outfile, false, REF_READ_FORWARD);
				} catch(bad_alloc& e) {
					if(autoMem) {
						cerr << "Switching to a packed string representation." << endl;
						packed = true;
					} else {
						throw e;
					}
				}
			}
			if(packed) {
				driver<S2bDnaString>(infile, infiles, outfile, true, REF_READ_FORWARD);
			}
		}
		int reverseType = reverseEach ? REF_READ_REVERSE_EACH : REF_READ_REVERSE;
		srand(seed);
		Timer timer(cout, "Total time for backward call to driver() for mirror index: ", verbose);
		if(!packed) {
			try {
				driver<SString<char> >(infile, infiles, outfile + ".rev", false, reverseType);
			} catch(bad_alloc& e) {
				if(autoMem) {
					cerr << "Switching to a packed string representation." << endl;
					packed = true;
				} else {
					throw e;
				}
			}
		}
		if(packed) {
			driver<S2bDnaString>(infile, infiles, outfile + ".rev", true, reverseType);
		}
		return 0;
	} catch(std::exception& e) {
		cerr << "Error: Encountered exception: '" << e.what() << "'" << endl;
		cerr << "Command: ";
		for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
		cerr << endl;
		deleteIdxFiles(outfile, writeRef || justRef, justRef);
		return 1;
	} catch(int e) {
		if(e != 0) {
			cerr << "Error: Encountered internal HISAT2 exception (#" << e << ")" << endl;
			cerr << "Command: ";
			for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
			cerr << endl;
		}
		deleteIdxFiles(outfile, writeRef || justRef, justRef);
		return e;
	}
}
}
