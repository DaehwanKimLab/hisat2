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

#include "aligner_cache.h"
#include "aligner_seed.h"
#include "search_globals.h"
#include "gfm.h"

using namespace std;

/**
 * Construct a constraint with no edits of any kind allowed.
 */
Constraint Constraint::exact() {
	Constraint c;
	c.edits = c.mms = c.ins = c.dels = c.penalty = 0;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint.
 */
Constraint Constraint::penaltyBased(int pen) {
	Constraint c;
	c.penalty = pen;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint related to the length of the read.
 */
Constraint Constraint::penaltyFuncBased(const SimpleFunc& f) {
	Constraint c;
	c.penFunc = f;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint.
 */
Constraint Constraint::mmBased(int mms) {
	Constraint c;
	c.mms = mms;
	c.edits = c.dels = c.ins = 0;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint.
 */
Constraint Constraint::editBased(int edits) {
	Constraint c;
	c.edits = edits;
	c.dels = c.ins = c.mms = 0;
	return c;
}

//
// Some static methods for constructing some standard SeedPolicies
//

/**
 * Given a read, depth and orientation, extract a seed data structure
 * from the read and fill in the steps & zones arrays.  The Seed
 * contains the sequence and quality values.
 */
bool
Seed::instantiate(
	const Read& read,
	const BTDnaString& seq, // seed read sequence
	const BTString& qual,   // seed quality sequence
	const Scoring& pens,
	int depth,
	int seedoffidx,
	int seedtypeidx,
	bool fw,
	InstantiatedSeed& is) const
{
	assert(overall != NULL);
	int seedlen = len;
	if((int)read.length() < seedlen) {
		// Shrink seed length to fit read if necessary
		seedlen = (int)read.length();
	}
	assert_gt(seedlen, 0);
	is.steps.resize(seedlen);
	is.zones.resize(seedlen);
	// Fill in 'steps' and 'zones'
	//
	// The 'steps' list indicates which read character should be
	// incorporated at each step of the search process.  Often we will
	// simply proceed from one end to the other, in which case the
	// 'steps' list is ascending or descending.  In some cases (e.g.
	// the 2mm case), we might want to switch directions at least once
	// during the search, in which case 'steps' will jump in the
	// middle.  When an element of the 'steps' list is negative, this
	// indicates that the next
	//
	// The 'zones' list indicates which zone constraint is active at
	// each step.  Each element of the 'zones' list is a pair; the
	// first pair element indicates the applicable zone when
	// considering either mismatch or delete (ref gap) events, while
	// the second pair element indicates the applicable zone when
	// considering insertion (read gap) events.  When either pair
	// element is a negative number, that indicates that we are about
	// to leave the zone for good, at which point we may need to
	// evaluate whether we have reached the zone's budget.
	//
	switch(type) {
		case SEED_TYPE_EXACT: {
			for(int k = 0; k < seedlen; k++) {
				is.steps[k] = -(seedlen - k);
				// Zone 0 all the way
				is.zones[k].first = is.zones[k].second = 0;
			}
			break;
		}
		case SEED_TYPE_LEFT_TO_RIGHT: {
			for(int k = 0; k < seedlen; k++) {
				is.steps[k] = k+1;
				// Zone 0 from 0 up to ceil(len/2), then 1
				is.zones[k].first = is.zones[k].second = ((k < (seedlen+1)/2) ? 0 : 1);
			}
			// Zone 1 ends at the RHS
			is.zones[seedlen-1].first = is.zones[seedlen-1].second = -1;
			break;
		}
		case SEED_TYPE_RIGHT_TO_LEFT: {
			for(int k = 0; k < seedlen; k++) {
				is.steps[k] = -(seedlen - k);
				// Zone 0 from 0 up to floor(len/2), then 1
				is.zones[k].first  = ((k < seedlen/2) ? 0 : 1);
				// Inserts: Zone 0 from 0 up to ceil(len/2)-1, then 1
				is.zones[k].second = ((k < (seedlen+1)/2+1) ? 0 : 1);
			}
			is.zones[seedlen-1].first = is.zones[seedlen-1].second = -1;
			break;
		}
		case SEED_TYPE_INSIDE_OUT: {
			// Zone 0 from ceil(N/4) up to N-floor(N/4)
			int step = 0;
			for(int k = (seedlen+3)/4; k < seedlen - (seedlen/4); k++) {
				is.zones[step].first = is.zones[step].second = 0;
				is.steps[step++] = k+1;
			}
			// Zone 1 from N-floor(N/4) up
			for(int k = seedlen - (seedlen/4); k < seedlen; k++) {
				is.zones[step].first = is.zones[step].second = 1;
				is.steps[step++] = k+1;
			}
			// No Zone 1 if seedlen is short (like 2)
			//assert_eq(1, is.zones[step-1].first);
			is.zones[step-1].first = is.zones[step-1].second = -1;
			// Zone 2 from ((seedlen+3)/4)-1 down to 0
			for(int k = ((seedlen+3)/4)-1; k >= 0; k--) {
				is.zones[step].first = is.zones[step].second = 2;
				is.steps[step++] = -(k+1);
			}
			assert_eq(2, is.zones[step-1].first);
			is.zones[step-1].first = is.zones[step-1].second = -2;
			assert_eq(seedlen, step);
			break;
		}
		default:
			throw 1;
	}
	// Instantiate constraints
	for(int i = 0; i < 3; i++) {
		is.cons[i] = zones[i];
		is.cons[i].instantiate(read.length());
	}
	is.overall = *overall;
	is.overall.instantiate(read.length());
	// Take a sweep through the seed sequence.  Consider where the Ns
	// occur and how zones are laid out.  Calculate the maximum number
	// of positions we can jump over initially (e.g. with the ftab) and
	// perhaps set this function's return value to false, indicating
	// that the arrangements of Ns prevents the seed from aligning.
	bool streak = true;
	is.maxjump = 0;
	bool ret = true;
	bool ltr = (is.steps[0] > 0); // true -> left-to-right
	for(size_t i = 0; i < is.steps.size(); i++) {
		assert_neq(0, is.steps[i]);
		int off = is.steps[i];
		off = abs(off)-1;
		Constraint& cons = is.cons[abs(is.zones[i].first)];
		int c = seq[off];  assert_range(0, 4, c);
		int q = qual[off];
		if(ltr != (is.steps[i] > 0) || // changed direction
		   is.zones[i].first != 0 ||   // changed zone
		   is.zones[i].second != 0)    // changed zone
		{
			streak = false;
		}
		if(c == 4) {
			// Induced mismatch
			if(cons.canN(q, pens)) {
				cons.chargeN(q, pens);
			} else {
				// Seed disqualified due to arrangement of Ns
				return false;
			}
		}
		if(streak) is.maxjump++;
	}
	is.seedoff = depth;
	is.seedoffidx = seedoffidx;
	is.fw = fw;
	is.s = *this;
	return ret;
}

/**
 * Return a set consisting of 1 seed encapsulating an exact matching
 * strategy.
 */
void
Seed::zeroMmSeeds(int ln, EList<Seed>& pols, Constraint& oall) {
	oall.init();
	// Seed policy 1: left-to-right search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_EXACT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::exact();
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
}

/**
 * Return a set of 2 seeds encapsulating a half-and-half 1mm strategy.
 */
void
Seed::oneMmSeeds(int ln, EList<Seed>& pols, Constraint& oall) {
	oall.init();
	// Seed policy 1: left-to-right search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_LEFT_TO_RIGHT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(1);
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
	// Seed policy 2: right-to-left search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_RIGHT_TO_LEFT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(1);
	pols.back().zones[1].mmsCeil = 0;
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
}

/**
 * Return a set of 3 seeds encapsulating search roots for:
 *
 * 1. Starting from the left-hand side and searching toward the
 *    right-hand side allowing 2 mismatches in the right half.
 * 2. Starting from the right-hand side and searching toward the
 *    left-hand side allowing 2 mismatches in the left half.
 * 3. Starting (effectively) from the center and searching out toward
 *    both the left and right-hand sides, allowing one mismatch on
 *    either side.
 *
 * This is not exhaustive.  There are 2 mismatch cases mised; if you
 * imagine the seed as divided into four successive quarters A, B, C
 * and D, the cases we miss are when mismatches occur in A and C or B
 * and D.
 */
void
Seed::twoMmSeeds(int ln, EList<Seed>& pols, Constraint& oall) {
	oall.init();
	// Seed policy 1: left-to-right search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_LEFT_TO_RIGHT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(2);
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
	// Seed policy 2: right-to-left search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_RIGHT_TO_LEFT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(2);
	pols.back().zones[1].mmsCeil = 1; // Must have used at least 1 mismatch
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
	// Seed policy 3: inside-out search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_INSIDE_OUT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(1);
	pols.back().zones[1].mmsCeil = 0; // Must have used at least 1 mismatch
	pols.back().zones[2] = Constraint::mmBased(1);
	pols.back().zones[2].mmsCeil = 0; // Must have used at least 1 mismatch
	pols.back().overall = &oall;
}

/**
 * Types of actions that can be taken by the SeedAligner.
 */
enum {
	SA_ACTION_TYPE_RESET = 1,
	SA_ACTION_TYPE_SEARCH_SEED, // 2
	SA_ACTION_TYPE_FTAB,        // 3
	SA_ACTION_TYPE_FCHR,        // 4
	SA_ACTION_TYPE_MATCH,       // 5
	SA_ACTION_TYPE_EDIT         // 6
};

#define MIN(x, y) ((x < y) ? x : y)

#ifdef ALIGNER_SEED_MAIN

#include <getopt.h>
#include <string>

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static int parseInt(const char *errmsg, const char *arg) {
	long l;
	char *endPtr = NULL;
	l = strtol(arg, &endPtr, 10);
	if (endPtr != NULL) {
		return (int32_t)l;
	}
	cerr << errmsg << endl;
	throw 1;
	return -1;
}

enum {
	ARG_NOFW = 256,
	ARG_NORC,
	ARG_MM,
	ARG_SHMEM,
	ARG_TESTS,
	ARG_RANDOM_TESTS,
	ARG_SEED
};

static const char *short_opts = "vCt";
static struct option long_opts[] = {
	{(char*)"verbose",  no_argument,       0, 'v'},
	{(char*)"color",    no_argument,       0, 'C'},
	{(char*)"timing",   no_argument,       0, 't'},
	{(char*)"nofw",     no_argument,       0, ARG_NOFW},
	{(char*)"norc",     no_argument,       0, ARG_NORC},
	{(char*)"mm",       no_argument,       0, ARG_MM},
	{(char*)"shmem",    no_argument,       0, ARG_SHMEM},
	{(char*)"tests",    no_argument,       0, ARG_TESTS},
	{(char*)"random",   required_argument, 0, ARG_RANDOM_TESTS},
	{(char*)"seed",     required_argument, 0, ARG_SEED},
};

static void printUsage(ostream& os) {
	os << "Usage: ac [options]* <index> <patterns>" << endl;
	os << "Options:" << endl;
	os << "  --mm                memory-mapped mode" << endl;
	os << "  --shmem             shared memory mode" << endl;
	os << "  --nofw              don't align forward-oriented read" << endl;
	os << "  --norc              don't align reverse-complemented read" << endl;
	os << "  -t/--timing         show timing information" << endl;
	os << "  -C/--color          colorspace mode" << endl;
	os << "  -v/--verbose        talkative mode" << endl;
}

bool gNorc = false;
bool gNofw = false;
bool gColor = false;
int gVerbose = 0;
int gGapBarrier = 1;
bool gColorExEnds = true;
int gSnpPhred = 30;
bool gReportOverhangs = true;

extern void aligner_seed_tests();
extern void aligner_random_seed_tests(
	int num_tests,
	uint32_t qslo,
	uint32_t qshi,
	bool color,
	uint32_t seed);

/**
 * A way of feeding simply tests to the seed alignment infrastructure.
 */
int main(int argc, char **argv) {
	bool useMm = false;
	bool useShmem = false;
	bool mmSweep = false;
	bool noRefNames = false;
	bool sanity = false;
	bool timing = false;
	int option_index = 0;
	int seed = 777;
	int next_option;
	do {
		next_option = getopt_long(
			argc, argv, short_opts, long_opts, &option_index);
		switch (next_option) {
			case 'v':       gVerbose = true; break;
			case 'C':       gColor   = true; break;
			case 't':       timing   = true; break;
			case ARG_NOFW:  gNofw    = true; break;
			case ARG_NORC:  gNorc    = true; break;
			case ARG_MM:    useMm    = true; break;
			case ARG_SHMEM: useShmem = true; break;
			case ARG_SEED:  seed = parseInt("", optarg); break;
			case ARG_TESTS: {
				aligner_seed_tests();
				aligner_random_seed_tests(
					100,     // num references
					100,   // queries per reference lo
					400,   // queries per reference hi
					false, // true -> generate colorspace reference/reads
					18);   // pseudo-random seed
				return 0;
			}
			case ARG_RANDOM_TESTS: {
				seed = parseInt("", optarg);
				aligner_random_seed_tests(
					100,   // num references
					100,   // queries per reference lo
					400,   // queries per reference hi
					false, // true -> generate colorspace reference/reads
					seed); // pseudo-random seed
				return 0;
			}
			case -1: break;
			default: {
				cerr << "Unknown option: " << (char)next_option << endl;
				printUsage(cerr);
				exit(1);
			}
		}
	} while(next_option != -1);
	char *reffn;
	if(optind >= argc) {
		cerr << "No reference; quitting..." << endl;
		return 1;
	}
	reffn = argv[optind++];
	if(optind >= argc) {
		cerr << "No reads; quitting..." << endl;
		return 1;
	}
	string gfmBase(reffn);
	BitPairReference ref(
		gfmBase,     // base path
		gColor,      // whether we expect it to be colorspace
		sanity,      // whether to sanity-check reference as it's loaded
		NULL,        // fasta files to sanity check reference against
		NULL,        // another way of specifying original sequences
		false,       // true -> infiles (2 args ago) contains raw seqs
		useMm,       // use memory mapping to load index?
		useShmem,    // use shared memory (not memory mapping)
		mmSweep,     // touch all the pages after memory-mapping the index
		gVerbose,    // verbose
		gVerbose);   // verbose but just for startup messages
	Timer *t = new Timer(cerr, "Time loading fw index: ", timing);
	GFM gfmFw(
		gfmBase,
		0,           // don't need entireReverse for fw index
		true,        // index is for the forward direction
		-1,          // offrate (irrelevant)
		useMm,       // whether to use memory-mapped files
		useShmem,    // whether to use shared memory
		mmSweep,     // sweep memory-mapped files
		!noRefNames, // load names?
		false,       // load SA sample?
		true,        // load ftab?
		true,        // load rstarts?
		NULL,        // reference map, or NULL if none is needed
		gVerbose,    // whether to be talkative
		gVerbose,    // talkative during initialization
		false,       // handle memory exceptions, don't pass them up
		sanity);
	delete t;
	t = new Timer(cerr, "Time loading bw index: ", timing);
	GFM gfmBw(
		gfmBase + ".rev",
		1,           // need entireReverse
		false,       // index is for the backward direction
		-1,          // offrate (irrelevant)
		useMm,       // whether to use memory-mapped files
		useShmem,    // whether to use shared memory
		mmSweep,     // sweep memory-mapped files
		!noRefNames, // load names?
		false,       // load SA sample?
		true,        // load ftab?
		false,       // load rstarts?
		NULL,        // reference map, or NULL if none is needed
		gVerbose,    // whether to be talkative
		gVerbose,    // talkative during initialization
		false,       // handle memory exceptions, don't pass them up
		sanity);
	delete t;
	for(int i = optind; i < argc; i++) {
	}
}
#endif
