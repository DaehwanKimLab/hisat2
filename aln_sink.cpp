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

#include <iomanip>
#include <limits>
#include "aln_sink.h"
#include "aligner_seed.h"
#include "util.h"

using namespace std;


/**
 * Initialize state machine with a new read.  The state we start in depends
 * on whether it's paired-end or unpaired.
 */
void ReportingState::nextRead(bool paired) {
	paired_ = paired;
	if(paired) {
		state_ = CONCORDANT_PAIRS;
		doneConcord_ = false;
		doneDiscord_ = p_.discord ? false : true;
		doneUnpair1_ = p_.mixed   ? false : true;
		doneUnpair2_ = p_.mixed   ? false : true;
		exitConcord_ = ReportingState::EXIT_DID_NOT_EXIT;
		exitDiscord_ = p_.discord ?
			ReportingState::EXIT_DID_NOT_EXIT :
			ReportingState::EXIT_DID_NOT_ENTER;
		exitUnpair1_ = p_.mixed   ?
			ReportingState::EXIT_DID_NOT_EXIT :
			ReportingState::EXIT_DID_NOT_ENTER;
		exitUnpair2_ = p_.mixed   ?
			ReportingState::EXIT_DID_NOT_EXIT :
			ReportingState::EXIT_DID_NOT_ENTER;
	} else {
		// Unpaired
		state_ = UNPAIRED;
		doneConcord_ = true;
		doneDiscord_ = true;
		doneUnpair1_ = false;
		doneUnpair2_ = true;
		exitConcord_ = ReportingState::EXIT_DID_NOT_ENTER; // not relevant
		exitDiscord_ = ReportingState::EXIT_DID_NOT_ENTER; // not relevant
		exitUnpair1_ = ReportingState::EXIT_DID_NOT_EXIT;
		exitUnpair2_ = ReportingState::EXIT_DID_NOT_ENTER; // not relevant
	}
	doneUnpair_ = doneUnpair1_ && doneUnpair2_;
	done_ = false;
	nconcord_ = ndiscord_ = nunpair1_ = nunpair2_ = 0;
}

/**
 * Caller uses this member function to indicate that one additional
 * concordant alignment has been found.
 */
bool ReportingState::foundConcordant() {
	assert(paired_);
	assert_geq(state_, ReportingState::CONCORDANT_PAIRS);
	assert(!doneConcord_);
	nconcord_++;
	areDone(nconcord_, doneConcord_, exitConcord_);
	// No need to search for discordant alignments if there are one or more
	// concordant alignments.
	doneDiscord_ = true;
	exitDiscord_ = ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED;
	if(doneConcord_) {
		// If we're finished looking for concordant alignments, do we have to
		// continue on to search for unpaired alignments?  Only if our exit
		// from the concordant stage is EXIT_SHORT_CIRCUIT_M.  If it's
		// EXIT_SHORT_CIRCUIT_k or EXIT_WITH_ALIGNMENTS, we can skip unpaired.
		assert_neq(ReportingState::EXIT_NO_ALIGNMENTS, exitConcord_);
		if(exitConcord_ != ReportingState::EXIT_SHORT_CIRCUIT_M) {
			if(!doneUnpair1_) {
				doneUnpair1_ = true;
				exitUnpair1_ = ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED;
			}
			if(!doneUnpair2_) {
				doneUnpair2_ = true;
				exitUnpair2_ = ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED;
			}
		}
	}
	updateDone();
	return done();
}

/**
 * Caller uses this member function to indicate that one additional unpaired
 * mate alignment has been found for the specified mate.
 */
bool ReportingState::foundUnpaired(bool mate1) {
	assert_gt(state_, ReportingState::NO_READ);
	// Note: it's not right to assert !doneUnpair1_/!doneUnpair2_ here.
	// Even if we're done with finding 
	if(mate1) {
		nunpair1_++;
		// Did we just finish with this mate?
		if(!doneUnpair1_) {
			areDone(nunpair1_, doneUnpair1_, exitUnpair1_);
			if(doneUnpair1_) {
				doneUnpair_ = doneUnpair1_ && doneUnpair2_;
				updateDone();
			}
		}
		if(nunpair1_ > 1) {
			doneDiscord_ = true;
			exitDiscord_ = ReportingState::EXIT_NO_ALIGNMENTS;
		}
	} else {
		nunpair2_++;
		// Did we just finish with this mate?
		if(!doneUnpair2_) {
			areDone(nunpair2_, doneUnpair2_, exitUnpair2_);
			if(doneUnpair2_) {
				doneUnpair_ = doneUnpair1_ && doneUnpair2_;
				updateDone();
			}
		}
		if(nunpair2_ > 1) {
			doneDiscord_ = true;
			exitDiscord_ = ReportingState::EXIT_NO_ALIGNMENTS;
		}
	}
	return done();
}
	
/**
 * Called to indicate that the aligner has finished searching for
 * alignments.  This gives us a chance to finalize our state.
 *
 * TODO: Keep track of short-circuiting information.
 */
void ReportingState::finish() {
	if(!doneConcord_) {
		doneConcord_ = true;
		exitConcord_ =
			((nconcord_ > 0) ?
				ReportingState::EXIT_WITH_ALIGNMENTS :
				ReportingState::EXIT_NO_ALIGNMENTS);
	}
	assert_gt(exitConcord_, EXIT_DID_NOT_EXIT);
	if(!doneUnpair1_) {
		doneUnpair1_ = true;
		exitUnpair1_ =
			((nunpair1_ > 0) ?
				ReportingState::EXIT_WITH_ALIGNMENTS :
				ReportingState::EXIT_NO_ALIGNMENTS);
	}
	assert_gt(exitUnpair1_, EXIT_DID_NOT_EXIT);
	if(!doneUnpair2_) {
		doneUnpair2_ = true;
		exitUnpair2_ =
			((nunpair2_ > 0) ?
				ReportingState::EXIT_WITH_ALIGNMENTS :
				ReportingState::EXIT_NO_ALIGNMENTS);
	}
	assert_gt(exitUnpair2_, EXIT_DID_NOT_EXIT);
	if(!doneDiscord_) {
		// Check if the unpaired alignments should be converted to a single
		// discordant paired-end alignment.
		assert_eq(0, ndiscord_);
		if(nconcord_ == 0 && nunpair1_ == 1 && nunpair2_ == 1) {
			convertUnpairedToDiscordant();
		}
		doneDiscord_ = true;
		exitDiscord_ =
			((ndiscord_ > 0) ?
				ReportingState::EXIT_WITH_ALIGNMENTS :
				ReportingState::EXIT_NO_ALIGNMENTS);
	}
	assert(!paired_ || exitDiscord_ > ReportingState::EXIT_DID_NOT_EXIT);
	doneUnpair_ = done_ = true;
	assert(done());
}
	
/**
 * Populate given counters with the number of various kinds of alignments
 * to report for this read.  Concordant alignments are preferable to (and
 * mutually exclusive with) discordant alignments, and paired-end
 * alignments are preferable to unpaired alignments.
 *
 * The caller also needs some additional information for the case where a
 * pair or unpaired read aligns repetitively.  If the read is paired-end
 * and the paired-end has repetitive concordant alignments, that should be
 * reported, and 'pairMax' is set to true to indicate this.  If the read is
 * paired-end, does not have any conordant alignments, but does have
 * repetitive alignments for one or both mates, then that should be
 * reported, and 'unpair1Max' and 'unpair2Max' are set accordingly.
 *
 * Note that it's possible in the case of a paired-end read for the read to
 * have repetitive concordant alignments, but for one mate to have a unique
 * unpaired alignment.
 */
void ReportingState::getReport(
	uint64_t& nconcordAln, // # concordant alignments to report
	uint64_t& ndiscordAln, // # discordant alignments to report
	uint64_t& nunpair1Aln, // # unpaired alignments for mate #1 to report
	uint64_t& nunpair2Aln, // # unpaired alignments for mate #2 to report
	bool& pairMax,         // repetitive concordant alignments
	bool& unpair1Max,      // repetitive alignments for mate #1
	bool& unpair2Max)      // repetitive alignments for mate #2
	const
{
	nconcordAln = ndiscordAln = nunpair1Aln = nunpair2Aln = 0;
	pairMax = unpair1Max = unpair2Max = false;
	assert_gt(p_.khits, 0);
	assert_gt(p_.mhits, 0);
	if(paired_) {
		// Do we have 1 or more concordant alignments to report?
		if(exitConcord_ == ReportingState::EXIT_SHORT_CIRCUIT_k) {
			// k at random
			assert_geq(nconcord_, (uint64_t)p_.khits);
			nconcordAln = p_.khits;
			return;
		} else if(exitConcord_ == ReportingState::EXIT_SHORT_CIRCUIT_M) {
			assert(p_.msample);
			assert_gt(nconcord_, 0);
			pairMax = true;  // repetitive concordant alignments
			if(p_.mixed) {
				unpair1Max = nunpair1_ > (uint64_t)p_.mhits;
				unpair2Max = nunpair2_ > (uint64_t)p_.mhits;
			}
			// Not sure if this is OK
			nconcordAln = 1; // 1 at random
			return;
		} else if(exitConcord_ == ReportingState::EXIT_WITH_ALIGNMENTS) {
			assert_gt(nconcord_, 0);
			// <= k at random
			nconcordAln = min<uint64_t>(nconcord_, p_.khits);
			return;
		}
		assert(!p_.mhitsSet() || nconcord_ <= (uint64_t)p_.mhits+1);
		
		// Do we have a discordant alignment to report?
		if(exitDiscord_ == ReportingState::EXIT_WITH_ALIGNMENTS) {
			// Report discordant
			assert(p_.discord);
			ndiscordAln = 1;
			return;
		}
	}
	
	assert_neq(ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED, exitUnpair1_);
	assert_neq(ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED, exitUnpair2_);

	if((paired_ && !p_.mixed) || nunpair1_ + nunpair2_ == 0) {
		// Unpaired alignments either not reportable or non-existant
		return;
	}

	// Do we have 1 or more alignments for mate #1 to report?
	if(exitUnpair1_ == ReportingState::EXIT_SHORT_CIRCUIT_k) {
		// k at random
		assert_geq(nunpair1_, (uint64_t)p_.khits);
		nunpair1Aln = p_.khits;
	} else if(exitUnpair1_ == ReportingState::EXIT_SHORT_CIRCUIT_M) {
		assert(p_.msample);
		assert_gt(nunpair1_, 0);
		unpair1Max = true;  // repetitive alignments for mate #1
		nunpair1Aln = 1; // 1 at random
	} else if(exitUnpair1_ == ReportingState::EXIT_WITH_ALIGNMENTS) {
		assert_gt(nunpair1_, 0);
		// <= k at random
		nunpair1Aln = min<uint64_t>(nunpair1_, (uint64_t)p_.khits);
	}
	assert(!p_.mhitsSet() || paired_ || nunpair1_ <= (uint64_t)p_.mhits+1);

	// Do we have 2 or more alignments for mate #2 to report?
	if(exitUnpair2_ == ReportingState::EXIT_SHORT_CIRCUIT_k) {
		// k at random
		nunpair2Aln = p_.khits;
	} else if(exitUnpair2_ == ReportingState::EXIT_SHORT_CIRCUIT_M) {
		assert(p_.msample);
		assert_gt(nunpair2_, 0);
		unpair2Max = true;  // repetitive alignments for mate #1
		nunpair2Aln = 1; // 1 at random
	} else if(exitUnpair2_ == ReportingState::EXIT_WITH_ALIGNMENTS) {
		assert_gt(nunpair2_, 0);
		// <= k at random
		nunpair2Aln = min<uint64_t>(nunpair2_, (uint64_t)p_.khits);
	}
	assert(!p_.mhitsSet() || paired_ || nunpair2_ <= (uint64_t)p_.mhits+1);
}

/**
 * Given the number of alignments in a category, check whether we
 * short-circuited out of the category.  Set the done and exit arguments to
 * indicate whether and how we short-circuited.
 */
inline void ReportingState::areDone(
	uint64_t cnt,    // # alignments in category
	bool& done,      // out: whether we short-circuited out of category
	int& exit) const // out: if done, how we short-circuited (-k? -m? etc)
{
	assert(!done);
	// Have we exceeded the -k limit?
	assert_gt(p_.khits, 0);
	assert_gt(p_.mhits, 0);
	if(cnt >= (uint64_t)p_.khits && !p_.mhitsSet()) {
		done = true;
		exit = ReportingState::EXIT_SHORT_CIRCUIT_k;
	}
	// Have we exceeded the -m or -M limit?
	else if(p_.mhitsSet() && cnt > (uint64_t)p_.mhits) {
		done = true;
		assert(p_.msample);
		exit = ReportingState::EXIT_SHORT_CIRCUIT_M;
	}
}

#ifdef ALN_SINK_MAIN

#include <iostream>

bool testDones(
	const ReportingState& st,
	bool done1,
	bool done2,
	bool done3,
	bool done4,
	bool done5,
	bool done6)
{
	assert(st.doneConcordant()    == done1);
	assert(st.doneDiscordant()    == done2);
	assert(st.doneUnpaired(true)  == done3);
	assert(st.doneUnpaired(false) == done4);
	assert(st.doneUnpaired()      == done5);
	assert(st.done()              == done6);
	assert(st.repOk());
	return true;
}

int main(void) {
	cerr << "Case 1 (simple unpaired 1) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,      // khits
			0,      // mhits
			0,      // pengap
			false,  // msample
			false,  // discord
			false); // mixed
		ReportingState st(rp);
		st.nextRead(false); // unpaired read
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, true, true, true, true));
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(0, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(2, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(2, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(!unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;

	cerr << "Case 2 (simple unpaired 1) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,      // khits
			3,      // mhits
			0,      // pengap
			false,  // msample
			false,  // discord
			false); // mixed
		ReportingState st(rp);
		st.nextRead(false); // unpaired read
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(0, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(0, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;

	cerr << "Case 3 (simple paired 1) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,      // khits
			3,      // mhits
			0,      // pengap
			false,  // msample
			false,  // discord
			false); // mixed
		ReportingState st(rp);
		st.nextRead(true); // unpaired read
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(4, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(4, st.numUnpaired2());
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(4, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(4, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(pairMax);
		assert(!unpair1Max); // because !mixed
		assert(!unpair2Max); // because !mixed
	}
	cerr << "PASSED" << endl;

	cerr << "Case 4 (simple paired 2) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,      // khits
			3,      // mhits
			0,      // pengap
			false,  // msample
			true,   // discord
			true);  // mixed
		ReportingState st(rp);
		st.nextRead(true); // unpaired read
		assert(testDones(st, false, false, false, false, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, false, false, false, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, false, false, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, false, false, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, false, false, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, false, false, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, false, false, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, false, false, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(4, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(4, st.numUnpaired2());
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(4, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(4, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(pairMax);
		assert(unpair1Max);
		assert(unpair2Max);
	}
	cerr << "PASSED" << endl;

	cerr << "Case 5 (potential discordant after concordant) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,      // khits
			3,      // mhits
			0,      // pengap
			false,  // msample
			true,   // discord
			true);  // mixed
		ReportingState st(rp);
		st.nextRead(true);
		assert(testDones(st, false, false, false, false, false, false));
		st.foundUnpaired(true);
		st.foundUnpaired(false);
		st.foundConcordant();
		assert(testDones(st, false, true, false, false, false, false));
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(1, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(1, st.numUnpaired1());
		assert_eq(1, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(1, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(!unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;

	cerr << "Case 6 (true discordant) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,      // khits
			3,      // mhits
			0,      // pengap
			false,  // msample
			true,   // discord
			true);  // mixed
		ReportingState st(rp);
		st.nextRead(true);
		assert(testDones(st, false, false, false, false, false, false));
		st.foundUnpaired(true);
		st.foundUnpaired(false);
		assert(testDones(st, false, false, false, false, false, false));
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(0, st.numConcordant());
		assert_eq(1, st.numDiscordant());
		assert_eq(0, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(1, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(!unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;

	cerr << "Case 7 (unaligned pair & uniquely aligned mate, mixed-mode) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			1,      // khits
			1,      // mhits
			0,      // pengap
			false,  // msample
			true,   // discord
			true);  // mixed
		ReportingState st(rp);
		st.nextRead(true); // unpaired read
		// assert(st.doneConcordant()    == done1);
		// assert(st.doneDiscordant()    == done2);
		// assert(st.doneUnpaired(true)  == done3);
		// assert(st.doneUnpaired(false) == done4);
		// assert(st.doneUnpaired()      == done5);
		// assert(st.done()              == done6);
		st.foundUnpaired(true);
		assert(testDones(st, false, false, false, false, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, false, false, false));
		assert_eq(0, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(2, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		st.finish();
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;

	cerr << "Case 8 (unaligned pair & uniquely aligned mate, NOT mixed-mode) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			1,      // khits
			1,      // mhits
			0,      // pengap
			false,  // msample
			true,   // discord
			false); // mixed
		ReportingState st(rp);
		st.nextRead(true); // unpaired read
		// assert(st.doneConcordant()    == done1);
		// assert(st.doneDiscordant()    == done2);
		// assert(st.doneUnpaired(true)  == done3);
		// assert(st.doneUnpaired(false) == done4);
		// assert(st.doneUnpaired()      == done5);
		// assert(st.done()              == done6);
		st.foundUnpaired(true);
		assert(testDones(st, false, false, true, true, true, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, true, true, false));
		assert_eq(0, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(2, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		st.finish();
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(!unpair1Max); // not really relevant
		assert(!unpair2Max); // not really relevant
	}
	cerr << "PASSED" << endl;

	cerr << "Case 9 (repetitive pair, only one mate repetitive) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			1,      // khits
			1,      // mhits
			0,      // pengap
			true,   // msample
			true,   // discord
			true);  // mixed
		ReportingState st(rp);
		st.nextRead(true); // unpaired read
		// assert(st.doneConcordant()    == done1);
		// assert(st.doneDiscordant()    == done2);
		// assert(st.doneUnpaired(true)  == done3);
		// assert(st.doneUnpaired(false) == done4);
		// assert(st.doneUnpaired()      == done5);
		// assert(st.done()              == done6);
		st.foundConcordant();
		assert(st.repOk());
		st.foundUnpaired(true);
		assert(st.repOk());
		st.foundUnpaired(false);
		assert(st.repOk());
		assert(testDones(st, false, true, false, false, false, false));
		assert(st.repOk());
		st.foundConcordant();
		assert(st.repOk());
		st.foundUnpaired(true);
		assert(st.repOk());
		assert(testDones(st, true, true, true, false, false, false));
		assert_eq(2, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(2, st.numUnpaired1());
		assert_eq(1, st.numUnpaired2());
		st.foundUnpaired(false);
		assert(st.repOk());
		assert(testDones(st, true, true, true, true, true, true));		
		assert_eq(2, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(2, st.numUnpaired1());
		assert_eq(2, st.numUnpaired2());
		st.finish();
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(1, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(pairMax);
		assert(unpair1Max); // not really relevant
		assert(unpair2Max); // not really relevant
	}
	cerr << "PASSED" << endl;
}

#endif /*def ALN_SINK_MAIN*/
