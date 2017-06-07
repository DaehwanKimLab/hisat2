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

#ifndef ALN_SINK_H_
#define ALN_SINK_H_

#include <limits>
#include "read.h"
#include "unique.h"
#include "sam.h"
#include "ds.h"
#include "simple_func.h"
#include "outq.h"
#include <utility>
#include "alt.h"
#include "splice_site.h"

// Forward decl
template <typename index_t>
class SeedResults;

enum {
	OUTPUT_SAM = 1
};

/**
 * Metrics summarizing the work done by the reporter and summarizing
 * the number of reads that align, that fail to align, and that align
 * non-uniquely.
 */
struct ReportingMetrics {

	ReportingMetrics():mutex_m() {
	    reset();
	}

	void reset() {
		init(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	}

	void init(
		uint64_t nread_,
		uint64_t npaired_,
		uint64_t nunpaired_,
		uint64_t nconcord_uni_,
		uint64_t nconcord_uni1_,
		uint64_t nconcord_uni2_,
		uint64_t nconcord_rep_,
		uint64_t nconcord_0_,
		uint64_t ndiscord_,
		uint64_t nunp_0_uni_,
		uint64_t nunp_0_uni1_,
		uint64_t nunp_0_uni2_,
		uint64_t nunp_0_rep_,
		uint64_t nunp_0_0_,
		uint64_t nunp_rep_uni_,
		uint64_t nunp_rep_uni1_,
		uint64_t nunp_rep_uni2_,
		uint64_t nunp_rep_rep_,
		uint64_t nunp_rep_0_,
		uint64_t nunp_uni_,
		uint64_t nunp_uni1_,
		uint64_t nunp_uni2_,
		uint64_t nunp_rep_,
		uint64_t nunp_0_,
		uint64_t sum_best1_,
		uint64_t sum_best2_,
		uint64_t sum_best_)
	{
		nread         = nread_;
		
		npaired       = npaired_;
		nunpaired     = nunpaired_;
		
		nconcord_uni  = nconcord_uni_;
		nconcord_uni1 = nconcord_uni1_;
		nconcord_uni2 = nconcord_uni2_;
		nconcord_rep  = nconcord_rep_;
		nconcord_0    = nconcord_0_;
		
		ndiscord      = ndiscord_;
		
		nunp_0_uni    = nunp_0_uni_;
		nunp_0_uni1   = nunp_0_uni1_;
		nunp_0_uni2   = nunp_0_uni2_;
		nunp_0_rep    = nunp_0_rep_;
		nunp_0_0      = nunp_0_0_;

		nunp_rep_uni  = nunp_rep_uni_;
		nunp_rep_uni1 = nunp_rep_uni1_;
		nunp_rep_uni2 = nunp_rep_uni2_;
		nunp_rep_rep  = nunp_rep_rep_;
		nunp_rep_0    = nunp_rep_0_;

		nunp_uni      = nunp_uni_;
		nunp_uni1     = nunp_uni1_;
		nunp_uni2     = nunp_uni2_;
		nunp_rep      = nunp_rep_;
		nunp_0        = nunp_0_;

		sum_best1     = sum_best1_;
		sum_best2     = sum_best2_;
		sum_best      = sum_best_;
	}
	
	/**
	 * Merge (add) the counters in the given ReportingMetrics object
	 * into this object.  This is the only safe way to update a
	 * ReportingMetrics shared by multiple threads.
	 */
	void merge(const ReportingMetrics& met, bool getLock = false) {
        ThreadSafe ts(&mutex_m, getLock);
		nread         += met.nread;

		npaired       += met.npaired;
		nunpaired     += met.nunpaired;

		nconcord_uni  += met.nconcord_uni;
		nconcord_uni1 += met.nconcord_uni1;
		nconcord_uni2 += met.nconcord_uni2;
		nconcord_rep  += met.nconcord_rep;
		nconcord_0    += met.nconcord_0;

		ndiscord      += met.ndiscord;

		nunp_0_uni    += met.nunp_0_uni;
		nunp_0_uni1   += met.nunp_0_uni1;
		nunp_0_uni2   += met.nunp_0_uni2;
		nunp_0_rep    += met.nunp_0_rep;
		nunp_0_0      += met.nunp_0_0;

		nunp_rep_uni  += met.nunp_rep_uni;
		nunp_rep_uni1 += met.nunp_rep_uni1;
		nunp_rep_uni2 += met.nunp_rep_uni2;
		nunp_rep_rep  += met.nunp_rep_rep;
		nunp_rep_0    += met.nunp_rep_0;

		nunp_uni      += met.nunp_uni;
		nunp_uni1     += met.nunp_uni1;
		nunp_uni2     += met.nunp_uni2;
		nunp_rep      += met.nunp_rep;
		nunp_0        += met.nunp_0;

		sum_best1     += met.sum_best1;
		sum_best2     += met.sum_best2;
		sum_best      += met.sum_best;
	}

	uint64_t  nread;         // # reads
	uint64_t  npaired;       // # pairs
	uint64_t  nunpaired;     // # unpaired reads
	
	// Paired
	
	//  Concordant
	uint64_t  nconcord_uni;  // # pairs with unique concordant alns
	uint64_t  nconcord_uni1; // # pairs with exactly 1 concordant alns
	uint64_t  nconcord_uni2; // # pairs with >1 concordant aln, still unique
	uint64_t  nconcord_rep;  // # pairs with repetitive concordant alns
	uint64_t  nconcord_0;    // # pairs with 0 concordant alns
	//  Discordant
	uint64_t  ndiscord;      // # pairs with 1 discordant aln
	
	//  Unpaired from failed pairs
	uint64_t  nunp_0_uni;    // # unique from nconcord_0_ - ndiscord_
	uint64_t  nunp_0_uni1;   // # pairs with exactly 1 concordant alns
	uint64_t  nunp_0_uni2;   // # pairs with >1 concordant aln, still unique
	uint64_t  nunp_0_rep;    // # repetitive from 
	uint64_t  nunp_0_0;      // # with 0 alignments

	//  Unpaired from repetitive pairs
	uint64_t  nunp_rep_uni;  // # pairs with unique concordant alns
	uint64_t  nunp_rep_uni1; // # pairs with exactly 1 concordant alns
	uint64_t  nunp_rep_uni2; // # pairs with >1 concordant aln, still unique
	uint64_t  nunp_rep_rep;  // # pairs with repetitive concordant alns
	uint64_t  nunp_rep_0;    // # pairs with 0 concordant alns
	
	// Unpaired
	
	uint64_t  nunp_uni;      // # unique from nconcord_0_ - ndiscord_
	uint64_t  nunp_uni1;     // # pairs with exactly 1 concordant alns
	uint64_t  nunp_uni2;     // # pairs with >1 concordant aln, still unique
	uint64_t  nunp_rep;      // # repetitive from 
	uint64_t  nunp_0;        // # with 0 alignments

	
	uint64_t  sum_best1;     // Sum of all the best alignment scores
	uint64_t  sum_best2;     // Sum of all the second-best alignment scores
	uint64_t  sum_best;      // Sum of all the best and second-best

	MUTEX_T mutex_m;
};

// Type for expression numbers of hits
typedef int64_t THitInt;

/**
 * Parameters affecting reporting of alignments, specifically -k & -a,
 * -m & -M.
 */
struct ReportingParams {

	explicit ReportingParams(
                             THitInt khits_,
                             THitInt kseeds_,
                             THitInt mhits_,
                             THitInt pengap_,
                             bool msample_,
                             bool discord_,
                             bool mixed_)
	{
		init(khits_, kseeds_, mhits_, pengap_, msample_, discord_, mixed_);
	}

	void init(
              THitInt khits_,
              THitInt kseeds_,
              THitInt mhits_,
              THitInt pengap_,
              bool msample_,
              bool discord_,
              bool mixed_)
	{
		khits   = khits_;     // -k (or high if -a)
        kseeds  = kseeds_;
		mhits   = ((mhits_ == 0) ? std::numeric_limits<THitInt>::max() : mhits_);
		pengap  = pengap_;
		msample = msample_;
		discord = discord_;
		mixed   = mixed_;
	}
	
#ifndef NDEBUG
	/**
	 * Check that reporting parameters are internally consistent.
	 */
	bool repOk() const {
		assert_geq(khits, 1);
		assert_geq(mhits, 1);
		return true;
	}
#endif
	
	/**
	 * Return true iff a -m or -M limit was set by the user.
	 */
	inline bool mhitsSet() const {
		return mhits < std::numeric_limits<THitInt>::max();
	}
	
	/**
	 * Return a multiplier that indicates how many alignments we might look for
	 * (max).  We can use this to boost parameters like ROWM and POSF
	 * appropriately.
	 */
	inline THitInt mult() const {
		if(mhitsSet()) {
			return mhits+1;
		}
		return khits;
	}

	/**
	 * Given ROWM, POSF thresholds, boost them according to mult().
	 */
	void boostThreshold(SimpleFunc& func) {
		THitInt mul = mult();
		assert_gt(mul, 0);
		if(mul == std::numeric_limits<THitInt>::max()) {
			func.setMin(std::numeric_limits<double>::max());
		} else if(mul > 1) {
			func.mult(mul);
		}
	}
	
	/**
	 * Return true iff we are reporting all hits.
	 */
	bool allHits() const {
		return khits == std::numeric_limits<THitInt>::max();
	}

	// Number of alignments to report
	THitInt khits;
    
    // Number of seeds allowed to extend
    THitInt kseeds;
	
	// Read is non-unique if mhits-1 next-best alignments are within
	// pengap of the best alignment
	THitInt mhits, pengap;
	
	// true if -M is specified, meaning that if the -M ceiling is
	// exceeded, we should report 'khits' alignments chosen at random
	// from those found
	bool msample;
	
	// true iff we should seek and report discordant paired-end alignments for
	// paired-end reads.
	bool discord;

	// true iff we should seek and report unpaired mate alignments when there
	// are paired-end alignments for a paired-end read, or if the number of
	// paired-end alignments exceeds the -m ceiling.
	bool mixed;
};

/**
 * A state machine keeping track of the number and type of alignments found so
 * far.  Its purpose is to inform the caller as to what stage the alignment is
 * in and what categories of alignment are still of interest.  This information
 * should allow the caller to short-circuit some alignment work.  Another
 * purpose is to tell the AlnSinkWrap how many and what type of alignment to
 * report.
 *
 * TODO: This class does not keep accurate information about what
 * short-circuiting took place.  If a read is identical to a previous read,
 * there should be a way to query this object to determine what work, if any,
 * has to be re-done for the new read.
 */
class ReportingState {

public:

	enum {
		NO_READ = 1,        // haven't got a read yet
		CONCORDANT_PAIRS,   // looking for concordant pairs
		DISCORDANT_PAIRS,   // looking for discordant pairs
		UNPAIRED,           // looking for unpaired
		DONE                // finished looking
	};

	// Flags for different ways we can finish out a category of potential
	// alignments.
	
	enum {
		EXIT_DID_NOT_EXIT = 1,        // haven't finished
		EXIT_DID_NOT_ENTER,           // never tried search	
		EXIT_SHORT_CIRCUIT_k,         // -k exceeded
		EXIT_SHORT_CIRCUIT_M,         // -M exceeded
		EXIT_SHORT_CIRCUIT_TRUMPED,   // made irrelevant
		EXIT_CONVERTED_TO_DISCORDANT, // unpair became discord
		EXIT_NO_ALIGNMENTS,           // none found
		EXIT_WITH_ALIGNMENTS          // some found
	};
	
	ReportingState(const ReportingParams& p) : p_(p) { reset(); }
	
	/**
	 * Set all state to uninitialized defaults.
	 */
	void reset() {
		state_ = ReportingState::NO_READ;
		paired_ = false;
		nconcord_ = 0;
		ndiscord_ = 0;
		nunpair1_ = 0;
		nunpair2_ = 0;
		doneConcord_ = false;
		doneDiscord_ = false;
		doneUnpair_  = false;
		doneUnpair1_ = false;
		doneUnpair2_ = false;
		exitConcord_ = ReportingState::EXIT_DID_NOT_ENTER;
		exitDiscord_ = ReportingState::EXIT_DID_NOT_ENTER;
		exitUnpair1_ = ReportingState::EXIT_DID_NOT_ENTER;
		exitUnpair2_ = ReportingState::EXIT_DID_NOT_ENTER;
		done_ = false;
	}
	
	/**
	 * Return true iff this ReportingState has been initialized with a call to
	 * nextRead() since the last time reset() was called.
	 */
	bool inited() const { return state_ != ReportingState::NO_READ; }

	/**
	 * Initialize state machine with a new read.  The state we start in depends
	 * on whether it's paired-end or unpaired.
	 */
	void nextRead(bool paired);

	/**
	 * Caller uses this member function to indicate that one additional
	 * concordant alignment has been found.
	 */
	bool foundConcordant();

	/**
	 * Caller uses this member function to indicate that one additional
	 * discordant alignment has been found.
	 */
	bool foundUnpaired(bool mate1);
	
	/**
	 * Called to indicate that the aligner has finished searching for
	 * alignments.  This gives us a chance to finalize our state.
	 *
	 * TODO: Keep track of short-circuiting information.
	 */
	void finish();
	
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
	void getReport(
		uint64_t& nconcordAln, // # concordant alignments to report
		uint64_t& ndiscordAln, // # discordant alignments to report
		uint64_t& nunpair1Aln, // # unpaired alignments for mate #1 to report
		uint64_t& nunpair2Aln, // # unpaired alignments for mate #2 to report
		bool& pairMax,         // repetitive concordant alignments
		bool& unpair1Max,      // repetitive alignments for mate #1
		bool& unpair2Max)      // repetitive alignments for mate #2
		const;

	/**
	 * Return an integer representing the alignment state we're in.
	 */
	inline int state() const { return state_; }
	
	/**
	 * If false, there's no need to solve any more dynamic programming problems
	 * for finding opposite mates.
	 */
	inline bool doneConcordant() const { return doneConcord_; }
	
	/**
	 * If false, there's no need to seek any more discordant alignment.
	 */
	inline bool doneDiscordant() const { return doneDiscord_; }
	
	/**
	 * If false, there's no need to seek any more unpaired alignments for the
	 * specified mate.  Note: this doesn't necessarily mean we can stop looking
	 * for alignments for the mate, since this might be necessary for finding
	 * concordant and discordant alignments.
	 */
	inline bool doneUnpaired(bool mate1) const {
		return mate1 ? doneUnpair1_ : doneUnpair2_;
	}
	
	/**
	 * If false, no further consideration of the given mate is necessary.  It's
	 * not needed for *any* class of alignment: concordant, discordant or
	 * unpaired.
	 */
	inline bool doneWithMate(bool mate1) const {
		bool doneUnpair = mate1 ? doneUnpair1_ : doneUnpair2_;
		uint64_t nun = mate1 ? nunpair1_ : nunpair2_;
		if(!doneUnpair || !doneConcord_) {
			return false; // still needed for future concordant/unpaired alns
		}
		if(!doneDiscord_ && nun == 0) {
			return false; // still needed for future discordant alignments
		}
		return true; // done
	}

	/**
	 * Return true iff there's no need to seek any more unpaired alignments.
	 */
	inline bool doneUnpaired() const { return doneUnpair_; }
	
	/**
	 * Return true iff all alignment stages have been exited.
	 */
	inline bool done() const { return done_; }

	inline uint64_t numConcordant() const { return nconcord_; }
	inline uint64_t numDiscordant() const { return ndiscord_; }
	inline uint64_t numUnpaired1()  const { return nunpair1_; }
	inline uint64_t numUnpaired2()  const { return nunpair2_; }

	inline int exitConcordant() const { return exitConcord_; }
	inline int exitDiscordant() const { return exitDiscord_; }
	inline int exitUnpaired1()  const { return exitUnpair1_; }
	inline int exitUnpaired2()  const { return exitUnpair2_; }

#ifndef NDEBUG
	/**
	 * Check that ReportingState is internally consistent.
	 */
	bool repOk() const {
		assert(p_.discord || doneDiscord_);
		assert(p_.mixed   || !paired_ || doneUnpair_);
		assert(doneUnpair_ || !doneUnpair1_ || !doneUnpair2_);
		if(p_.mhitsSet()) {
			assert_leq(numConcordant(), (uint64_t)p_.mhits+1);
			assert_leq(numDiscordant(), (uint64_t)p_.mhits+1);
			assert(paired_ || numUnpaired1() <= (uint64_t)p_.mhits+1);
			assert(paired_ || numUnpaired2() <= (uint64_t)p_.mhits+1);
		}
		assert(done() || !doneWithMate(true) || !doneWithMate(false));
		return true;
	}
#endif

	/**
	 * Return ReportingParams object governing this ReportingState.
	 */
	const ReportingParams& params() const {
		return p_;
	}

protected:

	/**
	 * Update state to reflect situation after converting two unique unpaired
	 * alignments, one for mate 1 and one for mate 2, into a single discordant
	 * alignment.
	 */
	void convertUnpairedToDiscordant() {
		assert_eq(1, numUnpaired1());
		assert_eq(1, numUnpaired2());
		assert_eq(0, numDiscordant());
		exitUnpair1_ = exitUnpair2_ = ReportingState::EXIT_CONVERTED_TO_DISCORDANT;
		nunpair1_ = nunpair2_ = 0;
		ndiscord_ = 1;
		assert_eq(1, numDiscordant());
	}

	/**
	 * Given the number of alignments in a category, check whether we
	 * short-circuited out of the category.  Set the done and exit arguments to
	 * indicate whether and how we short-circuited.
	 */
	inline void areDone(
		uint64_t cnt,     // # alignments in category
		bool& done,       // out: whether we short-circuited out of category
		int& exit) const; // out: if done, how we short-circuited (-k? -m? etc)
	
	/**
	 * Update done_ field to reflect whether we're totally done now.
	 */
	inline void updateDone() {
		doneUnpair_ = doneUnpair1_ && doneUnpair2_;
		done_ = doneUnpair_ && doneDiscord_ && doneConcord_;
	}

	const ReportingParams& p_;  // reporting parameters
	int state_;          // state we're currently in
	bool paired_;        // true iff read we're currently handling is paired
	uint64_t nconcord_;  // # concordants found so far
	uint64_t ndiscord_;  // # discordants found so far
	uint64_t nunpair1_;  // # unpaired alignments found so far for mate 1
	uint64_t nunpair2_;  // # unpaired alignments found so far for mate 2
	bool doneConcord_;   // true iff we're no longner interested in concordants
	bool doneDiscord_;   // true iff we're no longner interested in discordants
	bool doneUnpair_;    // no longner interested in unpaired alns
	bool doneUnpair1_;   // no longner interested in unpaired alns for mate 1
	bool doneUnpair2_;   // no longner interested in unpaired alns for mate 2
	int exitConcord_;    // flag indicating how we exited concordant state
	int exitDiscord_;    // flag indicating how we exited discordant state
	int exitUnpair1_;    // flag indicating how we exited unpaired 1 state
	int exitUnpair2_;    // flag indicating how we exited unpaired 2 state
	bool done_;          // done with all alignments
};

/**
 * Global hit sink for hits from the MultiSeed aligner.  Encapsulates
 * all aspects of the MultiSeed aligner hitsink that are global to all
 * threads.  This includes aspects relating to:
 *
 * (a) synchronized access to the output stream
 * (b) the policy to be enforced by the per-thread wrapper
 *
 * TODO: Implement splitting up of alignments into separate files
 * according to genomic coordinate.
 */
template <typename index_t>
class AlnSink {

	typedef EList<std::string> StrList;

public:

	explicit AlnSink(
                     OutputQueue& oq,
                     const StrList& refnames,
                     bool quiet,
                     ALTDB<index_t>* altdb = NULL,
                     SpliceSiteDB* ssdb = NULL) :
    oq_(oq),
    refnames_(refnames),
    quiet_(quiet),
    altdb_(altdb),
    spliceSiteDB_(ssdb)
	{ }

	/**
	 * Destroy HitSinkobject;
	 */
	virtual ~AlnSink() { }

	/**
	 * Called when the AlnSink is wrapped by a new AlnSinkWrap.  This helps us
	 * keep track of whether the main lock or any of the per-stream locks will
	 * be contended by multiple threads.
	 */
	void addWrapper() { numWrappers_++; }

	/**
	 * Append a single hit to the given output stream.  If
	 * synchronization is required, append() assumes the caller has
	 * already grabbed the appropriate lock.
	 */
	virtual void append(
		BTString&             o,
		StackedAln&           staln,
		size_t                threadId,
		const Read           *rd1,
		const Read           *rd2,
		const TReadId         rdid,
		AlnRes               *rs1,
		AlnRes               *rs2,
		const AlnSetSumm&     summ,
		const SeedAlSumm&     ssm1,
		const SeedAlSumm&     ssm2,
		const AlnFlags*       flags1,
		const AlnFlags*       flags2,
		const PerReadMetrics& prm,
		const Mapq&           mapq,
		const Scoring&        sc,
		bool                  report2) = 0;

	/**
	 * Report a given batch of hits for the given read or read pair.
	 * Should be called just once per read pair.  Assumes all the
	 * alignments are paired, split between rs1 and rs2.
	 *
	 * The caller hasn't decided which alignments get reported as primary
	 * or secondary; that's up to the routine.  Because the caller might
	 * want to know this, we use the pri1 and pri2 out arguments to
	 * convey this.
	 */
	virtual void reportHits(
		BTString&             o,              // write to this buffer
		StackedAln&           staln,       // StackedAln to write stacked alignment
		size_t                threadId,       // which thread am I?
		const Read           *rd1,            // mate #1
		const Read           *rd2,            // mate #2
		const TReadId         rdid,           // read ID
		const EList<size_t>&  select1,        // random subset of rd1s
		const EList<size_t>*  select2,        // random subset of rd2s
		EList<AlnRes>        *rs1,            // alignments for mate #1
		EList<AlnRes>        *rs2,            // alignments for mate #2
		bool                  maxed,          // true iff -m/-M exceeded
		const AlnSetSumm&     summ,           // summary
		const SeedAlSumm&     ssm1,           // seed alignment summ
		const SeedAlSumm&     ssm2,           // seed alignment summ
		const AlnFlags*       flags1,         // flags for mate #1
		const AlnFlags*       flags2,         // flags for mate #2
		const PerReadMetrics& prm,            // per-read metrics
		const Mapq&           mapq,           // MAPQ generator
		const Scoring&        sc,             // scoring scheme
		bool                  getLock = true) // true iff lock held by caller
	{
		// There are a few scenarios:
		// 1. Read is unpaired, in which case rd2 is NULL
		// 2. Read is paired-end and we're reporting concordant alignments
		// 3. Read is paired-end and we're reporting discordant alignments
		// 4. Read is paired-end and we're reporting unpaired alignments for
		//    both mates
		// 5. Read is paired-end and we're reporting an unpaired alignments for
		//    just one mate or the other
		assert(rd1 != NULL || rd2 != NULL);
		assert(rs1 != NULL || rs2 != NULL);
		AlnFlags flagscp1, flagscp2;
		if(flags1 != NULL) {
			flagscp1 = *flags1;
			flags1 = &flagscp1;
			flagscp1.setPrimary(true);
		}
		if(flags2 != NULL) {
			flagscp2 = *flags2;
			flags2 = &flagscp2;
			flagscp2.setPrimary(true);
		}
		if(select2 != NULL) {
			// Handle case 5
			assert(rd1 != NULL); assert(flags1 != NULL);
			assert(rd2 != NULL); assert(flags2 != NULL);
			assert_gt(select1.size(), 0);
			assert_gt(select2->size(), 0);
			AlnRes* r1pri = ((rs1 != NULL) ? &rs1->get(select1[0]) : NULL);
			AlnRes* r2pri = ((rs2 != NULL) ? &rs2->get((*select2)[0]) : NULL);
			append(o, staln, threadId, rd1, rd2, rdid, r1pri, r2pri, summ,
			       ssm1, ssm2, flags1, flags2, prm, mapq, sc, true);
			flagscp1.setPrimary(false);
			flagscp2.setPrimary(false);
			for(size_t i = 1; i < select1.size(); i++) {
				AlnRes* r1 = ((rs1 != NULL) ? &rs1->get(select1[i]) : NULL);
				append(o, staln, threadId, rd1, rd2, rdid, r1, r2pri, summ,
				       ssm1, ssm2, flags1, flags2, prm, mapq, sc, false);
			}
			for(size_t i = 1; i < select2->size(); i++) {
				AlnRes* r2 = ((rs2 != NULL) ? &rs2->get((*select2)[i]) : NULL);
				append(o, staln, threadId, rd2, rd1, rdid, r2, r1pri, summ,
				       ssm2, ssm1, flags2, flags1, prm, mapq, sc, false);
			}
		} else {
			// Handle cases 1-4
			for(size_t i = 0; i < select1.size(); i++) {
				AlnRes* r1 = ((rs1 != NULL) ? &rs1->get(select1[i]) : NULL);
				AlnRes* r2 = ((rs2 != NULL) ? &rs2->get(select1[i]) : NULL);
				append(o, staln, threadId, rd1, rd2, rdid, r1, r2, summ,
				       ssm1, ssm2, flags1, flags2, prm, mapq, sc, true);
				if(flags1 != NULL) {
					flagscp1.setPrimary(false);
				}
				if(flags2 != NULL) {
					flagscp2.setPrimary(false);
				}
			}
		}
	}

	/**
	 * Report an unaligned read.  Typically we do nothing, but we might
	 * want to print a placeholder when output is chained.
	 */
	virtual void reportUnaligned(
		BTString&             o,              // write to this string
		StackedAln&           staln,          // StackedAln to write stacked alignment
		size_t                threadId,       // which thread am I?
		const Read           *rd1,            // mate #1
		const Read           *rd2,            // mate #2
		const TReadId         rdid,           // read ID
		const AlnSetSumm&     summ,           // summary
		const SeedAlSumm&     ssm1,           // seed alignment summary
		const SeedAlSumm&     ssm2,           // seed alignment summary
		const AlnFlags*       flags1,         // flags for mate #1
		const AlnFlags*       flags2,         // flags for mate #2
		const PerReadMetrics& prm,            // per-read metrics
		const Mapq&           mapq,           // MAPQ calculator
		const Scoring&        sc,             // scoring scheme
		bool                  report2,        // report alns for both mates?
		bool                  getLock = true) // true iff lock held by caller
	{
		append(o, staln, threadId, rd1, rd2, rdid, NULL, NULL, summ,
		       ssm1, ssm2, flags1, flags2, prm, mapq, sc, report2);
	}

	/**
	 * Print summary of how many reads aligned, failed to align and aligned
	 * repetitively.  Write it to stderr.  Optionally write Hadoop counter
	 * updates.
	 */
	void printAlSumm(
        ostream& out,
		const ReportingMetrics& met,
		size_t repThresh, // threshold for uniqueness, or max if no thresh
		bool discord,     // looked for discordant alignments
		bool mixed,       // looked for unpaired alignments where paired failed?
        bool newSummary,  // alignment summary in a new style
		bool hadoopOut);  // output Hadoop counters?

	/**
	 * Called when all alignments are complete.  It is assumed that no
	 * synchronization is necessary.
	 */
	void finish(
                ostream& out,
                size_t repThresh,
                bool discord,
                bool mixed,
                bool newSummary,
                bool hadoopOut)
	{
		// Close output streams
		if(!quiet_) {
			printAlSumm(
                        out,
                        met_,
                        repThresh,
                        discord,
                        mixed,
                        newSummary,
                        hadoopOut);
		}
	}

#ifndef NDEBUG
	/**
	 * Check that hit sink is internally consistent.
	 */
	bool repOk() const { return true; }
#endif
	
	//
	// Related to reporting seed hits
	//

	/**
	 * Given a Read and associated, filled-in SeedResults objects,
	 * print a record summarizing the seed hits.
	 */
	void reportSeedSummary(
		BTString&          o,
		const Read&        rd,
		TReadId            rdid,
		size_t             threadId,
		const SeedResults<index_t>& rs,
		bool               getLock = true);

	/**
	 * Given a Read, print an empty record (all 0s).
	 */
	void reportEmptySeedSummary(
		BTString&          o,
		const Read&        rd,
		TReadId            rdid,
		size_t             threadId,
		bool               getLock = true);

	/**
	 * Append a batch of unresolved seed alignment results (i.e. seed
	 * alignments where all we know is the reference sequence aligned
	 * to and its SA range, not where it falls in the reference
	 * sequence) to the given output stream in Bowtie's seed-alignment
	 * verbose-mode format.
	 */
	virtual void appendSeedSummary(
		BTString&     o,
		const Read&   rd,
		const TReadId rdid,
		size_t        seedsTried,
		size_t        nonzero,
		size_t        ranges,
		size_t        elts,
		size_t        seedsTriedFw,
		size_t        nonzeroFw,
		size_t        rangesFw,
		size_t        eltsFw,
		size_t        seedsTriedRc,
		size_t        nonzeroRc,
		size_t        rangesRc,
		size_t        eltsRc);

	/**
	 * Merge given metrics in with ours by summing all individual metrics.
	 */
	void mergeMetrics(const ReportingMetrics& met, bool getLock = true) {
		met_.merge(met, getLock);
	}

	/**
	 * Return mutable reference to the shared OutputQueue.
	 */
	OutputQueue& outq() {
		return oq_;
	}

protected:

	OutputQueue&       oq_;           // output queue
	int                numWrappers_;  // # threads owning a wrapper for this HitSink
	const StrList&     refnames_;     // reference names
	bool               quiet_;        // true -> don't print alignment stats at the end
	ReportingMetrics   met_;          // global repository of reporting metrics
    ALTDB<index_t>*    altdb_;
    SpliceSiteDB*      spliceSiteDB_; //
};

/**
 * Per-thread hit sink "wrapper" for the MultiSeed aligner.  Encapsulates
 * aspects of the MultiSeed aligner hit sink that are per-thread.  This
 * includes aspects relating to:
 *
 * (a) Enforcement of the reporting policy
 * (b) Tallying of results
 * (c) Storing of results for the previous read in case this allows us to
 *     short-circuit some work for the next read (i.e. if it's identical)
 *
 * PHASED ALIGNMENT ASSUMPTION
 *
 * We make some assumptions about how alignment proceeds when we try to
 * short-circuit work for identical reads.  Specifically, we assume that for
 * each read the aligner proceeds in a series of stages (or perhaps just one
 * stage).  In each stage, the aligner either:
 *
 * (a)  Finds no alignments, or
 * (b)  Finds some alignments and short circuits out of the stage with some
 *      random reporting involved (e.g. in -k and/or -M modes), or
 * (c)  Finds all of the alignments in the stage
 *
 * In the event of (a), the aligner proceeds to the next stage and keeps
 * trying; we can skip the stage entirely for the next read if it's identical.
 * In the event of (b), or (c), the aligner stops and does not proceed to
 * further stages.  In the event of (b1), if the next read is identical we
 * would like to tell the aligner to start again at the beginning of the stage
 * that was short-circuited.
 *
 * In any event, the rs1_/rs2_/rs1u_/rs2u_ fields contain the alignments found
 * in the last alignment stage attempted.
 *
 * HANDLING REPORTING LIMITS
 *
 * The user can specify reporting limits, like -k (specifies number of
 * alignments to report out of those found) and -M (specifies a ceiling s.t. if
 * there are more alignments than the ceiling, read is called repetitive and
 * best found is reported).  Enforcing these limits is straightforward for
 * unpaired alignments: if a new alignment causes us to exceed the -M ceiling,
 * we can stop looking.
 *
 * The case where both paired-end and unpaired alignments are possible is
 * trickier.  Once we have a number of unpaired alignments that exceeds the
 * ceiling, we can stop looking *for unpaired alignments* - but we can't
 * necessarily stop looking for paired-end alignments, since there may yet be
 * more to find.  However, if the input read is not a pair, then we can stop at
 * this point.  If the input read is a pair and we have a number of paired
 * aligments that exceeds the -M ceiling, we can stop looking.
 *
 * CONCORDANT & DISCORDANT, PAIRED & UNPAIRED
 *
 * A note on paired-end alignment: Clearly, if an input read is
 * paired-end and we find either concordant or discordant paired-end
 * alignments for the read, then we would like to tally and report
 * those alignments as such (and not as groups of 2 unpaired
 * alignments).  And if we fail to find any paired-end alignments, but
 * we do find some unpaired alignments for one mate or the other, then
 * we should clearly tally and report those alignments as unpaired
 * alignments (if the user so desires).
 *
 * The situation is murkier when there are no paired-end alignments,
 * but there are unpaired alignments for *both* mates.  In this case,
 * we might want to pick out zero or more pairs of mates and classify
 * those pairs as discordant paired-end alignments.  And we might want
 * to classify the remaining alignments as unpaired.  But how do we
 * pick which pairs if any to call discordant?
 *
 * Because the most obvious use for discordant pairs is for identifying
 * large-scale variation, like rearrangements or large indels, we would
 * usually like to be conservative about what we call a discordant
 * alignment.  If there's a good chance that one or the other of the
 * two mates has a good alignment to another place on the genome, this
 * compromises the evidence for the large-scale variant.  For this
 * reason, Bowtie 2's policy is: if there are no paired-end alignments
 * and there is *exactly one alignment each* for both mates, then the
 * two alignments are paired and treated as a discordant paired-end
 * alignment.  Otherwise, all alignments are treated as unpaired
 * alignments.
 *
 * When both paired and unpaired alignments are discovered by the
 * aligner, only the paired alignments are reported by default.  This
 * is sensible considering relative likelihoods: if a good paired-end
 * alignment is found, it is much more likely that the placement of
 * the two mates implied by that paired alignment is correct than any
 * placement implied by an unpaired alignment.
 *
 * 
 */
template <typename index_t>
class AlnSinkWrap {
public:

	AlnSinkWrap(
                AlnSink<index_t>& g,       // AlnSink being wrapped
                const ReportingParams& rp, // Parameters governing reporting
                Mapq& mapq,                // Mapq calculator
                size_t threadId,           // Thread ID
                bool secondary = false,    // Secondary alignments
                const SpliceSiteDB* ssdb = NULL, // splice sites
                uint64_t threads_rids_mindist = 0) : // synchronization
		g_(g),
		rp_(rp),
        threadid_(threadId),
        mapq_(mapq),
    	secondary_(secondary),
        ssdb_(ssdb),
        threads_rids_mindist_(threads_rids_mindist),
		init_(false),   
		maxed1_(false),       // read is pair and we maxed out mate 1 unp alns
		maxed2_(false),       // read is pair and we maxed out mate 2 unp alns
		maxedOverall_(false), // alignments found so far exceed -m/-M ceiling
		bestPair_(std::numeric_limits<TAlScore>::min()),
		best2Pair_(std::numeric_limits<TAlScore>::min()),
		bestUnp1_(std::numeric_limits<TAlScore>::min()),
		best2Unp1_(std::numeric_limits<TAlScore>::min()),
		bestUnp2_(std::numeric_limits<TAlScore>::min()),
		best2Unp2_(std::numeric_limits<TAlScore>::min()),
        bestSplicedPair_(0),
        best2SplicedPair_(0),
        bestSplicedUnp1_(0),
        best2SplicedUnp1_(0),
        bestSplicedUnp2_(0),
        best2SplicedUnp2_(0),
		rd1_(NULL),    // mate 1
		rd2_(NULL),    // mate 2
		rdid_(std::numeric_limits<TReadId>::max()), // read id
		rs1_(),        // mate 1 alignments for paired-end alignments
		rs2_(),        // mate 2 alignments for paired-end alignments
		rs1u_(),       // mate 1 unpaired alignments
		rs2u_(),       // mate 2 unpaired alignments
		select1_(),    // for selecting random subsets for mate 1
		select2_(),    // for selecting random subsets for mate 2
		st_(rp)        // reporting state - what's left to do?
	{
		assert(rp_.repOk());
	}

	/**
	 * Initialize the wrapper with a new read pair and return an
	 * integer >= -1 indicating which stage the aligner should start
	 * at.  If -1 is returned, the aligner can skip the read entirely.
	 * at.  If .  Checks if the new read pair is identical to the
	 * previous pair.  If it is, then we return the id of the first
	 * stage to run.
	 */
	int nextRead(
		// One of the other of rd1, rd2 will = NULL if read is unpaired
		const Read* rd1,      // new mate #1
		const Read* rd2,      // new mate #2
		TReadId rdid,         // read ID for new pair
		bool qualitiesMatter);// aln policy distinguishes b/t quals?

	/**
	 * Inform global, shared AlnSink object that we're finished with
	 * this read.  The global AlnSink is responsible for updating
	 * counters, creating the output record, and delivering the record
	 * to the appropriate output stream.
	 */
	void finishRead(
		const SeedResults<index_t> *sr1, // seed alignment results for mate 1
		const SeedResults<index_t> *sr2, // seed alignment results for mate 2
		bool               exhaust1,     // mate 1 exhausted?
		bool               exhaust2,     // mate 2 exhausted?
		bool               nfilt1,       // mate 1 N-filtered?
		bool               nfilt2,       // mate 2 N-filtered?
		bool               scfilt1,      // mate 1 score-filtered?
		bool               scfilt2,      // mate 2 score-filtered?
		bool               lenfilt1,     // mate 1 length-filtered?
		bool               lenfilt2,     // mate 2 length-filtered?
		bool               qcfilt1,      // mate 1 qc-filtered?
		bool               qcfilt2,      // mate 2 qc-filtered?
		bool               sortByScore,  // prioritize alignments by score
		RandomSource&      rnd,          // pseudo-random generator
		ReportingMetrics&  met,          // reporting metrics
		const PerReadMetrics& prm,       // per-read metrics
		const Scoring& sc,               // scoring scheme
		bool suppressSeedSummary = true,
		bool suppressAlignments = false,
        bool templateLenAdjustment = true);
	
	/**
	 * Called by the aligner when a new unpaired or paired alignment is
	 * discovered in the given stage.  This function checks whether the
	 * addition of this alignment causes the reporting policy to be
	 * violated (by meeting or exceeding the limits set by -k, -m, -M),
	 * in which case true is returned immediately and the aligner is
	 * short circuited.  Otherwise, the alignment is tallied and false
	 * is returned.
	 */
	bool report(
		int stage,
		const AlnRes* rs1,
		const AlnRes* rs2);

#ifndef NDEBUG
	/**
	 * Check that hit sink wrapper is internally consistent.
	 */
	bool repOk() const {
		assert_eq(rs2_.size(), rs1_.size());
		if(rp_.mhitsSet()) {
			assert_gt(rp_.mhits, 0);
			assert_leq((int)rs1_.size(), rp_.mhits+1);
			assert_leq((int)rs2_.size(), rp_.mhits+1);
			assert(readIsPair() || (int)rs1u_.size() <= rp_.mhits+1);
			assert(readIsPair() || (int)rs2u_.size() <= rp_.mhits+1);
		}
		if(init_) {
			assert(rd1_ != NULL);
			assert_neq(std::numeric_limits<TReadId>::max(), rdid_);
		}
		assert_eq(st_.numConcordant() + st_.numDiscordant(), rs1_.size());
		//assert_eq(st_.numUnpaired1(), rs1u_.size());
		//assert_eq(st_.numUnpaired2(), rs2u_.size());
		assert(st_.repOk());
		return true;
	}
#endif
	
	/**
	 * Return true iff no alignments have been reported to this wrapper
	 * since the last call to nextRead().
	 */
	bool empty() const {
		return rs1_.empty() && rs1u_.empty() && rs2u_.empty();
	}
	
	/**
	 * Return true iff we have already encountered a number of alignments that
	 * exceeds the -m/-M ceiling.  TODO: how does this distinguish between
	 * pairs and mates?
	 */
	bool maxed() const {
		return maxedOverall_;
	}
	
	/**
	 * Return true if the current read is paired.
	 */
	bool readIsPair() const {
		return rd1_ != NULL && rd2_ != NULL;
	}
	
	/**
	 * Return true iff nextRead() has been called since the last time
	 * finishRead() was called.
	 */
	bool inited() const { return init_; }

	/**
	 * Return a const ref to the ReportingState object associated with the
	 * AlnSinkWrap.
	 */
	const ReportingState& state() const { return st_; }
    
    const ReportingParams& reportingParams() { return rp_;}
	
	/**
	 * Return true iff we're in -M mode.
	 */
	bool Mmode() const {
		return rp_.mhitsSet();
	}
	
	/**
	 * Return true iff the policy is to report all hits.
	 */
	bool allHits() const {
		return rp_.allHits();
	}
	
	/**
	 * Return true iff at least two alignments have been reported so far for an
	 * unpaired read or mate 1.
	 */
	bool hasSecondBestUnp1() const {
		return best2Unp1_ != std::numeric_limits<TAlScore>::min();
	}

	/**
	 * Return true iff at least two alignments have been reported so far for
	 * mate 2.
	 */
	bool hasSecondBestUnp2() const {
		return best2Unp2_ != std::numeric_limits<TAlScore>::min();
	}

	/**
	 * Return true iff at least two paired-end alignments have been reported so
	 * far.
	 */
	bool hasSecondBestPair() const {
		return best2Pair_ != std::numeric_limits<TAlScore>::min();
	}
	
	/**
	 * Get best score observed so far for an unpaired read or mate 1.
	 */
	TAlScore bestUnp1() const {
		return bestUnp1_;
	}

	/**
	 * Get second-best score observed so far for an unpaired read or mate 1.
	 */
	TAlScore secondBestUnp1() const {
		return best2Unp1_;
	}

	/**
	 * Get best score observed so far for mate 2.
	 */
	TAlScore bestUnp2() const {
		return bestUnp2_;
	}

	/**
	 * Get second-best score observed so far for mate 2.
	 */
	TAlScore secondBestUnp2() const {
		return best2Unp2_;
	}

	/**
	 * Get best score observed so far for paired-end read.
	 */
	TAlScore bestPair() const {
		return bestPair_;
	}

	/**
	 * Get second-best score observed so far for paired-end read.
	 */
	TAlScore secondBestPair() const {
		return best2Pair_;
	}
    
    index_t bestSplicedPair() const {
        return bestSplicedPair_;
    }
    
    index_t best2SplicedPair() const {
        return best2SplicedPair_;
    }
    
    index_t bestSplicedUnp1() const {
        return bestSplicedUnp1_;
    }
    
    index_t best2SplicedUnp1() const {
        return best2SplicedUnp1_;
    }
    
    index_t bestSplicedUnp2() const {
        return bestSplicedUnp2_;
    }
    
    index_t best2SplicedUnp2() const {
        return best2SplicedUnp2_;
    }
    
    bool secondary() const {
        return secondary_;
    }
    
    /**
     *
     */
    void getUnp1(const EList<AlnRes>*& rs) const { rs = &rs1u_; }
    void getUnp2(const EList<AlnRes>*& rs) const { rs = &rs2u_; }
    void getPair(const EList<AlnRes>*& rs1, const EList<AlnRes>*& rs2) const { rs1 = &rs1_; rs2 = &rs2_; }

protected:

	/**
	 * Return true iff the read in rd1/rd2 matches the last read handled, which
	 * should still be in rd1_/rd2_.
	 */
	bool sameRead(
		const Read* rd1,
		const Read* rd2,
		bool qualitiesMatter);

	/**
	 * If there is a configuration of unpaired alignments that fits our
	 * criteria for there being one or more discordant alignments, then
	 * shift the discordant alignments over to the rs1_/rs2_ lists, clear the
	 * rs1u_/rs2u_ lists and return true.  Otherwise, return false.
	 */
	bool prepareDiscordants();

	/**
	 * Given that rs is already populated with alignments, consider the
	 * alignment policy and make random selections where necessary.  E.g. if we
	 * found 10 alignments and the policy is -k 2 -m 20, select 2 alignments at
	 * random.  We "select" an alignment by setting the parallel entry in the
	 * 'select' list to true.
	 */
	size_t selectAlnsToReport(
		const EList<AlnRes>& rs,     // alignments to select from
		uint64_t             num,    // number of alignments to select
		EList<size_t>&       select, // list to put results in
		RandomSource&        rnd)
		const;

	/**
	 * rs1 (possibly together with rs2 if reads are paired) are populated with
	 * alignments.  Here we prioritize them according to alignment score, and
	 * some randomness to break ties.  Priorities are returned in the 'select'
	 * list.
	 */
	size_t selectByScore(
		const EList<AlnRes>* rs1,    // alignments to select from (mate 1)
		const EList<AlnRes>* rs2,    // alignments to select from (mate 2, or NULL)
		uint64_t             num,    // number of alignments to select
		EList<size_t>&       select, // prioritized list to put results in
		RandomSource&        rnd)
		const;

	AlnSink<index_t>& g_;     // global alignment sink
	ReportingParams   rp_;    // reporting parameters: khits, mhits etc
	size_t            threadid_; // thread ID
	Mapq&             mapq_;  // mapq calculator
    bool              secondary_; // allow for secondary alignments
    const SpliceSiteDB* ssdb_; // splice sites
    uint64_t threads_rids_mindist_; // synchronization
	bool              init_;  // whether we're initialized w/ read pair
	bool              maxed1_; // true iff # unpaired mate-1 alns reported so far exceeded -m/-M
	bool              maxed2_; // true iff # unpaired mate-2 alns reported so far exceeded -m/-M
	bool              maxedOverall_; // true iff # paired-end alns reported so far exceeded -m/-M
	TAlScore          bestPair_;     // greatest score so far for paired-end
	TAlScore          best2Pair_;    // second-greatest score so far for paired-end
	TAlScore          bestUnp1_;     // greatest score so far for unpaired/mate1
	TAlScore          best2Unp1_;    // second-greatest score so far for unpaired/mate1
	TAlScore          bestUnp2_;     // greatest score so far for mate 2
	TAlScore          best2Unp2_;    // second-greatest score so far for mate 2
    index_t           bestSplicedPair_;
    index_t           best2SplicedPair_;
    index_t           bestSplicedUnp1_;
    index_t           best2SplicedUnp1_;
    index_t           bestSplicedUnp2_;
    index_t           best2SplicedUnp2_;
	const Read*       rd1_;   // mate #1
	const Read*       rd2_;   // mate #2
	TReadId           rdid_;  // read ID (potentially used for ordering)
	EList<AlnRes>     rs1_;   // paired alignments for mate #1
	EList<AlnRes>     rs2_;   // paired alignments for mate #2
	EList<AlnRes>     rs1u_;  // unpaired alignments for mate #1
	EList<AlnRes>     rs2u_;  // unpaired alignments for mate #2
	EList<size_t>     select1_; // parallel to rs1_/rs2_ - which to report
	EList<size_t>     select2_; // parallel to rs1_/rs2_ - which to report
	ReportingState    st_;      // reporting state - what's left to do?
	
	EList<std::pair<TAlScore, size_t> > selectBuf_;
	BTString obuf_;
	StackedAln staln_;
    
    EList<SpliceSite> spliceSites_;
};

/**
 * An AlnSink concrete subclass for printing SAM alignments.  The user might
 * want to customize SAM output in various ways.  We encapsulate all these
 * customizations, and some of the key printing routines, in the SamConfig
 * class in sam.h/sam.cpp.
 */
template <typename index_t>
class AlnSinkSam : public AlnSink<index_t> {

	typedef EList<std::string> StrList;

public:

	AlnSinkSam(
               OutputQueue&     oq,            // output queue
               const SamConfig<index_t>& samc, // settings & routines for SAM output
               const StrList&   refnames,      // reference names
               bool             quiet,         // don't print alignment summary at end
               ALTDB<index_t>*  altdb = NULL,
               SpliceSiteDB*    ssdb  = NULL) :
		AlnSink<index_t>(
                         oq,
                         refnames,
                         quiet,
                         altdb,
                         ssdb),
    samc_(samc)
	{ }
	
	virtual ~AlnSinkSam() { }

	/**
	 * Append a single alignment result, which might be paired or
	 * unpaired, to the given output stream in Bowtie's verbose-mode
	 * format.  If the alignment is paired-end, print mate1's alignment
	 * then mate2's alignment.
	 */
	virtual void append(
		BTString&     o,           // write output to this string
		StackedAln&   staln,       // StackedAln to write stacked alignment
		size_t        threadId,    // which thread am I?
		const Read*   rd1,         // mate #1
		const Read*   rd2,         // mate #2
		const TReadId rdid,        // read ID
		AlnRes* rs1,               // alignments for mate #1
		AlnRes* rs2,               // alignments for mate #2
		const AlnSetSumm& summ,    // summary
		const SeedAlSumm& ssm1,    // seed alignment summary
		const SeedAlSumm& ssm2,    // seed alignment summary
		const AlnFlags* flags1,    // flags for mate #1
		const AlnFlags* flags2,    // flags for mate #2
		const PerReadMetrics& prm, // per-read metrics
		const Mapq& mapq,          // MAPQ calculator
		const Scoring& sc,         // scoring scheme
		bool report2)              // report alns for both mates
	{
		assert(rd1 != NULL || rd2 != NULL);
		if(rd1 != NULL) {
			assert(flags1 != NULL);
			appendMate(o, staln, *rd1, rd2, rdid, rs1, rs2, summ, ssm1, ssm2,
			           *flags1, prm, mapq, sc);
            if(rs1 != NULL && rs1->spliced() && this->spliceSiteDB_ != NULL) {
                this->spliceSiteDB_->addSpliceSite(*rd1, *rs1);
            }
		}
		if(rd2 != NULL && report2) {
			assert(flags2 != NULL);
			appendMate(o, staln, *rd2, rd1, rdid, rs2, rs1, summ, ssm2, ssm1,
			           *flags2, prm, mapq, sc);
            if(rs2 != NULL && rs2->spliced() && this->spliceSiteDB_ != NULL) {
                this->spliceSiteDB_->addSpliceSite(*rd2, *rs2);
            }
		}
	}

protected:

	/**
	 * Append a single per-mate alignment result to the given output
	 * stream.  If the alignment is part of a pair, information about
	 * the opposite mate and its alignment are given in rdo/rso.
	 */
	void appendMate(
		BTString&     o,
		StackedAln&   staln,
		const Read&   rd,
		const Read*   rdo,
		const TReadId rdid,
		AlnRes* rs,
		AlnRes* rso,
		const AlnSetSumm& summ,
		const SeedAlSumm& ssm,
		const SeedAlSumm& ssmo,
		const AlnFlags& flags,
		const PerReadMetrics& prm, // per-read metrics
		const Mapq& mapq,          // MAPQ calculator
		const Scoring& sc);        // scoring scheme

	const SamConfig<index_t>& samc_;    // settings & routines for SAM output
	BTDnaString               dseq_;    // buffer for decoded read sequence
	BTString                  dqual_;   // buffer for decoded quality sequence
};

static inline std::ostream& printPct(
							  std::ostream& os,
							  uint64_t num,
							  uint64_t denom)
{
	double pct = 0.0f;
	if(denom != 0) { pct = 100.0 * (double)num / (double)denom; }
	os << fixed << setprecision(2) << pct << '%';
	return os;
}

/**
 * Print a friendly summary of:
 *
 *  1. How many reads were aligned and had one or more alignments
 *     reported
 *  2. How many reads exceeded the -m or -M ceiling and therefore had
 *     their alignments suppressed or sampled
 *  3. How many reads failed to align entirely
 *
 * Optionally print a series of Hadoop streaming-style counter updates
 * with similar information.
 */
template <typename index_t>
void AlnSink<index_t>::printAlSumm(
                                   ostream& out,
								   const ReportingMetrics& met,
								   size_t repThresh,   // threshold for uniqueness, or max if no thresh
								   bool discord,       // looked for discordant alignments
								   bool mixed,         // looked for unpaired alignments where paired failed?
                                   bool newSummary,    // alignment summary in a new style
								   bool hadoopOut)     // output Hadoop counters?
{
	// NOTE: there's a filtering step at the very beginning, so everything
	// being reported here is post filtering
	
	bool canRep = repThresh != MAX_SIZE_T;
    if(hadoopOut) {
        out << "reporter:counter:HISAT2,Reads processed," << met.nread << endl;
    }
    uint64_t totread = met.nread;
    uint64_t totpair = met.npaired;
    uint64_t totunpair = met.nunpaired;
    uint64_t tot_al_cand = totunpair + totpair*2;
    uint64_t tot_al = (met.nconcord_uni + met.nconcord_rep) * 2 + (met.ndiscord) * 2 + met.nunp_0_uni + met.nunp_0_rep + met.nunp_uni + met.nunp_rep;
    assert_leq(tot_al, tot_al_cand);
    if(newSummary) {
        out << "HISAT2 summary stats:" << endl;
        if(totpair > 0) {
            uint64_t ncondiscord_0 = met.nconcord_0 - met.ndiscord;
            out << "\tTotal pairs: " << totpair << endl;
            out << "\t\tAligned concordantly or discordantly 0 time: " << ncondiscord_0 << " ("; printPct(out, ncondiscord_0, met.npaired); out << ")" << endl;
            out << "\t\tAligned concordantly 1 time: " << met.nconcord_uni1 << " ("; printPct(out, met.nconcord_uni1, met.npaired); out << ")" << endl;
            out << "\t\tAligned concordantly >1 times: " << met.nconcord_uni2 << " ("; printPct(out, met.nconcord_uni2, met.npaired); out << ")" << endl;
            out << "\t\tAligned discordantly 1 time: " << met.ndiscord << " ("; printPct(out, met.ndiscord, met.npaired); out << ")" << endl;
            
            out << "\tTotal unpaired reads: " << ncondiscord_0 * 2 << endl;
            out << "\t\tAligned 0 time: " << met.nunp_0_0 << " ("; printPct(out, met.nunp_0_0, ncondiscord_0 * 2); out << ")" << endl;
            out << "\t\tAligned 1 time: " << met.nunp_0_uni1 << " ("; printPct(out, met.nunp_0_uni1, ncondiscord_0 * 2); out << ")" << endl;
            out << "\t\tAligned >1 times: " << met.nunp_0_uni2 << " ("; printPct(out, met.nunp_0_uni2, ncondiscord_0 * 2); out << ")" << endl;
        } else {
            out << "\tTotal reads: " << totread << endl;
            out << "\t\tAligned 0 time: " << met.nunp_0 << " ("; printPct(out, met.nunp_0, met.nunpaired); out << ")" << endl;
            out << "\t\tAligned 1 time: " << met.nunp_uni1 << " ("; printPct(out, met.nunp_uni1, met.nunpaired); out << ")" << endl;
            out << "\t\tAligned >1 times: " << met.nunp_uni2 << " ("; printPct(out, met.nunp_uni2, met.nunpaired); out << ")" << endl;
        }
        out << "\tOverall alignment rate: "; printPct(out, tot_al, tot_al_cand); out << endl;
        
    } else {
        if(totread > 0) {
            out << "" << totread << " reads; of these:" << endl;
        } else {
            assert_eq(0, met.npaired);
            assert_eq(0, met.nunpaired);
            out << "" << totread << " reads" << endl;
        }
        if(totpair > 0) {
            // Paired output
            out << "  " << totpair << " (";
            printPct(out, totpair, totread);
            out << ") were paired; of these:" << endl;
            
            // Concordants
            out << "    " << met.nconcord_0 << " (";
            printPct(out, met.nconcord_0, met.npaired);
            out << ") aligned concordantly 0 times" << endl;
            if(canRep) {
                // Print the number that aligned concordantly exactly once
                assert_eq(met.nconcord_uni, met.nconcord_uni1+met.nconcord_uni2);
                out << "    " << met.nconcord_uni1 << " (";
                printPct(out, met.nconcord_uni1, met.npaired);
                out << ") aligned concordantly exactly 1 time" << endl;
                
                // Print the number that aligned concordantly more than once but
                // fewer times than the limit
                
                out << "    " << met.nconcord_uni2+met.nconcord_rep << " (";
                printPct(out, met.nconcord_uni2+met.nconcord_rep, met.npaired);
                out << ") aligned concordantly >1 times" << endl;
            } else {
                // Print the number that aligned concordantly exactly once
                assert_eq(met.nconcord_uni, met.nconcord_uni1+met.nconcord_uni2);
                out << "    " << met.nconcord_uni1 << " (";
                printPct(out, met.nconcord_uni1, met.npaired);
                out << ") aligned concordantly exactly 1 time" << endl;
                
                // Print the number that aligned concordantly more than once
                out << "    " << met.nconcord_uni2 << " (";
                printPct(out, met.nconcord_uni2, met.npaired);
                out << ") aligned concordantly >1 times" << endl;
            }
            if(discord) {
                // TODO: what about discoardant and on separate chromosomes?
                
                // Bring out the unaligned pair total so we can subtract discordants
                out << "    ----" << endl;
                out << "    " << met.nconcord_0
                << " pairs aligned concordantly 0 times; of these:" << endl;
                // Discordants
                out << "      " << met.ndiscord << " (";
                printPct(out, met.ndiscord, met.nconcord_0);
                out << ") aligned discordantly 1 time" << endl;
            }
            uint64_t ncondiscord_0 = met.nconcord_0 - met.ndiscord;
            if(mixed) {
                // Bring out the unaligned pair total so we can subtract discordants
                out << "    ----" << endl;
                out << "    " << ncondiscord_0
                    << " pairs aligned 0 times concordantly or discordantly; of these:" << endl;
                out << "      " << (ncondiscord_0 * 2) << " mates make up the pairs; of these:" << endl;
                out << "        " << met.nunp_0_0 << " " << "(";
                printPct(out, met.nunp_0_0, ncondiscord_0 * 2);
                out << ") aligned 0 times" << endl;
                if(canRep) {
                    // Print the number that aligned exactly once
                    assert_eq(met.nunp_0_uni, met.nunp_0_uni1+met.nunp_0_uni2);
                    out << "        " << met.nunp_0_uni1 << " (";
                    printPct(out, met.nunp_0_uni1, ncondiscord_0 * 2);
                    out << ") aligned exactly 1 time" << endl;
                    
                    // Print the number that aligned more than once but fewer times
                    // than the limit
                    out << "        " << met.nunp_0_uni2+met.nunp_0_rep << " (";
                    printPct(out, met.nunp_0_uni2+met.nunp_0_rep, ncondiscord_0 * 2);
                    out << ") aligned >1 times" << endl;
                } else {
                    // Print the number that aligned exactly once
                    assert_eq(met.nunp_0_uni, met.nunp_0_uni1+met.nunp_0_uni2);
                    out << "        " << met.nunp_0_uni1 << " (";
                    printPct(out, met.nunp_0_uni1, ncondiscord_0 * 2);
                    out << ") aligned exactly 1 time" << endl;
                    
                    // Print the number that aligned more than once but fewer times
                    // than the limit
                    out << "        " << met.nunp_0_uni2 << " (";
                    printPct(out, met.nunp_0_uni2, ncondiscord_0 * 2);
                    out << ") aligned >1 times" << endl;
                }
            }
        }
        if(totunpair > 0) {
            // Unpaired output
            out << "  " << totunpair << " (";
            printPct(out, totunpair, totread);
            out << ") were unpaired; of these:" << endl;
            
            out << "    " << met.nunp_0 << " (";
            printPct(out, met.nunp_0, met.nunpaired);
            out << ") aligned 0 times" << endl;
            if(hadoopOut) {
                out << "reporter:counter:HISAT 2,Unpaired reads with 0 alignments,"
                << met.nunpaired << endl;
            }
            
            if(canRep) {
                // Print the number that aligned exactly once
                assert_eq(met.nunp_uni, met.nunp_uni1+met.nunp_uni2);
                out << "    " << met.nunp_uni1 << " (";
                printPct(out, met.nunp_uni1, met.nunpaired);
                out << ") aligned exactly 1 time" << endl;
                
                // Print the number that aligned more than once but fewer times
                // than the limit
                out << "    " << met.nunp_uni2+met.nunp_rep << " (";
                printPct(out, met.nunp_uni2+met.nunp_rep, met.nunpaired);
                out << ") aligned >1 times" << endl;
            } else {
                // Print the number that aligned exactly once
                assert_eq(met.nunp_uni, met.nunp_uni1+met.nunp_uni2);
                out << "    " << met.nunp_uni1 << " (";
                printPct(out, met.nunp_uni1, met.nunpaired);
                out << ") aligned exactly 1 time" << endl;
                
                // Print the number that aligned more than once
                out << "    " << met.nunp_uni2 << " (";
                printPct(out, met.nunp_uni2, met.nunpaired);
                out << ") aligned >1 times" << endl;
            }
        }
        
        printPct(out, tot_al, tot_al_cand);
        out << " overall alignment rate" << endl;
    }
}

/**
 * Return true iff the read in rd1/rd2 matches the last read handled, which
 * should still be in rd1_/rd2_.
 */
template <typename index_t>
bool AlnSinkWrap<index_t>::sameRead(
									// One of the other of rd1, rd2 will = NULL if read is unpaired
									const Read* rd1,      // new mate #1
									const Read* rd2,      // new mate #2
									bool qualitiesMatter) // aln policy distinguishes b/t quals?
{
	bool same = false;
	if(rd1_ != NULL || rd2_ != NULL) {
		// This is not the first time the sink was initialized with
		// a read.  Check if new read/pair is identical to previous
		// read/pair
		if((rd1_ == NULL) == (rd1 == NULL) &&
		   (rd2_ == NULL) == (rd2 == NULL))
		{
			bool m1same = (rd1 == NULL && rd1_ == NULL);
			if(!m1same) {
				assert(rd1 != NULL);
				assert(rd1_ != NULL);
				m1same = Read::same(
									rd1->patFw,  // new seq
									rd1->qual,   // new quals
									rd1_->patFw, // old seq
									rd1_->qual,  // old quals
									qualitiesMatter);
			}
			if(m1same) {
				bool m2same = (rd2 == NULL && rd2_ == NULL);
				if(!m2same) {
					m2same = Read::same(
										rd2->patFw,  // new seq
										rd2->qual,   // new quals
										rd2_->patFw, // old seq
										rd2_->qual,  // old quals
										qualitiesMatter);
				}
				same = m2same;
			}
		}
	}
	return same;
}

/**
 * Initialize the wrapper with a new read pair and return an integer >= -1
 * indicating which stage the aligner should start at.  If -1 is returned, the
 * aligner can skip the read entirely.  Checks if the new read pair is
 * identical to the previous pair.  If it is, then we return the id of the
 * first stage to run.
 */
template <typename index_t>
int AlnSinkWrap<index_t>::nextRead(
								   // One of the other of rd1, rd2 will = NULL if read is unpaired
								   const Read* rd1,      // new mate #1
								   const Read* rd2,      // new mate #2
								   TReadId rdid,         // read ID for new pair
								   bool qualitiesMatter) // aln policy distinguishes b/t quals?
{
	assert(!init_);
	assert(rd1 != NULL || rd2 != NULL);
	init_ = true;
	// Keep copy of new read, so that we can compare it with the
	// next one
	if(rd1 != NULL) {
		rd1_ = rd1;
	} else rd1_ = NULL;
	if(rd2 != NULL) {
		rd2_ = rd2;
	} else rd2_ = NULL;
	rdid_ = rdid;
	// Caller must now align the read
	maxed1_ = false;
	maxed2_ = false;
	maxedOverall_ = false;
	bestPair_ = best2Pair_ =
	bestUnp1_ = best2Unp1_ =
	bestUnp2_ = best2Unp2_ = std::numeric_limits<THitInt>::min();
    bestSplicedPair_ = best2SplicedPair_ =
    bestSplicedUnp1_ = best2SplicedUnp1_ =
    bestSplicedUnp2_ = best2SplicedUnp2_ = 0;
	rs1_.clear();     // clear out paired-end alignments
	rs2_.clear();     // clear out paired-end alignments
	rs1u_.clear();    // clear out unpaired alignments for mate #1
	rs2u_.clear();    // clear out unpaired alignments for mate #2
	st_.nextRead(readIsPair()); // reset state
	assert(empty());
	assert(!maxed());
	// Start from the first stage
	return 0;
}

/**
 * Inform global, shared AlnSink object that we're finished with this read.
 * The global AlnSink is responsible for updating counters, creating the output
 * record, and delivering the record to the appropriate output stream.
 *
 * What gets reported for a paired-end alignment?
 *
 * 1. If there are reportable concordant alignments, report those and stop
 * 2. If there are reportable discordant alignments, report those and stop
 * 3. If unpaired alignments can be reported:
 *    3a. Report 
 #
 * Update metrics.  Only ambiguity is: what if a pair aligns repetitively and
 * one of its mates aligns uniquely?
 *
 * 	uint64_t al;   // # mates w/ >= 1 reported alignment
 *  uint64_t unal; // # mates w/ 0 alignments
 *  uint64_t max;  // # mates withheld for exceeding -M/-m ceiling
 *  uint64_t al_concord;  // # pairs w/ >= 1 concordant alignment
 *  uint64_t al_discord;  // # pairs w/ >= 1 discordant alignment
 *  uint64_t max_concord; // # pairs maxed out
 *  uint64_t unal_pair;   // # pairs where neither mate aligned
 */
template <typename index_t>
void AlnSinkWrap<index_t>::finishRead(
									  const SeedResults<index_t> *sr1, // seed alignment results for mate 1
									  const SeedResults<index_t> *sr2, // seed alignment results for mate 2
									  bool               exhaust1,     // mate 1 exhausted?
									  bool               exhaust2,     // mate 2 exhausted?
									  bool               nfilt1,       // mate 1 N-filtered?
									  bool               nfilt2,       // mate 2 N-filtered?
									  bool               scfilt1,      // mate 1 score-filtered?
									  bool               scfilt2,      // mate 2 score-filtered?
									  bool               lenfilt1,     // mate 1 length-filtered?
									  bool               lenfilt2,     // mate 2 length-filtered?
									  bool               qcfilt1,      // mate 1 qc-filtered?
									  bool               qcfilt2,      // mate 2 qc-filtered?
									  bool               sortByScore,  // prioritize alignments by score
									  RandomSource&      rnd,          // pseudo-random generator
									  ReportingMetrics&  met,          // reporting metrics
									  const PerReadMetrics& prm,       // per-read metrics
									  const Scoring& sc,               // scoring scheme
									  bool suppressSeedSummary,        // = true
                                      bool suppressAlignments,         // = false
                                      bool templateLenAdjustment)      // = true
{
	obuf_.clear();
	OutputQueueMark qqm(g_.outq(), obuf_, rdid_, threadid_);
	assert(init_);
	if(!suppressSeedSummary) {
		if(sr1 != NULL) {
			assert(rd1_ != NULL);
			// Mate exists and has non-empty SeedResults
			g_.reportSeedSummary(obuf_, *rd1_, rdid_, threadid_, *sr1, true);
		} else if(rd1_ != NULL) {
			// Mate exists but has NULL SeedResults
			g_.reportEmptySeedSummary(obuf_, *rd1_, rdid_, true);
		}
		if(sr2 != NULL) {
			assert(rd2_ != NULL);
			// Mate exists and has non-empty SeedResults
			g_.reportSeedSummary(obuf_, *rd2_, rdid_, threadid_, *sr2, true);
		} else if(rd2_ != NULL) {
			// Mate exists but has NULL SeedResults
			g_.reportEmptySeedSummary(obuf_, *rd2_, rdid_, true);
		}
	}
	if(!suppressAlignments) {
		// Ask the ReportingState what to report
		st_.finish();
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		st_.getReport(
					  nconcord,
					  ndiscord,
					  nunpair1,
					  nunpair2,
					  pairMax,
					  unpair1Max,
					  unpair2Max);
		assert_leq(nconcord, rs1_.size());
		assert_leq(nunpair1, rs1u_.size());
		assert_leq(nunpair2, rs2u_.size());
		assert_leq(ndiscord, 1);
		assert_gt(rp_.khits, 0);
		assert_gt(rp_.mhits, 0);
		assert(!pairMax    || rs1_.size()  >= (uint64_t)rp_.mhits);
		assert(!unpair1Max || rs1u_.size() >= (uint64_t)rp_.mhits);
		assert(!unpair2Max || rs2u_.size() >= (uint64_t)rp_.mhits);
		met.nread++;
		if(readIsPair()) {
			met.npaired++;
		} else {
			met.nunpaired++;
		}
		// Report concordant paired-end alignments if possible
		if(nconcord > 0) {
            AlnSetSumm concordSumm(
                                   rd1_, rd2_, &rs1_, &rs2_, &rs1u_, &rs2u_,
                                   exhaust1, exhaust2, -1, -1);
            
			// Possibly select a random subset
			size_t off;
			if(sortByScore) {
				// Sort by score then pick from low to high
				off = selectByScore(&rs1_, &rs2_, nconcord, select1_, rnd);
			} else {
				// Select subset randomly
				off = selectAlnsToReport(rs1_, nconcord, select1_, rnd);
			}
            
            concordSumm.numAlnsPaired(select1_.size());
            
			assert_lt(off, rs1_.size());
			const AlnRes *rs1 = &rs1_[off];
			const AlnRes *rs2 = &rs2_[off];
			AlnFlags flags1(
							ALN_FLAG_PAIR_CONCORD_MATE1,
							st_.params().mhitsSet(),
							unpair1Max,
							pairMax,
							nfilt1,
							scfilt1,
							lenfilt1,
							qcfilt1,
							st_.params().mixed,
							true,       // primary
							true,       // opp aligned
							rs2->fw()); // opp fw
			AlnFlags flags2(
							ALN_FLAG_PAIR_CONCORD_MATE2,
							st_.params().mhitsSet(),
							unpair2Max,
							pairMax,
							nfilt2,
							scfilt2,
							lenfilt2,
							qcfilt2,
							st_.params().mixed,
							false,      // primary
							true,       // opp aligned
							rs1->fw()); // opp fw
			// Issue: we only set the flags once, but some of the flags might
			// vary from pair to pair among the pairs we're reporting.  For
			// instance, whether a given mate aligns to the forward strand.
			SeedAlSumm ssm1, ssm2;
            if(sr1 != NULL && sr2 != NULL) {
                sr1->toSeedAlSumm(ssm1);
                sr2->toSeedAlSumm(ssm2);
            }
			for(size_t i = 0; i < rs1_.size(); i++) {
                spliceSites_.clear();
                if(templateLenAdjustment) {
                    rs1_[i].setMateParams(ALN_RES_TYPE_MATE1, &rs2_[i], flags1, ssdb_, threads_rids_mindist_, &spliceSites_);
                    rs2_[i].setMateParams(ALN_RES_TYPE_MATE2, &rs1_[i], flags2, ssdb_, threads_rids_mindist_, &spliceSites_);
                } else {
                    rs1_[i].setMateParams(ALN_RES_TYPE_MATE1, &rs2_[i], flags1);
                    rs2_[i].setMateParams(ALN_RES_TYPE_MATE2, &rs1_[i], flags2);
                }
				assert_eq(abs(rs1_[i].fragmentLength()), abs(rs2_[i].fragmentLength()));
			}
			assert(!select1_.empty());
			g_.reportHits(
						  obuf_,
						  staln_,
						  threadid_,
						  rd1_,
						  rd2_,
						  rdid_,
						  select1_,
						  NULL,
						  &rs1_,
						  &rs2_,
						  pairMax,
						  concordSumm,
						  ssm1,
						  ssm2,
						  &flags1,
						  &flags2,
						  prm,
						  mapq_,
						  sc);
			if(pairMax) {
				met.nconcord_rep++;
			} else {
				met.nconcord_uni++;
				assert(!rs1_.empty());
				if(select1_.size() == 1) {
					met.nconcord_uni1++;
				} else {
					met.nconcord_uni2++;
				}
			}
			init_ = false;
			//g_.outq().finishRead(obuf_, rdid_, threadid_);
			return;
		}
		// Report concordant paired-end alignments if possible
		else if(ndiscord > 0) {
			ASSERT_ONLY(bool ret =) prepareDiscordants();
			assert(ret);
			assert_eq(1, rs1_.size());
			assert_eq(1, rs2_.size());
			AlnSetSumm discordSumm(
								   rd1_, rd2_, &rs1_, &rs2_, &rs1u_, &rs2u_,
								   exhaust1, exhaust2, -1, -1);
			const AlnRes *rs1 = &rs1_[0];
			const AlnRes *rs2 = &rs2_[0];
			AlnFlags flags1(
							ALN_FLAG_PAIR_DISCORD_MATE1,
							st_.params().mhitsSet(),
							false,
							pairMax,
							nfilt1,
							scfilt1,
							lenfilt1,
							qcfilt1,
							st_.params().mixed,
							true,       // primary
							true,       // opp aligned
							rs2->fw()); // opp fw
			AlnFlags flags2(
							ALN_FLAG_PAIR_DISCORD_MATE2,
							st_.params().mhitsSet(),
							false,
							pairMax,
							nfilt2,
							scfilt2,
							lenfilt2,
							qcfilt2,
							st_.params().mixed,
							false,      // primary
							true,       // opp aligned
							rs1->fw()); // opp fw
			SeedAlSumm ssm1, ssm2;
            if(sr1 != NULL) sr1->toSeedAlSumm(ssm1);
			if(sr2 != NULL) sr2->toSeedAlSumm(ssm2);
			for(size_t i = 0; i < rs1_.size(); i++) {
				rs1_[i].setMateParams(ALN_RES_TYPE_MATE1, &rs2_[i], flags1);
				rs2_[i].setMateParams(ALN_RES_TYPE_MATE2, &rs1_[i], flags2);
				assert(rs1_[i].isFraglenSet() == rs2_[i].isFraglenSet());
				assert(!rs1_[i].isFraglenSet() || abs(rs1_[i].fragmentLength()) == abs(rs2_[i].fragmentLength()));
			}
			ASSERT_ONLY(size_t off);
			if(sortByScore) {
				// Sort by score then pick from low to high
				ASSERT_ONLY(off =) selectByScore(&rs1_, &rs2_, ndiscord, select1_, rnd);
			} else {
				// Select subset randomly
				ASSERT_ONLY(off =) selectAlnsToReport(rs1_, ndiscord, select1_, rnd);
			}
			assert_eq(0, off);
			assert(!select1_.empty());
			g_.reportHits(
						  obuf_,
						  staln_,
						  threadid_,
						  rd1_,
						  rd2_,
						  rdid_,
						  select1_,
						  NULL,
						  &rs1_,
						  &rs2_,
						  pairMax,
						  discordSumm,
						  ssm1,
						  ssm2,
						  &flags1,
						  &flags2,
						  prm,
						  mapq_,
						  sc);
			met.nconcord_0++;
			met.ndiscord++;
			init_ = false;
			//g_.outq().finishRead(obuf_, rdid_, threadid_);
			return;
		}
		// If we're at this point, at least one mate failed to align.
		// BTL: That's not true.  It could be that there are no concordant
		// alignments but both mates have unpaired alignments, with one of
		// the mates having more than one.
		//assert(nunpair1 == 0 || nunpair2 == 0);
		assert(!pairMax);
			
		const AlnRes *repRs1 = NULL, *repRs2 = NULL;
		AlnSetSumm summ1, summ2;
		AlnFlags flags1, flags2;
		TRefId refid = -1; TRefOff refoff = -1;
		bool rep1 = rd1_ != NULL && nunpair1 > 0;
		bool rep2 = rd2_ != NULL && nunpair2 > 0;
		
		// This is the preliminary if statement for mate 1 - here we're
		// gathering some preliminary information, making it possible to call
		// g_.reportHits(...) with information about both mates potentially
		if(rep1) {
			// Mate 1 aligned at least once
            if(rep2) {
                summ1.init(
					   rd1_, rd2_, NULL, NULL, &rs1u_, &rs2u_,
					   exhaust1, exhaust2, -1, -1);
            } else {
                summ1.init(
                           rd1_, NULL, NULL, NULL, &rs1u_, NULL,
                           exhaust1, exhaust2, -1, -1);
            }
			size_t off;
			if(sortByScore) {
				// Sort by score then pick from low to high
				off = selectByScore(&rs1u_, NULL, nunpair1, select1_, rnd);
			} else {
				// Select subset randomly
				off = selectAlnsToReport(rs1u_, nunpair1, select1_, rnd);
			}
            summ1.numAlns1(select1_.size());
            summ2.numAlns1(select1_.size());
			repRs1 = &rs1u_[off];
		} else if(rd1_ != NULL) {
			// Mate 1 failed to align - don't do anything yet.  First we want
			// to collect information on mate 2 in case that factors into the
			// summary
			assert(!unpair1Max);
		}
		
		if(rep2) {
            if(rep1) {
                summ2.init(
                           rd1_, rd2_, NULL, NULL, &rs1u_, &rs2u_,
                           exhaust1, exhaust2, -1, -1);
            } else {
                summ2.init(
                           NULL, rd2_, NULL, NULL, NULL, &rs2u_,
                           exhaust1, exhaust2, -1, -1);
            }
			size_t off;
			if(sortByScore) {
				// Sort by score then pick from low to high
				off = selectByScore(&rs2u_, NULL, nunpair2, select2_, rnd);
			} else {
				// Select subset randomly
				off = selectAlnsToReport(rs2u_, nunpair2, select2_, rnd);
			}
			repRs2 = &rs2u_[off];
            summ1.numAlns2(select2_.size());
            summ2.numAlns2(select2_.size());
		} else if(rd2_ != NULL) {
			// Mate 2 failed to align - don't do anything yet.  First we want
			// to collect information on mate 1 in case that factors into the
			// summary
			assert(!unpair2Max);
		}
        
        // Update counters given that one mate didn't align
        if(readIsPair()) {
            met.nconcord_0++;
        }
        if(rd1_ != NULL) {
            if(nunpair1 > 0) {
                // Update counters
                if(readIsPair()) {
                    if(unpair1Max) met.nunp_0_rep++;
                    else {
                        met.nunp_0_uni++;
                        assert(!rs1u_.empty());
                        if(select1_.size() == 1) {
                            met.nunp_0_uni1++;
                        } else {
                            met.nunp_0_uni2++;
                        }
                    }
                } else {
                    if(unpair1Max) met.nunp_rep++;
                    else {
                        met.nunp_uni++;
                        assert(!rs1u_.empty());
                        if(select1_.size() == 1) {
                            met.nunp_uni1++;
                        } else {
                            met.nunp_uni2++;
                        }
                    }
                }
            } else if(unpair1Max) {
                // Update counters
                if(readIsPair())   met.nunp_0_rep++;
                else               met.nunp_rep++;
            } else {
                // Update counters
                if(readIsPair())   met.nunp_0_0++;
                else               met.nunp_0++;
            }
        }
        if(rd2_ != NULL) {
            if(nunpair2 > 0) {
                // Update counters
                if(readIsPair()) {
                    if(unpair2Max) met.nunp_0_rep++;
                    else {
                        assert(!rs2u_.empty());
                        met.nunp_0_uni++;
                        if(select2_.size() == 1) {
                            met.nunp_0_uni1++;
                        } else {
                            met.nunp_0_uni2++;
                        }
                    }
                } else {
                    if(unpair2Max) met.nunp_rep++;
                    else {
                        assert(!rs2u_.empty());
                        met.nunp_uni++;
                        if(select2_.size() == 1) {
                            met.nunp_uni1++;
                        } else {
                            met.nunp_uni2++;
                        }
                    }
                }
            } else if(unpair2Max) {
                // Update counters
                if(readIsPair())   met.nunp_0_rep++;
                else               met.nunp_rep++;
            } else {
                // Update counters
                if(readIsPair())   met.nunp_0_0++;
                else               met.nunp_0++;
            }
        }
		
		// Now set up flags
		if(rep1) {
			// Initialize flags.  Note: We want to have information about how
			// the other mate aligned (if it did) at this point
			flags1.init(
						readIsPair() ?
						ALN_FLAG_PAIR_UNPAIRED_MATE1 :
						ALN_FLAG_PAIR_UNPAIRED,
						st_.params().mhitsSet(),
						unpair1Max,
						pairMax,
						nfilt1,
						scfilt1,
						lenfilt1,
						qcfilt1,
						st_.params().mixed,
						true,   // primary
						repRs2 != NULL,                    // opp aligned
						repRs2 == NULL || repRs2->fw());   // opp fw
			for(size_t i = 0; i < rs1u_.size(); i++) {
				rs1u_[i].setMateParams(ALN_RES_TYPE_UNPAIRED_MATE1, NULL, flags1);
			}
		}
		if(rep2) {
			// Initialize flags.  Note: We want to have information about how
			// the other mate aligned (if it did) at this point
			flags2.init(
						readIsPair() ?
						ALN_FLAG_PAIR_UNPAIRED_MATE2 :
						ALN_FLAG_PAIR_UNPAIRED,
						st_.params().mhitsSet(),
						unpair2Max,
						pairMax,
						nfilt2,
						scfilt2,
						lenfilt2,
						qcfilt2,
						st_.params().mixed,
						true,   // primary
						repRs1 != NULL,                  // opp aligned
						repRs1 == NULL || repRs1->fw()); // opp fw
			for(size_t i = 0; i < rs2u_.size(); i++) {
				rs2u_[i].setMateParams(ALN_RES_TYPE_UNPAIRED_MATE2, NULL, flags2);
			}
		}
		
		// Now report mate 1
		if(rep1) {
			SeedAlSumm ssm1, ssm2;
			if(sr1 != NULL) sr1->toSeedAlSumm(ssm1);
			if(sr2 != NULL) sr2->toSeedAlSumm(ssm2);
			assert(!select1_.empty());
			g_.reportHits(
						  obuf_,
						  staln_,
						  threadid_,
						  rd1_,
						  repRs2 != NULL ? rd2_ : NULL,
						  rdid_,
						  select1_,
						  repRs2 != NULL ? &select2_ : NULL,
						  &rs1u_,
						  repRs2 != NULL ? &rs2u_ : NULL,
						  unpair1Max,
						  summ1,
						  ssm1,
						  ssm2,
						  &flags1,
						  repRs2 != NULL ? &flags2 : NULL,
						  prm,
						  mapq_,
						  sc);
			assert_lt(select1_[0], rs1u_.size());
			refid = rs1u_[select1_[0]].refid();
			refoff = rs1u_[select1_[0]].refoff();
		}
		
		// Now report mate 2
		if(rep2 && !rep1) {
			SeedAlSumm ssm1, ssm2;
			if(sr1 != NULL) sr1->toSeedAlSumm(ssm1);
			if(sr2 != NULL) sr2->toSeedAlSumm(ssm2);
			assert(!select2_.empty());
			g_.reportHits(
						  obuf_,
						  staln_,
						  threadid_,
						  rd2_,
						  repRs1 != NULL ? rd1_ : NULL,
						  rdid_,
						  select2_,
						  repRs1 != NULL ? &select1_ : NULL,
						  &rs2u_,
						  repRs1 != NULL ? &rs1u_ : NULL,
						  unpair2Max,
						  summ2,
						  ssm1,
						  ssm2,
						  &flags2,
						  repRs1 != NULL ? &flags1 : NULL,
						  prm,
						  mapq_,
						  sc);
			assert_lt(select2_[0], rs2u_.size());
			refid = rs2u_[select2_[0]].refid();
			refoff = rs2u_[select2_[0]].refoff();
		}
		
		if(rd1_ != NULL && nunpair1 == 0) {
			if(nunpair2 > 0) {
				assert_neq(-1, refid);
				summ1.init(
						   rd1_, NULL, NULL, NULL, NULL, NULL,
						   exhaust1, exhaust2, refid, refoff);
			} else {
				summ1.init(
						   rd1_, NULL, NULL, NULL, NULL, NULL,
						   exhaust1, exhaust2, -1, -1);
			}
			SeedAlSumm ssm1, ssm2;
			if(sr1 != NULL) sr1->toSeedAlSumm(ssm1);
			if(sr2 != NULL) sr2->toSeedAlSumm(ssm2);
			flags1.init(
						readIsPair() ?
						ALN_FLAG_PAIR_UNPAIRED_MATE1 :
						ALN_FLAG_PAIR_UNPAIRED,
						st_.params().mhitsSet(),
						false,
						false,
						nfilt1,
						scfilt1,
						lenfilt1,
						qcfilt1,
						st_.params().mixed,
						true,           // primary
						repRs2 != NULL, // opp aligned
						(repRs2 != NULL) ? repRs2->fw() : false); // opp fw
			g_.reportUnaligned(
							   obuf_,      // string to write output to
							   staln_,
							   threadid_,
							   rd1_,    // read 1
							   NULL,    // read 2
							   rdid_,   // read id
							   summ1,   // summ
							   ssm1,    // 
							   ssm2,
							   &flags1, // flags 1
							   NULL,    // flags 2
							   prm,     // per-read metrics
							   mapq_,   // MAPQ calculator
							   sc,      // scoring scheme
							   true);   // get lock?
		}
		if(rd2_ != NULL && nunpair2 == 0) {
			if(nunpair1 > 0) {
				assert_neq(-1, refid);
				summ2.init(
						   NULL, rd2_, NULL, NULL, NULL, NULL,
						   exhaust1, exhaust2, refid, refoff);
			} else {
				summ2.init(
						   NULL, rd2_, NULL, NULL, NULL, NULL,
						   exhaust1, exhaust2, -1, -1);
			}
			SeedAlSumm ssm1, ssm2;
			if(sr1 != NULL) sr1->toSeedAlSumm(ssm1);
			if(sr2 != NULL) sr2->toSeedAlSumm(ssm2);
			flags2.init(
						readIsPair() ?
						ALN_FLAG_PAIR_UNPAIRED_MATE2 :
						ALN_FLAG_PAIR_UNPAIRED,
						st_.params().mhitsSet(),
						false,
						false,
						nfilt2,
						scfilt2,
						lenfilt2,
						qcfilt2,
						st_.params().mixed,
						true,           // primary
						repRs1 != NULL, // opp aligned
						(repRs1 != NULL) ? repRs1->fw() : false); // opp fw
			g_.reportUnaligned(
							   obuf_,      // string to write output to
							   staln_,
							   threadid_,
							   rd2_,    // read 1
							   NULL,    // read 2
							   rdid_,   // read id
							   summ2,   // summ
							   ssm1,
							   ssm2,
							   &flags2, // flags 1
							   NULL,    // flags 2
							   prm,     // per-read metrics
							   mapq_,   // MAPQ calculator
							   sc,      // scoring scheme
							   true);   // get lock?
		}
	} // if(suppress alignments)
	init_ = false;
	return;
}

/**
 * Called by the aligner when a new unpaired or paired alignment is
 * discovered in the given stage.  This function checks whether the
 * addition of this alignment causes the reporting policy to be
 * violated (by meeting or exceeding the limits set by -k, -m, -M),
 * in which case true is returned immediately and the aligner is
 * short circuited.  Otherwise, the alignment is tallied and false
 * is returned.
 */
template <typename index_t>
bool AlnSinkWrap<index_t>::report(
								  int stage,
								  const AlnRes* rs1,
								  const AlnRes* rs2)
{
	assert(init_);
	assert(rs1 != NULL || rs2 != NULL);
	assert(rs1 == NULL || !rs1->empty());
	assert(rs2 == NULL || !rs2->empty());
	assert(rs1 == NULL || rs1->repOk());
	assert(rs2 == NULL || rs2->repOk());
	bool paired = (rs1 != NULL && rs2 != NULL);
	bool one = (rs1 != NULL);
	const AlnRes* rsa = one ? rs1 : rs2;
	const AlnRes* rsb = one ? rs2 : rs1;
	if(paired) {
		assert(readIsPair());
		st_.foundConcordant();
		rs1_.push_back(*rs1);
		rs2_.push_back(*rs2);
	} else {
        st_.foundUnpaired(one);
		if(one) {
			rs1u_.push_back(*rs1);
  		} else {
			rs2u_.push_back(*rs2);
		}
	}
	// Tally overall alignment score
	TAlScore score = rsa->score().score();
	if(rsb != NULL) score += rsb->score().score();
    index_t num_spliced = (index_t)rsa->num_spliced();
    if(rsb != NULL) num_spliced += (index_t)rsb->num_spliced();
	// Update best score so far
	if(paired) {
		if(score > bestPair_) {
			best2Pair_ = bestPair_;
			bestPair_ = score;
            best2SplicedPair_ = bestSplicedPair_;
            bestSplicedPair_ = num_spliced;
		} else if(score > best2Pair_) {
			best2Pair_ = score;
            best2SplicedPair_ = num_spliced;
		}
	} else {
		if(one) {
			if(score > bestUnp1_) {
				best2Unp1_ = bestUnp1_;
				bestUnp1_ = score;
                best2SplicedUnp1_ = bestSplicedUnp1_;
                bestSplicedUnp1_ = num_spliced;
			} else if(score > best2Unp1_) {
				best2Unp1_ = score;
                best2SplicedUnp1_ = num_spliced;
			}
		} else {
			if(score > bestUnp2_) {
				best2Unp2_ = bestUnp2_;
				bestUnp2_ = score;
                best2SplicedUnp2_ = bestSplicedUnp2_;
                bestSplicedUnp2_ = num_spliced;
			} else if(score > best2Unp2_) {
				best2Unp2_ = score;
                best2SplicedUnp1_ = num_spliced;
			}
		}
	}
	return st_.done();
}

/**
 * If there is a configuration of unpaired alignments that fits our
 * criteria for there being one or more discordant alignments, then
 * shift the discordant alignments over to the rs1_/rs2_ lists, clear the
 * rs1u_/rs2u_ lists and return true.  Otherwise, return false.
 */
template <typename index_t>
bool AlnSinkWrap<index_t>::prepareDiscordants() {
	if(rs1u_.size() == 1 && rs2u_.size() == 1) {
		assert(rs1_.empty());
		assert(rs2_.empty());
		rs1_.push_back(rs1u_[0]);
		rs2_.push_back(rs2u_[0]);
		return true;
	}
	return false;
}

/**
 * rs1 (possibly together with rs2 if reads are paired) are populated with
 * alignments.  Here we prioritize them according to alignment score, and
 * some randomness to break ties.  Priorities are returned in the 'select'
 * list.
 */
template <typename index_t>
size_t AlnSinkWrap<index_t>::selectByScore(
										   const EList<AlnRes>* rs1,    // alignments to select from (mate 1)
										   const EList<AlnRes>* rs2,    // alignments to select from (mate 2, or NULL)
										   uint64_t             num,    // number of alignments to select
										   EList<size_t>&       select, // prioritized list to put results in
										   RandomSource&        rnd)
const
{
	assert(init_);
	assert(repOk());
	assert_gt(num, 0);
	assert(rs1 != NULL);
	size_t sz = rs1->size(); // sz = # alignments found
	assert_leq(num, sz);
	if(sz < num) {
		num = sz;
	}
	// num = # to select
	if(sz < 1) {
		return 0;
	}
	select.resize((size_t)num);
	// Use 'selectBuf_' as a temporary list for sorting purposes
	EList<std::pair<TAlScore, size_t> >& buf =
	const_cast<EList<std::pair<TAlScore, size_t> >& >(selectBuf_);
	buf.resize(sz);
	// Sort by score.  If reads are pairs, sort by sum of mate scores.
	for(size_t i = 0; i < sz; i++) {
		buf[i].first = (*rs1)[i].score().hisat2_score();
		if(rs2 != NULL) {
			buf[i].first += (*rs2)[i].score().hisat2_score();
		}
		buf[i].second = i; // original offset
	}
	buf.sort(); buf.reverse(); // sort in descending order by score
	
	// Randomize streaks of alignments that are equal by score
	size_t streak = 0;
	for(size_t i = 1; i < buf.size(); i++) {
		if(buf[i].first == buf[i-1].first) {
			if(streak == 0) { streak = 1; }
			streak++;
		} else {
			if(streak > 1) {
				assert_geq(i, streak);
				buf.shufflePortion(i-streak, streak, rnd);
			}
			streak = 0;
		}
	}
	if(streak > 1) {
		buf.shufflePortion(buf.size() - streak, streak, rnd);
	}
	
	for(size_t i = 0; i < num; i++) { select[i] = buf[i].second; }
    
    if(!secondary_) {
        assert_geq(buf.size(), select.size());
        for(size_t i = 0; i + 1 < select.size(); i++) {
            if(buf[i].first != buf[i+1].first) {
                select.resize(i+1);
                break;
            }
        }
    }
    
	// Returns index of the representative alignment, but in 'select' also
	// returns the indexes of the next best selected alignments in order by
	// score.
	return selectBuf_[0].second;
}

/**
 * Given that rs is already populated with alignments, consider the
 * alignment policy and make random selections where necessary.  E.g. if we
 * found 10 alignments and the policy is -k 2 -m 20, select 2 alignments at
 * random.  We "select" an alignment by setting the parallel entry in the
 * 'select' list to true.
 *
 * Return the "representative" alignment.  This is simply the first one
 * selected.  That will also be what SAM calls the "primary" alignment.
 */
template <typename index_t>
size_t AlnSinkWrap<index_t>::selectAlnsToReport(
												const EList<AlnRes>& rs,     // alignments to select from
												uint64_t             num,    // number of alignments to select
												EList<size_t>&       select, // list to put results in
												RandomSource&        rnd)
const
{
	assert(init_);
	assert(repOk());
	assert_gt(num, 0);
	size_t sz = rs.size();
	if(sz < num) {
		num = sz;
	}
	if(sz < 1) {
		return 0;
	}
	select.resize((size_t)num);
	if(sz == 1) {
		assert_eq(1, num);
		select[0] = 0;
		return 0;
	}
	// Select a random offset into the list of alignments
	uint32_t off = rnd.nextU32() % (uint32_t)sz;
	uint32_t offOrig = off;
	// Now take elements starting at that offset, wrapping around to 0 if
	// necessary.  Leave the rest.
	for(size_t i = 0; i < num; i++) {
		select[i] = off;
		off++;
		if(off == sz) {
			off = 0;
		}
	}
	return offOrig;
}

#define NOT_SUPPRESSED !suppress_[field++]
#define BEGIN_FIELD { \
if(firstfield) firstfield = false; \
else o.append('\t'); \
}
#define WRITE_TAB { \
if(firstfield) firstfield = false; \
else o.append('\t'); \
}
#define WRITE_NUM(o, x) { \
itoa10(x, buf); \
o.append(buf); \
}

/**
 * Print a seed summary to the first output stream in the outs_ list.
 */
template <typename index_t>
void AlnSink<index_t>::reportSeedSummary(
										 BTString&          o,
										 const Read&        rd,
										 TReadId            rdid,
										 size_t             threadId,
										 const SeedResults<index_t>& rs,
										 bool               getLock)
{
	appendSeedSummary(
					  o,                     // string to write to
					  rd,                    // read
					  rdid,                  // read id
					  rs.numOffs()*2,        // # seeds tried
					  rs.nonzeroOffsets(),   // # seeds with non-empty results
					  rs.numRanges(),        // # ranges for all seed hits
					  rs.numElts(),          // # elements for all seed hits
					  rs.numOffs(),          // # seeds tried from fw read
					  rs.nonzeroOffsetsFw(), // # seeds with non-empty results from fw read
					  rs.numRangesFw(),      // # ranges for seed hits from fw read
					  rs.numEltsFw(),        // # elements for seed hits from fw read
					  rs.numOffs(),          // # seeds tried from rc read
					  rs.nonzeroOffsetsRc(), // # seeds with non-empty results from fw read
					  rs.numRangesRc(),      // # ranges for seed hits from fw read
					  rs.numEltsRc());       // # elements for seed hits from fw read
}

/**
 * Print an empty seed summary to the first output stream in the outs_ list.
 */
template <typename index_t>
void AlnSink<index_t>::reportEmptySeedSummary(
											  BTString&          o,
											  const Read&        rd,
											  TReadId            rdid,
											  size_t             threadId,
											  bool               getLock)
{
	appendSeedSummary(
					  o,                     // string to append to
					  rd,                    // read
					  rdid,                  // read id
					  0,                     // # seeds tried
					  0,                     // # seeds with non-empty results
					  0,                     // # ranges for all seed hits
					  0,                     // # elements for all seed hits
					  0,                     // # seeds tried from fw read
					  0,                     // # seeds with non-empty results from fw read
					  0,                     // # ranges for seed hits from fw read
					  0,                     // # elements for seed hits from fw read
					  0,                     // # seeds tried from rc read
					  0,                     // # seeds with non-empty results from fw read
					  0,                     // # ranges for seed hits from fw read
					  0);                    // # elements for seed hits from fw read
}

/**
 * Print the given string.  If ws = true, print only up to and not
 * including the first space or tab.  Useful for printing reference
 * names.
 */
template<typename T>
static inline void printUptoWs(
							   BTString& s,
							   const T& str,
							   bool chopws)
{
	size_t len = str.length();
	for(size_t i = 0; i < len; i++) {
		if(!chopws || (str[i] != ' ' && str[i] != '\t')) {
			s.append(str[i]);
		} else {
			break;
		}
	}
}

/**
 * Append a batch of unresolved seed alignment summary results (i.e.
 * seed alignments where all we know is the reference sequence aligned
 * to and its SA range, not where it falls in the reference
 * sequence) to the given output stream in Bowtie's seed-sumamry
 * verbose-mode format.
 *
 * The seed summary format is:
 *
 *  - One line per read
 *  - A typical line consists of a set of tab-delimited fields:
 *
 *    1. Read name
 *    2. Total number of seeds extracted from the read
 *    3. Total number of seeds that aligned to the reference at
 *       least once (always <= field 2)
 *    4. Total number of distinct BW ranges found in all seed hits
 *       (always >= field 3)
 *    5. Total number of distinct BW elements found in all seed
 *       hits (always >= field 4)
 *    6-9.:   Like 2-5. but just for seeds extracted from the
 *            forward representation of the read
 *    10-13.: Like 2-5. but just for seeds extracted from the
 *            reverse-complement representation of the read
 *
 *    Note that fields 6 and 10 should add to field 2, 7 and 11
 *    should add to 3, etc.
 *
 *  - Lines for reads that are filtered out for any reason (e.g. too
 *    many Ns) have columns 2 through 13 set to 0.
 */
template <typename index_t>
void AlnSink<index_t>::appendSeedSummary(
										 BTString&     o,
										 const Read&   rd,
										 const TReadId rdid,
										 size_t        seedsTried,
										 size_t        nonzero,
										 size_t        ranges,
										 size_t        elts,
										 size_t        seedsTriedFw,
										 size_t        nonzeroFw,
										 size_t        rangesFw,
										 size_t        eltsFw,
										 size_t        seedsTriedRc,
										 size_t        nonzeroRc,
										 size_t        rangesRc,
										 size_t        eltsRc)
{
	char buf[1024];
	bool firstfield = true;
	//
	// Read name
	//
	BEGIN_FIELD;
	printUptoWs(o, rd.name, true);
	
	//
	// Total number of seeds tried
	//
	BEGIN_FIELD;
	WRITE_NUM(o, seedsTried);
	
	//
	// Total number of seeds tried where at least one range was found.
	//
	BEGIN_FIELD;
	WRITE_NUM(o, nonzero);
	
	//
	// Total number of ranges found
	//
	BEGIN_FIELD;
	WRITE_NUM(o, ranges);
	
	//
	// Total number of elements found
	//
	BEGIN_FIELD;
	WRITE_NUM(o, elts);
	
	//
	// The same four numbers, but only for seeds extracted from the
	// forward read representation.
	//
	BEGIN_FIELD;
	WRITE_NUM(o, seedsTriedFw);
	
	BEGIN_FIELD;
	WRITE_NUM(o, nonzeroFw);
	
	BEGIN_FIELD;
	WRITE_NUM(o, rangesFw);
	
	BEGIN_FIELD;
	WRITE_NUM(o, eltsFw);
	
	//
	// The same four numbers, but only for seeds extracted from the
	// reverse complement read representation.
	//
	BEGIN_FIELD;
	WRITE_NUM(o, seedsTriedRc);
	
	BEGIN_FIELD;
	WRITE_NUM(o, nonzeroRc);
	
	BEGIN_FIELD;
	WRITE_NUM(o, rangesRc);
	
	BEGIN_FIELD;
	WRITE_NUM(o, eltsRc);
	
	o.append('\n');
}

/**
 * Append a single hit to the given output stream in Bowtie's
 * verbose-mode format.
 */
template <typename index_t>
void AlnSinkSam<index_t>::appendMate(
									 BTString&     o,           // append to this string
									 StackedAln&   staln,       // store stacked alignment struct here
									 const Read&   rd,
									 const Read*   rdo,
									 const TReadId rdid,
									 AlnRes* rs,
									 AlnRes* rso,
									 const AlnSetSumm& summ,
									 const SeedAlSumm& ssm,
									 const SeedAlSumm& ssmo,
									 const AlnFlags& flags,
									 const PerReadMetrics& prm,
									 const Mapq& mapqCalc,
									 const Scoring& sc)
{
	if(rs == NULL && samc_.omitUnalignedReads()) {
		return;
	}
	char buf[1024];
	char mapqInps[1024];
	if(rs != NULL) {
		staln.reset();
		rs->initStacked(rd, staln);
		staln.leftAlign(false /* not past MMs */);
	}
	int offAdj = 0;
	// QNAME
	samc_.printReadName(o, rd.name, flags.partOfPair());
	o.append('\t');
	// FLAG
	int fl = 0;
	if(flags.partOfPair()) {
		fl |= SAM_FLAG_PAIRED;
		if(flags.alignedConcordant()) {
			fl |= SAM_FLAG_MAPPED_PAIRED;
 		}
		if(!flags.mateAligned()) {
			// Other fragment is unmapped
			fl |= SAM_FLAG_MATE_UNMAPPED;
		}
		fl |= (flags.readMate1() ?
			   SAM_FLAG_FIRST_IN_PAIR : SAM_FLAG_SECOND_IN_PAIR);
		if(flags.mateAligned() && rso != NULL) {
			if(!rso->fw()) {
				fl |= SAM_FLAG_MATE_STRAND;
			}
		}
	}
	if(!flags.isPrimary()) {
		fl |= SAM_FLAG_NOT_PRIMARY;
	}
	if(rs != NULL && !rs->fw()) {
		fl |= SAM_FLAG_QUERY_STRAND;
	}
	if(rs == NULL) {
		// Failed to align
		fl |= SAM_FLAG_UNMAPPED;
	}
	itoa10<int>(fl, buf);
	o.append(buf);
	o.append('\t');
	// RNAME
	if(rs != NULL) {
		samc_.printRefNameFromIndex(o, (size_t)rs->refid());
		o.append('\t');
	} else {
		if(summ.orefid() != -1) {
			// Opposite mate aligned but this one didn't - print the opposite
			// mate's RNAME and POS as is customary
			assert(flags.partOfPair());
			samc_.printRefNameFromIndex(o, (size_t)summ.orefid());
		} else {		
			// No alignment
			o.append('*');
		}
		o.append('\t');
	}
	// POS
	// Note: POS is *after* soft clipping.  I.e. POS points to the
	// upstream-most character *involved in the clipped alignment*.
	if(rs != NULL) {
		itoa10<int64_t>(rs->refoff()+1+offAdj, buf);
		o.append(buf);
		o.append('\t');
	} else {
		if(summ.orefid() != -1) {
			// Opposite mate aligned but this one didn't - print the opposite
			// mate's RNAME and POS as is customary
			assert(flags.partOfPair());
			itoa10<int64_t>(summ.orefoff()+1+offAdj, buf);
			o.append(buf);
		} else {
			// No alignment
			o.append('0');
		}
		o.append('\t');
	}
	// MAPQ
	mapqInps[0] = '\0';
	if(rs != NULL) {
		itoa10<TMapq>(mapqCalc.mapq(
									summ, flags, rd.mate < 2, rd.length(),
									rdo == NULL ? 0 : rdo->length(), mapqInps), buf);
		o.append(buf);
		o.append('\t');
	} else {
		// No alignment
		o.append("0\t");
	}
	// CIGAR
	if(rs != NULL) {
		staln.buildCigar(false);
		staln.writeCigar(&o, NULL);
		o.append('\t');
	} else {
		// No alignment
		o.append("*\t");
	}
	// RNEXT
	if(rs != NULL && flags.partOfPair()) {
		if(rso != NULL && rs->refid() != rso->refid()) {
			samc_.printRefNameFromIndex(o, (size_t)rso->refid());
			o.append('\t');
		} else {
			o.append("=\t");
		}
	} else if(summ.orefid() != -1) {
		// The convention if this mate fails to align but the other doesn't is
		// to copy the mate's details into here
		o.append("=\t");
	} else {
		o.append("*\t");
	}
	// PNEXT
	if(rs != NULL && flags.partOfPair()) {
		if(rso != NULL) {
			itoa10<int64_t>(rso->refoff()+1, buf);
			o.append(buf);
			o.append('\t');
		} else {
			// The convenstion is that if this mate aligns but the opposite
			// doesn't, we print this mate's offset here
			itoa10<int64_t>(rs->refoff()+1, buf);
			o.append(buf);
			o.append('\t');
		}
	} else if(summ.orefid() != -1) {
		// The convention if this mate fails to align but the other doesn't is
		// to copy the mate's details into here
		itoa10<int64_t>(summ.orefoff()+1, buf);
		o.append(buf);
		o.append('\t');
	} else {
		o.append("0\t");
	}
	// ISIZE
	if(rs != NULL && rs->isFraglenSet()) {
		itoa10<int64_t>(rs->fragmentLength(), buf);
		o.append(buf);
		o.append('\t');
	} else {
		// No fragment
		o.append("0\t");
	}
	// SEQ
	if(!flags.isPrimary() && samc_.omitSecondarySeqQual()) {
		o.append('*');
	} else {
		// Print the read
		if(rd.patFw.length() == 0) {
			o.append('*');
		} else {
			if(rs == NULL || rs->fw()) {
				o.append(rd.patFw.toZBuf());
			} else {
				o.append(rd.patRc.toZBuf());
			}
		}
	}
	o.append('\t');
	// QUAL
	if(!flags.isPrimary() && samc_.omitSecondarySeqQual()) {
		o.append('*');
	} else {
		// Print the quals
		if(rd.qual.length() == 0) {
			o.append('*');
		} else {
			if(rs == NULL || rs->fw()) {
				o.append(rd.qual.toZBuf());
			} else {
				o.append(rd.qualRev.toZBuf());
			}
		}
	}
	o.append('\t');
	//
	// Optional fields
	//
	if(rs != NULL) {
		samc_.printAlignedOptFlags(
								   o,           // output buffer
								   true,        // first opt flag printed is first overall?
								   rd,          // read
								   *rs,         // individual alignment result
								   staln,       // stacked alignment
								   flags,       // alignment flags
								   summ,        // summary of alignments for this read
								   ssm,         // seed alignment summary
								   prm,         // per-read metrics
								   sc,          // scoring scheme
								   mapqInps,    // inputs to MAPQ calculation
                                   this->altdb_);
	} else {
		samc_.printEmptyOptFlags(
								 o,           // output buffer
								 true,        // first opt flag printed is first overall?
								 rd,          // read
								 flags,       // alignment flags
								 summ,        // summary of alignments for this read
								 ssm,         // seed alignment summary
								 prm,         // per-read metrics
								 sc);         // scoring scheme
	}
	o.append('\n');
}

#endif /*ndef ALN_SINK_H_*/
