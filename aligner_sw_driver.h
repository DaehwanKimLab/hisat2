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

/*
 *  aligner_sw_driver.h
 *
 * REDUNDANT SEED HITS
 *
 * We say that two seed hits are redundant if they trigger identical
 * seed-extend dynamic programming problems.  Put another way, they both lie on
 * the same diagonal of the overall read/reference dynamic programming matrix.
 * Detecting redundant seed hits is simple when the seed hits are ungapped.  We
 * do this after offset resolution but before the offset is converted to genome
 * coordinates (see uses of the seenDiags1_/seenDiags2_ fields for examples).
 *
 * REDUNDANT ALIGNMENTS
 *
 * In an unpaired context, we say that two alignments are redundant if they
 * share any cells in the global DP table.  Roughly speaking, this is like
 * saying that two alignments are redundant if any read character aligns to the
 * same reference character (same reference sequence, same strand, same offset)
 * in both alignments.
 *
 * In a paired-end context, we say that two paired-end alignments are redundant
 * if the mate #1s are redundant and the mate #2s are redundant.
 *
 * How do we enforce this?  In the unpaired context, this is relatively simple:
 * the cells from each alignment are checked against a set containing all cells
 * from all previous alignments.  Given a new alignment, for each cell in the
 * new alignment we check whether it is in the set.  If there is any overlap,
 * the new alignment is rejected as redundant.  Otherwise, the new alignment is
 * accepted and its cells are added to the set.
 *
 * Enforcement in a paired context is a little trickier.  Consider the
 * following approaches:
 *
 * 1. Skip anchors that are redundant with any previous anchor or opposite
 *    alignment.  This is sufficient to ensure no two concordant alignments
 *    found are redundant.
 *
 * 2. Same as scheme 1, but with a "transitive closure" scheme for finding all
 *    concordant pairs in the vicinity of an anchor.  Consider the AB/AC
 *    scenario from the previous paragraph.  If B is the anchor alignment, we
 *    will find AB but not AC.  But under this scheme, once we find AB we then
 *    let B be a new anchor and immediately look for its opposites.  Likewise,
 *    if we find any opposite, we make them anchors and continue searching.  We
 *    don't stop searching until every opposite is used as an anchor.
 *
 * 3. Skip anchors that are redundant with any previous anchor alignment (but
 *    allow anchors that are redundant with previous opposite alignments).
 *    This isn't sufficient to avoid redundant concordant alignments.  To avoid
 *    redundant concordants, we need an additional procedure that checks each
 *    new concordant alignment one-by-one against a list of previous concordant
 *    alignments to see if it is redundant.
 *
 * We take approach 1.
 */

#ifndef ALIGNER_SW_DRIVER_H_
#define ALIGNER_SW_DRIVER_H_

#include <iostream>
// -- BTL remove --
#include <stdlib.h>
#include <sys/time.h>
// -- --
#include <utility>
#include "ds.h"
#include "aligner_seed.h"
#include "aligner_sw.h"
#include "aligner_cache.h"
#include "reference.h"
#include "group_walk.h"
#include "gfm.h"
#include "mem_ids.h"
#include "aln_sink.h"
#include "pe.h"
#include "ival_list.h"
#include "simple_func.h"
#include "random_util.h"
#include "dp_framer.h"

using namespace std;

template <typename index_t>
struct SeedPos {

	SeedPos() : fw(false), offidx(0), rdoff(0), seedlen(0) { }

	SeedPos(
		bool fw_,
		index_t offidx_,
		index_t rdoff_,
		index_t seedlen_)
	{
		init(fw_, offidx_, rdoff_, seedlen_);
	}
	
	void init(
		bool fw_,
		index_t offidx_,
		index_t rdoff_,
		index_t seedlen_)
	{
		fw      = fw_;
		offidx  = offidx_;
		rdoff   = rdoff_;
		seedlen = seedlen_;
	}
	
	bool operator<(const SeedPos& o) const {
		if(offidx < o.offidx)   return true;
		if(offidx > o.offidx)   return false;
		if(rdoff < o.rdoff)     return true;
		if(rdoff > o.rdoff)     return false;
		if(seedlen < o.seedlen) return true;
		if(seedlen > o.seedlen) return false;
		if(fw && !o.fw)         return true;
		if(!fw && o.fw)         return false;
		return false;
	}
	
	bool operator>(const SeedPos& o) const {
		if(offidx < o.offidx)   return false;
		if(offidx > o.offidx)   return true;
		if(rdoff < o.rdoff)     return false;
		if(rdoff > o.rdoff)     return true;
		if(seedlen < o.seedlen) return false;
		if(seedlen > o.seedlen) return true;
		if(fw && !o.fw)         return false;
		if(!fw && o.fw)         return true;
		return false;
	}
	
	bool operator==(const SeedPos& o) const {
		return fw == o.fw && offidx == o.offidx &&
		       rdoff == o.rdoff && seedlen == o.seedlen;
	}

	bool fw;
	index_t offidx;
	index_t rdoff;
	index_t seedlen;
};

/**
 * An SATuple along with the associated seed position.
 */
template <typename index_t>
struct SATupleAndPos {
	
	SATuple<index_t> sat;    // result for this seed hit
	SeedPos<index_t> pos;    // seed position that yielded the range this was taken from
	index_t          origSz; // size of range this was taken from
	index_t          nlex;   // # position we can extend seed hit to left w/o edit
	index_t          nrex;   // # position we can extend seed hit to right w/o edit
	
	bool operator<(const SATupleAndPos& o) const {
		if(sat < o.sat) return true;
		if(sat > o.sat) return false;
		return pos < o.pos;
	}

	bool operator==(const SATupleAndPos& o) const {
		return sat == o.sat && pos == o.pos;
	}
};

/**
 * Encapsulates the weighted random sampling scheme we want to use to pick
 * which seed hit range to sample a row from.
 */
template <typename index_t>
class RowSampler {

public:

	RowSampler(int cat = 0) : elim_(cat), masses_(cat) { 
		mass_ = 0.0f;
	}
	
	/**
	 * Initialze sampler with respect to a range of elements in a list of
	 * SATupleAndPos's.
	 */
	void init(
		const EList<SATupleAndPos<index_t>, 16>& salist,
		index_t sai,
		index_t saf,
		bool lensq, // whether to square the numerator, which = extended length
		bool szsq)  // whether to square denominator, which = 
	{
		assert_gt(saf, sai);
		elim_.resize(saf - sai);
		elim_.fill(false);
		// Initialize mass
		mass_ = 0.0f;
		masses_.resize(saf - sai);
		for(index_t i = sai; i < saf; i++) {
			index_t len = salist[i].nlex + salist[i].nrex + 1; // + salist[i].sat.key.len;
			double num = (double)len;
			if(lensq) {
				num *= num;
			}
			double denom = (double)salist[i].sat.size();
			if(szsq) {
				denom *= denom;
			}
			masses_[i - sai] = num / denom;
			mass_ += masses_[i - sai];
		}
	}
	
	/**
	 * Caller is indicating that the bin at index i is exhausted and we should
	 * exclude it from our sampling from now on.
	 */
	void finishedRange(index_t i) {
		assert_lt(i, masses_.size());
		elim_[i] = true;
		mass_ -= masses_[i];
	}
	
	/**
	 * Sample randomly from the mass.
	 */
	size_t next(RandomSource& rnd) {
		// Throw the dart
		double rd = rnd.nextFloat() * mass_;
		double mass_sofar = 0.0f;
		size_t sz = masses_.size();
		size_t last_unelim = std::numeric_limits<size_t>::max();
		for(size_t i = 0; i < sz; i++) {
			if(!elim_[i]) {
				last_unelim = i;
				mass_sofar += masses_[i];
				if(rd < mass_sofar) {
					// This is the one we hit
					return i;
				}
			}
		}
		assert_neq(std::numeric_limits<size_t>::max(), last_unelim);
		return last_unelim;
	}

protected:
	double        mass_;    // total probability mass to throw darts at
	EList<bool>   elim_;    // whether the range is eliminated
	EList<double> masses_;  // mass of each range
};

/**
 * Return values from extendSeeds and extendSeedsPaired.
 */
enum {
	// All end-to-end and seed hits were examined
	// The policy does not need us to look any further
	EXTEND_EXHAUSTED_CANDIDATES = 1,
	EXTEND_POLICY_FULFILLED,
	// We stopped because we reached a point where the only remaining
	// alignments of interest have perfect scores, but we already investigated
	// perfect alignments
	EXTEND_PERFECT_SCORE,
	// We stopped because we ran up against a limit on how much work we should
	// do for one set of seed ranges, e.g. the limit on number of consecutive
	// unproductive DP extensions
	EXTEND_EXCEEDED_SOFT_LIMIT,
	// We stopped because we ran up against a limit on how much work we should
	// do for overall before giving up on a mate
	EXTEND_EXCEEDED_HARD_LIMIT
};

/**
 * Data structure encapsulating a range that's been extended out in two
 * directions.
 */
struct ExtendRange {

	void init(size_t off_, size_t len_, size_t sz_) {
		off = off_; len = len_; sz = sz_;
	}

	size_t off; // offset of extended region
	size_t len; // length between extremes of extended region
	size_t sz;  // # of elements in SA range
};

template <typename index_t>
class SwDriver {

	typedef PList<index_t, CACHE_PAGE_SZ> TSAList;

public:

	SwDriver(size_t bytes) :
		satups_(DP_CAT),
		gws_(DP_CAT),
		seenDiags1_(DP_CAT),
		seenDiags2_(DP_CAT),
		redAnchor_(DP_CAT),
		redMate1_(DP_CAT),
		redMate2_(DP_CAT),
		pool_(bytes, CACHE_PAGE_SZ, DP_CAT),
		salistEe_(DP_CAT),
		gwstate_(GW_CAT) { }

	/**
	 * Given a collection of SeedHits for a single read, extend seed alignments
	 * into full alignments.  Where possible, try to avoid redundant offset
	 * lookups and dynamic programming problems.  Optionally report alignments
	 * to a AlnSinkWrap object as they are discovered.
	 *
	 * If 'reportImmediately' is true, returns true iff a call to
	 * mhs->report() returned true (indicating that the reporting
	 * policy is satisfied and we can stop).  Otherwise, returns false.
	 */
	int extendSeeds(
		Read& rd,                    // read to align
		bool mate1,                  // true iff rd is mate #1
		SeedResults<index_t>& sh,    // seed hits to extend into full alignments
		const GFM<index_t>& gfmFw,   // BWT
		const GFM<index_t>* gfmBw,   // BWT'
		const BitPairReference& ref, // Reference strings
		SwAligner& swa,              // dynamic programming aligner
		const Scoring& sc,           // scoring scheme
		int seedmms,                 // # mismatches allowed in seed
		int seedlen,                 // length of seed
		int seedival,                // interval between seeds
		TAlScore& minsc,             // minimum score for anchor
		int nceil,                   // maximum # Ns permitted in ref portion
		size_t maxhalf,              // maximum width on one side of DP table
		bool doUngapped,             // do ungapped alignment
		size_t maxIters,             // stop after this many seed-extend loop iters
		size_t maxUg,                // max # ungapped extends
		size_t maxDp,                // max # DPs
		size_t maxUgStreak,          // stop after streak of this many ungap fails
		size_t maxDpStreak,          // stop after streak of this many dp fails
		bool doExtend,               // do seed extension
		bool enable8,                // use 8-bit SSE where possible
		size_t cminlen,              // use checkpointer if read longer than this
		size_t cpow2,                // interval between diagonals to checkpoint
		bool doTri,                  // triangular mini-fills
		int tighten,                 // -M score tightening mode
		AlignmentCacheIface<index_t>& ca,     // alignment cache for seed hits
		RandomSource& rnd,           // pseudo-random source
		WalkMetrics& wlm,            // group walk left metrics
		SwMetrics& swmSeed,          // DP metrics for seed-extend
		PerReadMetrics& prm,         // per-read metrics
		AlnSinkWrap<index_t>* mhs,   // HitSink for multiseed-style aligner
		bool reportImmediately,      // whether to report hits immediately to mhs
		bool& exhaustive);

	/**
	 * Given a collection of SeedHits for a read pair, extend seed
	 * alignments into full alignments and then look for the opposite
	 * mate using dynamic programming.  Where possible, try to avoid
	 * redundant offset lookups.  Optionally report alignments to a
	 * AlnSinkWrap object as they are discovered.
	 *
	 * If 'reportImmediately' is true, returns true iff a call to
	 * mhs->report() returned true (indicating that the reporting
	 * policy is satisfied and we can stop).  Otherwise, returns false.
	 */
	int extendSeedsPaired(
		Read& rd,                    // mate to align as anchor
		Read& ord,                   // mate to align as opposite
		bool anchor1,                // true iff anchor mate is mate1
		bool oppFilt,                // true iff opposite mate was filtered out
		SeedResults<index_t>& sh,    // seed hits for anchor
		const GFM<index_t>& gfmFw,   // BWT
		const GFM<index_t>* gfmBw,   // BWT'
		const BitPairReference& ref, // Reference strings
		SwAligner& swa,              // dyn programming aligner for anchor
		SwAligner& swao,             // dyn programming aligner for opposite
		const Scoring& sc,           // scoring scheme
		const PairedEndPolicy& pepol,// paired-end policy
		int seedmms,                 // # mismatches allowed in seed
		int seedlen,                 // length of seed
		int seedival,                // interval between seeds
		TAlScore& minsc,             // minimum score for anchor
		TAlScore& ominsc,            // minimum score for opposite
		int nceil,                   // max # Ns permitted in ref for anchor
		int onceil,                  // max # Ns permitted in ref for opposite
		bool nofw,                   // don't align forward read
		bool norc,                   // don't align revcomp read
		size_t maxhalf,              // maximum width on one side of DP table
		bool doUngapped,             // do ungapped alignment
		size_t maxIters,             // stop after this many seed-extend loop iters
		size_t maxUg,                // max # ungapped extends
		size_t maxDp,                // max # DPs
		size_t maxEeStreak,          // stop after streak of this many end-to-end fails
		size_t maxUgStreak,          // stop after streak of this many ungap fails
		size_t maxDpStreak,          // stop after streak of this many dp fails
		size_t maxMateStreak,        // stop seed range after N mate-find fails
		bool doExtend,               // do seed extension
		bool enable8,                // use 8-bit SSE where possible
		size_t cminlen,              // use checkpointer if read longer than this
		size_t cpow2,                // interval between diagonals to checkpoint
		bool doTri,                  // triangular mini-fills
		int tighten,                 // -M score tightening mode
		AlignmentCacheIface<index_t>& cs,     // alignment cache for seed hits
		RandomSource& rnd,           // pseudo-random source
		WalkMetrics& wlm,            // group walk left metrics
		SwMetrics& swmSeed,          // DP metrics for seed-extend
		SwMetrics& swmMate,          // DP metrics for mate finidng
		PerReadMetrics& prm,         // per-read metrics for anchor
		AlnSinkWrap<index_t>* msink, // AlnSink wrapper for multiseed-style aligner
		bool swMateImmediately,      // whether to look for mate immediately
		bool reportImmediately,      // whether to report hits immediately to msink
		bool discord,                // look for discordant alignments?
		bool mixed,                  // look for unpaired as well as paired alns?
		bool& exhaustive);
        
	/**
	 * Prepare for a new read.
	 */
	void nextRead(bool paired, size_t mate1len, size_t mate2len) {
		redAnchor_.reset();
		seenDiags1_.reset();
		seenDiags2_.reset();
		seedExRangeFw_[0].clear(); // mate 1 fw
		seedExRangeFw_[1].clear(); // mate 2 fw
		seedExRangeRc_[0].clear(); // mate 1 rc
		seedExRangeRc_[1].clear(); // mate 2 rc
		size_t maxlen = mate1len;
		if(paired) {
			redMate1_.reset();
			redMate1_.init(mate1len);
			redMate2_.reset();
			redMate2_.init(mate2len);
			if(mate2len > maxlen) {
				maxlen = mate2len;
			}
		}
		redAnchor_.init(maxlen);
	}

protected:

	bool eeSaTups(
		const Read& rd,              // read
		SeedResults<index_t>& sh,    // seed hits to extend into full alignments
		const GFM<index_t>& gfm,     // BWT
		const BitPairReference& ref, // Reference strings
		RandomSource& rnd,           // pseudo-random generator
		WalkMetrics& wlm,            // group walk left metrics
		SwMetrics& swmSeed,          // metrics for seed extensions
		index_t& nelt_out,           // out: # elements total
        index_t maxelts,             // max # elts to report
		bool all);                   // report all hits?
    
    void extend(
		const Read& rd,       // read
		const GFM<index_t>& gfmFw,   // Forward Bowtie index
		const GFM<index_t>* gfmBw,   // Backward Bowtie index
		index_t topf,        // top in fw index
		index_t botf,        // bot in fw index
		index_t topb,        // top in bw index
		index_t botb,        // bot in bw index
		bool fw,             // seed orientation
		index_t off,         // seed offset from 5' end
		index_t len,         // seed length
		PerReadMetrics& prm, // per-read metrics
		index_t& nlex,       // # positions we can extend to left w/o edit
		index_t& nrex);      // # positions we can extend to right w/o edit

	void prioritizeSATups(
		const Read& rd,              // read
		SeedResults<index_t>& sh,    // seed hits to extend into full alignments
		const GFM<index_t>& gfmFw, // BWT
		const GFM<index_t>* gfmBw, // BWT'
		const BitPairReference& ref, // Reference strings
		int seedmms,                 // # seed mismatches allowed
		index_t maxelt,              // max elts we'll consider
		bool doExtend,               // extend out seeds
		bool lensq,                  // square extended length
		bool szsq,                   // square SA range size
		index_t nsm,                 // if range as <= nsm elts, it's "small"
		AlignmentCacheIface<index_t>& ca,     // alignment cache for seed hits
		RandomSource& rnd,           // pseudo-random generator
		WalkMetrics& wlm,            // group walk left metrics
		PerReadMetrics& prm,         // per-read metrics
		index_t& nelt_out,           // out: # elements total
		bool all);                   // report all hits?

	Random1toN               rand_;    // random number generators
	EList<Random1toN, 16>    rands_;   // random number generators
	EList<Random1toN, 16>    rands2_;  // random number generators
	EList<EEHit<index_t>, 16>         eehits_;  // holds end-to-end hits
	EList<SATupleAndPos<index_t>, 16> satpos_;  // holds SATuple, SeedPos pairs
	EList<SATupleAndPos<index_t>, 16> satpos2_; // holds SATuple, SeedPos pairs
	EList<SATuple<index_t>, 16>       satups_;  // holds SATuples to explore elements from
	EList<GroupWalk2S<index_t, TSlice, 16> > gws_;   // list of GroupWalks; no particular order
	EList<size_t>            mateStreaks_; // mate-find fail streaks
	RowSampler<index_t>      rowsamp_;     // row sampler
	
	// Ranges that we've extended through when extending seed hits
	EList<ExtendRange> seedExRangeFw_[2];
	EList<ExtendRange> seedExRangeRc_[2];

	// Data structures encapsulating the diagonals that have already been used
	// to seed alignment for mate 1 and mate 2.
	EIvalMergeListBinned seenDiags1_;
	EIvalMergeListBinned seenDiags2_;

	// For weeding out redundant alignments
	RedundantAlns  redAnchor_;  // database of cells used for anchor alignments
	RedundantAlns  redMate1_;   // database of cells used for mate 1 alignments
	RedundantAlns  redMate2_;   // database of cells used for mate 2 alignments

	// For holding results for anchor (res_) and opposite (ores_) mates
	SwResult       resGap_;    // temp holder for alignment result
	SwResult       oresGap_;   // temp holder for alignment result, opp mate
	SwResult       resUngap_;  // temp holder for ungapped alignment result
	SwResult       oresUngap_; // temp holder for ungap. aln. opp mate
	SwResult       resEe_;     // temp holder for ungapped alignment result
	SwResult       oresEe_;    // temp holder for ungap. aln. opp mate
	
	Pool           pool_;      // memory pages for salistExact_
	TSAList        salistEe_;  // PList for offsets for end-to-end hits
	GroupWalkState<index_t> gwstate_;   // some per-thread state shared by all GroupWalks
	
	// For AlnRes::matchesRef:
	ASSERT_ONLY(SStringExpandable<char>     raw_refbuf_);
	ASSERT_ONLY(SStringExpandable<uint32_t> raw_destU32_);
	ASSERT_ONLY(EList<bool>                 raw_matches_);
	ASSERT_ONLY(BTDnaString                 tmp_rf_);
	ASSERT_ONLY(BTDnaString                 tmp_rdseq_);
	ASSERT_ONLY(BTString                    tmp_qseq_);
	ASSERT_ONLY(EList<index_t>              tmp_reflens_);
	ASSERT_ONLY(EList<index_t>              tmp_refoffs_);
};

#define TIMER_START() \
struct timeval tv_i, tv_f; \
struct timezone tz_i, tz_f; \
size_t total_usecs; \
gettimeofday(&tv_i, &tz_i)

#define IF_TIMER_END() \
gettimeofday(&tv_f, &tz_f); \
total_usecs = \
(tv_f.tv_sec - tv_i.tv_sec) * 1000000 + (tv_f.tv_usec - tv_i.tv_usec); \
if(total_usecs > 300000)

/*
 * aligner_sw_driver.cpp
 *
 * Routines that drive the alignment process given a collection of seed hits.
 * This is generally done in a few stages: extendSeeds visits the set of
 * seed-hit BW elements in some order; for each element visited it resolves its
 * reference offset; once the reference offset is known, bounds for a dynamic
 * programming subproblem are established; if these bounds are distinct from
 * the bounds we've already tried, we solve the dynamic programming subproblem
 * and report the hit; if the AlnSinkWrap indicates that we can stop, we
 * return, otherwise we continue on to the next BW element.
 */


/**
 * Given end-to-end alignment results stored in the SeedResults structure, set
 * up all of our state for resolving and keeping track of reference offsets for
 * hits.  Order the list of ranges to examine such that all exact end-to-end
 * alignments are examined before any 1mm end-to-end alignments.
 *
 * Note: there might be a lot of hits and a lot of wide ranges to look for
 * here.  We use 'maxelt'.
 */
template <typename index_t>
bool SwDriver<index_t>::eeSaTups(
								 const Read& rd,              // read
								 SeedResults<index_t>& sh,    // seed hits to extend into full alignments
                                 const GFM<index_t>& gfm,     // BWT
								 const BitPairReference& ref, // Reference strings
								 RandomSource& rnd,           // pseudo-random generator
								 WalkMetrics& wlm,            // group walk left metrics
								 SwMetrics& swmSeed,          // metrics for seed extensions
								 index_t& nelt_out,           // out: # elements total
								 index_t maxelt,              // max elts we'll consider
								 bool all)                    // report all hits?
{
    assert_eq(0, nelt_out);
	gws_.clear();
	rands_.clear();
	satpos_.clear();
	eehits_.clear();
	// First, count up the total number of satpos_, rands_, eehits_, and gws_
	// we're going to tuse
	index_t nobj = 0;
	if(!sh.exactFwEEHit().empty()) nobj++;
	if(!sh.exactRcEEHit().empty()) nobj++;
	nobj += sh.mm1EEHits().size();
    nobj = min(nobj, maxelt);
	gws_.ensure(nobj);
	rands_.ensure(nobj);
	satpos_.ensure(nobj);
	eehits_.ensure(nobj);
	index_t tot = sh.exactFwEEHit().size() + sh.exactRcEEHit().size();
	bool succ = false;
	bool firstEe = true;
    bool done = false;
	if(tot > 0) {
		bool fwFirst = true;
        // Pick fw / rc to go first in a weighted random fashion
#ifdef BOWTIE_64BIT_INDEX
		index_t rn64 = rnd.nextU64();
		index_t rn = rn64 % (uint64_t)tot;
#else
		index_t rn32 = rnd.nextU32();
		index_t rn = rn32 % (uint32_t)tot;
#endif
		if(rn >= sh.exactFwEEHit().size()) {
			fwFirst = false;
		}
		for(int fwi = 0; fwi < 2 && !done; fwi++) {
			bool fw = ((fwi == 0) == fwFirst);
			EEHit<index_t> hit = fw ? sh.exactFwEEHit() : sh.exactRcEEHit();
			if(hit.empty()) {
				continue;
			}
			assert(hit.fw == fw);
			if(hit.bot > hit.top) {
                // Possibly adjust bot and width if we would have exceeded maxelt
                index_t tops[2] = { hit.top, 0 };
                index_t bots[2] = { hit.bot, 0 };
                index_t width = hit.bot - hit.top;
                if(nelt_out + width > maxelt) {
                    index_t trim = (index_t)((nelt_out + width) - maxelt);
#ifdef BOWTIE_64BIT_INDEX
                    index_t rn = rnd.nextU64() % width;
#else
                    index_t rn = rnd.nextU32() % width;
#endif
                    index_t newwidth = width - trim;
                    if(hit.top + rn + newwidth > hit.bot) {
                        // Two pieces
                        tops[0] = hit.top + rn;
                        bots[0] = hit.bot;
                        tops[1] = hit.top;
                        bots[1] = hit.top + newwidth - (bots[0] - tops[0]);
                    } else {
                        // One piece
                        tops[0] = hit.top + rn;
                        bots[0] = tops[0] + newwidth;
                    }
                    assert_leq(bots[0], hit.bot);
                    assert_leq(bots[1], hit.bot);
                    assert_geq(bots[0], tops[0]);
                    assert_geq(bots[1], tops[1]);
                    assert_eq(newwidth, (bots[0] - tops[0]) + (bots[1] - tops[1]));
                }
                for(int i = 0; i < 2 && !done; i++) {
                    if(bots[i] <= tops[i]) break;
                    index_t width = bots[i] - tops[i];
                    index_t top = tops[i];
                    // Clear list where resolved offsets are stored
                    swmSeed.exranges++;
                    swmSeed.exrows += width;
                    if(!succ) {
                        swmSeed.exsucc++;
                        succ = true;
                    }
                    if(firstEe) {
                        salistEe_.clear();
                        pool_.clear();
                        firstEe = false;
                    }
                    // We have to be careful not to allocate excessive amounts of memory here
                    TSlice o(salistEe_, (index_t)salistEe_.size(), width);
                    for(index_t i = 0; i < width; i++) {
                        if(!salistEe_.add(pool_, (index_t)OFF_MASK)) {
                            swmSeed.exooms++;
                            return false;
                        }
                    }
                    assert(!done);
                    eehits_.push_back(hit);
                    satpos_.expand();
                    satpos_.back().sat.init(SAKey(), top, (index_t)OFF_MASK, o);
                    satpos_.back().sat.key.seq = MAX_U64;
                    satpos_.back().sat.key.len = (index_t)rd.length();
                    satpos_.back().pos.init(fw, 0, 0, (index_t)rd.length());
                    satpos_.back().origSz = width;
                    rands_.expand();
                    rands_.back().init(width, all);
                    gws_.expand();
					SARangeWithOffs<TSlice, index_t> sa;
					sa.topf = satpos_.back().sat.topf;
					sa.len = satpos_.back().sat.key.len;
					sa.offs = satpos_.back().sat.offs;
                    gws_.back().init(
									 gfm,                // forward Bowtie index
									 ref,                // reference sequences
									 sa,                 // SATuple
									 rnd,                // pseudo-random generator
									 wlm);               // metrics
                    assert(gws_.back().repOk(sa));
                    nelt_out += width;
                    if(nelt_out >= maxelt) {
                        done = true;
                    }
                }
			}
		}
	}
	succ = false;
	if(!done && !sh.mm1EEHits().empty()) {
		sh.sort1mmEe(rnd);
		index_t sz = sh.mm1EEHits().size();
		for(index_t i = 0; i < sz && !done; i++) {
			EEHit<index_t> hit = sh.mm1EEHits()[i];
			assert(hit.repOk(rd));
			assert(!hit.empty());
            // Possibly adjust bot and width if we would have exceeded maxelt
            index_t tops[2] = { hit.top, 0 };
            index_t bots[2] = { hit.bot, 0 };
            index_t width = hit.bot - hit.top;
            if(nelt_out + width > maxelt) {
                index_t trim = (index_t)((nelt_out + width) - maxelt);
#ifdef BOWTIE_64BIT_INDEX
                index_t rn = rnd.nextU64() % width;
#else
                index_t rn = rnd.nextU32() % width;
#endif
                index_t newwidth = width - trim;
                if(hit.top + rn + newwidth > hit.bot) {
                    // Two pieces
                    tops[0] = hit.top + rn;
                    bots[0] = hit.bot;
                    tops[1] = hit.top;
                    bots[1] = hit.top + newwidth - (bots[0] - tops[0]);
                } else {
                    // One piece
                    tops[0] = hit.top + rn;
                    bots[0] = tops[0] + newwidth;
                }
                assert_leq(bots[0], hit.bot);
                assert_leq(bots[1], hit.bot);
                assert_geq(bots[0], tops[0]);
                assert_geq(bots[1], tops[1]);
                assert_eq(newwidth, (bots[0] - tops[0]) + (bots[1] - tops[1]));
            }
            for(int i = 0; i < 2 && !done; i++) {
                if(bots[i] <= tops[i]) break;
                index_t width = bots[i] - tops[i];
                index_t top = tops[i];
                // Clear list where resolved offsets are stored
                swmSeed.mm1ranges++;
                swmSeed.mm1rows += width;
                if(!succ) {
                    swmSeed.mm1succ++;
                    succ = true;
                }
                if(firstEe) {
                    salistEe_.clear();
                    pool_.clear();
                    firstEe = false;
                }
                TSlice o(salistEe_, (index_t)salistEe_.size(), width);
                for(size_t i = 0; i < width; i++) {
                    if(!salistEe_.add(pool_, (index_t)OFF_MASK)) {
                        swmSeed.mm1ooms++;
                        return false;
                    }
                }
                eehits_.push_back(hit);
                satpos_.expand();
                satpos_.back().sat.init(SAKey(), top, (index_t)OFF_MASK, o);
                satpos_.back().sat.key.seq = MAX_U64;
                satpos_.back().sat.key.len = (index_t)rd.length();
                satpos_.back().pos.init(hit.fw, 0, 0, (index_t)rd.length());
                satpos_.back().origSz = width;
                rands_.expand();
                rands_.back().init(width, all);
                gws_.expand();
				SARangeWithOffs<TSlice, index_t> sa;
				sa.topf = satpos_.back().sat.topf;
				sa.len = satpos_.back().sat.key.len;
				sa.offs = satpos_.back().sat.offs;
                gws_.back().init(
								 gfm,  // forward Bowtie index
								 ref,  // reference sequences
								 sa,   // SATuple
								 rnd,  // pseudo-random generator
								 wlm); // metrics
                assert(gws_.back().repOk(sa));
                nelt_out += width;
                if(nelt_out >= maxelt) {
                    done = true;
                }
            }
		}
	}
	return true;
}

/**
 * Extend a seed hit out on either side.  Requires that we know the seed hit's
 * offset into the read and orientation.  Also requires that we know top/bot
 * for the seed hit in both the forward and (if we want to extend to the right)
 * reverse index.
 */
template <typename index_t>
void SwDriver<index_t>::extend(
							   const Read& rd,       // read
							   const GFM<index_t>& gfmFw,   // Forward Bowtie index
							   const GFM<index_t>* gfmBw,   // Backward Bowtie index
							   index_t topf,         // top in fw index
							   index_t botf,         // bot in fw index
							   index_t topb,         // top in bw index
							   index_t botb,         // bot in bw index
							   bool fw,              // seed orientation
							   index_t off,          // seed offset from 5' end
							   index_t len,          // seed length
							   PerReadMetrics& prm,  // per-read metrics
							   index_t& nlex,        // # positions we can extend to left w/o edit
							   index_t& nrex)        // # positions we can extend to right w/o edit
{
	index_t t[4], b[4];
	index_t tp[4], bp[4];
	SideLocus<index_t> tloc, bloc;
	index_t rdlen = (index_t)rd.length();
	index_t lim = fw ? off : rdlen - len - off;
	// We're about to add onto the beginning, so reverse it
#ifndef NDEBUG
	if(false) {
		// TODO: This will sometimes fail even when the extension is legitimate
		// This is because contains() comes in from one extreme or the other,
		// whereas we started from the inside and worked outwards.  This
		// affects which Ns are OK and which are not OK.
		
		// Have to do both because whether we can get through an N depends on
		// which direction we're coming in
		bool fwContains = gfmFw.contains(tmp_rdseq_);
		tmp_rdseq_.reverse();
		bool bwContains = gfmBw != NULL && gfmBw->contains(tmp_rdseq_);
		tmp_rdseq_.reverse();
		assert(fwContains || bwContains);
	}
#endif
	ASSERT_ONLY(tmp_rdseq_.reverse());
	if(lim > 0) {
		const GFM<index_t> *gfm = &gfmFw;
		assert(gfm != NULL);
		// Extend left using forward index
		const BTDnaString& seq = fw ? rd.patFw : rd.patRc;
		// See what we get by extending 
		index_t top = topf, bot = botf;
		t[0] = t[1] = t[2] = t[3] = 0;
		b[0] = b[1] = b[2] = b[3] = 0;
		tp[0] = tp[1] = tp[2] = tp[3] = topb;
		bp[0] = bp[1] = bp[2] = bp[3] = botb;
		SideLocus<index_t> tloc, bloc;
		INIT_LOCS(top, bot, tloc, bloc, *gfm);
		for(index_t ii = 0; ii < lim; ii++) {
			// Starting to left of seed (<off) and moving left
			index_t i = 0;
			if(fw) {
				i = off - ii - 1;
			} else {
				i = rdlen - off - len - 1 - ii;
			}
			// Get char from read
			int rdc = seq.get(i);
			// See what we get by extending 
			if(bloc.valid()) {
				prm.nSdFmops++;
				t[0] = t[1] = t[2] = t[3] =
				b[0] = b[1] = b[2] = b[3] = 0;
				gfm->mapBiLFEx(tloc, bloc, t, b, tp, bp);
				SANITY_CHECK_4TUP(t, b, tp, bp);
				int nonz = -1;
				bool abort = false;
				size_t origSz = bot - top;
				for(int j = 0; j < 4; j++) {
					if(b[j] > t[j]) {
						if(nonz >= 0) {
							abort = true;
							break;
						}
						nonz = j;
						top = t[j]; bot = b[j];
					}
				}
				assert_leq(bot - top, origSz);
				if(abort || (nonz != rdc && rdc <= 3) || bot - top < origSz) {
					break;
				}
			} else {
				assert_eq(bot, top+1);
				prm.nSdFmops++;
				int c = gfm->mapLF1(top, tloc);
				if(c != rdc && rdc <= 3) {
					break;
				}
				bot = top + 1;
			}
			ASSERT_ONLY(tmp_rdseq_.append(rdc));
			if(++nlex == 255) {
				break;
			}
			INIT_LOCS(top, bot, tloc, bloc, *gfm);
		}
	}
	// We're about to add onto the end, so re-reverse
	ASSERT_ONLY(tmp_rdseq_.reverse());
	lim = fw ? rdlen - len - off : off;
	if(lim > 0 && gfmBw != NULL) {
		const GFM<index_t> *gfm = gfmBw;
		assert(gfm != NULL);
		// Extend right using backward index
		const BTDnaString& seq = fw ? rd.patFw : rd.patRc;
		// See what we get by extending 
		index_t top = topb, bot = botb;
		t[0] = t[1] = t[2] = t[3] = 0;
		b[0] = b[1] = b[2] = b[3] = 0;
		tp[0] = tp[1] = tp[2] = tp[3] = topf;
		bp[0] = bp[1] = bp[2] = bp[3] = botf;
		INIT_LOCS(top, bot, tloc, bloc, *gfm);
		for(index_t ii = 0; ii < lim; ii++) {
			// Starting to right of seed (<off) and moving right
			index_t i;
			if(fw) {
				i = ii + len + off;
			} else {
				i = rdlen - off + ii;
			}
			// Get char from read
			int rdc = seq.get(i);
			// See what we get by extending 
			if(bloc.valid()) {
				prm.nSdFmops++;
				t[0] = t[1] = t[2] = t[3] =
				b[0] = b[1] = b[2] = b[3] = 0;
				gfm->mapBiLFEx(tloc, bloc, t, b, tp, bp);
				SANITY_CHECK_4TUP(t, b, tp, bp);
				int nonz = -1;
				bool abort = false;
				size_t origSz = bot - top;
				for(int j = 0; j < 4; j++) {
					if(b[j] > t[j]) {
						if(nonz >= 0) {
							abort = true;
							break;
						}
						nonz = j;
						top = t[j]; bot = b[j];
					}
				}
				assert_leq(bot - top, origSz);
				if(abort || (nonz != rdc && rdc <= 3) || bot - top < origSz) {
					break;
				}
			} else {
				assert_eq(bot, top+1);
				prm.nSdFmops++;
				int c = gfm->mapLF1(top, tloc);
				if(c != rdc && rdc <= 3) {
					break;
				}
				bot = top + 1;
			}
			ASSERT_ONLY(tmp_rdseq_.append(rdc));
			if(++nrex == 255) {
				break;
			}
			INIT_LOCS(top, bot, tloc, bloc, *gfm);
		}
	}
#ifndef NDEBUG
	if(false) {
		// TODO: This will sometimes fail even when the extension is legitimate
		// This is because contains() comes in from one extreme or the other,
		// whereas we started from the inside and worked outwards.  This
		// affects which Ns are OK and which are not OK.
		
		// Have to do both because whether we can get through an N depends on
		// which direction we're coming in
		bool fwContains = gfmFw.contains(tmp_rdseq_);
		tmp_rdseq_.reverse();
		bool bwContains = gfmBw != NULL && gfmBw->contains(tmp_rdseq_);
		tmp_rdseq_.reverse();
		assert(fwContains || bwContains);
	}
#endif
	assert_lt(nlex, rdlen);
	assert_lt(nrex, rdlen);
	return;
}

/**
 * Given seed results, set up all of our state for resolving and keeping
 * track of reference offsets for hits.
 */
template <typename index_t>
void SwDriver<index_t>::prioritizeSATups(
										 const Read& read,            // read
										 SeedResults<index_t>& sh,    // seed hits to extend into full alignments
										 const GFM<index_t>& gfmFw,   // BWT
										 const GFM<index_t>* gfmBw,   // BWT
										 const BitPairReference& ref, // Reference strings
										 int seedmms,                 // # mismatches allowed in seed
										 index_t maxelt,              // max elts we'll consider
										 bool doExtend,               // do extension of seed hits?
										 bool lensq,                  // square length in weight calculation
										 bool szsq,                   // square range size in weight calculation
										 index_t nsm,                 // if range as <= nsm elts, it's "small"
										 AlignmentCacheIface<index_t>& ca,     // alignment cache for seed hits
										 RandomSource& rnd,           // pseudo-random generator
										 WalkMetrics& wlm,            // group walk left metrics
										 PerReadMetrics& prm,         // per-read metrics
										 index_t& nelt_out,           // out: # elements total
										 bool all)                    // report all hits?
{
	const index_t nonz = sh.nonzeroOffsets(); // non-zero positions
	const int matei = (read.mate <= 1 ? 0 : 1);
	satups_.clear();
	gws_.clear();
	rands_.clear();
	rands2_.clear();
	satpos_.clear();
	satpos2_.clear();
	index_t nrange = 0, nelt = 0, nsmall = 0, nsmall_elts = 0;
	bool keepWhole = false;
	EList<SATupleAndPos<index_t>, 16>& satpos = keepWhole ? satpos_ : satpos2_;
	for(index_t i = 0; i < nonz; i++) {
		bool fw = true;
		index_t offidx = 0, rdoff = 0, seedlen = 0;
		QVal<index_t> qv = sh.hitsByRank(i, offidx, rdoff, fw, seedlen);
		assert(qv.valid());
		assert(!qv.empty());
		assert(qv.repOk(ca.current()));
		ca.queryQval(qv, satups_, nrange, nelt);
		for(size_t j = 0; j < satups_.size(); j++) {
			const index_t sz = satups_[j].size();
			// Check whether this hit occurs inside the extended boundaries of
			// another hit we already processed for this read.
			if(seedmms == 0) {
				// See if we're covered by a previous extended seed hit
				EList<ExtendRange>& range =
				fw ? seedExRangeFw_[matei] : seedExRangeRc_[matei];
				bool skip = false;
				for(index_t k = 0; k < range.size(); k++) {
					index_t p5 = range[k].off;
					index_t len = range[k].len;
					if(p5 <= rdoff && p5 + len >= (rdoff + seedlen)) {
						if(sz <= range[k].sz) {
							skip = true;
							break;
						}
					}
				}
				if(skip) {
					assert_gt(nrange, 0);
					nrange--;
					assert_geq(nelt, sz);
					nelt -= sz;
					continue; // Skip this seed
				}
			}
			satpos.expand();
			satpos.back().sat = satups_[j];
			satpos.back().origSz = sz;
			satpos.back().pos.init(fw, offidx, rdoff, seedlen);
			if(sz <= nsm) {
				nsmall++;
				nsmall_elts += sz;
			}
			satpos.back().nlex = satpos.back().nrex = 0;
#ifndef NDEBUG
			tmp_rdseq_.clear();
			uint64_t key = satpos.back().sat.key.seq;
			for(size_t k = 0; k < seedlen; k++) {
				int c = (int)(key & 3);
				tmp_rdseq_.append(c);
				key >>= 2;
			}
			tmp_rdseq_.reverse();
#endif
			index_t nlex = 0, nrex = 0;
			if(doExtend) {
				extend(
					   read,
					   gfmFw,
					   gfmBw,
					   satpos.back().sat.topf,
					   (index_t)(satpos.back().sat.topf + sz),
					   satpos.back().sat.topb,
					   (index_t)(satpos.back().sat.topb + sz),
					   fw,
					   rdoff,
					   seedlen,
					   prm,
					   nlex,
					   nrex);
			}
			satpos.back().nlex = nlex;
			satpos.back().nrex = nrex;
			if(seedmms == 0 && (nlex > 0 || nrex > 0)) {
				assert_geq(rdoff, (fw ? nlex : nrex));
				index_t p5 = rdoff - (fw ? nlex : nrex);
				EList<ExtendRange>& range =
				fw ? seedExRangeFw_[matei] : seedExRangeRc_[matei];
				range.expand();
				range.back().off = p5;
				range.back().len = seedlen + nlex + nrex;
				range.back().sz = sz;
			}
		}
		satups_.clear();
	}
	assert_leq(nsmall, nrange);
	nelt_out = nelt; // return the total number of elements
	assert_eq(nrange, satpos.size());
	satpos.sort();
	if(keepWhole) {
		gws_.ensure(nrange);
		rands_.ensure(nrange);
		for(index_t i = 0; i < nrange; i++) {
			gws_.expand();
			SARangeWithOffs<TSlice, index_t> sa;
			sa.topf = satpos_.back().sat.topf;
			sa.len = satpos_.back().sat.key.len;
			sa.offs = satpos_.back().sat.offs;
			gws_.back().init(
							 gfmFw,  // forward Bowtie index
							 ref,    // reference sequences
							 sa,     // SA tuples: ref hit, salist range
							 rnd,    // pseudo-random generator
							 wlm);   // metrics
			assert(gws_.back().initialized());
			rands_.expand();
			rands_.back().init(satpos_[i].sat.size(), all);
		}
		return;
	}
	// Resize satups_ list so that ranges having elements that we might
	// possibly explore are present
	satpos_.ensure(min(maxelt, nelt));
	gws_.ensure(min(maxelt, nelt));
	rands_.ensure(min(maxelt, nelt));
	rands2_.ensure(min(maxelt, nelt));
	size_t nlarge_elts = nelt - nsmall_elts;
	if(maxelt < nelt) {
		size_t diff = nelt - maxelt;
		if(diff >= nlarge_elts) {
			nlarge_elts = 0;
		} else {
			nlarge_elts -= diff;
		}
	}
	index_t nelt_added = 0;
	// Now we have a collection of ranges in satpos2_.  Now we want to decide
	// how we explore elements from them.  The basic idea is that: for very
	// small guys, where "very small" means that the size of the range is less
	// than or equal to the parameter 'nsz', we explore them in their entirety
	// right away.  For the rest, we want to select in a way that is (a)
	// random, and (b) weighted toward examining elements from the smaller
	// ranges more frequently (and first).
	//
	// 1. do the smalls
	for(index_t j = 0; j < nsmall && nelt_added < maxelt; j++) {
		satpos_.expand();
		satpos_.back() = satpos2_[j];
		gws_.expand();
		SARangeWithOffs<TSlice, index_t> sa;
		sa.topf = satpos_.back().sat.topf;
		sa.len = satpos_.back().sat.key.len;
		sa.offs = satpos_.back().sat.offs;
		gws_.back().init(
						 gfmFw,  // forward Bowtie index
						 ref,    // reference sequences
						 sa,     // SA tuples: ref hit, salist range
						 rnd,    // pseudo-random generator
						 wlm);   // metrics
		assert(gws_.back().initialized());
		rands_.expand();
		rands_.back().init(satpos_.back().sat.size(), all);
		nelt_added += satpos_.back().sat.size();
#ifndef NDEBUG
		for(size_t k = 0; k < satpos_.size()-1; k++) {
			assert(!(satpos_[k] == satpos_.back()));
		}
#endif
	}
	if(nelt_added >= maxelt || nsmall == satpos2_.size()) {
		nelt_out = nelt_added;
		return;
	}
	// 2. do the non-smalls
	// Initialize the row sampler
	rowsamp_.init(satpos2_, nsmall, satpos2_.size(), lensq, szsq);
	// Initialize the random choosers
	rands2_.resize(satpos2_.size());
	for(index_t j = 0; j < satpos2_.size(); j++) {
		rands2_[j].reset();
	}
	while(nelt_added < maxelt && nelt_added < nelt) {
		// Pick a non-small range to sample from
		index_t ri = rowsamp_.next(rnd) + nsmall;
		assert_geq(ri, nsmall);
		assert_lt(ri, satpos2_.size());
		// Initialize random element chooser for that range
		if(!rands2_[ri].inited()) {
			rands2_[ri].init(satpos2_[ri].sat.size(), all);
			assert(!rands2_[ri].done());
		}
		assert(!rands2_[ri].done());
		// Choose an element from the range
		uint32_t r = rands2_[ri].next(rnd);
		if(rands2_[ri].done()) {
			// Tell the row sampler this range is done
			rowsamp_.finishedRange(ri - nsmall);
		}
		// Add the element to the satpos_ list
		SATuple<index_t> sat;
		TSlice o;
		o.init(satpos2_[ri].sat.offs, r, r+1);
		sat.init(satpos2_[ri].sat.key, (index_t)(satpos2_[ri].sat.topf + r), (index_t)OFF_MASK, o);
		satpos_.expand();
		satpos_.back().sat = sat;
		satpos_.back().origSz = satpos2_[ri].origSz;
		satpos_.back().pos = satpos2_[ri].pos;
		// Initialize GroupWalk object
		gws_.expand();
		SARangeWithOffs<TSlice, index_t> sa;
		sa.topf = sat.topf;
		sa.len = sat.key.len;
		sa.offs = sat.offs;
		gws_.back().init(
						 gfmFw,  // forward Bowtie index
						 ref,    // reference sequences
						 sa,     // SA tuples: ref hit, salist range
						 rnd,    // pseudo-random generator
						 wlm);   // metrics
		assert(gws_.back().initialized());
		// Initialize random selector
		rands_.expand();
		rands_.back().init(1, all);
		nelt_added++;
	}
	nelt_out = nelt_added;
	return;
}

enum {
	FOUND_NONE = 0,
	FOUND_EE,
	FOUND_UNGAPPED,
};

/**
 * Given a collection of SeedHits for a single read, extend seed alignments
 * into full alignments.  Where possible, try to avoid redundant offset lookups
 * and dynamic programming wherever possible.  Optionally report alignments to
 * a AlnSinkWrap object as they are discovered.
 *
 * If 'reportImmediately' is true, returns true iff a call to msink->report()
 * returned true (indicating that the reporting policy is satisfied and we can
 * stop).  Otherwise, returns false.
 */
template <typename index_t>
int SwDriver<index_t>::extendSeeds(
								   Read& rd,                    // read to align
								   bool mate1,                  // true iff rd is mate #1
								   SeedResults<index_t>& sh,    // seed hits to extend into full alignments
								   const GFM<index_t>& gfmFw,   // BWT
								   const GFM<index_t>* gfmBw,   // BWT'
								   const BitPairReference& ref, // Reference strings
								   SwAligner& swa,              // dynamic programming aligner
								   const Scoring& sc,           // scoring scheme
								   int seedmms,                 // # mismatches allowed in seed
								   int seedlen,                 // length of seed
								   int seedival,                // interval between seeds
								   TAlScore& minsc,             // minimum score for anchor
								   int nceil,                   // maximum # Ns permitted in reference portion
								   size_t maxhalf,  	         // max width in either direction for DP tables
								   bool doUngapped,             // do ungapped alignment
								   size_t maxIters,             // stop after this many seed-extend loop iters
								   size_t maxUg,                // stop after this many ungaps
								   size_t maxDp,                // stop after this many dps
								   size_t maxUgStreak,          // stop after streak of this many ungap fails
								   size_t maxDpStreak,          // stop after streak of this many dp fails
								   bool doExtend,               // do seed extension
								   bool enable8,                // use 8-bit SSE where possible
								   size_t cminlen,              // use checkpointer if read longer than this
								   size_t cpow2,                // interval between diagonals to checkpoint
								   bool doTri,                  // triangular mini-fills?
								   int tighten,                 // -M score tightening mode
								   AlignmentCacheIface<index_t>& ca,     // alignment cache for seed hits
								   RandomSource& rnd,           // pseudo-random source
								   WalkMetrics& wlm,            // group walk left metrics
								   SwMetrics& swmSeed,          // DP metrics for seed-extend
								   PerReadMetrics& prm,         // per-read metrics
								   AlnSinkWrap<index_t>* msink, // AlnSink wrapper for multiseed-style aligner
								   bool reportImmediately,      // whether to report hits immediately to msink
								   bool& exhaustive)            // set to true iff we searched all seeds exhaustively
{
	bool all = msink->allHits();
	// typedef std::pair<index_t, index_t> UPair;
	
	assert(!reportImmediately || msink != NULL);
	assert(!reportImmediately || !msink->maxed());
	
	assert_geq(nceil, 0);
	assert_leq((size_t)nceil, rd.length());
	
	// Calculate the largest possible number of read and reference gaps
	const index_t rdlen = (index_t)rd.length();
	TAlScore perfectScore = sc.perfectScore(rdlen);
	
	DynProgFramer dpframe(!gReportOverhangs);
	swa.reset();
	
	// Initialize a set of GroupWalks, one for each seed.  Also, intialize the
	// accompanying lists of reference seed hits (satups*)
	const index_t nsm = 5;
	const index_t nonz = sh.nonzeroOffsets(); // non-zero positions
	index_t eeHits = sh.numE2eHits();
	bool eeMode = eeHits > 0;
	bool firstEe = true;
	bool firstExtend = true;
	
	// Reset all the counters related to streaks
	prm.nEeFail = 0;
	prm.nUgFail = 0;
	prm.nDpFail = 0;
	
	index_t nelt = 0, neltLeft = 0;
	index_t rows = rdlen;
	index_t eltsDone = 0;
	// cerr << "===" << endl;
	while(true) {
		if(eeMode) {
			if(firstEe) {
				firstEe = false;
				eeMode = eeSaTups(
								  rd,           // read
								  sh,           // seed hits to extend into full alignments
								  gfmFw,        // BWT
								  ref,          // Reference strings
								  rnd,          // pseudo-random generator
								  wlm,          // group walk left metrics
								  swmSeed,      // seed-extend metrics
								  nelt,         // out: # elements total
								  maxIters,     // max # to report
								  all);         // report all hits?
				assert_eq(gws_.size(), rands_.size());
				assert_eq(gws_.size(), satpos_.size());
			} else {
				eeMode = false;
			}
		}
		if(!eeMode) {
			if(nonz == 0) {
				return EXTEND_EXHAUSTED_CANDIDATES; // No seed hits!  Bail.
			}
			if(minsc == perfectScore) {
				return EXTEND_PERFECT_SCORE; // Already found all perfect hits!
			}
			if(firstExtend) {
				nelt = 0;
				prioritizeSATups(
								 rd,            // read
								 sh,            // seed hits to extend into full alignments
								 gfmFw,         // BWT
								 gfmBw,         // BWT'
								 ref,           // Reference strings
								 seedmms,       // # seed mismatches allowed
								 maxIters,      // max rows to consider per position
								 doExtend,      // extend out seeds
								 true,          // square extended length
								 true,          // square SA range size
								 nsm,           // smallness threshold
								 ca,            // alignment cache for seed hits
								 rnd,           // pseudo-random generator
								 wlm,           // group walk left metrics
								 prm,           // per-read metrics
								 nelt,          // out: # elements total
								 all);          // report all hits?
				assert_eq(gws_.size(), rands_.size());
				assert_eq(gws_.size(), satpos_.size());
				neltLeft = nelt;
				firstExtend = false;
			}
			if(neltLeft == 0) {
				// Finished examining gapped candidates
				break;
			}
		}
		for(size_t i = 0; i < gws_.size(); i++) {
			if(eeMode && eehits_[i].score < minsc) {
				return EXTEND_PERFECT_SCORE;
			}
			bool is_small      = satpos_[i].sat.size() < nsm;
			bool fw            = satpos_[i].pos.fw;
			index_t rdoff      = satpos_[i].pos.rdoff;
			index_t seedhitlen = satpos_[i].pos.seedlen;
			if(!fw) {
				// 'rdoff' and 'offidx' are with respect to the 5' end of
				// the read.  Here we convert rdoff to be with respect to
				// the upstream (3') end of ther read.
				rdoff = (index_t)(rdlen - rdoff - seedhitlen);
			}
			bool first = true;
			// If the range is small, investigate all elements now.  If the
			// range is large, just investigate one and move on - we might come
			// back to this range later.
			index_t riter = 0;
			while(!rands_[i].done() && (first || is_small || eeMode)) {
				assert(!gws_[i].done());
				riter++;
				if(minsc == perfectScore) {
					if(!eeMode || eehits_[i].score < perfectScore) {
						return EXTEND_PERFECT_SCORE;
					}
				} else if(eeMode && eehits_[i].score < minsc) {
					break;
				}
				if(prm.nExDps >= maxDp || prm.nMateDps >= maxDp) {
					return EXTEND_EXCEEDED_HARD_LIMIT;
				}
				if(prm.nExUgs >= maxUg || prm.nMateUgs >= maxUg) {
					return EXTEND_EXCEEDED_HARD_LIMIT;
				}
				if(prm.nExIters >= maxIters) {
					return EXTEND_EXCEEDED_HARD_LIMIT;
				}
				prm.nExIters++;
				first = false;
				// Resolve next element offset
				WalkResult<index_t> wr;
				uint32_t elt = rands_[i].next(rnd);
				//cerr << "elt=" << elt << endl;
				SARangeWithOffs<TSlice, index_t> sa;
				sa.topf = satpos_[i].sat.topf;
				sa.len = satpos_[i].sat.key.len;
				sa.offs = satpos_[i].sat.offs;
				gws_[i].advanceElement((index_t)elt, gfmFw, ref, sa, gwstate_, wr, wlm, prm);
				eltsDone++;
				if(!eeMode) {
					assert_gt(neltLeft, 0);
					neltLeft--;
				}
				assert_neq((index_t)OFF_MASK, wr.toff);
				index_t tidx = 0, toff = 0, tlen = 0;
				bool straddled = false;
				gfmFw.joinedToTextOff(
                                      wr.elt.len,
                                      wr.toff,
                                      tidx,
                                      toff,
                                      tlen,
                                      eeMode,     // reject straddlers?
                                      straddled); // did it straddle?
				if(tidx == (index_t)OFF_MASK) {
					// The seed hit straddled a reference boundary so the seed hit
					// isn't valid
					continue;
				}
#ifndef NDEBUG
				if(!eeMode && !straddled) { // Check that seed hit matches reference
					uint64_t key = satpos_[i].sat.key.seq;
					for(index_t k = 0; k < wr.elt.len; k++) {
						int c = ref.getBase(tidx, toff + wr.elt.len - k - 1);
						assert_leq(c, 3);
						int ck = (int)(key & 3);
						key >>= 2;
						assert_eq(c, ck);
					}
				}
#endif
				// Find offset of alignment's upstream base assuming net gaps=0
				// between beginning of read and beginning of seed hit
				int64_t refoff = (int64_t)toff - rdoff;
				// Coordinate of the seed hit w/r/t the pasted reference string
				Coord refcoord(tidx, refoff, fw);
				if(seenDiags1_.locusPresent(refcoord)) {
					// Already handled alignments seeded on this diagonal
					prm.nRedundants++;
					swmSeed.rshit++;
					continue;
				}
				// Now that we have a seed hit, there are many issues to solve
				// before we have a completely framed dynamic programming problem.
				// They include:
				//
				// 1. Setting reference offsets on either side of the seed hit,
				//    accounting for where the seed occurs in the read
				// 2. Adjusting the width of the banded dynamic programming problem
				//    and adjusting reference bounds to allow for gaps in the
				//    alignment
				// 3. Accounting for the edges of the reference, which can impact
				//    the width of the DP problem and reference bounds.
				// 4. Perhaps filtering the problem down to a smaller problem based
				//    on what DPs we've already solved for this read
				//
				// We do #1 here, since it is simple and we have all the seed-hit
				// information here.  #2 and #3 are handled in the DynProgFramer.
				int readGaps = 0, refGaps = 0;
				bool ungapped = false;
				if(!eeMode) {
					readGaps = sc.maxReadGaps(minsc, rdlen);
					refGaps  = sc.maxRefGaps(minsc, rdlen);
					ungapped = (readGaps == 0 && refGaps == 0);
				}
				int state = FOUND_NONE;
				bool found = false;
				if(eeMode) {
					resEe_.reset();
					resEe_.alres.reset();
					const EEHit<index_t>& h = eehits_[i];
					assert_leq(h.score, perfectScore);
					resEe_.alres.setScore(AlnScore(h.score, h.ns(), 0));
					resEe_.alres.setShape(
										  refcoord.ref(),  // ref id
										  refcoord.off(),  // 0-based ref offset
										  tlen,            // length of reference
										  fw,              // aligned to Watson?
										  rdlen,           // read length
										  true,            // pretrim soft?
										  0,               // pretrim 5' end
										  0,               // pretrim 3' end
										  true,            // alignment trim soft?
										  0,               // alignment trim 5' end
										  0);              // alignment trim 3' end
					resEe_.alres.setRefNs(h.refns());
					if(h.mms() > 0) {
						assert_eq(1, h.mms());
						assert_lt(h.e1.pos, rd.length());
						resEe_.alres.ned().push_back(h.e1);
					}
					assert(resEe_.repOk(rd));
					state = FOUND_EE;
					found = true;
					Interval refival(refcoord, 1);
					seenDiags1_.add(refival);
				} else if(doUngapped && ungapped) {
					resUngap_.reset();
					int al = swa.ungappedAlign(
											   fw ? rd.patFw : rd.patRc,
											   fw ? rd.qual  : rd.qualRev,
											   refcoord,
											   ref,
											   tlen,
											   sc,
											   gReportOverhangs,
											   minsc,
											   resUngap_);
					Interval refival(refcoord, 1);
					seenDiags1_.add(refival);
					prm.nExUgs++;
					if(al == 0) {
						prm.nExUgFails++;
						prm.nUgFail++;
						if(prm.nUgFail >= maxUgStreak) {
							return EXTEND_EXCEEDED_SOFT_LIMIT;
						}
						swmSeed.ungapfail++;
						continue;
					} else if(al == -1) {
						prm.nExUgFails++;
						prm.nUgFail++; // count this as failure
						if(prm.nUgFail >= maxUgStreak) {
							return EXTEND_EXCEEDED_SOFT_LIMIT;
						}
						swmSeed.ungapnodec++;
					} else {
						prm.nExUgSuccs++;
						prm.nUgLastSucc = prm.nExUgs-1;
						if(prm.nUgFail > prm.nUgFailStreak) {
							prm.nUgFailStreak = prm.nUgFail;
						}
						prm.nUgFail = 0;
						found = true;
						state = FOUND_UNGAPPED;
						swmSeed.ungapsucc++;
					}
				}
				int64_t pastedRefoff = (int64_t)wr.toff - rdoff;
				DPRect rect;
				if(state == FOUND_NONE) {
					found = dpframe.frameSeedExtensionRect(
														   refoff,   // ref offset implied by seed hit assuming no gaps
														   rows,     // length of read sequence used in DP table
														   tlen,     // length of reference
														   readGaps, // max # of read gaps permitted in opp mate alignment
														   refGaps,  // max # of ref gaps permitted in opp mate alignment
														   (size_t)nceil, // # Ns permitted
														   maxhalf,  // max width in either direction
														   rect);    // DP rectangle
					assert(rect.repOk());
					// Add the seed diagonal at least
					seenDiags1_.add(Interval(refcoord, 1));
					if(!found) {
						continue;
					}
				}
				int64_t leftShift = refoff - rect.refl;
				size_t nwindow = 0;
				if(toff >= rect.refl) {
					nwindow = (size_t)(toff - rect.refl);
				}
				// NOTE: We might be taking off more than we should because the
				// pasted string omits non-A/C/G/T characters, but we included them
				// when calculating leftShift.  We'll account for this later.
				pastedRefoff -= leftShift;
				size_t nsInLeftShift = 0;
				if(state == FOUND_NONE) {
					if(!swa.initedRead()) {
						// Initialize the aligner with a new read
						swa.initRead(
									 rd.patFw,  // fw version of query
									 rd.patRc,  // rc version of query
									 rd.qual,   // fw version of qualities
									 rd.qualRev,// rc version of qualities
									 0,         // off of first char in 'rd' to consider
									 rdlen,     // off of last char (excl) in 'rd' to consider
									 sc);       // scoring scheme
					}
					swa.initRef(
								fw,        // whether to align forward or revcomp read
								tidx,      // reference aligned against
								rect,      // DP rectangle
								ref,       // Reference strings
								tlen,      // length of reference sequence
								sc,        // scoring scheme
								minsc,     // minimum score permitted
								enable8,   // use 8-bit SSE if possible?
								cminlen,   // minimum length for using checkpointing scheme
								cpow2,     // interval b/t checkpointed diags; 1 << this
								doTri,     // triangular mini-fills?
								true,      // this is a seed extension - not finding a mate
								nwindow,
								nsInLeftShift);
					// Because of how we framed the problem, we can say that we've
					// exhaustively scored the seed diagonal as well as maxgaps
					// diagonals on either side
					Interval refival(tidx, 0, fw, 0);
					rect.initIval(refival);
					seenDiags1_.add(refival);
					// Now fill the dynamic programming matrix and return true iff
					// there is at least one valid alignment
					TAlScore bestCell = std::numeric_limits<TAlScore>::min();
					found = swa.align(rnd, bestCell);
					swmSeed.tallyGappedDp(readGaps, refGaps);
					prm.nExDps++;
					if(!found) {
						prm.nExDpFails++;
						prm.nDpFail++;
						if(prm.nDpFail >= maxDpStreak) {
							return EXTEND_EXCEEDED_SOFT_LIMIT;
						}
						if(bestCell > std::numeric_limits<TAlScore>::min() && bestCell > prm.bestLtMinscMate1) {
							prm.bestLtMinscMate1 = bestCell;
						}
						continue; // Look for more anchor alignments
					} else {
						prm.nExDpSuccs++;
						prm.nDpLastSucc = prm.nExDps-1;
						if(prm.nDpFail > prm.nDpFailStreak) {
							prm.nDpFailStreak = prm.nDpFail;
						}
						prm.nDpFail = 0;
					}
				}
				bool firstInner = true;
				while(true) {
					assert(found);
					SwResult *res = NULL;
					if(state == FOUND_EE) {
						if(!firstInner) {
							break;
						}
						res = &resEe_;
					} else if(state == FOUND_UNGAPPED) {
						if(!firstInner) {
							break;
						}
						res = &resUngap_;
					} else {
						resGap_.reset();
						assert(resGap_.empty());
						if(swa.done()) {
							break;
						}
						swa.nextAlignment(resGap_, minsc, rnd);
						found = !resGap_.empty();
						if(!found) {
							break;
						}
						res = &resGap_;
					}
					assert(res != NULL);
					firstInner = false;
					assert(res->alres.matchesRef(
												 rd,
												 ref,
												 tmp_rf_,
												 tmp_rdseq_,
												 tmp_qseq_,
												 raw_refbuf_,
												 raw_destU32_,
												 raw_matches_,
												 tmp_reflens_,
												 tmp_refoffs_));
					Interval refival(tidx, 0, fw, tlen);
					assert_gt(res->alres.refExtent(), 0);
					if(gReportOverhangs &&
					   !refival.containsIgnoreOrient(res->alres.refival()))
					{
						res->alres.clipOutside(true, 0, tlen);
						if(res->alres.refExtent() == 0) {
							continue;
						}
					}
					assert(gReportOverhangs ||
					       refival.containsIgnoreOrient(res->alres.refival()));
					// Did the alignment fall entirely outside the reference?
					if(!refival.overlapsIgnoreOrient(res->alres.refival())) {
						continue;
					}
					// Is this alignment redundant with one we've seen previously?
					if(redAnchor_.overlap(res->alres)) {
						// Redundant with an alignment we found already
						continue;
					}
					redAnchor_.add(res->alres);
					// Annotate the AlnRes object with some key parameters
					// that were used to obtain the alignment.
					res->alres.setParams(
										 seedmms,   // # mismatches allowed in seed
										 seedlen,   // length of seed
										 seedival,  // interval between seeds
										 minsc);    // minimum score for valid alignment
					
					if(reportImmediately) {
						assert(msink != NULL);
						assert(res->repOk());
						// Check that alignment accurately reflects the
						// reference characters aligned to
						assert(res->alres.matchesRef(
													 rd,
													 ref,
													 tmp_rf_,
													 tmp_rdseq_,
													 tmp_qseq_,
													 raw_refbuf_,
													 raw_destU32_,
													 raw_matches_,
													 tmp_reflens_,
													 tmp_refoffs_));
						// Report an unpaired alignment
						assert(!msink->maxed());
						if(msink->report(
										 0,
										 mate1 ? &res->alres : NULL,
										 mate1 ? NULL : &res->alres))
						{
							// Short-circuited because a limit, e.g. -k, -m or
							// -M, was exceeded
							return EXTEND_POLICY_FULFILLED;
						}
						if(tighten > 0 &&
						   msink->Mmode() &&
						   msink->hasSecondBestUnp1())
						{
							if(tighten == 1) {
								if(msink->bestUnp1() >= minsc) {
									minsc = msink->bestUnp1();
									if(minsc < perfectScore &&
									   msink->bestUnp1() == msink->secondBestUnp1())
									{
										minsc++;
									}
								}
							} else if(tighten == 2) {
								if(msink->secondBestUnp1() >= minsc) {
									minsc = msink->secondBestUnp1();
									if(minsc < perfectScore) {
										minsc++;
									}
								}
							} else {
								TAlScore diff = msink->bestUnp1() - msink->secondBestUnp1();
								TAlScore bot = msink->secondBestUnp1() + ((diff*3)/4);
								if(bot >= minsc) {
									minsc = bot;
									if(minsc < perfectScore) {
										minsc++;
									}
								}
							}
							assert_leq(minsc, perfectScore);
						}
					}
				}
				
				// At this point we know that we aren't bailing, and will
				// continue to resolve seed hits.  
				
			} // while(!gws_[i].done())
		}
	}
	// Short-circuited because a limit, e.g. -k, -m or -M, was exceeded
	return EXTEND_EXHAUSTED_CANDIDATES;
}

/**
 * Given a collection of SeedHits for both mates in a read pair, extend seed
 * alignments into full alignments and then look for the opposite mate using
 * dynamic programming.  Where possible, try to avoid redundant offset lookups.
 * Optionally report alignments to a AlnSinkWrap object as they are discovered.
 *
 * If 'reportImmediately' is true, returns true iff a call to
 * msink->report() returned true (indicating that the reporting
 * policy is satisfied and we can stop).  Otherwise, returns false.
 *
 * REDUNDANT SEED HITS
 *
 * See notes at top of aligner_sw_driver.h.
 *
 * REDUNDANT ALIGNMENTS
 *
 * See notes at top of aligner_sw_driver.h.
 *
 * MIXING PAIRED AND UNPAIRED ALIGNMENTS
 *
 * There are distinct paired-end alignment modes for the cases where (a) the
 * user does or does not want to see unpaired alignments for individual mates
 * when there are no reportable paired-end alignments involving both mates, and
 * (b) the user does or does not want to see discordant paired-end alignments.
 * The modes have implications for this function and for the AlnSinkWrap, since
 * it affects when we're "done."  Also, whether the user has asked us to report
 * discordant alignments affects whether and how much searching for unpaired
 * alignments we must do (i.e. if there are no paired-end alignments, we must
 * at least do -m 1 for both mates).
 *
 * Mode 1: Just concordant paired-end.  Print only concordant paired-end
 * alignments.  As soon as any limits (-k/-m/-M) are reached, stop.
 *
 * Mode 2: Concordant and discordant paired-end.  If -k/-m/-M limits are
 * reached for paired-end alignments, stop.  Otherwise, if no paired-end
 * alignments are found, align both mates in an unpaired -m 1 fashion.  If
 * there is exactly one unpaired alignment for each mate, report the
 * combination as a discordant alignment.
 *
 * Mode 3: Concordant paired-end if possible, otherwise unpaired.  If -k/-M
 * limit is reached for paired-end alignmnts, stop.  If -m limit is reached for
 * paired-end alignments or no paired-end alignments are found, align both
 * mates in an unpaired fashion.  All the same settings governing validity and
 * reportability in paired-end mode apply here too (-k/-m/-M/etc).
 *
 * Mode 4: Concordant or discordant paired-end if possible, otherwise unpaired.
 * If -k/-M limit is reached for paired-end alignmnts, stop.  If -m limit is
 * reached for paired-end alignments or no paired-end alignments are found,
 * align both mates in an unpaired fashion.  If the -m limit was reached, there
 * is no need to search for a discordant alignment, and unapired alignment can
 * proceed as in Mode 3.  If no paired-end alignments were found, then unpaired
 * alignment proceeds as in Mode 3 but with this caveat: alignment must be at
 * least as thorough as dictated by -m 1 up until the point where
 *
 * Print paired-end alignments when there are reportable paired-end
 * alignments, otherwise report reportable unpaired alignments.  If -k limit is
 * reached for paired-end alignments, stop.  If -m/-M limit is reached for
 * paired-end alignments, stop searching for paired-end alignments and look
 * only for unpaired alignments.  If searching only for unpaired alignments,
 * respect -k/-m/-M limits separately for both mates.
 *
 * The return value from the AlnSinkWrap's report member function must be
 * specific enough to distinguish between:
 *
 * 1. Stop searching for paired-end alignments
 * 2. Stop searching for alignments for unpaired alignments for mate #1
 * 3. Stop searching for alignments for unpaired alignments for mate #2
 * 4. Stop searching for any alignments
 *
 * Note that in Mode 2, options affecting validity and reportability of
 * alignments apply .  E.g. if -m 1 is specified
 *
 * WORKFLOW
 *
 * Our general approach to finding paired and unpaired alignments here
 * is as follows:
 *
 * - For mate in mate1, mate2:
 *   - For each seed hit in mate:
 *     - Try to extend it into a full alignment; if we can't, continue
 *       to the next seed hit
 *     - Look for alignment for opposite mate; if we can't find one,
 *     - 
 *     - 
 *
 */
template <typename index_t>
int SwDriver<index_t>::extendSeedsPaired(
										 Read& rd,                    // mate to align as anchor
										 Read& ord,                   // mate to align as opposite
										 bool anchor1,                // true iff anchor mate is mate1
										 bool oppFilt,                // true iff opposite mate was filtered out
										 SeedResults<index_t>& sh,    // seed hits for anchor
										 const GFM<index_t>& gfmFw,   // BWT
										 const GFM<index_t>* gfmBw,   // BWT'
										 const BitPairReference& ref, // Reference strings
										 SwAligner& swa,              // dynamic programming aligner for anchor
										 SwAligner& oswa,             // dynamic programming aligner for opposite
										 const Scoring& sc,           // scoring scheme
										 const PairedEndPolicy& pepol,// paired-end policy
										 int seedmms,                 // # mismatches allowed in seed
										 int seedlen,                 // length of seed
										 int seedival,                // interval between seeds
										 TAlScore& minsc,             // minimum score for valid anchor aln
										 TAlScore& ominsc,            // minimum score for valid opposite aln
										 int nceil,                   // max # Ns permitted in ref for anchor
										 int onceil,                  // max # Ns permitted in ref for opposite
										 bool nofw,                   // don't align forward read
										 bool norc,                   // don't align revcomp read
										 size_t maxhalf,              // max width in either direction for DP tables
										 bool doUngapped,             // do ungapped alignment
										 size_t maxIters,             // stop after this many seed-extend loop iters
										 size_t maxUg,                // stop after this many ungaps
										 size_t maxDp,                // stop after this many dps
										 size_t maxEeStreak,          // stop after streak of this many end-to-end fails
										 size_t maxUgStreak,          // stop after streak of this many ungap fails
										 size_t maxDpStreak,          // stop after streak of this many dp fails
										 size_t maxMateStreak,        // stop seed range after N mate-find fails
										 bool doExtend,               // do seed extension
										 bool enable8,                // use 8-bit SSE where possible
										 size_t cminlen,              // use checkpointer if read longer than this
										 size_t cpow2,                // interval between diagonals to checkpoint
										 bool doTri,                  // triangular mini-fills?
										 int tighten,                 // -M score tightening mode
										 AlignmentCacheIface<index_t>& ca,     // alignment cache for seed hits
										 RandomSource& rnd,           // pseudo-random source
										 WalkMetrics& wlm,            // group walk left metrics
										 SwMetrics& swmSeed,          // DP metrics for seed-extend
										 SwMetrics& swmMate,          // DP metrics for mate finidng
										 PerReadMetrics& prm,         // per-read metrics
										 AlnSinkWrap<index_t>* msink, // AlnSink wrapper for multiseed-style aligner
										 bool swMateImmediately,      // whether to look for mate immediately
										 bool reportImmediately,      // whether to report hits immediately to msink
										 bool discord,                // look for discordant alignments?
										 bool mixed,                  // look for unpaired as well as paired alns?
										 bool& exhaustive)
{
	bool all = msink->allHits();
	// typedef std::pair<uint32_t, uint32_t> U32Pair;
	
	assert(!reportImmediately || msink != NULL);
	assert(!reportImmediately || !msink->maxed());
	assert(!msink->state().doneWithMate(anchor1));
	
	assert_geq(nceil, 0);
	assert_geq(onceil, 0);
	assert_leq((size_t)nceil,  rd.length());
	assert_leq((size_t)onceil, ord.length());
	
	const index_t rdlen  = rd.length();
	const index_t ordlen = ord.length();
	const TAlScore perfectScore = sc.perfectScore(rdlen);
	const TAlScore operfectScore = sc.perfectScore(ordlen);
	
	assert_leq(minsc, perfectScore);
	assert(oppFilt || ominsc <= operfectScore);
	
	TAlScore bestPairScore = perfectScore + operfectScore;
	if(tighten > 0 && msink->Mmode() && msink->hasSecondBestPair()) {
		// Paired-end alignments should have at least this score from now
		TAlScore ps;
		if(tighten == 1) {
			ps = msink->bestPair();
		} else if(tighten == 2) {
			ps = msink->secondBestPair();
		} else {
			TAlScore diff = msink->bestPair() - msink->secondBestPair();
			ps = msink->secondBestPair() + (diff * 3)/4;
		}
		if(tighten == 1 && ps < bestPairScore &&
		   msink->bestPair() == msink->secondBestPair())
		{
			ps++;
		}
		if(tighten >= 2 && ps < bestPairScore) {
			ps++;
		}
		// Anchor mate must have score at least 'ps' minus the best possible
		// score for the opposite mate.
		TAlScore nc = ps - operfectScore;
		if(nc > minsc) {
			minsc = nc;
		}
		assert_leq(minsc, perfectScore);
	}
	
	DynProgFramer dpframe(!gReportOverhangs);
	swa.reset();
	oswa.reset();
	
	// Initialize a set of GroupWalks, one for each seed.  Also, intialize the
	// accompanying lists of reference seed hits (satups*)
	const index_t nsm = 5;
	const index_t nonz = sh.nonzeroOffsets(); // non-zero positions
	index_t eeHits = sh.numE2eHits();
	bool eeMode = eeHits > 0;
	bool firstEe = true;
	bool firstExtend = true;
	
	// Reset all the counters related to streaks
	prm.nEeFail = 0;
	prm.nUgFail = 0;
	prm.nDpFail = 0;
	
	index_t nelt = 0, neltLeft = 0;
	const index_t rows = rdlen;
	const index_t orows  = ordlen;
	index_t eltsDone = 0;
	while(true) {
		if(eeMode) {
			if(firstEe) {
				firstEe = false;
				eeMode = eeSaTups(
								  rd,           // read
								  sh,           // seed hits to extend into full alignments
								  gfmFw,        // BWT
								  ref,          // Reference strings
								  rnd,          // pseudo-random generator
								  wlm,          // group walk left metrics
								  swmSeed,      // seed-extend metrics
								  nelt,         // out: # elements total
								  maxIters,     // max elts to report
								  all);         // report all hits
				assert_eq(gws_.size(), rands_.size());
				assert_eq(gws_.size(), satpos_.size());
				neltLeft = nelt;
				// Initialize list that contains the mate-finding failure
				// streak for each range
				mateStreaks_.resize(gws_.size());
				mateStreaks_.fill(0);
			} else {
				eeMode = false;
			}
		}
		if(!eeMode) {
			if(nonz == 0) {
				// No seed hits!  Bail.
				return EXTEND_EXHAUSTED_CANDIDATES;
			}
			if(msink->Mmode() && minsc == perfectScore) {
				// Already found all perfect hits!
				return EXTEND_PERFECT_SCORE;
			}
			if(firstExtend) {
				nelt = 0;
				prioritizeSATups(
								 rd,            // read
								 sh,            // seed hits to extend into full alignments
								 gfmFw,         // BWT
								 gfmBw,         // BWT'
								 ref,           // Reference strings
								 seedmms,       // # seed mismatches allowed
								 maxIters,      // max rows to consider per position
								 doExtend,      // extend out seeds
								 true,          // square extended length
								 true,          // square SA range size
								 nsm,           // smallness threshold
								 ca,            // alignment cache for seed hits
								 rnd,           // pseudo-random generator
								 wlm,           // group walk left metrics
								 prm,           // per-read metrics
								 nelt,          // out: # elements total
								 all);          // report all hits?
				assert_eq(gws_.size(), rands_.size());
				assert_eq(gws_.size(), satpos_.size());
				neltLeft = nelt;
				firstExtend = false;
				mateStreaks_.resize(gws_.size());
				mateStreaks_.fill(0);
			}
			if(neltLeft == 0) {
				// Finished examining gapped candidates
				break;
			}
		}
		for(index_t i = 0; i < gws_.size(); i++) {
			if(eeMode && eehits_[i].score < minsc) {
				return EXTEND_PERFECT_SCORE;
			}
			bool is_small      = satpos_[i].sat.size() < nsm;
			bool fw            = satpos_[i].pos.fw;
			index_t rdoff      = satpos_[i].pos.rdoff;
			index_t seedhitlen = satpos_[i].pos.seedlen;
			if(!fw) {
				// 'rdoff' and 'offidx' are with respect to the 5' end of
				// the read.  Here we convert rdoff to be with respect to
				// the upstream (3') end of ther read.
				rdoff = (index_t)(rdlen - rdoff - seedhitlen);
			}
			bool first = true;
			// If the range is small, investigate all elements now.  If the
			// range is large, just investigate one and move on - we might come
			// back to this range later.
			while(!rands_[i].done() && (first || is_small || eeMode)) {
				if(minsc == perfectScore) {
					if(!eeMode || eehits_[i].score < perfectScore) {
						return EXTEND_PERFECT_SCORE;
					}
				} else if(eeMode && eehits_[i].score < minsc) {
					break;
				}
				if(prm.nExDps >= maxDp || prm.nMateDps >= maxDp) {
					return EXTEND_EXCEEDED_HARD_LIMIT;
				}
				if(prm.nExUgs >= maxUg || prm.nMateUgs >= maxUg) {
					return EXTEND_EXCEEDED_HARD_LIMIT;
				}
				if(prm.nExIters >= maxIters) {
					return EXTEND_EXCEEDED_HARD_LIMIT;
				}
				if(eeMode && prm.nEeFail >= maxEeStreak) {
					return EXTEND_EXCEEDED_SOFT_LIMIT;
				}
				if(!eeMode && prm.nDpFail >= maxDpStreak) {
					return EXTEND_EXCEEDED_SOFT_LIMIT;
				}
				if(!eeMode && prm.nUgFail >= maxUgStreak) {
					return EXTEND_EXCEEDED_SOFT_LIMIT;
				}
				if(mateStreaks_[i] >= maxMateStreak) {
					// Don't try this seed range anymore
					rands_[i].setDone();
					assert(rands_[i].done());
					break;
				}
				prm.nExIters++;
				first = false;
				assert(!gws_[i].done());
				// Resolve next element offset
				WalkResult<index_t> wr;
				uint32_t elt = rands_[i].next(rnd);
				SARangeWithOffs<TSlice, index_t> sa;
				sa.topf = satpos_[i].sat.topf;
				sa.len = satpos_[i].sat.key.len;
				sa.offs = satpos_[i].sat.offs;
				gws_[i].advanceElement((index_t)elt, gfmFw, ref, sa, gwstate_, wr, wlm, prm);
				eltsDone++;
				assert_gt(neltLeft, 0);
				neltLeft--;
				assert_neq((index_t)OFF_MASK, wr.toff);
				index_t tidx = 0, toff = 0, tlen = 0;
				bool straddled = false;
				gfmFw.joinedToTextOff(
                                      wr.elt.len,
                                      wr.toff,
                                      tidx,
                                      toff,
                                      tlen,
                                      eeMode,       // reject straddlers?
                                      straddled);   // straddled?
				if(tidx == (index_t)OFF_MASK) {
					// The seed hit straddled a reference boundary so the seed hit
					// isn't valid
					continue;
				}
#ifndef NDEBUG
				if(!eeMode && !straddled) { // Check that seed hit matches reference
					uint64_t key = satpos_[i].sat.key.seq;
					for(index_t k = 0; k < wr.elt.len; k++) {
						int c = ref.getBase(tidx, toff + wr.elt.len - k - 1);
						assert_leq(c, 3);
						int ck = (int)(key & 3);
						key >>= 2;
						assert_eq(c, ck);
					}
				}
#endif
				// Find offset of alignment's upstream base assuming net gaps=0
				// between beginning of read and beginning of seed hit
				int64_t refoff = (int64_t)toff - rdoff;
				EIvalMergeListBinned& seenDiags  = anchor1 ? seenDiags1_ : seenDiags2_;
				// Coordinate of the seed hit w/r/t the pasted reference string
				Coord refcoord(tidx, refoff, fw);
				if(seenDiags.locusPresent(refcoord)) {
					// Already handled alignments seeded on this diagonal
					prm.nRedundants++;
					swmSeed.rshit++;
					continue;
				}
				// Now that we have a seed hit, there are many issues to solve
				// before we have a completely framed dynamic programming problem.
				// They include:
				//
				// 1. Setting reference offsets on either side of the seed hit,
				//    accounting for where the seed occurs in the read
				// 2. Adjusting the width of the banded dynamic programming problem
				//    and adjusting reference bounds to allow for gaps in the
				//    alignment
				// 3. Accounting for the edges of the reference, which can impact
				//    the width of the DP problem and reference bounds.
				// 4. Perhaps filtering the problem down to a smaller problem based
				//    on what DPs we've already solved for this read
				//
				// We do #1 here, since it is simple and we have all the seed-hit
				// information here.  #2 and #3 are handled in the DynProgFramer.
				int readGaps = 0, refGaps = 0;
				bool ungapped = false;
				if(!eeMode) {
					readGaps = sc.maxReadGaps(minsc, rdlen);
					refGaps  = sc.maxRefGaps(minsc, rdlen);
					ungapped = (readGaps == 0 && refGaps == 0);
				}
				int state = FOUND_NONE;
				bool found = false;
				// In unpaired mode, a seed extension is successful if it
				// results in a full alignment that meets the minimum score
				// threshold.  In paired-end mode, a seed extension is
				// successful if it results in a *full paired-end* alignment
				// that meets the minimum score threshold.
				if(eeMode) {
					resEe_.reset();
					resEe_.alres.reset();
					const EEHit<index_t>& h = eehits_[i];
					assert_leq(h.score, perfectScore);
					resEe_.alres.setScore(AlnScore(h.score, h.ns(), 0));
					resEe_.alres.setShape(
										  refcoord.ref(),  // ref id
										  refcoord.off(),  // 0-based ref offset
										  tlen,            // reference length
										  fw,              // aligned to Watson?
										  rdlen,           // read length
										  true,            // pretrim soft?
										  0,               // pretrim 5' end
										  0,               // pretrim 3' end
										  true,            // alignment trim soft?
										  0,               // alignment trim 5' end
										  0);              // alignment trim 3' end
					resEe_.alres.setRefNs(h.refns());
					if(h.mms() > 0) {
						assert_eq(1, h.mms());
						assert_lt(h.e1.pos, rd.length());
						resEe_.alres.ned().push_back(h.e1);
					}
					assert(resEe_.repOk(rd));
					state = FOUND_EE;
					found = true;
					Interval refival(refcoord, 1);
					seenDiags.add(refival);
					prm.nExEes++;
					prm.nEeFail++; // say it's failed until proven successful
					prm.nExEeFails++;
				} else if(doUngapped && ungapped) {
					resUngap_.reset();
					int al = swa.ungappedAlign(
											   fw ? rd.patFw : rd.patRc,
											   fw ? rd.qual  : rd.qualRev,
											   refcoord,
											   ref,
											   tlen,
											   sc,
											   gReportOverhangs,
											   minsc, // minimum
											   resUngap_);
					Interval refival(refcoord, 1);
					seenDiags.add(refival);
					prm.nExUgs++;
					prm.nUgFail++; // say it's failed until proven successful
					prm.nExUgFails++;
					if(al == 0) {
						swmSeed.ungapfail++;
						continue;
					} else if(al == -1) {
						swmSeed.ungapnodec++;
					} else {
						found = true;
						state = FOUND_UNGAPPED;
						swmSeed.ungapsucc++;
					}
				}
				int64_t pastedRefoff = (int64_t)wr.toff - rdoff;
				DPRect rect;
				if(state == FOUND_NONE) {
					found = dpframe.frameSeedExtensionRect(
														   refoff,   // ref offset implied by seed hit assuming no gaps
														   rows,     // length of read sequence used in DP table
														   tlen,     // length of reference
														   readGaps, // max # of read gaps permitted in opp mate alignment
														   refGaps,  // max # of ref gaps permitted in opp mate alignment
														   (size_t)nceil, // # Ns permitted
														   maxhalf,  // max width in either direction
														   rect);    // DP rectangle
					assert(rect.repOk());
					// Add the seed diagonal at least
					seenDiags.add(Interval(refcoord, 1));
					if(!found) {
						continue;
					}
				}
				int64_t leftShift = refoff - rect.refl;
				size_t nwindow = 0;
				if(toff >= rect.refl) {
					nwindow = (size_t)(toff - rect.refl);
				}
				// NOTE: We might be taking off more than we should because the
				// pasted string omits non-A/C/G/T characters, but we included them
				// when calculating leftShift.  We'll account for this later.
				pastedRefoff -= leftShift;
				size_t nsInLeftShift = 0;
				if(state == FOUND_NONE) {
					if(!swa.initedRead()) {
						// Initialize the aligner with a new read
						swa.initRead(
									 rd.patFw,  // fw version of query
									 rd.patRc,  // rc version of query
									 rd.qual,   // fw version of qualities
									 rd.qualRev,// rc version of qualities
									 0,         // off of first char in 'rd' to consider
									 rdlen,     // off of last char (excl) in 'rd' to consider
									 sc);       // scoring scheme
					}
					swa.initRef(
								fw,        // whether to align forward or revcomp read
								tidx,      // reference aligned against
								rect,      // DP rectangle
								ref,       // Reference strings
								tlen,      // length of reference sequence
								sc,        // scoring scheme
								minsc,     // minimum score permitted
								enable8,   // use 8-bit SSE if possible?
								cminlen,   // minimum length for using checkpointing scheme
								cpow2,     // interval b/t checkpointed diags; 1 << this
								doTri,     // triangular mini-fills?
								true,      // this is a seed extension - not finding a mate
								nwindow,
								nsInLeftShift);
					// Because of how we framed the problem, we can say that we've
					// exhaustively scored the seed diagonal as well as maxgaps
					// diagonals on either side
					Interval refival(tidx, 0, fw, 0);
					rect.initIval(refival);
					seenDiags.add(refival);
					// Now fill the dynamic programming matrix and return true iff
					// there is at least one valid alignment
					TAlScore bestCell = std::numeric_limits<TAlScore>::min();
					found = swa.align(rnd, bestCell);
					swmSeed.tallyGappedDp(readGaps, refGaps);
					prm.nExDps++;
					prm.nDpFail++;    // failed until proven successful
					prm.nExDpFails++; // failed until proven successful
					if(!found) {
						TAlScore bestLast = anchor1 ? prm.bestLtMinscMate1 : prm.bestLtMinscMate2;
						if(bestCell > std::numeric_limits<TAlScore>::min() && bestCell > bestLast) {
							if(anchor1) {
								prm.bestLtMinscMate1 = bestCell;
							} else {
								prm.bestLtMinscMate2 = bestCell;
							}
						}
						continue; // Look for more anchor alignments
					}
				}
				bool firstInner = true;
				bool foundConcordant = false;
				while(true) {
					assert(found);
					SwResult *res = NULL;
					if(state == FOUND_EE) {
						if(!firstInner) {
							break;
						}
						res = &resEe_;
						assert(res->repOk(rd));
					} else if(state == FOUND_UNGAPPED) {
						if(!firstInner) {
							break;
						}
						res = &resUngap_;
						assert(res->repOk(rd));
					} else {
						resGap_.reset();
						assert(resGap_.empty());
						if(swa.done()) {
							break;
						}
						swa.nextAlignment(resGap_, minsc, rnd);
						found = !resGap_.empty();
						if(!found) {
							break;
						}
						res = &resGap_;
						assert(res->repOk(rd));
					}
					// TODO: If we're just taking anchor alignments out of the
					// same rectangle, aren't we getting very similar
					// rectangles for the opposite mate each time?  Seems like
					// we could save some work by detecting this.
					assert(res != NULL);
					firstInner = false;
					assert(res->alres.matchesRef(
												 rd,
												 ref,
												 tmp_rf_,
												 tmp_rdseq_,
												 tmp_qseq_,
												 raw_refbuf_,
												 raw_destU32_,
												 raw_matches_,
												 tmp_reflens_,
												 tmp_refoffs_));
					Interval refival(tidx, 0, fw, tlen);
					assert_gt(res->alres.refExtent(), 0);
					if(gReportOverhangs &&
					   !refival.containsIgnoreOrient(res->alres.refival()))
					{
						res->alres.clipOutside(true, 0, tlen);
						if(res->alres.refExtent() == 0) {
							continue;
						}
					}
					assert(gReportOverhangs ||
					       refival.containsIgnoreOrient(res->alres.refival()));
					// Did the alignment fall entirely outside the reference?
					if(!refival.overlapsIgnoreOrient(res->alres.refival())) {
						continue;
					}
					// Is this alignment redundant with one we've seen previously?
					if(redAnchor_.overlap(res->alres)) {
						continue;
					}
					redAnchor_.add(res->alres);
					// Annotate the AlnRes object with some key parameters
					// that were used to obtain the alignment.
					res->alres.setParams(
										 seedmms,   // # mismatches allowed in seed
										 seedlen,   // length of seed
										 seedival,  // interval between seeds
										 minsc);    // minimum score for valid alignment
					bool foundMate = false;
					TRefOff off = res->alres.refoff();
					if( msink->state().doneWithMate(!anchor1) &&
					   !msink->state().doneWithMate( anchor1))
					{
						// We're done with the opposite mate but not with the
						// anchor mate; don't try to mate up the anchor.
						swMateImmediately = false;
					}
					if(found && swMateImmediately) {
						assert(!msink->state().doneWithMate(!anchor1));
						bool oleft = false, ofw = false;
						int64_t oll = 0, olr = 0, orl = 0, orr = 0;
						assert(!msink->state().done());
						foundMate = !oppFilt;
						TAlScore ominsc_cur = ominsc;
						//bool oungapped = false;
						int oreadGaps = 0, orefGaps = 0;
						//int oungappedAlign = -1; // defer
						if(foundMate) {
							// Adjust ominsc given the alignment score of the
							// anchor mate
							ominsc_cur = ominsc;
							if(tighten > 0 && msink->Mmode() && msink->hasSecondBestPair()) {
								// Paired-end alignments should have at least this score from now
								TAlScore ps;
								if(tighten == 1) {
									ps = msink->bestPair();
								} else if(tighten == 2) {
									ps = msink->secondBestPair();
								} else {
									TAlScore diff = msink->bestPair() - msink->secondBestPair();
									ps = msink->secondBestPair() + (diff * 3)/4;
								}
								if(tighten == 1 && ps < bestPairScore &&
								   msink->bestPair() == msink->secondBestPair())
								{
									ps++;
								}
								if(tighten >= 2 && ps < bestPairScore) {
									ps++;
								}
								// Anchor mate must have score at least 'ps' minus the best possible
								// score for the opposite mate.
								TAlScore nc = ps - res->alres.score().score();
								if(nc > ominsc_cur) {
									ominsc_cur = nc;
									assert_leq(ominsc_cur, operfectScore);
								}
							}
							oreadGaps = sc.maxReadGaps(ominsc_cur, ordlen);
							orefGaps  = sc.maxRefGaps (ominsc_cur, ordlen);
							//oungapped = (oreadGaps == 0 && orefGaps == 0);
							// TODO: Something lighter-weight than DP to scan
							// for other mate??
							//if(oungapped) {
							//	oresUngap_.reset();
							//	oungappedAlign = oswa.ungappedAlign(
							//		ofw ? ord.patFw : ord.patRc,
							//		ofw ? ord.qual  : ord.qualRev,
							//		orefcoord,
							//		ref,
							//		otlen,
							//		sc,
							//		gReportOverhangs,
							//		ominsc_cur,
							//		0,
							//		oresUngap_);
							//}
							foundMate = pepol.otherMate(
														anchor1,             // anchor mate is mate #1?
														fw,                  // anchor aligned to Watson?
														off,                 // offset of anchor mate
														orows + oreadGaps,   // max # columns spanned by alignment
														tlen,                // reference length
														anchor1 ? rd.length() : ord.length(), // mate 1 len
														anchor1 ? ord.length() : rd.length(), // mate 2 len
														oleft,               // out: look left for opposite mate?
														oll,
														olr,
														orl,
														orr,
														ofw);
						}
						DPRect orect;
						if(foundMate) {
							foundMate = dpframe.frameFindMateRect(
																  !oleft,      // true iff anchor alignment is to the left
																  oll,         // leftmost Watson off for LHS of opp aln
																  olr,         // rightmost Watson off for LHS of opp aln
																  orl,         // leftmost Watson off for RHS of opp aln
																  orr,         // rightmost Watson off for RHS of opp aln
																  orows,       // length of opposite mate
																  tlen,        // length of reference sequence aligned to
																  oreadGaps,   // max # of read gaps in opp mate aln
																  orefGaps,    // max # of ref gaps in opp mate aln
																  (size_t)onceil, // max # Ns on opp mate
																  maxhalf,     // max width in either direction
																  orect);      // DP rectangle
							assert(!foundMate || orect.refr >= orect.refl);
						}
						if(foundMate) {
							oresGap_.reset();
							assert(oresGap_.empty());
							if(!oswa.initedRead()) {
								oswa.initRead(
											  ord.patFw,  // read to align
											  ord.patRc,  // qualities
											  ord.qual,   // read to align
											  ord.qualRev,// qualities
											  0,          // off of first char to consider
											  ordlen,     // off of last char (ex) to consider
											  sc);        // scoring scheme
							}
							// Given the boundaries defined by refi and reff, initilize
							// the SwAligner with the dynamic programming problem that
							// aligns the read to this reference stretch.
							size_t onsInLeftShift = 0;
							assert_geq(orect.refr, orect.refl);
							oswa.initRef(
										 ofw,       // align forward or revcomp read?
										 tidx,      // reference aligned against
										 orect,     // DP rectangle
										 ref,       // Reference strings
										 tlen,      // length of reference sequence
										 sc,        // scoring scheme
										 ominsc_cur,// min score for valid alignments
										 enable8,   // use 8-bit SSE if possible?
										 cminlen,   // minimum length for using checkpointing scheme
										 cpow2,     // interval b/t checkpointed diags; 1 << this
										 doTri,     // triangular mini-fills?
										 false,     // this is finding a mate - not seed ext
										 0,         // nwindow?
										 onsInLeftShift);
							// TODO: Can't we add some diagonals to the
							// opposite mate's seenDiags when we fill in the
							// opposite mate's DP?  Or can we?  We might want
							// to use this again as an anchor - will that still
							// happen?  Also, isn't there a problem with
							// consistency of the minimum score?  Minimum score
							// here depends in part on the score of the anchor
							// alignment here, but it won't when the current
							// opposite becomes the anchor.
							
							// Because of how we framed the problem, we can say
							// that we've exhaustively explored the "core"
							// diagonals
							//Interval orefival(tidx, 0, ofw, 0);
							//orect.initIval(orefival);
							//oseenDiags.add(orefival);
							
							// Now fill the dynamic programming matrix, return true
							// iff there is at least one valid alignment
							TAlScore bestCell = std::numeric_limits<TAlScore>::min();
							foundMate = oswa.align(rnd, bestCell);
							prm.nMateDps++;
							swmMate.tallyGappedDp(oreadGaps, orefGaps);
							if(!foundMate) {
								TAlScore bestLast = anchor1 ? prm.bestLtMinscMate2 : prm.bestLtMinscMate1;
								if(bestCell > std::numeric_limits<TAlScore>::min() && bestCell > bestLast) {
									if(anchor1) {
										prm.bestLtMinscMate2 = bestCell;
									} else {
										prm.bestLtMinscMate1 = bestCell;
									}
								}
							}
						}
						bool didAnchor = false;
						do {
							oresGap_.reset();
							assert(oresGap_.empty());
							if(foundMate && oswa.done()) {
								foundMate = false;
							} else if(foundMate) {
								oswa.nextAlignment(oresGap_, ominsc_cur, rnd);
								foundMate = !oresGap_.empty();
								assert(!foundMate || oresGap_.alres.matchesRef(
																			   ord,
																			   ref,
																			   tmp_rf_,
																			   tmp_rdseq_,
																			   tmp_qseq_,
																			   raw_refbuf_,
																			   raw_destU32_,
																			   raw_matches_,
																			   tmp_reflens_,
																			   tmp_refoffs_));
							}
							if(foundMate) {
								// Redundant with one we've seen previously?
								if(!redAnchor_.overlap(oresGap_.alres)) {
									redAnchor_.add(oresGap_.alres);
								}
								assert_eq(ofw, oresGap_.alres.fw());
								// Annotate the AlnRes object with some key parameters
								// that were used to obtain the alignment.
								oresGap_.alres.setParams(
														 seedmms,    // # mismatches allowed in seed
														 seedlen,    // length of seed
														 seedival,   // interval between seeds
														 ominsc);    // minimum score for valid alignment
								assert_gt(oresGap_.alres.refExtent(), 0);
								if(gReportOverhangs &&
								   !refival.containsIgnoreOrient(oresGap_.alres.refival()))
								{
									oresGap_.alres.clipOutside(true, 0, tlen);
									foundMate = oresGap_.alres.refExtent() > 0;
								}
								if(foundMate && 
								   ((!gReportOverhangs &&
									 !refival.containsIgnoreOrient(oresGap_.alres.refival())) ||
									!refival.overlapsIgnoreOrient(oresGap_.alres.refival())))
								{
									foundMate = false;
								}
							}
							ASSERT_ONLY(TRefId refid);
							TRefOff off1, off2;
							size_t len1, len2;
							bool fw1, fw2;
							int pairCl = PE_ALS_DISCORD;
							if(foundMate) {
								ASSERT_ONLY(refid =) res->alres.refid();
								assert_eq(refid, oresGap_.alres.refid());
								off1 = anchor1 ? off : oresGap_.alres.refoff();
								off2 = anchor1 ? oresGap_.alres.refoff() : off;
								len1 = anchor1 ?
								res->alres.refExtent() : oresGap_.alres.refExtent();
								len2 = anchor1 ?
								oresGap_.alres.refExtent() : res->alres.refExtent();
								fw1  = anchor1 ? res->alres.fw() : oresGap_.alres.fw();
								fw2  = anchor1 ? oresGap_.alres.fw() : res->alres.fw();
								// Check that final mate alignments are consistent with
								// paired-end fragment constraints
								pairCl = pepol.peClassifyPair(
															  off1,
															  len1,
															  fw1,
															  off2,
															  len2,
															  fw2);
								// Instead of trying
								//foundMate = pairCl != PE_ALS_DISCORD;
							}
							if(msink->state().doneConcordant()) {
								foundMate = false;
							}
							if(reportImmediately) {
								if(foundMate) {
									// Report pair to the AlnSinkWrap
									assert(!msink->state().doneConcordant());
									assert(msink != NULL);
									assert(res->repOk());
									assert(oresGap_.repOk());
									// Report an unpaired alignment
									assert(!msink->maxed());
									assert(!msink->state().done());
									bool doneUnpaired = false;
									//if(mixed || discord) {
									// Report alignment for mate #1 as an
									// unpaired alignment.
									if(!anchor1 || !didAnchor) {
										if(anchor1) {
											didAnchor = true;
										}
										const AlnRes& r1 = anchor1 ?
										res->alres : oresGap_.alres;
										if(!redMate1_.overlap(r1)) {
											redMate1_.add(r1);
											if(msink->report(0, &r1, NULL)) {
												doneUnpaired = true; // Short-circuited
											}
										}
									}
									// Report alignment for mate #2 as an
									// unpaired alignment.
									if(anchor1 || !didAnchor) {
										if(!anchor1) {
											didAnchor = true;
										}
										const AlnRes& r2 = anchor1 ?
										oresGap_.alres : res->alres;
										if(!redMate2_.overlap(r2)) {
											redMate2_.add(r2);
											if(msink->report(0, NULL, &r2)) {
												doneUnpaired = true; // Short-circuited
											}
										}
									}
									//} // if(mixed || discord)
									bool donePaired = false;
									if(pairCl != PE_ALS_DISCORD) {
										foundConcordant = true;
										if(msink->report(
														 0,
														 anchor1 ? &res->alres : &oresGap_.alres,
														 anchor1 ? &oresGap_.alres : &res->alres))
										{
											// Short-circuited because a limit, e.g.
											// -k, -m or -M, was exceeded
											donePaired = true;
										} else {
											if(tighten > 0 && msink->Mmode() && msink->hasSecondBestPair()) {
												// Paired-end alignments should have at least this score from now
												TAlScore ps;
												if(tighten == 1) {
													ps = msink->bestPair();
												} else if(tighten == 2) {
													ps = msink->secondBestPair();
												} else {
													TAlScore diff = msink->bestPair() - msink->secondBestPair();
													ps = msink->secondBestPair() + (diff * 3)/4;
												}
												if(tighten == 1 && ps < bestPairScore &&
												   msink->bestPair() == msink->secondBestPair())
												{
													ps++;
												}
												if(tighten >= 2 && ps < bestPairScore) {
													ps++;
												}
												// Anchor mate must have score at least 'ps' minus the best possible
												// score for the opposite mate.
												TAlScore nc = ps - operfectScore;
												if(nc > minsc) {
													minsc = nc;
													assert_leq(minsc, perfectScore);
													if(minsc > res->alres.score().score()) {
														// We're done with this anchor
														break;
													}
												}
												assert_leq(minsc, perfectScore);
											}
										}
									} // if(pairCl != PE_ALS_DISCORD)
									if(donePaired || doneUnpaired) {
										return EXTEND_POLICY_FULFILLED;
									}
									if(msink->state().doneWithMate(anchor1)) {
										// We're now done with the mate that we're
										// currently using as our anchor.  We're not
										// with the read overall.
										return EXTEND_POLICY_FULFILLED;
									}
								} else if((mixed || discord) && !didAnchor) {
									didAnchor = true;
									// Report unpaired hit for anchor
									assert(msink != NULL);
									assert(res->repOk());
									// Check that alignment accurately reflects the
									// reference characters aligned to
									assert(res->alres.matchesRef(
																 rd,
																 ref,
																 tmp_rf_,
																 tmp_rdseq_,
																 tmp_qseq_,
																 raw_refbuf_,
																 raw_destU32_,
																 raw_matches_,
																 tmp_reflens_,
																 tmp_refoffs_));
									// Report an unpaired alignment
									assert(!msink->maxed());
									assert(!msink->state().done());
									// Report alignment for mate #1 as an
									// unpaired alignment.
									if(!msink->state().doneUnpaired(anchor1)) {
										const AlnRes& r = res->alres;
										RedundantAlns& red = anchor1 ? redMate1_ : redMate2_;
										const AlnRes* r1 = anchor1 ? &res->alres : NULL;
										const AlnRes* r2 = anchor1 ? NULL : &res->alres;
										if(!red.overlap(r)) {
											red.add(r);
											if(msink->report(0, r1, r2)) {
												return EXTEND_POLICY_FULFILLED; // Short-circuited
											}
										}
									}
									if(msink->state().doneWithMate(anchor1)) {
										// Done with mate, but not read overall
										return EXTEND_POLICY_FULFILLED;
									}
								}
							}
						} while(!oresGap_.empty());
					} // if(found && swMateImmediately)
					else if(found) {
						assert(!msink->state().doneWithMate(anchor1));
						// We found an anchor alignment but did not attempt to find
						// an alignment for the opposite mate (probably because
						// we're done with it)
						if(reportImmediately && (mixed || discord)) {
							// Report unpaired hit for anchor
							assert(msink != NULL);
							assert(res->repOk());
							// Check that alignment accurately reflects the
							// reference characters aligned to
							assert(res->alres.matchesRef(
														 rd,
														 ref,
														 tmp_rf_,
														 tmp_rdseq_,
														 tmp_qseq_,
														 raw_refbuf_,
														 raw_destU32_,
														 raw_matches_,
														 tmp_reflens_,
														 tmp_refoffs_));
							// Report an unpaired alignment
							assert(!msink->maxed());
							assert(!msink->state().done());
							// Report alignment for mate #1 as an
							// unpaired alignment.
							if(!msink->state().doneUnpaired(anchor1)) {
								const AlnRes& r = res->alres;
								RedundantAlns& red = anchor1 ? redMate1_ : redMate2_;
								const AlnRes* r1 = anchor1 ? &res->alres : NULL;
								const AlnRes* r2 = anchor1 ? NULL : &res->alres;
								if(!red.overlap(r)) {
									red.add(r);
									if(msink->report(0, r1, r2)) {
										return EXTEND_POLICY_FULFILLED; // Short-circuited
									}
								}
							}
							if(msink->state().doneWithMate(anchor1)) {
								// Done with mate, but not read overall
								return EXTEND_POLICY_FULFILLED;
							}
						}
					}
				} // while(true)
				
				if(foundConcordant) {
					prm.nMateDpSuccs++;
					mateStreaks_[i] = 0;
					// Register this as a success.  Now we need to
					// make the streak variables reflect the
					// success.
					if(state == FOUND_UNGAPPED) {
						assert_gt(prm.nUgFail, 0);
						assert_gt(prm.nExUgFails, 0);
						prm.nExUgFails--;
						prm.nExUgSuccs++;
						prm.nUgLastSucc = prm.nExUgs-1;
						if(prm.nUgFail > prm.nUgFailStreak) {
							prm.nUgFailStreak = prm.nUgFail;
						}
						prm.nUgFail = 0;
					} else if(state == FOUND_EE) {
						assert_gt(prm.nEeFail, 0);
						assert_gt(prm.nExEeFails, 0);
						prm.nExEeFails--;
						prm.nExEeSuccs++;
						prm.nEeLastSucc = prm.nExEes-1;
						if(prm.nEeFail > prm.nEeFailStreak) {
							prm.nEeFailStreak = prm.nEeFail;
						}
						prm.nEeFail = 0;
					} else {
						assert_gt(prm.nDpFail, 0);
						assert_gt(prm.nExDpFails, 0);
						prm.nExDpFails--;
						prm.nExDpSuccs++;
						prm.nDpLastSucc = prm.nExDps-1;
						if(prm.nDpFail > prm.nDpFailStreak) {
							prm.nDpFailStreak = prm.nDpFail;
						}
						prm.nDpFail = 0;
					}
				} else {
					prm.nMateDpFails++;
					mateStreaks_[i]++;
				}
				// At this point we know that we aren't bailing, and will continue to resolve seed hits.  
				
			} // while(!gw.done())
		} // for(size_t i = 0; i < gws_.size(); i++)
	}
	return EXTEND_EXHAUSTED_CANDIDATES;
}

#endif /*ALIGNER_SW_DRIVER_H_*/
