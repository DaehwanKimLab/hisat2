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

#ifndef ALIGNER_SEED2_H_
#define ALIGNER_SEED2_H_

/**
 * The user of the DescentDriver class specifies a collection of search roots.
 * Logic for picking these search roots is located elsewhere, not in this
 * module.  The search roots are annotated with a priority score, which 
 *
 * The heap is a min-heap over pairs, where the first element of each pair is
 * the score associated with a descent and the second element of each pair is
 * the descent ID.
 *
 * Weeding out redundant descents is key; otherwise we end up reporting slight
 * variations on the same alignment repeatedly, including variations with poor
 * scores.  What criteria do we use to determine whether two paths are
 * redundant?
 *
 * Here's an example where the same set of read characters have been aligned in
 * all three cases:
 *
 * Alignment 1 (sc = 0):
 * Rd: GCTATATAGCGCGCTCGCATCATTTTGTGT
 *     ||||||||||||||||||||||||||||||
 * Rf: GCTATATAGCGCGCTCGCATCATTTTGTGT
 *
 * Alignment 2 (sc = -22):
 * Rd: GCTATATAGCGCGCTCGCATCATTTTGTGT
 *     |||||||||||||||||||||||  | |||
 * Rf: GCTATATAGCGCGCTCGCATCAT--TTTGT
 *
 * Alignment 3 (sc = -22):
 * Rd: GCTATATAGCGCGCTCGCATCATT--TTGTGT
 *     ||||||||||||||||||||||||   |||||
 * Rf: GCTATATAGCGCGCTCGCATCATTTTGTGTGT
 *
 * Rf from aln 1: GCTATATAGCGCGCTCGCATCATTTTGTGT
 * Rf from aln 2: GCTATATAGCGCGCTCGCATCATTTTGT
 * Rf from aln 3: GCTATATAGCGCGCTCGCATCATTTTGTGTGT
 *
 * Are alignments 2 and 3 redundant with alignment 1?  We can't totally say
 * without knowing the associated SA ranges.  Take alignments 1 and 2.  Either
 * the SA ranges are the same or the SA range for 2 contains the SA range for
 * 1.  If they're the same, then alignment 2 is redundant with alignment 1.
 * Otherwise, *some* of the elements in the SA range for alignment 2 are not
 * redundant.
 *
 * In that example, the same read characters are aligned in all three
 * alignments.  Is it possible and profitable to consider scenarios where an
 * alignment might be redundant with another alignment 
 *
 * Another question is *when* do we try to detect the redundancy?  Before we
 * try to extend through the matches, or after.  After is easier, but less work
 * has been avoided.
 *
 * What data structure do we query to determine whether there's redundancy?
 * The situation is harder when we try to detect overlaps between SA ranges
 * rather than identical SA ranges.  Maybe: read intervals -> intersection tree -> penalties.
 *
 * 1. If we're introducing a gap and we could have introduced it deeper in the
 *    descent with the same effect w/r/t homopolymer length.
 * 2. If we have Descent A with penalty B and Descent a with penalty b, and A
 *    aligns read characters [X, Y] to SA range [Z, W], and B aligns read
 *    characters [x, y] to SA range [z, w], then A is redundant with B if
 *    [x, y] is within [X, Y].
 *
 * Found an alignment with total penalty = 3
 * GCAATATAGCGCGCTCGCATCATTTTGTGT
 * || |||||||||||||||||||||||||||
 * GCTATATAGCGCGCTCGCATCATTTTGTGT
 *
 * Found an alignment with total penalty = 27
 * gCAATATAGCGCGCTCGCATCATTTTGTGT
 *   |   ||||||||||||||||||||||||
 *  TATA-TAGCGCGCTCGCATCATTTTGTGT
 */

#include <stdint.h>
#include <math.h>
#include <utility>
#include <limits>
#include "assert_helpers.h"
#include "random_util.h"
#include "aligner_result.h"
#include "gfm.h"
#include "simple_func.h"
#include "scoring.h"
#include "edit.h"
#include "read.h"
#include "ds.h"
#include "group_walk.h"
#include "btypes.h"

typedef size_t   TReadOff;
typedef int64_t  TScore;
typedef float    TRootPri;
typedef size_t   TDescentId;
typedef size_t   TRootId;

/**
 * enum encapsulating a few different policies for how we might extend descents
 * in the direction opposite from their primary direction.  
 */
enum {
	// Never extened in the direction opposite from the primary.  Just go in
	// the primary direction until the bounce.
	DESC_EX_NONE = 1,
	
	// When we're finished extending out the matches for a descent, try to
	// extend in the opposite direction in a way that extends all branches
	// simultaneously.  The Descent.nex_ field contains the number of positions
	// we were able to extend through in this way.
	DESC_EX_FROM_1ST_BRANCH = 2,
	
	// Each time we add an edge to the summary, extend it in the opposite
	// direction.  The DescentEdge.nex field contains the number of positions
	// we were able to extend through, and this in turn gets propagated to
	// Descent.nex_ if and when we branch from the DescentEdge.
	DESC_EX_EACH_EDGE = 3
};

/**
 * Counters to keep track of how much work is being done.
 */
struct DescentMetrics {

	DescentMetrics() { reset(); }

	void reset() {
		bwops = bwops_1 = bwops_bi = recalc = branch = branch_mm =
		branch_del = branch_ins = heap_max = descent_max = descentpos_max = 
		nex = 0;
	}

	uint64_t bwops;          // # FM Index opbs
	uint64_t bwops_1;        // # LF1 FM Index opbs
	uint64_t bwops_bi;       // # BiEx FM Index opbs
	uint64_t recalc;         // # times outgoing edge summary was recalculated
	uint64_t branch;         // # times we descended from another descent
	uint64_t branch_mm;      // # times branch was on a mismatch
	uint64_t branch_del;     // # times branch was on a deletion
	uint64_t branch_ins;     // # times branch was on a insertion
	uint64_t heap_max;       // maximum size of Descent heap
	uint64_t descent_max;    // maximum size of Descent factory
	uint64_t descentpos_max; // maximum size of DescentPos factory
	uint64_t nex;            // # extensions
};

/**
 * Priority used to rank which descent we should branch from next.  Right now,
 * priority is governed by a 4-tuple.  From higher to lower priority:
 *
 *  1. Penalty accumulated so far
 *  2. Depth into the search space, including extensions
 *  3. Width of the SA range (i.e. uniqueness)
 *  4. Root priority
 */
struct DescentPriority {

	DescentPriority() { reset(); }

	DescentPriority(
		TScore pen_,
		size_t depth_,
		TIndexOffU width_,
		float rootpri_)
	{
		pen = pen_;
		depth = depth_;
		width = width_;
		rootpri = rootpri_;
	}
	
	/**
	 * Initialize new DescentPriority.
	 */
	void init(TScore pen_, size_t depth_, TIndexOffU width_, float rootpri_) {
		pen = pen_;
		depth = depth_;
		width = width_;
		rootpri = rootpri_;
	}
	
	/**
	 * Reset to uninitialized state.
	 */
	void reset() {
		width = 0;
	}
	
	/**
	 * Return true iff DescentPriority is initialized.
	 */
	bool inited() const {
		return width > 0;
	}
	
	/**
	 * Return true iff this priority is prior to given priority.
	 */
	bool operator<(const DescentPriority& o) const {
		assert(inited());
		assert(o.inited());
		// 1st priority: penalty accumulated so far
		if(pen < o.pen) return true;
		if(pen > o.pen) return false;
		// 2nd priority: depth into the search space, including extensions
		if(depth > o.depth) return true;
		if(depth < o.depth) return false;
		// 3rd priority: width of the SA range (i.e. uniqueness)
		if(width < o.width) return true;
		if(width > o.width) return false;
		// 4th priority: root priority
		if(rootpri > o.rootpri) return true;
		return false;
	}

	/**
	 * Return true iff this priority is prior to or equal to given priority.
	 */
	bool operator<=(const DescentPriority& o) const {
		assert(inited());
		assert(o.inited());
		// 1st priority: penalty accumulated so far
		if(pen < o.pen) return true;
		if(pen > o.pen) return false;
		// 2nd priority: depth into the search space, including extensions
		if(depth > o.depth) return true;
		if(depth < o.depth) return false;
		// 3rd priority: width of the SA range (i.e. uniqueness)
		if(width < o.depth) return true;
		if(width > o.width) return false;
		// 4th priority: root priority
		if(rootpri > o.rootpri) return true;
		return true;
	}

	/**
	 * Return true iff this priority is prior to or equal to given priority.
	 */
	bool operator==(const DescentPriority& o) const {
		assert(inited());
		assert(o.inited());
		return pen == o.pen && depth == o.depth && width == o.width && rootpri == o.rootpri;
	}

	TScore pen;      // total penalty accumulated so far
	size_t depth;    // depth from root of descent
	TIndexOffU width; // width of the SA range
	float  rootpri;  // priority of the root
};

static inline std::ostream& operator<<(
	std::ostream& os,
	const DescentPriority& o)
{
	os << "[" << o.pen << ", " << o.depth << ", " << o.width << ", " << o.rootpri << "]";
	return os;
}

static inline std::ostream& operator<<(
	std::ostream& os,
	const std::pair<DescentPriority, TDescentId>& o)
{
	os << "{[" << o.first.pen << ", " << o.first.depth << ", "
	   << o.first.width << ", " << o.first.rootpri << "], " << o.second << "}";
	return os;
}

typedef std::pair<DescentPriority, TDescentId> TDescentPair;

/**
 * Encapsulates the constraints limiting which outgoing edges are permitted.
 * Specifically, we constrain the total penalty accumulated so far so that some
 * outgoing edges will exceed the limit and be pruned.  The limit is set
 * according to our "depth" into the search, as measured by the number of read
 * characters aligned so far.  We divide the depth domain into two pieces, a
 * piece close to the root, where the penty is constrained to be 0, and the
 * remainder, where the maximum penalty is an interpolation between 0 and the
 * maximum penalty
 */
struct DescentConstraints {

	DescentConstraints() { reset(); }
	
	/**
	 * Initialize with new constraint function.
	 */
	DescentConstraints(size_t nzero, double exp) {
        init(nzero, exp);
	}
    
    /**
     * Initialize with given function.
     */
    void init(size_t nzero_, double exp_) {
		nzero = nzero_ > 0 ? nzero_ : 1;
		exp = exp_;
#ifndef NDEBUG
		for(size_t i = 1; i < nzero_ + 5; i++) {
			assert_geq(get(i, nzero_ + 10, 100), get(i-1, nzero_ + 10, 100));
		}
#endif
    }
    
    /**
     * Reset to uninitialized state.
     */
    void reset() {
        nzero = 0;
		exp = -1.0f;
    }
    
    /**
     * Return true iff the DescentConstraints has been initialized.
     */
    bool inited() const {
        return exp >= 0.0f;
    }
	
	/**
	 * Get the maximum penalty total for depth 'off'.
	 */
	inline TScore get(TReadOff off, TReadOff rdlen, TAlScore maxpen) const {
		if(off < nzero || nzero >= rdlen) {
			return 0;
		}
		double frac = (double)(off - nzero) / (rdlen - nzero);
		if(fabs(exp - 1.0f) > 0.00001) {
			if(fabs(exp - 2.0f) < 0.00001) {
				frac *= frac;
			} else {
				frac = pow(frac, exp);
			}
		}
		return (TAlScore)(frac * maxpen + 0.5f);
	}

	size_t nzero;
	double exp;
};

/**
 * Encapsulates settings governing how we descent.
 */
struct DescentConfig {

    DescentConfig() { reset(); }
    
    /**
     * Reset the DescentConfig to an uninitialized state.
     */
    void reset() { expol = 0; }
    
    /**
     * Return true iff this DescentConfig is initialized.
     */
    bool inited() const { return expol != 0; }

    DescentConstraints cons; // constraints
    int expol; // extend policy
};

/**
 * Encapsulates the state of a Descent that allows us to determine whether it
 * is redundant with another Descent.  Two Descents are redundant if:
 *
 * 1. Both are aligning the same read orientation (fw or rc)
 * 2. Both are growing the alignment in the same direction (left-to-right or
 *    right-to-left)
 * 3. They have aligned exactly the same read characters (which are always
 *    consecutive in the read)
 * 4. The corresponding reference strings are identical
 */
struct DescentRedundancyKey {

	DescentRedundancyKey() { reset(); }
	
	DescentRedundancyKey(
		TReadOff  al5pf_,
		size_t    rflen_,
		TIndexOffU topf_,
		TIndexOffU botf_)
	{
		init(al5pf_, rflen_, topf_, botf_);
	}

	void reset() {
		al5pf = 0;
		rflen = 0;
		topf = botf = 0;
	}
	
	bool inited() const { return rflen > 0; }

	void init(
		TReadOff  al5pf_,
		size_t    rflen_,
		TIndexOffU topf_,
		TIndexOffU botf_)
	{
		al5pf = al5pf_;
		rflen = rflen_;
		topf = topf_;
		botf = botf_;
	}
	
	bool operator==(const DescentRedundancyKey& o) const {
		return al5pf == o.al5pf && rflen == o.rflen && topf == o.topf && botf == o.botf;
	}

	bool operator<(const DescentRedundancyKey& o) const {
		if(al5pf < o.al5pf) return true;
		if(al5pf > o.al5pf) return false;
		if(rflen < o.rflen) return true;
		if(rflen > o.rflen) return false;
		if(topf < o.topf) return true;
		if(topf > o.topf) return false;
		return botf < o.botf;
	}

	TReadOff al5pf; // 3'-most aligned char, as offset from 5' end
	size_t rflen;   // number of reference characters involved in alignment
	TIndexOffU topf; // top w/r/t forward index
	TIndexOffU botf; // bot w/r/t forward index
};

/**
 * Map from pairs to top, bot, penalty triples.
 */
class DescentRedundancyChecker {

public:

	DescentRedundancyChecker() { reset(); }

	void clear() { reset(); }
	
	/**
	 * Reset to uninitialized state.
	 */
	void reset() {
		bits_.reset();
		inited_ = false;
		totsz_ = 0;  // total size
		totcap_ = 0; // total capacity
	}
	
	const static int NPARTS = 8;
	const static int PART_MASK = 7;
	const static int NBITS = (1 << 16);

	/**
	 * Initialize using given read length.
	 */
	void init(TReadOff rdlen) {
		reset();
		// daehwan - for debugging purposes
#if 0
		bits_.resize(NBITS);
		maplist_fl_.resize(NPARTS);
		maplist_fr_.resize(NPARTS);
		maplist_rl_.resize(NPARTS);
		maplist_rr_.resize(NPARTS);
		for(int i = 0; i < NPARTS; i++) {
			maplist_fl_[i].resize(rdlen);
			maplist_fr_[i].resize(rdlen);
			maplist_rl_[i].resize(rdlen);
			maplist_rr_[i].resize(rdlen);
			totcap_ += maplist_fl_[i].totalCapacityBytes();
			totcap_ += maplist_fr_[i].totalCapacityBytes();
			totcap_ += maplist_rl_[i].totalCapacityBytes();
			totcap_ += maplist_rr_[i].totalCapacityBytes();
			for(size_t j = 0; j < rdlen; j++) {
				maplist_fl_[i][j].clear();
				maplist_fr_[i][j].clear();
				maplist_rl_[i][j].clear();
				maplist_rr_[i][j].clear();
				totcap_ += maplist_fl_[i][j].totalCapacityBytes();
				totcap_ += maplist_fr_[i][j].totalCapacityBytes();
				totcap_ += maplist_rl_[i][j].totalCapacityBytes();
				totcap_ += maplist_rr_[i][j].totalCapacityBytes();
			}
		}
#endif
		inited_ = true;
	}
	
	/**
	 * Return true iff the checker is initialized.
	 */
	bool inited() const {
		return inited_;
	}

	/**
	 * Check if this partial alignment is redundant with one that we've already
	 * explored.
	 */
	bool check(
		bool fw,
		bool l2r,
		TReadOff al5pi,
		TReadOff al5pf,
		size_t rflen,
		TIndexOffU topf,
		TIndexOffU botf,
		TScore pen)
	{
		// daehwan - for debugging purposes
		return true;
		
		assert(inited_);
		assert(topf > 0 || botf > 0);
		DescentRedundancyKey k(al5pf, rflen, topf, botf);
		size_t i = std::numeric_limits<size_t>::max();
		size_t mask = topf & PART_MASK;
		EMap<DescentRedundancyKey, TScore>& map =
			(fw ? (l2r ? maplist_fl_[mask][al5pi] : maplist_fr_[mask][al5pi]) :
			      (l2r ? maplist_rl_[mask][al5pi] : maplist_rr_[mask][al5pi]));
		size_t key = (topf & 255) | ((botf & 255) << 8);
		if(bits_.test(key) && map.containsEx(k, i)) {
			// Already contains the key
			assert_lt(i, map.size());
			assert_geq(pen, map[i].second);
			return false;
		}
		assert(!map.containsEx(k, i));
		size_t oldsz = map.totalSizeBytes();
		size_t oldcap = map.totalCapacityBytes();
		map.insert(make_pair(k, pen));
		bits_.set(key);
		totsz_ += (map.totalSizeBytes() - oldsz);
		totcap_ += (map.totalCapacityBytes() - oldcap);
		return true;
	}

	/**
	 * Check if this partial alignment is redundant with one that we've already
	 * explored using the Bw index SA range.
	 */
	bool contains(
		bool fw,
		bool l2r,
		TReadOff al5pi,
		TReadOff al5pf,
		size_t rflen,
		TIndexOffU topf,
		TIndexOffU botf,
		TScore pen)
	{
		// daehwan - for debugging purposes
		return false;
		
		assert(inited_);
		size_t key = (topf & 255) | ((botf & 255) << 8);
		if(!bits_.test(key)) {
			return false;
		}
		DescentRedundancyKey k(al5pf, rflen, topf, botf);
		size_t mask = topf & PART_MASK;
		EMap<DescentRedundancyKey, TScore>& map =
			(fw ? (l2r ? maplist_fl_[mask][al5pi] : maplist_fr_[mask][al5pi]) :
			      (l2r ? maplist_rl_[mask][al5pi] : maplist_rr_[mask][al5pi]));
		return map.contains(k);
	}
	
	/**
	 * Return the total size of the redundancy map.
	 */
	size_t totalSizeBytes() const {
		return totsz_;
	}

	/**
	 * Return the total capacity of the redundancy map.
	 */
	size_t totalCapacityBytes() const {
		return totcap_;
	}

protected:

	bool inited_;   // initialized?
	size_t totsz_;  // total size
	size_t totcap_; // total capacity
	
	// List of maps.  Each entry is a map for all the DescentRedundancyKeys
	// with al5pi equal to the offset into the list.
	ELList<EMap<DescentRedundancyKey, TScore>, NPARTS, 100> maplist_fl_; //  fw,  l2r
	ELList<EMap<DescentRedundancyKey, TScore>, NPARTS, 100> maplist_rl_; // !fw,  l2r
	ELList<EMap<DescentRedundancyKey, TScore>, NPARTS, 100> maplist_fr_; //  fw, !l2r
	ELList<EMap<DescentRedundancyKey, TScore>, NPARTS, 100> maplist_rr_; // !fw, !l2r
		
	EBitList<128> bits_;
};

/**
 * A search root.  Consists of an offset from the 5' end read and flags
 * indicating (a) whether we're initially heading left-to-right or
 * right-to-left, and (b) whether we're examining the read or its reverse
 * complement.
 *
 * A root also comes with a priority ("pri") score indicating how promising it
 * is as a root.  Promising roots have long stretches of high-quality,
 * non-repetitive nucleotides in the first several ply of the search tree.
 * Also, roots beginning at the 5' end of the read may receive a higher
 * priority.
 */
struct DescentRoot {

	DescentRoot() { reset(); }

	DescentRoot(size_t off5p_, bool l2r_, bool fw_, size_t len, float pri_) {
		init(off5p_, l2r_, fw_, len, pri_);
	}
	
	/**
	 * Reset this DescentRoot to uninitialized state.
	 */
	void reset() {
		off5p = std::numeric_limits<size_t>::max();
	}
	
	/**
	 * Return true iff this DescentRoot is uninitialized.
	 */
	bool inited() const {
		return off5p == std::numeric_limits<size_t>::max();
	}
	
	/**
	 * Initialize a new descent root.
	 */
	void init(size_t off5p_, bool l2r_, bool fw_, size_t len, float pri_) {
		off5p = off5p_;
		l2r = l2r_;
		fw = fw_;
		pri = pri_;
		assert_lt(off5p, len);
	}

	TReadOff off5p;   // root origin offset, expressed as offset from 5' end
	bool     l2r;     // true -> move in left-to-right direction
	bool     fw;      // true -> work with forward read, false -> revcomp
	float    pri;     // priority of seed
};

/**
 * Set of flags indicating outgoing edges we've tried from a DescentPos.
 */
struct DescentPosFlags {

	DescentPosFlags() { reset(); }
	
	/**
	 * Set all flags to 1, indicating all outgoing edges are yet to be
	 * explored.
	 */
	void reset() {
		mm_a = mm_c = mm_g = mm_t = rdg_a = rdg_c = rdg_g = rdg_t = rfg = 1;
		reserved = 0;
	}
	
	/**
	 * Return true iff all outgoing edges have already been explored.
	 */
	bool exhausted() const {
		return ((uint16_t*)this)[0] == 0;
	}
	
	/**
	 * Return false iff the specified mismatch has already been explored.
	 */
	bool mmExplore(int c) {
		assert_range(0, 3, c);
		if(c == 0) {
			return mm_a;
		} else if(c == 1) {
			return mm_c;
		} else if(c == 2) {
			return mm_g;
		} else {
			return mm_t;
		}
	}

	/**
	 * Try to explore a mismatch.  Return false iff it has already been
	 * explored.
	 */
	bool mmSet(int c) {
		assert_range(0, 3, c);
		if(c == 0) {
			bool ret = mm_a; mm_a = 0; return ret;
		} else if(c == 1) {
			bool ret = mm_c; mm_c = 0; return ret;
		} else if(c == 2) {
			bool ret = mm_g; mm_g = 0; return ret;
		} else {
			bool ret = mm_t; mm_t = 0; return ret;
		}
	}

	/**
	 * Return false iff specified read gap has already been explored.
	 */
	bool rdgExplore(int c) {
		assert_range(0, 3, c);
		if(c == 0) {
			return rdg_a;
		} else if(c == 1) {
			return rdg_c;
		} else if(c == 2) {
			return rdg_g;
		} else {
			return rdg_t;
		}
	}

	/**
	 * Try to explore a read gap.  Return false iff it has already been
	 * explored.
	 */
	bool rdgSet(int c) {
		assert_range(0, 3, c);
		if(c == 0) {
			bool ret = rdg_a; rdg_a = 0; return ret;
		} else if(c == 1) {
			bool ret = rdg_c; rdg_c = 0; return ret;
		} else if(c == 2) {
			bool ret = rdg_g; rdg_g = 0; return ret;
		} else {
			bool ret = rdg_t; rdg_t = 0; return ret;
		}
	}

	/**
	 * Return false iff the reference gap has already been explored.
	 */
	bool rfgExplore() {
		return rfg;
	}

	/**
	 * Try to explore a reference gap.  Return false iff it has already been
	 * explored.
	 */
	bool rfgSet() {
		bool ret = rfg; rfg = 0; return ret;
	}

	uint16_t mm_a     : 1;
	uint16_t mm_c     : 1;
	uint16_t mm_g     : 1;
	uint16_t mm_t     : 1;

	uint16_t rdg_a    : 1;
	uint16_t rdg_c    : 1;
	uint16_t rdg_g    : 1;
	uint16_t rdg_t    : 1;

	uint16_t rfg      : 1;
	
	uint16_t reserved : 7;
};

/**
 * FM Index state associated with a single position in a descent.  For both the
 * forward and backward indexes, it stores the four SA ranges corresponding to
 * the four nucleotides.
 */
struct DescentPos {

	/**
	 * Reset all tops and bots to 0.
	 */
	void reset() {
		topf[0] = topf[1] = topf[2] = topf[3] = 0;
		botf[0] = botf[1] = botf[2] = botf[3] = 0;
		topb[0] = topb[1] = topb[2] = topb[3] = 0;
		botb[0] = botb[1] = botb[2] = botb[3] = 0;
        c = -1;
		flags.reset();
	}
    
    /**
     * Return true iff DescentPos has been initialized.
     */
    bool inited() const {
        return c >= 0;
    }
	
#ifndef NDEBUG
	/**
	 * Check that DescentPos is internally consistent.
	 */
	bool repOk() const {
		assert_range(0, 3, (int)c);
		return true;
	}
#endif
	
	TIndexOffU       topf[4]; // SA range top indexes in fw index
	TIndexOffU       botf[4]; // SA range bottom indexes (exclusive) in fw index
	TIndexOffU       topb[4]; // SA range top indexes in bw index
	TIndexOffU       botb[4]; // SA range bottom indexes (exclusive) in bw index
    char            c;       // read char that would yield match
	DescentPosFlags flags;   // flags 
};

/**
 * Encapsulates an edge outgoing from a descent.
 */
struct DescentEdge {

	DescentEdge() { reset(); }

	DescentEdge(
		Edit e_,
		TReadOff off5p_,
		DescentPriority pri_,
        size_t posFlag_,
		TReadOff nex_
#ifndef NDEBUG
        ,
        size_t d_,
		TIndexOffU topf_,
		TIndexOffU botf_,
		TIndexOffU topb_,
		TIndexOffU botb_
#endif
        )
	{
		init(e_, off5p_, pri_, posFlag_
#ifndef NDEBUG
        , d_, topf_, botf_, topb_, botb_
#endif
        );
	}

	/**
	 * Return true iff edge is initialized.
	 */
	bool inited() const { return e.inited(); }

	/**
	 * Reset to uninitialized state.
	 */
	void reset() { e.reset(); }
	
	/**
	 * Initialize DescentEdge given 5' offset, nucleotide, and priority.
	 */
	void init(
		Edit e_,
		TReadOff off5p_,
		DescentPriority pri_,
        size_t posFlag_
#ifndef NDEBUG
        ,
        size_t d_,
		TIndexOffU topf_,
		TIndexOffU botf_,
		TIndexOffU topb_,
		TIndexOffU botb_
#endif
        )
	{
		e = e_;
		off5p = off5p_;
		pri = pri_;
        posFlag = posFlag_;
#ifndef NDEBUG
        d = d_;
		topf = topf_;
		botf = botf_;
		topb = topb_;
		botb = botb_;
#endif
	}

    /**
     * Update flags to show this edge as visited.
     */
    void updateFlags(EFactory<DescentPos>& pf) {
        if(inited()) {
            if(e.isReadGap()) {
                assert_neq('-', e.chr);
                pf[posFlag].flags.rdgSet(asc2dna[e.chr]);
            } else if(e.isRefGap()) {
                pf[posFlag].flags.rfgSet();
            } else {
                assert_neq('-', e.chr);
                pf[posFlag].flags.mmSet(asc2dna[e.chr]);
            }
        }
    }
	
	/**
	 * Return true iff this edge has higher priority than the given edge.
	 */
	bool operator<(const DescentEdge& o) const {
        if(inited() && !o.inited()) {
            return true;
        } else if(!inited()) {
            return false;
        }
		return pri < o.pri;
	}

	DescentPriority pri; // priority of the edge
	//TReadOff nex;        // # extends possible from this edge
    size_t posFlag;      // depth of DescentPos where flag should be set


#ifndef NDEBUG
    // This can be recreated by looking at the edit, the paren't descent's
    // len_, al5pi_, al5pf_.  I have it here so we can sanity check.
    size_t d;
	TIndexOffU topf, botf, topb, botb;
#endif

	Edit e;
	TReadOff off5p;
};

/**
 * Encapsulates an incomplete summary of the outgoing edges from a descent.  We
 * don't try to store information about all outgoing edges, because doing so
 * will generally be wasteful.  We'll typically only try a handful of them per
 * descent.
 */
class DescentOutgoing {

public:

	/**
	 * Return the best edge and rotate in preparation for next call.
	 */
	DescentEdge rotate() {
		DescentEdge tmp = best1;
        assert(!(best2 < tmp));
		best1 = best2;
        assert(!(best3 < best2));
		best2 = best3;
        assert(!(best4 < best3));
		best3 = best4;
        assert(!(best5 < best4));
		best4 = best5;
		best5.reset();
		return tmp;
	}
	
	/**
	 * Given a potental outgoing edge, place it where it belongs in the running
	 * list of best 5 outgoing edges from this descent.
	 */
	void update(DescentEdge e) {
		if(!best1.inited()) {
			best1 = e;
		} else if(e < best1) {
			best5 = best4;
			best4 = best3;
			best3 = best2;
			best2 = best1;
			best1 = e;
		} else if(!best2.inited()) {
			best2 = e;
		} else if(e < best2) {
			best5 = best4;
			best4 = best3;
			best3 = best2;
			best2 = e;
		} else if(!best3.inited()) {
			best3 = e;
		} else if(e < best3) {
			best5 = best4;
			best4 = best3;
			best3 = e;
		} else if(!best4.inited()) {
			best4 = e;
		} else if(e < best4) {
			best5 = best4;
			best4 = e;
		}  else if(!best5.inited() || e < best5) {
			best5 = e;
		}
	}
	
	/**
	 * Clear all the outgoing edges stored here.
	 */
	void clear() {
		best1.reset();
		best2.reset();
		best3.reset();
		best4.reset();
		best5.reset();
	}
	
	/**
	 * Return true iff there are no outgoing edges currently represented in
	 * this summary.  There may still be outgoing edges, they just haven't
	 * been added to the summary.
	 */
	bool empty() const {
		return !best1.inited();
	}
	
	/**
	 * Return the DescentPriority of the best outgoing edge.
	 */
	DescentPriority bestPri() const {
		assert(!empty());
		return best1.pri;
	}

	DescentEdge best1; // best
	DescentEdge best2; // 2nd-best
	DescentEdge best3; // 3rd-best
	DescentEdge best4; // 4th-best
	DescentEdge best5; // 5th-best
};

template <typename index_t>
class DescentAlignmentSink;

/**
 * Encapsulates a descent through a search tree, along a path of matches.
 * Descents that are part of the same alignment form a chain.  Two aligments
 * adjacent in the chain are connected either by an edit, or by a switch in
 * direction.  Because a descent might have a different direction from the
 * DescentRoot it ultimately came from, it has its own 'l2r' field, which might
 * differ from the root's.
 */
template <typename index_t>
class Descent {

public:

	Descent() { reset(); }

	/**
	 * Initialize a new descent branching from the given descent via the given
	 * edit.  Return false if the Descent has no outgoing edges (and can
     * therefore have its memory freed), true otherwise.
	 */
	bool init(
		const Read& q,                  // query
		TRootId rid,                    // root id
		const Scoring& sc,              // scoring scheme
		TAlScore minsc,                 // minimum score
		TAlScore maxpen,                // maximum penalty
		TReadOff al5pi,                 // offset from 5' of 1st aligned char
		TReadOff al5pf,                 // offset from 5' of last aligned char
		TIndexOffU topf,                 // SA range top in FW index
		TIndexOffU botf,                 // SA range bottom in FW index
		TIndexOffU topb,                 // SA range top in BW index
		TIndexOffU botb,                 // SA range bottom in BW index
		bool l2r,                       // direction this descent will go in
		size_t descid,                  // my ID
		TDescentId parent,              // parent ID
		TScore pen,                     // total penalties so far
		const Edit& e,                  // edit for incoming edge
		const GFM<index_t>& gfmFw,      // forward index
		const GFM<index_t>& gfmBw,      // mirror index
		DescentRedundancyChecker& re,   // redundancy checker
		EFactory<Descent>& df,          // Descent factory
		EFactory<DescentPos>& pf,       // DescentPos factory
        const EList<DescentRoot>& rs,   // roots
        const EList<DescentConfig>& cs, // configs
		EHeap<TDescentPair>& heap,      // heap
        DescentAlignmentSink<index_t>& alsink,   // alignment sink
		DescentMetrics& met,            // metrics
		PerReadMetrics& prm);           // per-read metrics

	/**
	 * Initialize a new descent beginning at the given root.  Return false if
     * the Descent has no outgoing edges (and can therefore have its memory
     * freed), true otherwise.
	 */
	bool init(
        const Read& q,                  // query
        TRootId rid,                    // root id
        const Scoring& sc,              // scoring scheme
		TAlScore minsc,                 // minimum score
		TAlScore maxpen,                // maximum penalty
        size_t descid,                  // id of this Descent
        const GFM<index_t>& gfmFw,      // forward index
        const GFM<index_t>& gfmBw,      // mirror index
		DescentRedundancyChecker& re,   // redundancy checker
        EFactory<Descent>& df,          // Descent factory
        EFactory<DescentPos>& pf,       // DescentPos factory
        const EList<DescentRoot>& rs,   // roots
        const EList<DescentConfig>& cs, // configs
        EHeap<TDescentPair>& heap,      // heap
        DescentAlignmentSink<index_t>& alsink,   // alignment sink
        DescentMetrics& met,            // metrics
		PerReadMetrics& prm);           // per-read metrics
	
	/**
	 * Return true iff this Descent has been initialized.
	 */
	bool inited() const {
		return descid_ != std::numeric_limits<size_t>::max();
	}
	
	/**
	 * Reset to uninitialized state.
	 */
	void reset() {
        lastRecalc_ = true;
		descid_ = std::numeric_limits<size_t>::max();
	}
	
	/**
	 * Return true iff this Descent is a search root.
	 */
	bool root() const {
		return parent_ == std::numeric_limits<TDescentId>::max();
	}
	
	/**
	 * Return the edit.
	 */
	const Edit& edit() const {
		return edit_;
	}
	
	/**
	 * Return id of parent.
	 */
	TDescentId parent() const {
		return parent_;
	}
	
	/**
	 * Take the best outgoing edge and follow it.
	 */
	void followBestOutgoing(
        const Read& q,                  // read
        const GFM<index_t>& gfmFw,      // forward index
        const GFM<index_t>& gfmBw,      // mirror index
        const Scoring& sc,              // scoring scheme
		TAlScore minsc,                 // minimum score
		TAlScore maxpen,                // maximum penalty
		DescentRedundancyChecker& re,   // redundancy checker
		EFactory<Descent>& df,          // factory with Descent
		EFactory<DescentPos>& pf,       // factory with DescentPoss
        const EList<DescentRoot>& rs,   // roots
        const EList<DescentConfig>& cs, // configs
        EHeap<TDescentPair>& heap,      // heap of descents
        DescentAlignmentSink<index_t>& alsink,   // alignment sink
        DescentMetrics& met,            // metrics
		PerReadMetrics& prm);           // per-read metrics
	
	/**
	 * Return true iff no outgoing edges from this descent remain unexplored.
	 */
	bool empty() const { return lastRecalc_ && out_.empty(); }
	
#ifndef NDEBUG
	/**
	 * Return true iff the Descent is internally consistent.
	 */
	bool repOk(const Read *q) const {
		// A non-root can have an uninitialized edit_ if it is from a bounce
		//assert( root() ||  edit_.inited());
		assert(!root() || !edit_.inited());
		assert_eq(botf_ - topf_, botb_ - topb_);
		if(q != NULL) {
			assert_leq(len_, q->length());
		}
		return true;
	}
#endif
	
	size_t al5pi() const { return al5pi_; }
	size_t al5pf() const { return al5pf_; }
	bool l2r() const { return l2r_; }

	/**
	 * Print a stacked representation of this descent and all its parents.  Assumes that
	 */
	void print(
		std::ostream* os,
		const char *prefix,
        const Read& q,
		size_t trimLf,
		size_t trimRg,
		bool fw,
		const EList<Edit>& edits,
		size_t ei,
		size_t en,
		BTDnaString& rf) const;
	
	/**
	 * Collect all the edits
	 */
	void collectEdits(
		EList<Edit>& edits,
		const Edit *e,
		EFactory<Descent>& df)
	{
		// Take just the portion of the read that has aligned up until this
		// point
		size_t nuninited = 0;
		size_t ei = edits.size();
		size_t en = 0;
		if(e != NULL && e->inited()) {
			edits.push_back(*e);
			en++;
		}
		size_t cur = descid_;
		while(cur != std::numeric_limits<TDescentId>::max()) {
			if(!df[cur].edit().inited()) {
				nuninited++;
				assert_leq(nuninited, 2);
			} else {
				edits.push_back(df[cur].edit());
				en++;
			}
			cur = df[cur].parent();
		}
		// Sort just the edits we just added
		edits.sortPortion(ei, en);
	}

protected:

	/**
	 *
	 */
    bool bounce(
        const Read& q,                  // query string
        TIndexOffU topf,                 // SA range top in fw index
        TIndexOffU botf,                 // SA range bottom in fw index
        TIndexOffU topb,                 // SA range top in bw index
        TIndexOffU botb,                 // SA range bottom in bw index
        const GFM<index_t>& gfmFw,      // forward index
        const GFM<index_t>& gfmBw,      // mirror index
        const Scoring& sc,              // scoring scheme
		TAlScore minsc,                 // minimum score
		TAlScore maxpen,                // maximum penalty
		DescentRedundancyChecker& re,   // redundancy checker
		EFactory<Descent>& df,          // factory with Descent
		EFactory<DescentPos>& pf,       // factory with DescentPoss
        const EList<DescentRoot>& rs,   // roots
        const EList<DescentConfig>& cs, // configs
        EHeap<TDescentPair>& heap,      // heap of descents
        DescentAlignmentSink<index_t>& alsink,   // alignment sink
        DescentMetrics& met,            // metrics
		PerReadMetrics& prm);           // per-read metrics

	/**
	 * Given the forward and backward indexes, and given topf/botf/topb/botb,
	 * get tloc, bloc ready for the next step.
	 */
	void nextLocsBi(
		const GFM<index_t>& gfmFw,    // forward index
		const GFM<index_t>& gfmBw,    // mirror index
		SideLocus<index_t>&  tloc,    // top locus
		SideLocus<index_t>&  bloc,    // bot locus
		index_t topf,                 // top in BWT
		index_t botf,                 // bot in BWT
		index_t topb,                 // top in BWT'
		index_t botb);                // bot in BWT'

	/**
	 * Advance this descent by following read matches as far as possible.
	 */
    bool followMatches(
        const Read& q,     // query string
		const Scoring& sc,         // scoring scheme
        const GFM<index_t>& gfmFw,    // forward index
        const GFM<index_t>& gfmBw,    // mirror index
		DescentRedundancyChecker& re, // redundancy checker
        EFactory<Descent>& df,     // Descent factory
        EFactory<DescentPos>& pf,  // DescentPos factory
        const EList<DescentRoot>& rs,   // roots
        const EList<DescentConfig>& cs, // configs
        EHeap<TDescentPair>& heap, // heap
        DescentAlignmentSink<index_t>& alsink, // alignment sink
        DescentMetrics& met,       // metrics
		PerReadMetrics& prm,       // per-read metrics
        bool& branches,            // out: true -> there are > 0 ways to branch
        bool& hitEnd,              // out: true -> hit read end with non-empty range
        bool& done,                // out: true -> we made a full alignment
        TReadOff& off5p_i,         // out: initial 5' offset
        TIndexOffU& topf_bounce,    // out: top of SA range for fw idx for bounce
        TIndexOffU& botf_bounce,    // out: bot of SA range for fw idx for bounce
        TIndexOffU& topb_bounce,    // out: top of SA range for bw idx for bounce
        TIndexOffU& botb_bounce);   // out: bot of SA range for bw idx for bounce

	/**
	 * Recalculate our summary of the outgoing edges from this descent.  When
	 * deciding what outgoing edges are legal, we abide by constraints.
	 * Typically, they limit the total of the penalties accumulated so far, as
	 * a function of distance from the search root.  E.g. a constraint might
	 * disallow any gaps or mismatches within 20 ply of the search root, then
	 * allow 1 mismatch within 30 ply, then allow up to 1 mismatch or 1 gap
	 * within 40 ply, etc.
	 */
	size_t recalcOutgoing(
		const Read& q,                   // query string
		const Scoring& sc,               // scoring scheme
		TAlScore minsc,                  // minimum score
		TAlScore maxpen,                 // maximum penalty
		DescentRedundancyChecker& re,    // redundancy checker
		EFactory<DescentPos>& pf,        // factory with DescentPoss
        const EList<DescentRoot>& rs,    // roots
        const EList<DescentConfig>& cs,  // configs
		PerReadMetrics& prm);            // per-read metrics

    TRootId         rid_;         // root id

	TReadOff        al5pi_;       // lo offset from 5' end of aligned read char
	TReadOff        al5pf_;       // hi offset from 5' end of aligned read char
	bool            l2r_;         // left-to-right?
	int             gapadd_;      // net ref characters additional
    TReadOff        off5p_i_;     // offset we started out at for this descent

	TIndexOffU       topf_, botf_; // incoming SA range w/r/t forward index
	TIndexOffU       topb_, botb_; // incoming SA range w/r/t forward index

	size_t          descid_;      // ID of this descent
	TDescentId      parent_;      // ID of parent descent
	TScore          pen_;         // total penalties accumulated so far
	size_t          posid_;       // ID of 1st elt of the DescentPos factory w/
	                              // descent pos info for this descent
	size_t          len_;         // length of stretch of matches
	DescentOutgoing out_;         // summary of outgoing edges
	Edit            edit_;        // edit joining this descent with parent
	bool            lastRecalc_;  // set by recalcOutgoing if out edges empty
};

/**
 * An alignment result from a Descent.
 */
struct DescentAlignment {

	DescentAlignment() { reset(); }

	/**
	 * Reset DescentAlignment to be uninitialized.
	 */
	void reset() {
		topf = botf = 0;
		pen = 0;
		fw = false;
		ei = en = 0;
	}

	/**
	 * Initialize this DescentAlignment.
	 */
	void init(
		TScore pen_,
		bool fw_,
		TIndexOffU topf_,
		TIndexOffU botf_,
		size_t ei_,
		size_t en_)
	{
		assert_gt(botf_, topf_);
		pen = pen_;
		fw = fw_;
		topf = topf_;
		botf = botf_;
		ei = ei_;
		en = en_;
	}
	
	/**
	 * Return true iff DescentAlignment is initialized.
	 */
	bool inited() const {
		return botf > topf;
	}
	
	/**
	 * Return true iff the alignment is perfect (has no edits)
	 */
	bool perfect() const {
		return pen == 0;
	}
	
	/**
	 * Return the number of elements in this range.
	 */
	size_t size() const {
		return botf - topf;
	}

	TScore pen; // score
	
	bool fw; // forward or revcomp aligned?

	TIndexOffU topf; // top in forward index
	TIndexOffU botf; // bot in forward index

	size_t ei; // First edit in DescentAlignmentSink::edits_ involved in aln
	size_t en; // # edits in DescentAlignmentSink::edits_ involved in aln
};

/**
 * A partial alignment result from a Descent where the reference offset has
 * been resolved.
 */
struct DescentPartialResolvedAlignment {

	DescentPartialResolvedAlignment() { reset(); }

	/**
	 * Reset DescentAlignment to be uninitialized.
	 */
	void reset() {
		topf = botf = 0;
		pen = 0;
		fw = false;
		ei = en = 0;
		refcoord.reset();
	}

	/**
	 * Initialize this DescentAlignment.
	 */
	void init(
		TScore pen_,
		bool fw_,
		TIndexOffU topf_,
		TIndexOffU botf_,
		size_t ei_,
		size_t en_,
		const Coord& refcoord_)
	{
		assert_gt(botf_, topf_);
		pen = pen_;
		fw = fw_;
		topf = topf_;
		botf = botf_;
		ei = ei_;
		en = en_;
		refcoord = refcoord_;
	}
	
	/**
	 * Return true iff DescentAlignment is initialized.
	 */
	bool inited() const {
		return botf > topf;
	}
	
	/**
	 * Return the number of elements in this range.
	 */
	size_t size() const {
		return botf - topf;
	}

	TScore pen;     // score
	
	bool fw;        // forward or revcomp aligned?

	TIndexOffU topf; // top in forward index
	TIndexOffU botf; // bot in forward index

	size_t ei;      // First edit in DescentAlignmentSink::edits_ involved in aln
	size_t en;      // # edits in DescentAlignmentSink::edits_ involved in aln
	
	Coord refcoord; // reference coord of leftmost ref char involved
};

/**
 * Class that accepts alignments found during descent and maintains the state
 * required to dispense them to consumers in an appropriate order.
 *
 * As for order in which they are dispensed, in order to maintain uniform
 * distribution over equal-scoring alignments, a good policy may be not to
 * dispense alignments at a given score stratum until *all* alignments at that
 * stratum have been accumulated (i.e. until our best-first search has moved on
 * to a worse stratum).  This also has the advantage that, for each alignment,
 * we can also report the number of other alignments in that cost stratum.
 *
 * A lazier alternative is to assume that the order in which alignments in a
 * given stratum arrive is already pseudo-random, which frees us from having to
 * wait until the entire stratum has been explored.  But there is reason to
 * think that this order is not truly pseudo-random, since our root placement
 * and root priorities will tend to first lead us to alignments with certain
 * patterns of edits.
 */
template <typename index_t>
class DescentAlignmentSink {

public:

    /**
     * If this is the final descent in a complete end-to-end alignment, report
     * the alignment.
     */
    bool reportAlignment(
        const Read& q,           // query string
		const GFM<index_t>& gfmFw,       // forward index
		const GFM<index_t>& gfmBw,       // mirror index
		TIndexOffU topf,                  // SA range top in forward index
		TIndexOffU botf,                  // SA range bottom in forward index
		TIndexOffU topb,                  // SA range top in backward index
		TIndexOffU botb,                  // SA range bottom in backward index
        TDescentId id,                   // id of leaf Descent
		TRootId rid,                     // id of search root
        const Edit& e,                   // final edit, if needed
        TScore pen,                      // total penalty
        EFactory<Descent<index_t> >& df, // factory with Descent
        EFactory<DescentPos>& pf,        // factory with DescentPoss
        const EList<DescentRoot>& rs,    // roots
        const EList<DescentConfig>& cs); // configs
    
    /**
     * Reset to uninitialized state.
     */
    void reset() {
		edits_.clear();
		als_.clear();
		lhs_.clear();
		rhs_.clear();
		nelt_ = 0;
		bestPen_ = worstPen_ = std::numeric_limits<TAlScore>::max();
    }

	/**
	 * Return the total size occupued by the Descent driver and all its
	 * constituent parts.
	 */
	size_t totalSizeBytes() const {
		return edits_.totalSizeBytes() +
		       als_.totalSizeBytes() +
			   lhs_.totalSizeBytes() +
			   rhs_.totalSizeBytes() +
			   sizeof(size_t);
	}

	/**
	 * Return the total capacity of the Descent driver and all its constituent
	 * parts.
	 */
	size_t totalCapacityBytes() const {
		return edits_.totalCapacityBytes() +
		       als_.totalCapacityBytes() +
			   lhs_.totalCapacityBytes() +
			   rhs_.totalCapacityBytes() +
			   sizeof(size_t);
	}
	
	/**
	 * Return the number of SA ranges involved in hits.
	 */
	size_t nrange() const {
		return als_.size();
	}

	/**
	 * Return the number of SA elements involved in hits.
	 */
	size_t nelt() const {
		return nelt_;
	}
	
	/**
	 * The caller provides 'i', which is an offset of a particular element in
	 * one of the SA ranges in the current stratum.  This function returns, in
	 * 'al' and 'off', information about the element in terms of the range it's
	 * part of and its offset into that range.
	 */
	void elt(size_t i, DescentAlignment& al, size_t& ri, size_t& off) const {
		assert_lt(i, nelt());
		for(size_t j = 0; j < als_.size(); j++) {
			if(i < als_[j].size()) {
				al = als_[j];
				ri = j;
				off = i;
				return;
			}
			i -= als_[j].size();
		}
		assert(false);
	}
	
	/**
	 * Get a particular alignment.
	 */
	const DescentAlignment& operator[](size_t i) const {
		return als_[i];
	}

	/**
	 * Return true iff (a) we found an alignment since the sink was initialized
	 * or since the last time advanceStratum() was called, and (b) the penalty
	 * associated with the current-best task on the heap ('best') is worse
	 * (higher) than the penalty associated with the alignments found most
	 * recently (worstPen_).
	 */
	bool stratumDone(TAlScore bestPen) const {
		if(nelt_ > 0 && bestPen > worstPen_) {
			return true;
		}
		return false;
	}
	
	/**
	 * The alignment consumer calls this to indicate that they are done with
	 * all the alignments in the current best non-empty stratum.  We can
	 * therefore mark all those alignments as "reported" and start collecting
	 * results for the next stratum.
	 */
	void advanceStratum() {
		assert_gt(nelt_, 0);
		edits_.clear();
		als_.clear();
		// Don't reset lhs_ or rhs_
		nelt_ = 0;
		bestPen_ = worstPen_ = std::numeric_limits<TAlScore>::max();
	}
	
#ifndef NDEBUG
	/**
	 * Check that alignment sink is internally consistent.
	 */
	bool repOk() const {
		assert_geq(nelt_, als_.size());
		for(size_t i = 1; i < als_.size(); i++) {
			assert_geq(als_[i].pen, als_[i-1].pen);
		}
		assert(bestPen_ == std::numeric_limits<TAlScore>::max() || worstPen_ >= bestPen_);
		return true;
	}
#endif
	
	TAlScore bestPenalty() const { return bestPen_; }
	TAlScore worstPenalty() const { return worstPen_; }

	size_t editsSize() const { return edits_.size(); }
	size_t alsSize() const { return als_.size(); }
	size_t lhsSize() const { return lhs_.size(); }
	size_t rhsSize() const { return rhs_.size(); }
	
	const EList<Edit>& edits() const { return edits_; }

protected:

	EList<Edit> edits_;
	EList<DescentAlignment> als_;
	ESet<Triple<TIndexOffU, TIndexOffU, size_t> > lhs_;
	ESet<Triple<TIndexOffU, TIndexOffU, size_t> > rhs_;
	size_t nelt_;
	TAlScore bestPen_;  // best (smallest) penalty among as-yet-unreported alns
	TAlScore worstPen_; // worst (greatest) penalty among as-yet-unreported alns
#ifndef NDEBUG
	BTDnaString tmprfdnastr_;
#endif

};

/**
 * Class that aggregates partial alignments taken from a snapshot of the
 * DescentDriver heap.
 */
class DescentPartialResolvedAlignmentSink {

public:
   
    /**
     * Reset to uninitialized state.
     */
    void reset() {
		edits_.clear();
		als_.clear();
		nelt_ = 0;
		bestPen_ = worstPen_ = std::numeric_limits<TAlScore>::max();
    }

	/**
	 * Return the total size occupued by the Descent driver and all its
	 * constituent parts.
	 */
	size_t totalSizeBytes() const {
		return edits_.totalSizeBytes() +
		       als_.totalSizeBytes() +
			   sizeof(size_t);
	}

	/**
	 * Return the total capacity of the Descent driver and all its constituent
	 * parts.
	 */
	size_t totalCapacityBytes() const {
		return edits_.totalCapacityBytes() +
		       als_.totalCapacityBytes() +
			   sizeof(size_t);
	}
	
	/**
	 * Return the number of SA ranges involved in hits.
	 */
	size_t nrange() const {
		return als_.size();
	}

	/**
	 * Return the number of SA elements involved in hits.
	 */
	size_t nelt() const {
		return nelt_;
	}
	
	/**
	 * The caller provides 'i', which is an offset of a particular element in
	 * one of the SA ranges in the current stratum.  This function returns, in
	 * 'al' and 'off', information about the element in terms of the range it's
	 * part of and its offset into that range.
	 */
	void elt(size_t i, DescentPartialResolvedAlignment& al, size_t& ri, size_t& off) const {
		assert_lt(i, nelt());
		for(size_t j = 0; j < als_.size(); j++) {
			if(i < als_[j].size()) {
				al = als_[j];
				ri = j;
				off = i;
				return;
			}
			i -= als_[j].size();
		}
		assert(false);
	}
	
	/**
	 * Get a particular alignment.
	 */
	const DescentPartialResolvedAlignment& operator[](size_t i) const {
		return als_[i];
	}

	/**
	 * Return true iff (a) we found an alignment since the sink was initialized
	 * or since the last time advanceStratum() was called, and (b) the penalty
	 * associated with the current-best task on the heap ('best') is worse
	 * (higher) than the penalty associated with the alignments found most
	 * recently (worstPen_).
	 */
	bool stratumDone(TAlScore bestPen) const {
		if(nelt_ > 0 && bestPen > worstPen_) {
			return true;
		}
		return false;
	}
	
	/**
	 * The alignment consumer calls this to indicate that they are done with
	 * all the alignments in the current best non-empty stratum.  We can
	 * therefore mark all those alignments as "reported" and start collecting
	 * results for the next stratum.
	 */
	void advanceStratum() {
		assert_gt(nelt_, 0);
		edits_.clear();
		als_.clear();
		nelt_ = 0;
		bestPen_ = worstPen_ = std::numeric_limits<TAlScore>::max();
	}
	
#ifndef NDEBUG
	/**
	 * Check that partial alignment sink is internally consistent.
	 */
	bool repOk() const {
		assert_geq(nelt_, als_.size());
		//for(size_t i = 1; i < als_.size(); i++) {
		//	assert_geq(als_[i].pen, als_[i-1].pen);
		//}
		assert(bestPen_ == std::numeric_limits<TAlScore>::max() || worstPen_ >= bestPen_);
		return true;
	}
#endif
	
	TAlScore bestPenalty() const { return bestPen_; }
	TAlScore worstPenalty() const { return worstPen_; }

	size_t editsSize() const { return edits_.size(); }
	size_t alsSize() const { return als_.size(); }
	
	const EList<Edit>& edits() const { return edits_; }

protected:

	EList<Edit> edits_;
	EList<DescentPartialResolvedAlignment> als_;
	size_t nelt_;
	TAlScore bestPen_;  // best (smallest) penalty among as-yet-unreported alns
	TAlScore worstPen_; // worst (greatest) penalty among as-yet-unreported alns
};

/**
 * Abstract parent for classes that select descent roots and descent
 * configurations given information about the read.
 */
class DescentRootSelector {

public:

	virtual ~DescentRootSelector() { }

	virtual void select(
		const Read& q,          // read that we're selecting roots for
		const Read* qo,         // opposite mate, if applicable
		bool nofw,              // don't add roots for fw read
		bool norc,              // don't add roots for rc read
		EList<DescentConfig>& confs,    // put DescentConfigs here
		EList<DescentRoot>& roots) = 0; // put DescentRoot here
};

/**
 * Encapsulates a set of conditions governing when the DescentDriver should
 * stop.
 */
struct DescentStoppingConditions {

	DescentStoppingConditions() { reset(); }

	DescentStoppingConditions(
		size_t totsz_,
		size_t nfound_,
		bool stra_,
		size_t nbwop_)
	{
		init(totsz_, nfound_, stra_, nbwop_);
	}
	
	/**
	 * Reset to uninitialized state.
	 */
	void reset() {
		totsz = nfound = nbwop = std::numeric_limits<size_t>::max();
		stra = false;
		assert(!inited());
	}

	/**
	 * Initialize this DescentStoppingConditions.
	 */
	void init(
		size_t totsz_,
		size_t nfound_,
		bool stra_,
		size_t nbwop_)
	{
		totsz = totsz_;
		nfound = nfound_;
		stra = stra_;
		nbwop = nbwop_;
		assert(inited());
	}
	
	/**
	 * Return true iff this instance is initialized.
	 */
	bool inited() const {
		return totsz != std::numeric_limits<size_t>::max();
	}

	size_t totsz;  // total size of all the expandable data structures in bytes
	size_t nfound; // # alignments found
	bool stra;     // stop after each non-empty stratum
	size_t nbwop;  // # Burrows-Wheeler (rank) operations performed
};

enum {
	DESCENT_DRIVER_ALN = 1,
	DESCENT_DRIVER_STRATA = 2,
	DESCENT_DRIVER_MEM = 4,
	DESCENT_DRIVER_BWOPS = 8,
	DESCENT_DRIVER_DONE = 16
};

/**
 * Class responsible for advancing all the descents.  The initial descents may
 * emanate from several different locations in the read.  Note that descents
 * may become redundant with each other, and should then be eliminated.
 */
template <typename index_t>
class DescentDriver {
public:

	DescentDriver(bool veryVerbose) :
		veryVerbose_(veryVerbose)
	{
		reset();
	}
	
	/**
	 * Initialize driver with respect to a new read.  If a DescentRootSelector
	 * is specified, then it is used to obtain roots as well.
	 */
	void initRead(
		const Read& q,
		bool nofw,
		bool norc,
		TAlScore minsc,
		TAlScore maxpen,
		const Read* qu = NULL,
		DescentRootSelector *sel = NULL)
	{
		reset();
		q_ = q;
		minsc_ = minsc;
		maxpen_ = maxpen;
		if(sel != NULL) {
			sel->select(q_, qu, nofw, norc, confs_, roots_);
		}
		re_.init(q.length());
	}
	
	/**
	 * Add a new search root, which might (a) prefer to move in a left-to-right
	 * direction, and might (b) be with respect to the read or its reverse
	 * complement.
	 */
	void addRoot(
        const DescentConfig& conf,
        TReadOff off,
        bool l2r,
        bool fw,
        float pri)
    {
        confs_.push_back(conf);
		assert_lt(off, q_.length());
		if(l2r && off == q_.length()-1) {
			l2r = !l2r;
		} else if(!l2r && off == 0) {
			l2r = !l2r;
		}
		roots_.push_back(DescentRoot(off, l2r, fw, q_.length(), pri));
	}
	
	/**
	 * Clear out the DescentRoots currently configured.
	 */
	void clearRoots() {
		confs_.clear();
		roots_.clear();
	}
	
	/**
	 * Clear the Descent driver so that we're ready to re-start seed alignment
	 * for the current read.
	 */
	void resetRead() {
		df_.clear();     // clear Descents
		assert_leq(df_.totalSizeBytes(), 100);
		pf_.clear();     // clear DescentPoss
		assert_leq(pf_.totalSizeBytes(), 100);
		heap_.clear();   // clear Heap
		assert_leq(heap_.totalSizeBytes(), 100);
		roots_.clear();  // clear roots
		assert_leq(roots_.totalSizeBytes(), 100);
		confs_.clear();  // clear confs
		assert_leq(confs_.totalSizeBytes(), 100);
        alsink_.reset(); // clear alignment sink
		assert_leq(alsink_.totalSizeBytes(), 100);
		re_.reset();
		assert_leq(re_.totalSizeBytes(), 100);
		rootsInited_ = 0; // haven't yet created initial descents
		curPen_ = 0;      //
	}
	
	/**
	 * Clear the Descent driver so that we're ready to re-start seed alignment
	 * for the current read.
	 */
	void reset() {
		resetRead();
	}

	/**
	 * Perform seed alignment.
	 */
	void go(
        const Scoring& sc,    // scoring scheme
		const GFM<index_t>& gfmFw,   // forward index
		const GFM<index_t>& gfmBw,   // mirror index
        DescentMetrics& met,  // metrics
		PerReadMetrics& prm); // per-read metrics

	/**
	 * Perform seed alignment until some stopping condition is satisfied.
	 */
	int advance(
		const DescentStoppingConditions& stopc, // stopping conditions
        const Scoring& sc,    // scoring scheme
		const GFM<index_t>& gfmFw,   // forward index
		const GFM<index_t>& gfmBw,   // mirror index
        DescentMetrics& met,  // metrics
		PerReadMetrics& prm); // per-read metrics

#ifndef NDEBUG
	/**
	 * Return true iff this DescentDriver is well formed.  Throw an assertion
	 * otherwise.
	 */
	bool repOk() const {
		return true;
	}
#endif

	/**
	 * Return the number of end-to-end alignments reported.
	 */
	size_t numAlignments() const {
		return alsink_.nelt();
	}
	
	/**
	 * Return the associated DescentAlignmentSink object.
	 */
	const DescentAlignmentSink<index_t>& sink() const {
		return alsink_;
	}

	/**
	 * Return the associated DescentAlignmentSink object.
	 */
	DescentAlignmentSink<index_t>& sink() {
		return alsink_;
	}
	
	/**
	 * Return the total size occupued by the Descent driver and all its
	 * constituent parts.
	 */
	size_t totalSizeBytes() const {
		return df_.totalSizeBytes() +
		       pf_.totalSizeBytes() +
			   heap_.totalSizeBytes() +
			   roots_.totalSizeBytes() +
			   confs_.totalSizeBytes() +
		       alsink_.totalSizeBytes() +
			   re_.totalSizeBytes();
	}

	/**
	 * Return the total capacity of the Descent driver and all its constituent
	 * parts.
	 */
	size_t totalCapacityBytes() const {
		return df_.totalCapacityBytes() +
		       pf_.totalCapacityBytes() +
			   heap_.totalCapacityBytes() +
			   roots_.totalCapacityBytes() +
			   confs_.totalCapacityBytes() +
		       alsink_.totalCapacityBytes() +
			   re_.totalCapacityBytes();
	}
	
	/**
	 * Return a const ref to the query.
	 */
	const Read& query() const {
		return q_;
	}
	
	/**
	 * Return the minimum score that must be achieved by an alignment in order
	 * for it to be considered "valid".
	 */
	TAlScore minScore() const {
		return minsc_;
	}

protected:

	Read                        q_;      // query nucleotide and quality strings
	TAlScore                    minsc_;  // minimum score
	TAlScore                    maxpen_; // maximum penalty
	EFactory<Descent<index_t> > df_;     // factory holding all the Descents, which
	                              // must be referred to by ID
	EFactory<DescentPos> pf_;     // factory holding all the DescentPoss, which
	                              // must be referred to by ID
	EList<DescentRoot>   roots_;  // search roots
    EList<DescentConfig> confs_;  // configuration params for each root
	size_t rootsInited_;          // # initial Descents already created
	EHeap<TDescentPair>  heap_;   // priority queue of Descents
    DescentAlignmentSink<index_t> alsink_; // alignment sink
	DescentRedundancyChecker re_; // redundancy checker
	TAlScore             curPen_; // current penalty
	bool veryVerbose_;            // print lots of partial alignments

	EList<Edit> tmpedit_;
	BTDnaString tmprfdnastr_;
};

/**
 * Selects alignments to report from a complete non-empty stratum of
 * alignments stored in the DescentAlignmentSink.
 */
template <typename index_t>
class DescentAlignmentSelector {

public:

	DescentAlignmentSelector() : gwstate_(GW_CAT) { reset(); }

	/**
	 * Initialize a new selector w/r/t a DescentAlignmentSink holding a
	 * non-empty alignment stratum.
	 */
	void init(
		const Read& q,
		const DescentAlignmentSink<index_t>& sink,
		const GFM<index_t>& gfmFw, // forward Bowtie index for walking left
		const BitPairReference& ref, // bitpair-encoded reference
		RandomSource& rnd,           // pseudo-random generator for sampling rows
		WalkMetrics& met)
	{
		// We're going to sample from space of *alignments*, not ranges.  So
		// when we extract a sample, we'll have to do a little extra work to
		// convert it to a <range, offset> coordinate.
		rnd_.init(
			sink.nelt(), // # elements to choose from
			true);       // without replacement
		offs_.resize(sink.nelt());
		offs_.fill(std::numeric_limits<TIndexOffU>::max());
		sas_.resize(sink.nrange());
		gws_.resize(sink.nrange());
		size_t ei = 0;
		for(size_t i = 0; i < sas_.size(); i++) {
			size_t en = sink[i].botf - sink[i].topf;
			sas_[i].init(sink[i].topf, q.length(), EListSlice<TIndexOffU, 16>(offs_, ei, en));
			gws_[i].init(gfmFw, ref, sas_[i], rnd, met);
			ei += en;
		}
	}
	
	/**
	 * Reset the selector.
	 */
	void reset() {
		rnd_.reset();
	}
	
	/**
	 * Return true iff the selector is currently initialized.
	 */
	bool inited() const {
		return rnd_.size() > 0;
	}
	
	/**
	 * Get next alignment and convert it to an AlnRes.
	 */
	bool next(
		const DescentDriver<index_t>& dr,
		const GFM<index_t>& gfmFw, // forward Bowtie index for walking left
		const BitPairReference& ref, // bitpair-encoded reference
		RandomSource& rnd,
		AlnRes& rs,
		WalkMetrics& met,
		PerReadMetrics& prm)
	{
		// Sample one alignment randomly from pool of remaining alignments
		size_t ri = (size_t)rnd_.next(rnd);
		size_t off = 0;
		DescentAlignment al;
		size_t rangei = 0;
		// Convert random alignment index into a <range, offset> coordinate
		dr.sink().elt(ri, al, rangei, off);
		assert_lt(off, al.size());
		Coord refcoord;
		WalkResult<index_t> wr;
		TIndexOffU tidx = 0, toff = 0, tlen = 0;
		gws_[rangei].advanceElement(
			(TIndexOffU)off,
			gfmFw,        // forward Bowtie index for walking left
			ref,          // bitpair-encoded reference
			sas_[rangei], // SA range with offsets
			gwstate_,     // GroupWalk state; scratch space
			wr,           // put the result here
			met,          // metrics
			prm);         // per-read metrics
		assert_neq(OFF_MASK, wr.toff);
		bool straddled = false;
		gfmFw.joinedToTextOff(
			wr.elt.len,
			wr.toff,
			tidx,
			toff,
			tlen,
			true,        // reject straddlers?
			straddled);  // straddled?
		if(tidx == OFF_MASK) {
			// The seed hit straddled a reference boundary so the seed
			// hit isn't valid
			return false;
		}
		// Coordinate of the seed hit w/r/t the pasted reference string
		refcoord.init(tidx, (int64_t)toff, dr.sink()[rangei].fw);
		const EList<Edit>& edits = dr.sink().edits();
		size_t ns = 0, ngap = 0, nrefn = 0;
		for(size_t i = al.ei; i < al.ei + al.en; i++) {
			if(edits[i].qchr == 'N' || edits[i].chr == 'N') ns++;
			if(edits[i].chr == 'N') nrefn++;
			if(edits[i].isGap()) ngap++;
		}
		AlnScore asc(
			-dr.sink().bestPenalty(),  // numeric score
			ns,                        // # Ns
			ngap);                     // # gaps
		rs.init(
			dr.query().length(),       // # chars after hard trimming
			asc,                       // alignment score
			&dr.sink().edits(),        // nucleotide edits array
			al.ei,                     // nucleotide edits first pos
			al.en,                     // nucleotide edits last pos
			NULL,                      // ambig base array
			0,                         // ambig base first pos
			0,                         // ambig base last pos
			refcoord,                  // coord of leftmost aligned char in ref
			tlen,                      // length of reference aligned to
			-1,                        // # seed mms allowed
			-1,                        // seed length
			-1,                        // seed interval
			dr.minScore(),             // minimum score for valid alignment
			-1,                        // nuc5p (for colorspace)
			-1,                        // nuc3p (for colorspace)
			false,                     // soft pre-trimming?
			0,                         // 5p pre-trimming
			0,                         // 3p pre-trimming
			false,                     // soft trimming?
			0,                         // 5p trimming
			0);                        // 3p trimming
		rs.setRefNs(nrefn);
		return true;
	}
	
	/**
	 * Return true iff all elements have been reported.
	 */
	bool done() const {
		return rnd_.done();
	}

	/**
	 * Return the total size occupued by the Descent driver and all its
	 * constituent parts.
	 */
	size_t totalSizeBytes() const {
		return rnd_.totalSizeBytes() +
		       offs_.totalSizeBytes() +
			   sas_.totalSizeBytes() +
			   gws_.totalSizeBytes();
	}

	/**
	 * Return the total capacity of the Descent driver and all its constituent
	 * parts.
	 */
	size_t totalCapacityBytes() const {
		return rnd_.totalCapacityBytes() +
		       offs_.totalCapacityBytes() +
			   sas_.totalCapacityBytes() +
			   gws_.totalCapacityBytes();
	}
	
protected:

	Random1toN rnd_;
	EList<TIndexOffU, 16> offs_;
	EList<SARangeWithOffs<EListSlice<TIndexOffU, 16>, index_t> > sas_;
	EList<GroupWalk2S<index_t, EListSlice<TIndexOffU, 16>, 16> > gws_;
	GroupWalkState<index_t> gwstate_;
};

/**
 * Selects and prioritizes partial alignments from the heap of the
 * DescentDriver.  We assume that the heap is no longer changing (i.e. that the
 * DescentDriver is done).  Usually, the user will then attempt to extend the
 * partial alignments into full alignments.  This can happen incrementally;
 * that is, the user might ask for the partial alignments one "batch" at a
 * time, and the selector will only do as much work is necessary to supply each
 * requesteded batch.
 *
 * The actual work done here includes: (a) scanning the heap for high-priority
 * partial alignments, (b) setting up the rnd_, offs_, sas_, gws_, and gwstate_
 * fields and resolving offsets of partial alignments, (c) packaging and
 * delivering batches of results to the caller.
 *
 * How to prioritize partial alignments?  One idea is to use the same
 * penalty-based prioritization used in the heap.  This has pros: (a) maintains
 * the guarantee that we're visiting alignments in best-to-worst order in
 * end-to-end alignment mode, (b) the heap is already prioritized this way, so
 * it's easier for us to compile high-priority partial alignments.  But the con
 * is that it doesn't take depth into account, which could mean that we're
 * extending a lot of very short partial alignments first.
 *
 * A problem we should keep in mind is that some 
 */
template <typename index_t>
class DescentPartialAlignmentSelector {

public:

	DescentPartialAlignmentSelector() : gwstate_(GW_CAT) { reset(); }

	/**
	 * Initialize a new selector w/r/t a read, index and heap of partial
	 * alignments.
	 */
	void init(
		const Read& q,                   // read
		const EHeap<TDescentPair>& heap, // the heap w/ the partial alns
		TAlScore depthBonus,             // use depth when prioritizing
		size_t nbatch,                   // # of alignments in a batch
		const GFM<index_t>& gfmFw,       // forward Bowtie index for walk-left
		const BitPairReference& ref,     // bitpair-encoded reference
		RandomSource& rnd,               // pseudo-randoms for sampling rows
		WalkMetrics& met)                // metrics re: offset resolution
	{
		// Make our internal heap
		if(depthBonus > 0) {
			heap_.clear();
			for(size_t i = 0; i < heap.size(); i++) {
				TDescentPair p = heap[i];
				p.first.pen += depthBonus * p.first.depth;
				heap_.insert(p);
			}
		} else {
			heap_ = heap;
		}
#if 0
		// We're going to sample from space of *alignments*, not ranges.  So
		// when we extract a sample, we'll have to do a little extra work to
		// convert it to a <range, offset> coordinate.
		rnd_.init(
			sink.nelt(), // # elements to choose from
			true);       // without replacement
		offs_.resize(sink.nelt());
		offs_.fill(std::numeric_limits<TIndexOff>::max());
		sas_.resize(sink.nrange());
		gws_.resize(sink.nrange());
		size_t ei = 0;
		for(size_t i = 0; i < sas_.size(); i++) {
			size_t en = sink[i].botf - sink[i].topf;
			sas_[i].init(sink[i].topf, q.length(), EListSlice<TIndexOff, 16>(offs_, ei, en));
			gws_[i].init(gfmFw, ref, sas_[i], rnd, met);
			ei += en;
		}
#endif
	}
	
	/**
	 *
	 */
	void compileBatch() {
	}
	
	/**
	 * Reset the selector.
	 */
	void reset() {
		heap_.clear();
	}
	
	/**
	 * Return true iff the selector is currently initialized.
	 */
	bool inited() const {
		return !heap_.empty();
	}
	
	/**
	 * Get next alignment and convert it to an AlnRes.
	 */
	bool next(
		const DescentDriver<index_t>& dr,
		const GFM<index_t>& gfmFw,          // forward Bowtie index for walking left
		const BitPairReference& ref, // bitpair-encoded reference
		RandomSource& rnd,
		AlnRes& rs,
		WalkMetrics& met,
		PerReadMetrics& prm)
	{
		// Sample one alignment randomly from pool of remaining alignments
		size_t ri = (size_t)rnd_.next(rnd);
		size_t off = 0;
		DescentAlignment al;
		size_t rangei = 0;
		// Convert random alignment index into a <range, offset> coordinate
		dr.sink().elt(ri, al, rangei, off);
		assert_lt(off, al.size());
		Coord refcoord;
		WalkResult<index_t> wr;
		uint32_t tidx = 0, toff = 0, tlen = 0;
		gws_[rangei].advanceElement(
			(uint32_t)off,
			gfmFw,        // forward Bowtie index for walking left
			ref,          // bitpair-encoded reference
			sas_[rangei], // SA range with offsets
			gwstate_,     // GroupWalk state; scratch space
			wr,           // put the result here
			met,          // metrics
			prm);         // per-read metrics
		assert_neq(0xffffffff, wr.toff);
		bool straddled = false;
		gfmFw.joinedToTextOff(
			wr.elt.len,
			wr.toff,
			tidx,
			toff,
			tlen,
			true,        // reject straddlers?
			straddled);  // straddled?
		if(tidx == 0xffffffff) {
			// The seed hit straddled a reference boundary so the seed
			// hit isn't valid
			return false;
		}
		// Coordinate of the seed hit w/r/t the pasted reference string
		refcoord.init(tidx, (int64_t)toff, dr.sink()[rangei].fw);
		const EList<Edit>& edits = dr.sink().edits();
		size_t ns = 0, ngap = 0, nrefn = 0;
		for(size_t i = al.ei; i < al.ei + al.en; i++) {
			if(edits[i].qchr == 'N' || edits[i].chr == 'N') ns++;
			if(edits[i].chr == 'N') nrefn++;
			if(edits[i].isGap()) ngap++;
		}
		return true;
	}
	
	/**
	 * Return true iff all elements have been reported.
	 */
	bool done() const {
		return rnd_.done();
	}

	/**
	 * Return the total size occupued by the Descent driver and all its
	 * constituent parts.
	 */
	size_t totalSizeBytes() const {
		return heap_.totalSizeBytes() +
		       rnd_.totalSizeBytes() +
		       offs_.totalSizeBytes() +
			   sas_.totalSizeBytes() +
			   gws_.totalSizeBytes();
	}

	/**
	 * Return the total capacity of the Descent driver and all its constituent
	 * parts.
	 */
	size_t totalCapacityBytes() const {
		return heap_.totalCapacityBytes() +
		       rnd_.totalCapacityBytes() +
		       offs_.totalCapacityBytes() +
			   sas_.totalCapacityBytes() +
			   gws_.totalCapacityBytes();
	}
	
protected:

	// This class's working heap.  This might simply be a copy of the original
	// heap, or it might be re-prioritized in some way.
	EHeap<TDescentPair> heap_;

	Random1toN rnd_;
	EList<TIndexOff, 16> offs_;
	EList<SARangeWithOffs<EListSlice<TIndexOff, 16>, index_t> > sas_;
	EList<GroupWalk2S<index_t, EListSlice<TIndexOff, 16>, 16> > gws_;
	GroupWalkState<index_t> gwstate_;
};

/**
 * Drive the process of descending from all search roots.
 */
template <typename index_t>
void DescentDriver<index_t>::go(
								const Scoring& sc,    // scoring scheme
								const GFM<index_t>& gfmFw,   // forward index
								const GFM<index_t>& gfmBw,   // mirror index
								DescentMetrics& met,  // metrics
								PerReadMetrics& prm)  // per-read metrics
{
	assert(q_.repOk());
    // Convert DescentRoots to the initial Descents
    for(size_t i = 0; i < roots_.size(); i++) {
        size_t dfsz = df_.size();
        size_t pfsz = pf_.size();
        TDescentId id = df_.alloc();
        Edit e_null;
        assert(!e_null.inited());
        bool succ = df_[id].init(
								 q_,        // read
								 i,         // root and conf id
								 sc,        // scoring scheme
								 minsc_,    // minimum score
								 maxpen_,   // maximum penalty
								 id,        // new Descent's id
								 gfmFw,     // forward index
								 gfmBw,     // mirror index
								 re_,       // redundancy checker
								 df_,       // Descent factory
								 pf_,       // DescentPos factory
								 roots_,    // DescentRoots
								 confs_,    // DescentConfs
								 heap_,     // heap
								 alsink_,   // alignment sink
								 met,       // metrics
								 prm);      // per-read metrics
		if(veryVerbose_) {
			bool fw = roots_[i].fw;
			tmpedit_.clear();
			df_[id].print(
						  &cerr,
						  "",
						  q_,
						  0,
						  0,
						  fw,
						  tmpedit_,
						  0,
						  tmpedit_.size(),
						  tmprfdnastr_);
		}
        if(!succ) {
            // Reclaim memory we had used for this descent and its DescentPos info
            df_.resize(dfsz);
            pf_.resize(pfsz);
        }
    }
    // Advance until some stopping condition
    bool stop = heap_.empty();
    while(!stop) {
		// Pop off the highest-priority descent.  Note that some outgoing edges
		// might have since been explored, which could reduce the priority of
		// the descent once we .
        TDescentPair p = heap_.pop();
		df_.alloc(); df_.pop();
        df_[p.second].followBestOutgoing(
										 q_,        // read
										 gfmFw,     // index over text
										 gfmBw,     // index over reverse text
										 sc,        // scoring scheme
										 minsc_,    // minimum score
										 maxpen_,   // maximum penalty
										 re_,       // redundancy checker
										 df_,       // Descent factory
										 pf_,       // DescentPos factory
										 roots_,    //
										 confs_,    //
										 heap_,     // priority queue for Descents
										 alsink_,   // alignment sink
										 met,       // metrics
										 prm);      // per-read metrics
        stop = heap_.empty();
    }
}

/**
 * Perform seed alignment until some stopping condition is satisfied.
 */
template <typename index_t>
int DescentDriver<index_t>::advance(
									const DescentStoppingConditions& stopc, // stopping conditions
									const Scoring& sc,    // scoring scheme
									const GFM<index_t>& gfmFw,   // forward index
									const GFM<index_t>& gfmBw,   // mirror index
									DescentMetrics& met,  // metrics
									PerReadMetrics& prm)  // per-read metrics
{
	size_t nbwop_i = met.bwops;
	while(rootsInited_ < roots_.size()) {
		size_t dfsz = df_.size();
		size_t pfsz = pf_.size();
		TDescentId id = df_.alloc();
		Edit e_null;
		assert(!e_null.inited());
        bool succ = df_[id].init(
								 q_,        // query
								 rootsInited_, // root and conf id
								 sc,        // scoring scheme
								 minsc_,    // minimum score
								 maxpen_,   // maximum penalty
								 id,        // new Descent's id
								 gfmFw,     // forward index
								 gfmBw,     // mirror index
								 re_,       // redundancy checker
								 df_,       // Descent factory
								 pf_,       // DescentPos factory
								 roots_,    // DescentRoots
								 confs_,    // DescentConfs
								 heap_,     // heap
								 alsink_,   // alignment sink
								 met,       // metrics
								 prm);      // per-read metrics
        if(!succ) {
            // Reclaim memory we had used for this descent and its DescentPos info
            df_.resize(dfsz);
            pf_.resize(pfsz);
        }
		rootsInited_++;
		TAlScore best = std::numeric_limits<TAlScore>::max();
		if(!heap_.empty()) {
			best = heap_.top().first.pen;
		}
		if(stopc.nfound > 0 && alsink_.nelt() > stopc.nfound) {
			return DESCENT_DRIVER_ALN;
		}
		if(alsink_.stratumDone(best)) {
			return DESCENT_DRIVER_STRATA;
		}
		if(stopc.nbwop > 0 && (met.bwops - nbwop_i) > stopc.nbwop) {
			return DESCENT_DRIVER_BWOPS;
		}
		if(stopc.totsz > 0 && totalSizeBytes() > stopc.totsz) {
			return DESCENT_DRIVER_MEM;
		}
    }
    // Advance until some stopping condition
    bool stop = heap_.empty();
    while(!stop) {
		// Pop off the highest-priority descent.  Note that some outgoing edges
		// might have since been explored, which could reduce the priority of
		// the descent once we .
        TDescentPair p = heap_.pop();
		df_.alloc(); df_.pop();
        df_[p.second].followBestOutgoing(
										 q_,
										 gfmFw,
										 gfmBw,
										 sc,
										 minsc_,    // minimum score
										 maxpen_,   // maximum penalty
										 re_,       // redundancy checker
										 df_,       // Descent factory
										 pf_,       // DescentPos factory
										 roots_,
										 confs_,
										 heap_,
										 alsink_,
										 met,
										 prm);      // per-read metrics
		TAlScore best = std::numeric_limits<TAlScore>::max();
		if(!heap_.empty()) {
			best = heap_.top().first.pen;
		}
		if(stopc.nfound > 0 && alsink_.nelt() > stopc.nfound) {
			return DESCENT_DRIVER_ALN;
		}
		if(alsink_.stratumDone(best)) {
			return DESCENT_DRIVER_STRATA;
		}
		if(stopc.nbwop > 0 && (met.bwops - nbwop_i) > stopc.nbwop) {
			return DESCENT_DRIVER_BWOPS;
		}
		if(stopc.totsz > 0 && totalSizeBytes() > stopc.totsz) {
			return DESCENT_DRIVER_MEM;
		}
        stop = heap_.empty();
    }
	return DESCENT_DRIVER_DONE;
}

/**
 * If this is the final descent in a complete end-to-end alignment, report
 * the alignment.
 */
template <typename index_t>
bool DescentAlignmentSink<index_t>::reportAlignment(
													const Read& q,                  // query string
													const GFM<index_t>& gfmFw,      // forward index
													const GFM<index_t>& gfmBw,      // mirror index
													TIndexOffU topf,                 // SA range top in forward index
													TIndexOffU botf,                 // SA range bottom in forward index
													TIndexOffU topb,                 // SA range top in backward index
													TIndexOffU botb,                 // SA range bottom in backward index
													TDescentId id,                  // id of leaf Descent
													TRootId rid,                    // id of search root
													const Edit& e,                  // final edit, if needed
													TScore pen,                     // total penalty
													EFactory<Descent<index_t> >& df,          // factory with Descent
													EFactory<DescentPos>& pf,       // factory with DescentPoss
													const EList<DescentRoot>& rs,   // roots
													const EList<DescentConfig>& cs) // configs
{
	TDescentId cur = id;
	ASSERT_ONLY(const Descent<index_t>& desc = df[id]);
	const bool fw = rs[rid].fw;
	ASSERT_ONLY(size_t len = q.length());
	assert(q.repOk());
	assert_lt(desc.al5pf(), len);
	// Adjust al5pi and al5pf to take the final edit into account (if
	// there is one)
	// Check if this is redundant with a previous reported alignment
	Triple<TIndexOffU, TIndexOffU, size_t> lhs(topf, botf, 0);
	Triple<TIndexOffU, TIndexOffU, size_t> rhs(topb, botb, q.length()-1);
	if(!lhs_.insert(lhs)) {
		rhs_.insert(rhs);
		return false; // Already there
	}
	if(!rhs_.insert(rhs)) {
		return false; // Already there
	}
	size_t ei = edits_.size();
	df[cur].collectEdits(edits_, &e, df);
	size_t en = edits_.size() - ei;
#ifndef NDEBUG
	{
		for(size_t i = 1; i < en; i++) {
			assert_geq(edits_[ei+i].pos, edits_[ei+i-1].pos);
		}
		// Now figure out how much we refrained from aligning on either
		// side.
		size_t trimLf = 0;
		size_t trimRg = 0;
		BTDnaString& rf = tmprfdnastr_;
		rf.clear();
		if(!fw) {
			// Edit offsets are w/r/t 5' end, but desc.print wants them w/r/t
			// the *left* end of the read sequence that aligned
			Edit::invertPoss(edits_, len, ei, en, true);
		}
		desc.print(NULL, "", q, trimLf, trimRg, fw, edits_, ei, en, rf);
		if(!fw) {
			// Invert them back to how they were before
			Edit::invertPoss(edits_, len, ei, en, true);
		}
		ASSERT_ONLY(TIndexOffU toptmp = 0);
		ASSERT_ONLY(TIndexOffU bottmp = 0);
		// Check that the edited string occurs in the reference
		if(!gfmFw.contains(rf, &toptmp, &bottmp)) {
			std::cerr << rf << std::endl;
			assert(false);
		}
	}
#endif
	als_.expand();
	als_.back().init(pen, fw, topf, botf, ei, en);
	nelt_ += (botf - topf);
	if(bestPen_ == std::numeric_limits<TAlScore>::max() || pen < bestPen_) {
		bestPen_ = pen;
	}
	if(worstPen_ == std::numeric_limits<TAlScore>::max() || pen > worstPen_) {
		worstPen_ = pen;
	}
	return true;
}

/**
 * Initialize a new descent branching from the given descent via the given
 * edit.  Return false if the Descent has no outgoing edges (and can
 * therefore have its memory freed), true otherwise.
 */
template <typename index_t>
bool Descent<index_t>::init(
							const Read& q,                  // query
							TRootId rid,                    // root id
							const Scoring& sc,              // scoring scheme
							TAlScore minsc,                 // minimum score
							TAlScore maxpen,                // maximum penalty
							TReadOff al5pi,                 // offset from 5' of 1st aligned char
							TReadOff al5pf,                 // offset from 5' of last aligned char
							TIndexOffU topf,                 // SA range top in FW index
							TIndexOffU botf,                 // SA range bottom in FW index
							TIndexOffU topb,                 // SA range top in BW index
							TIndexOffU botb,                 // SA range bottom in BW index
							bool l2r,                       // direction this descent will go in
							size_t descid,                  // my ID
							TDescentId parent,              // parent ID
							TScore pen,                     // total penalties so far
							const Edit& e,                  // edit for incoming edge
							const GFM<index_t>& gfmFw,      // forward index
							const GFM<index_t>& gfmBw,      // mirror index
							DescentRedundancyChecker& re,   // redundancy checker
							EFactory<Descent>& df,          // Descent factory
							EFactory<DescentPos>& pf,       // DescentPos factory
							const EList<DescentRoot>& rs,   // roots
							const EList<DescentConfig>& cs, // configs
							EHeap<TDescentPair>& heap,      // heap
							DescentAlignmentSink<index_t>& alsink,   // alignment sink
							DescentMetrics& met,            // metrics
							PerReadMetrics& prm)            // per-read metrics
{
	assert(q.repOk());
    rid_ = rid;
    al5pi_ = al5pi;
    al5pf_ = al5pf;
    l2r_ = l2r;
    topf_ = topf;
    botf_ = botf;
    topb_ = topb;
    botb_ = botb;
    descid_ = descid;
    parent_ = parent;
    pen_ = pen;
    posid_ = std::numeric_limits<size_t>::max();
    len_ = 0;
    out_.clear();
    edit_ = e;
    lastRecalc_ = true;
	gapadd_ = df[parent].gapadd_;
	if(e.inited()) {
		if(e.isReadGap()) {
			gapadd_++;
		} else if(e.isRefGap()) {
			gapadd_--;
		}
	}
    bool branches = false, hitEnd = false, done = false;
    TIndexOffU topf_new = 0, botf_new = 0, topb_new = 0, botb_new = 0;
    off5p_i_ = 0;
#ifndef NDEBUG
    size_t depth = al5pf_ - al5pi_ + 1;
    TAlScore maxpend = cs[rid_].cons.get(depth, q.length(), maxpen);
    assert_geq(maxpend, pen_);    // can't have already exceeded max penalty
#endif
    bool matchSucc = followMatches(
								   q,
								   sc,
								   gfmFw,
								   gfmBw,
								   re,
								   df,
								   pf,
								   rs,
								   cs,
								   heap,
								   alsink,
								   met,
								   prm,
								   branches,
								   hitEnd,
								   done,
								   off5p_i_,
								   topf_new,
								   botf_new,
								   topb_new,
								   botb_new);
    bool bounceSucc = false;
    if(matchSucc && hitEnd && !done) {
		assert(topf_new > 0 || botf_new > 0);
        bounceSucc = bounce(
							q,
							topf_new,
							botf_new,
							topb_new,
							botb_new,
							gfmFw,
							gfmBw,
							sc,
							minsc,    // minimum score
							maxpen,   // maximum penalty
							re,
							df,
							pf,
							rs,
							cs,
							heap,
							alsink,
							met,      // descent metrics
							prm);     // per-read metrics
    }
	if(matchSucc) {
		// Calculate info about outgoing edges
		recalcOutgoing(q, sc, minsc, maxpen, re, pf, rs, cs, prm);
		if(!empty()) {
			heap.insert(make_pair(out_.bestPri(), descid)); // Add to heap
		}
	}
    return !empty() || bounceSucc;
}

/**
 * Initialize a new descent beginning at the given root.  Return false if
 * the Descent has no outgoing edges (and can therefore have its memory
 * freed), true otherwise.
 */
template <typename index_t>
bool Descent<index_t>::init(
							const Read& q,                  // query
							TRootId rid,                    // root id
							const Scoring& sc,              // scoring scheme
							TAlScore minsc,                 // minimum score
							TAlScore maxpen,                // maximum penalty
							size_t descid,                  // id of this Descent
							const GFM<index_t>& gfmFw,      // forward index
							const GFM<index_t>& gfmBw,      // mirror index
							DescentRedundancyChecker& re,   // redundancy checker
							EFactory<Descent<index_t> >& df,          // Descent factory
							EFactory<DescentPos>& pf,       // DescentPos factory
							const EList<DescentRoot>& rs,   // roots
							const EList<DescentConfig>& cs, // configs
							EHeap<TDescentPair>& heap,      // heap
							DescentAlignmentSink<index_t>& alsink,   // alignment sink
							DescentMetrics& met,            // metrics
							PerReadMetrics& prm)            // per-read metrics
{
    rid_ = rid;
    al5pi_ = rs[rid].off5p;
    al5pf_ = rs[rid].off5p;
	assert_lt(al5pi_, q.length());
	assert_lt(al5pf_, q.length());
    l2r_ = rs[rid].l2r;
    topf_ = botf_ = topb_ = botb_ = 0;
    descid_ = descid;
    parent_ = std::numeric_limits<size_t>::max();
    pen_ = 0;
    posid_ = std::numeric_limits<size_t>::max();
    len_ = 0;
    out_.clear();
    edit_.reset();
    lastRecalc_ = true;
	gapadd_ = 0;
    bool branches = false, hitEnd = false, done = false;
    TIndexOffU topf_new = 0, botf_new = 0, topb_new = 0, botb_new = 0;
    off5p_i_ = 0;
    bool matchSucc = followMatches(
								   q,
								   sc,
								   gfmFw,
								   gfmBw,
								   re,
								   df,
								   pf,
								   rs,
								   cs,
								   heap,
								   alsink,
								   met,
								   prm,
								   branches,
								   hitEnd,
								   done,
								   off5p_i_,
								   topf_new,
								   botf_new,
								   topb_new,
								   botb_new);
    bool bounceSucc = false;
    if(matchSucc && hitEnd && !done) {
		assert(topf_new > 0 || botf_new > 0);
        bounceSucc = bounce(
							q,
							topf_new,
							botf_new,
							topb_new,
							botb_new,
							gfmFw,
							gfmBw,
							sc,
							minsc,    // minimum score
							maxpen,   // maximum penalty
							re,
							df,
							pf,
							rs,
							cs,
							heap,
							alsink,
							met,      // descent metrics
							prm);     // per-read metrics
    }
    // Calculate info about outgoing edges
    assert(empty());
	if(matchSucc) {
		recalcOutgoing(q, sc, minsc, maxpen, re, pf, rs, cs, prm);
		if(!empty()) {
			heap.insert(make_pair(out_.bestPri(), descid)); // Add to heap
		}
	}
    return !empty() || bounceSucc;
}

/**
 * Recalculate our summary of the outgoing edges from this descent.  When
 * deciding what outgoing edges are legal, we abide by constraints.
 * Typically, they limit the total of the penalties accumulated so far, as
 * a function of distance from the search root.  E.g. a constraint might
 * disallow any gaps or mismatches within 20 ply of the search root, then
 * allow 1 mismatch within 30 ply, then allow up to 1 mismatch or 1 gap
 * within 40 ply, etc.
 *
 * Return the total number of valid outgoing edges found.
 *
 * TODO: Eliminate outgoing gap edges that are redundant with others owing to
 *       the DNA sequence and the fact that we don't care to distinguish among
 *       "equivalent" homopolymer extensinos and retractions.
 */
template <typename index_t>
size_t Descent<index_t>::recalcOutgoing(
										const Read& q,                   // query string
										const Scoring& sc,               // scoring scheme
										TAlScore minsc,                  // minimum score
										TAlScore maxpen,                 // maximum penalty
										DescentRedundancyChecker& re,    // redundancy checker
										EFactory<DescentPos>& pf,        // factory with DescentPoss
										const EList<DescentRoot>& rs,    // roots
										const EList<DescentConfig>& cs,  // configs
										PerReadMetrics& prm)             // per-read metrics
{
    assert_eq(botf_ - topf_, botb_ - topb_);
	assert(out_.empty());
	assert(repOk(&q));
	// Get initial 5' and 3' offsets
    bool fw = rs[rid_].fw;
    float rootpri = rs[rid_].pri;
	bool toward3p = (l2r_ == fw);
	size_t off5p = off5p_i_;
	assert_geq(al5pf_, al5pi_);
	size_t off3p = q.length() - off5p - 1;
	// By "depth" we essentially mean the number of characters already aligned
	size_t depth, extrai = 0, extraf = 0;
	size_t cur5pi = al5pi_, cur5pf = al5pf_;
    if(toward3p) {
		// Toward 3'
		cur5pf = off5p;
        depth = off5p - al5pi_;
		// Failed to match out to the end?
		if(al5pf_ < q.length() - 1) {
			extraf = 1; // extra 
		}
    } else {
		// Toward 5'
		cur5pi = off5p;
        depth = al5pf_ - off5p;
		if(al5pi_ > 0) {
			extrai = 1;
		}
    }
	// Get gap penalties
	TScore pen_rdg_ex = sc.readGapExtend(), pen_rfg_ex = sc.refGapExtend();
	TScore pen_rdg_op = sc.readGapOpen(),   pen_rfg_op = sc.refGapOpen();
	// Top and bot in the direction of the descent
	TIndexOffU top  = l2r_ ? topb_ : topf_;
	TIndexOffU bot  = l2r_ ? botb_ : botf_;
	// Top and bot in the opposite direction
	TIndexOffU topp = l2r_ ? topf_ : topb_;
	TIndexOffU botp = l2r_ ? botf_ : botb_;
	assert_eq(botp - topp, bot - top);
	DescentEdge edge;
	size_t nout = 0;
	// Enumerate all outgoing edges, starting at the root and going out
    size_t d = posid_;
	// At first glance, we might think we should be bounded by al5pi_ and
	// al5pf_, but those delimit the positions that matched between reference
	// and read.  If we hit a position that failed to match as part of
	// followMatches, then we also want to evaluate ways of leaving that
	// position, which adds one more position to viist.
	while(off5p >= al5pi_ - extrai && off5p <= al5pf_ + extraf) {
        assert_lt(off5p, q.length());
        assert_lt(off3p, q.length());
		TScore maxpend = cs[rid_].cons.get(depth, q.length(), maxpen);
		assert(depth > 0 || maxpend == 0);
		assert_geq(maxpend, pen_);    // can't have already exceeded max penalty
		TScore diff = maxpend - pen_; // room we have left
		// Get pointer to SA ranges in the direction of descent
		const TIndexOffU *t  = l2r_ ? pf[d].topb : pf[d].topf;
		const TIndexOffU *b  = l2r_ ? pf[d].botb : pf[d].botf;
		const TIndexOffU *tp = l2r_ ? pf[d].topf : pf[d].topb;
		const TIndexOffU *bp = l2r_ ? pf[d].botf : pf[d].botb;
		assert_eq(pf[d].botf - pf[d].topf, pf[d].botb - pf[d].topb);
		// What are the read char / quality?
		std::pair<int, int> p = q.get(off5p, fw);
		int c = p.first;
		assert_range(0, 4, c);
		// Only entertain edits if there is at least one type of edit left and
		// there is some penalty budget left
		if(!pf[d].flags.exhausted() && diff > 0) {
			// What would the penalty be if we mismatched at this position?
			// This includes the case where the mismatch is for an N in the
			// read.
			int qq = p.second;
            assert_geq(qq, 0);
			TScore pen_mm = sc.mm(c, qq);
			if(pen_mm <= diff) {
				for(int j = 0; j < 4; j++) {
					if(j == c) continue; // Match, not mismatch
					if(b[j] <= t[j]) {
						continue; // No outgoing edge with this nucleotide
					}
					if(!pf[d].flags.mmExplore(j)) {
						continue; // Already been explored
					}
					TIndexOffU topf = pf[d].topf[j], botf = pf[d].botf[j];
					ASSERT_ONLY(TIndexOffU topb = pf[d].topb[j], botb = pf[d].botb[j]);
					if(re.contains(fw, l2r_, cur5pi, cur5pf, cur5pf - cur5pi + 1 + gapadd_, topf, botf, pen_ + pen_mm)) {
						prm.nRedSkip++;
						continue; // Redundant with a path already explored
					}
					prm.nRedFail++;
					TIndexOffU width = b[j] - t[j];
					Edit edit((uint32_t)off5p, (int)("ACGTN"[j]), (int)("ACGTN"[c]), EDIT_TYPE_MM);
					DescentPriority pri(pen_ + pen_mm, depth, width, rootpri);
                    assert(topf != 0 || botf != 0);
                    assert(topb != 0 || botb != 0);
					assert_eq(botb - topb, botf - topf);
					edge.init(edit, off5p, pri, d
#ifndef NDEBUG
							  , d, topf, botf, topb, botb
#endif
							  );
					out_.update(edge);
					nout++;
				}
			}
			bool gapsAllowed = (off5p >= (size_t)sc.gapbar && off3p >= (size_t)sc.gapbar);
			if(gapsAllowed) {
				assert_gt(depth, 0);
				// An easy redundancy check is: if all ways of proceeding are
				// matches, then there's no need to entertain gaps here.
				// Shifting the gap one position further downstream is
				// guarnteed not to be worse.
				size_t totwidth = (b[0] - t[0]) +
				(b[1] - t[1]) +
				(b[2] - t[2]) +
				(b[3] - t[3]);
				assert(c > 3 || b[c] - t[c] <= totwidth);
				bool allmatch = c < 4 && (totwidth == (b[c] - t[c]));
				bool rdex = false, rfex = false;
				size_t cur5pi_i = cur5pi, cur5pf_i = cur5pf;
				if(toward3p) {
					cur5pf_i--;
				} else {
					cur5pi_i++;
				}
				if(off5p == off5p_i_ && edit_.inited()) {
					// If we're at the root of the descent, and the descent
					// branched on a gap, then this could be scored as an
					// extension of that gap.
					if(pen_rdg_ex <= diff && edit_.isReadGap()) {
						// Extension of a read gap
						rdex = true;
						for(int j = 0; j < 4; j++) {
							if(b[j] <= t[j]) {
								continue; // No outgoing edge with this nucleotide
							}
							if(!pf[d].flags.rdgExplore(j)) {
								continue; // Already been explored
							}
							TIndexOffU topf = pf[d].topf[j], botf = pf[d].botf[j];
							ASSERT_ONLY(TIndexOffU topb = pf[d].topb[j], botb = pf[d].botb[j]);
							assert(topf != 0 || botf != 0);
							assert(topb != 0 || botb != 0);
							if(re.contains(fw, l2r_, cur5pi_i, cur5pf_i, cur5pf - cur5pi + 1 + gapadd_, topf, botf, pen_ + pen_rdg_ex)) {
								prm.nRedSkip++;
								continue; // Redundant with a path already explored
							}
							prm.nRedFail++;
							TIndexOffU width = b[j] - t[j];
							// off5p holds the offset from the 5' of the next
							// character we were trying to align when we decided to
							// introduce a read gap (before that character).  If we
							// were walking toward the 5' end, we need to increment
							// by 1.
							uint32_t off = (uint32_t)off5p + (toward3p ? 0 : 1);
							Edit edit(off, (int)("ACGTN"[j]), '-', EDIT_TYPE_READ_GAP);
							assert(edit.pos2 != std::numeric_limits<uint32_t>::max());
							edit.pos2 = edit_.pos2 + (toward3p ? 1 : -1);
							DescentPriority pri(pen_ + pen_rdg_ex, depth, width, rootpri);
                            assert(topf != 0 || botf != 0);
                            assert(topb != 0 || botb != 0);
							assert_eq(botb - topb, botf - topf);
							edge.init(edit, off5p, pri, d
#ifndef NDEBUG
									  , d,
									  topf, botf, topb, botb
#endif
									  );
							out_.update(edge);
							nout++;
						}
					}
					if(pen_rfg_ex <= diff && edit_.isRefGap()) {
						// Extension of a reference gap
						rfex = true;
						if(pf[d].flags.rfgExplore()) {
                            TIndexOffU topf = l2r_ ? topp : top;
                            TIndexOffU botf = l2r_ ? botp : bot;
							ASSERT_ONLY(TIndexOffU topb = l2r_ ? top : topp);
							ASSERT_ONLY(TIndexOffU botb = l2r_ ? bot : botp);
							assert(topf != 0 || botf != 0);
							assert(topb != 0 || botb != 0);
							size_t nrefal = cur5pf - cur5pi + gapadd_;
							if(!re.contains(fw, l2r_, cur5pi, cur5pf, nrefal, topf, botf, pen_ + pen_rfg_ex)) {
								TIndexOffU width = bot - top;
								Edit edit((uint32_t)off5p, '-', (int)("ACGTN"[c]), EDIT_TYPE_REF_GAP);
								DescentPriority pri(pen_ + pen_rfg_ex, depth, width, rootpri);
								assert(topf != 0 || botf != 0);
								assert(topb != 0 || botb != 0);
								edge.init(edit, off5p, pri, d
#ifndef NDEBUG
										  // It's a little unclear what the depth ought to be.
										  // Is it the depth we were at when we did the ref
										  // gap?  I.e. the depth of the flags where rfgExplore()
										  // returned true?  Or is it the depth where we can
										  // retrieve the appropriate top/bot?  We make it the
										  // latter, might wrap around, indicating we should get
										  // top/bot from the descent's topf_, ... fields.
										  , (d == posid_) ? std::numeric_limits<size_t>::max() : (d - 1),
										  topf, botf, topb, botb
#endif
										  );
								out_.update(edge);
								nout++;
								prm.nRedFail++;
							} else {
								prm.nRedSkip++;
							}
						}
					}
				}
				if(!allmatch && pen_rdg_op <= diff && !rdex) {
					// Opening a new read gap
					for(int j = 0; j < 4; j++) {
						if(b[j] <= t[j]) {
							continue; // No outgoing edge with this nucleotide
						}
						if(!pf[d].flags.rdgExplore(j)) {
							continue; // Already been explored
						}
						TIndexOffU topf = pf[d].topf[j], botf = pf[d].botf[j];
						ASSERT_ONLY(TIndexOffU topb = pf[d].topb[j], botb = pf[d].botb[j]);
						assert(topf != 0 || botf != 0);
						assert(topb != 0 || botb != 0);
						if(re.contains(fw, l2r_, cur5pi_i, cur5pf_i, cur5pf - cur5pi + 1 + gapadd_, topf, botf, pen_ + pen_rdg_op)) {
							prm.nRedSkip++;
							continue; // Redundant with a path already explored
						}
						prm.nRedFail++;
						TIndexOffU width = b[j] - t[j];
						// off5p holds the offset from the 5' of the next
						// character we were trying to align when we decided to
						// introduce a read gap (before that character).  If we
						// were walking toward the 5' end, we need to increment
						// by 1.
						uint32_t off = (uint32_t)off5p + (toward3p ? 0 : 1);
						Edit edit(off, (int)("ACGTN"[j]), '-', EDIT_TYPE_READ_GAP);
						assert(edit.pos2 != std::numeric_limits<uint32_t>::max());
						DescentPriority pri(pen_ + pen_rdg_op, depth, width, rootpri);
                        assert(topf != 0 || botf != 0);
                        assert(topb != 0 || botb != 0);
						assert_eq(botb - topb, botf - topf);
						edge.init(edit, off5p, pri, d
#ifndef NDEBUG
								  , d, topf, botf, topb, botb
#endif
								  );
						out_.update(edge);
						nout++;
					}
				}
				if(!allmatch && pen_rfg_op <= diff && !rfex) {
					// Opening a new reference gap
                    if(pf[d].flags.rfgExplore()) {
                        TIndexOffU topf = l2r_ ? topp : top;
                        TIndexOffU botf = l2r_ ? botp : bot;
						ASSERT_ONLY(TIndexOffU topb = l2r_ ? top : topp);
						ASSERT_ONLY(TIndexOffU botb = l2r_ ? bot : botp);
						assert(topf != 0 || botf != 0);
						assert(topb != 0 || botb != 0);
						size_t nrefal = cur5pf - cur5pi + gapadd_;
						if(!re.contains(fw, l2r_, cur5pi, cur5pf, nrefal, topf, botf, pen_ + pen_rfg_op)) {
							TIndexOffU width = bot - top;
							Edit edit((uint32_t)off5p, '-', (int)("ACGTN"[c]), EDIT_TYPE_REF_GAP);
							DescentPriority pri(pen_ + pen_rfg_op, depth, width, rootpri);
							assert(topf != 0 || botf != 0);
							assert(topb != 0 || botb != 0);
							edge.init(edit, off5p, pri, d
#ifndef NDEBUG
									  // It's a little unclear what the depth ought to be.
									  // Is it the depth we were at when we did the ref
									  // gap?  I.e. the depth of the flags where rfgExplore()
									  // returned true?  Or is it the depth where we can
									  // retrieve the appropriate top/bot?  We make it the
									  // latter, might wrap around, indicating we should get
									  // top/bot from the descent's topf_, ... fields.
									  , (d == posid_) ? std::numeric_limits<size_t>::max() : (d - 1),
									  topf, botf, topb, botb
#endif
									  );
							out_.update(edge);
							nout++;
							prm.nRedFail++;
						} else {
							prm.nRedSkip++;
						}
                    }
				}
			}
		}
		// Update off5p, off3p, depth
        d++;
		depth++;
        assert_leq(depth, al5pf_ - al5pi_ + 2);
        if(toward3p) {
            if(off3p == 0) {
                break;
            }
            off5p++;
            off3p--;
			cur5pf++;
        } else {
            if(off5p == 0) {
                break;
            }
            off3p++;
            off5p--;
			cur5pi--;
        }
		// Update top and bot
		if(off5p >= al5pi_ - extrai && off5p <= al5pf_ + extraf) {
			assert_range(0, 3, c);
			top = t[c], topp = tp[c];
			bot = b[c], botp = bp[c];
			assert_eq(bot-top, botp-topp);
		}
	}
	lastRecalc_ = (nout <= 5);
    out_.best1.updateFlags(pf);
    out_.best2.updateFlags(pf);
    out_.best3.updateFlags(pf);
    out_.best4.updateFlags(pf);
    out_.best5.updateFlags(pf);
	return nout;
}

template <typename index_t>
void Descent<index_t>::print(
							 std::ostream *os,
							 const char *prefix,
							 const Read& q,
							 size_t trimLf,
							 size_t trimRg,
							 bool fw,
							 const EList<Edit>& edits,
							 size_t ei,
							 size_t en,
							 BTDnaString& rf) const
{
	const BTDnaString& read = fw ? q.patFw : q.patRc;
	size_t eidx = ei;
	if(os != NULL) { *os << prefix; }
	// Print read
	for(size_t i = 0; i < read.length(); i++) {
		if(i < trimLf || i >= read.length() - trimRg) {
			if(os != NULL) { *os << (char)tolower(read.toChar(i)); }
			continue;
		}
		bool del = false, mm = false;
		while(eidx < ei + en && edits[eidx].pos == i) {
			if(edits[eidx].isReadGap()) {
				if(os != NULL) { *os << '-'; }
			} else if(edits[eidx].isRefGap()) {
				del = true;
				assert_eq((int)edits[eidx].qchr, read.toChar(i));
				if(os != NULL) { *os << read.toChar(i); }
			} else {
				mm = true;
				assert(edits[eidx].isMismatch());
				assert_eq((int)edits[eidx].qchr, read.toChar(i));
				if(os != NULL) { *os << (char)edits[eidx].qchr; }
			}
			eidx++;
		}
		if(!del && !mm) {
			// Print read character
			if(os != NULL) { *os << read.toChar(i); }
		}
	}
	if(os != NULL) {
		*os << endl;
		*os << prefix;
	}
	eidx = ei;
	// Print match bars
	for(size_t i = 0; i < read.length(); i++) {
		if(i < trimLf || i >= read.length() - trimRg) {
			if(os != NULL) { *os << ' '; }
			continue;
		}
		bool del = false, mm = false;
		while(eidx < ei + en && edits[eidx].pos == i) {
			if(edits[eidx].isReadGap()) {
				if(os != NULL) { *os << ' '; }
			} else if(edits[eidx].isRefGap()) {
				del = true;
				if(os != NULL) { *os << ' '; }
			} else {
				mm = true;
				assert(edits[eidx].isMismatch());
				if(os != NULL) { *os << ' '; }
			}
			eidx++;
		}
		if(!del && !mm && os != NULL) { *os << '|'; }
	}
	if(os != NULL) {
		*os << endl;
		*os << prefix;
	}
	eidx = ei;
	// Print reference
	for(size_t i = 0; i < read.length(); i++) {
		if(i < trimLf || i >= read.length() - trimRg) {
			if(os != NULL) { *os << ' '; }
			continue;
		}
		bool del = false, mm = false;
		while(eidx < ei + en && edits[eidx].pos == i) {
			if(edits[eidx].isReadGap()) {
				rf.appendChar((char)edits[eidx].chr);
				if(os != NULL) { *os << (char)edits[eidx].chr; }
			} else if(edits[eidx].isRefGap()) {
				del = true;
				if(os != NULL) { *os << '-'; }
			} else {
				mm = true;
				assert(edits[eidx].isMismatch());
				rf.appendChar((char)edits[eidx].chr);
				if(os != NULL) { *os << (char)edits[eidx].chr; }
			}
			eidx++;
		}
		if(!del && !mm) {
			rf.append(read[i]);
			if(os != NULL) { *os << read.toChar(i); }
		}
	}
	if(os != NULL) { *os << endl; }
}

/**
 * Create a new Descent 
 */
template <typename index_t>
bool Descent<index_t>::bounce(
							  const Read& q,                  // query string
							  TIndexOffU topf,                 // SA range top in fw index
							  TIndexOffU botf,                 // SA range bottom in fw index
							  TIndexOffU topb,                 // SA range top in bw index
							  TIndexOffU botb,                 // SA range bottom in bw index
							  const GFM<index_t>& gfmFw,             // forward index
							  const GFM<index_t>& gfmBw,             // mirror index
							  const Scoring& sc,              // scoring scheme
							  TAlScore minsc,                 // minimum score
							  TAlScore maxpen,                // maximum penalty
							  DescentRedundancyChecker& re,   // redundancy checker
							  EFactory<Descent<index_t> >& df,          // factory with Descent
							  EFactory<DescentPos>& pf,       // factory with DescentPoss
							  const EList<DescentRoot>& rs,   // roots
							  const EList<DescentConfig>& cs, // configs
							  EHeap<TDescentPair>& heap,      // heap of descents
							  DescentAlignmentSink<index_t>& alsink,   // alignment sink
							  DescentMetrics& met,            // metrics
							  PerReadMetrics& prm)            // per-read metrics
{
    assert_gt(botf, topf);
    assert(al5pi_ == 0 || al5pf_ == q.length()-1);
    assert(!(al5pi_ == 0 && al5pf_ == q.length()-1));
    size_t dfsz = df.size();
    size_t pfsz = pf.size();
	TDescentId id = df.alloc();
    Edit e_null;
    assert(!e_null.inited());
	// Follow matches 
	bool succ = df[id].init(
							q,         // query
							rid_,      // root id
							sc,        // scoring scheme
							minsc,     // minimum score
							maxpen,    // maximum penalty
							al5pi_,    // new near-5' extreme
							al5pf_,    // new far-5' extreme
							topf,      // SA range top in FW index
							botf,      // SA range bottom in FW index
							topb,      // SA range top in BW index
							botb,      // SA range bottom in BW index
							!l2r_,     // direction this descent will go in; opposite from parent
							id,        // my ID
							descid_,   // parent ID
							pen_,      // total penalties so far - same as parent
							e_null,    // edit for incoming edge; uninitialized if bounced
							gfmFw,     // forward index
							gfmBw,     // mirror index
							re,        // redundancy checker
							df,        // Descent factory
							pf,        // DescentPos factory
							rs,        // DescentRoot list
							cs,        // DescentConfig list
							heap,      // heap
							alsink,    // alignment sink
							met,       // metrics
							prm);      // per-read metrics
    if(!succ) {
        // Reclaim memory we had used for this descent and its DescentPos info
        df.resize(dfsz);
        pf.resize(pfsz);
    }
    return succ;
}

/**
 * Take the best outgoing edge and place it in the heap.  When deciding what
 * outgoing edges exist, abide by constraints in DescentConfig.  These
 * constraints limit total penalty accumulated so far versus distance from
 * search root.  E.g. a constraint might disallow any gaps or mismatches within
 * 20 ply of the root, then allow 1 mismatch within 30 ply, 1 mismatch or 1 gap
 * within 40 ply, etc.
 */
template <typename index_t>
void Descent<index_t>::followBestOutgoing(
										  const Read& q,                   // query string
										  const GFM<index_t>& gfmFw,       // forward index
										  const GFM<index_t>& gfmBw,       // mirror index
										  const Scoring& sc,               // scoring scheme
										  TAlScore minsc,                  // minimum score
										  TAlScore maxpen,                 // maximum penalty
										  DescentRedundancyChecker& re,    // redundancy checker
										  EFactory<Descent<index_t> >& df, // factory with Descent
										  EFactory<DescentPos>& pf,        // factory with DescentPoss
										  const EList<DescentRoot>& rs,    // roots
										  const EList<DescentConfig>& cs,  // configs
										  EHeap<TDescentPair>& heap,       // heap of descents
										  DescentAlignmentSink<index_t>& alsink,   // alignment sink
										  DescentMetrics& met,             // metrics
										  PerReadMetrics& prm)             // per-read metrics
{
	// We assume this descent has been popped off the heap.  We'll re-add it if
	// it hasn't been exhausted.
	assert(q.repOk());
	assert(!empty());
	assert(!out_.empty());
	while(!out_.empty()) {
		DescentPriority best = out_.bestPri();
		DescentEdge e = out_.rotate();
		TReadOff al5pi_new = al5pi_, al5pf_new = al5pf_;
		bool fw = rs[rid_].fw;
		bool toward3p = (l2r_ == fw);
		TReadOff edoff = e.off5p; // 5' offset of edit
		assert_leq(edoff, al5pf_ + 1);
		assert_geq(edoff + 1, al5pi_);
		if(out_.empty()) {
			if(!lastRecalc_) {
				// This might allocate new Descents
				recalcOutgoing(q, sc, minsc, maxpen, re, pf, rs, cs, prm);
				if(empty()) {
					// Could happen, since some outgoing edges may have become
					// redundant in the meantime.
					break;
				}
			} else {
				assert(empty());
			}
		}
		TReadOff doff; // edit's offset into this descent
		int chr = asc2dna[e.e.chr];
		// hitEnd is set to true iff this edit pushes us to the extreme 5' or 3'
		// end of the alignment
		bool hitEnd = false;
		// done is set to true iff this edit aligns the only remaining character of
		// the read
		bool done = false;
		if(toward3p) {
			// The 3' extreme of the new Descent is further in (away from the 3'
			// end) than the parent's.
			al5pf_new = doff = edoff;
			if(e.e.isReadGap()) {
				// We didn't actually consume the read character at 'edoff', so
				// retract al5pf_new by one position.  This doesn't effect the
				// "depth" (doff) of the SA range we took, though.
				assert_gt(al5pf_new, 0);
				al5pf_new--;
			}
			assert_lt(al5pf_new, q.length());
			hitEnd = (al5pf_new == q.length() - 1);
			done = (hitEnd && al5pi_new == 0);
			assert_geq(doff, off5p_i_);
			doff = doff - off5p_i_;
			assert_leq(doff, len_);
		} else {
			// The 5' extreme of the new Descent is further in (away from the 5'
			// end) than the parent's.
			al5pi_new = doff = edoff;
			if(e.e.isReadGap()) {
				// We didn't actually consume the read character at 'edoff', so
				// move al5pi_new closer to the 3' end by one position.  This
				// doesn't effect the "depth" (doff) of the SA range we took,
				// though.
				al5pi_new++;
			}
			hitEnd = (al5pi_new == 0);
			done = (hitEnd && al5pf_new == q.length() - 1);
			assert_geq(off5p_i_, doff);
			doff = off5p_i_ - doff;
			assert_leq(doff, len_);
		}
		// Check if this is redundant with an already-explored path
		bool l2r = l2r_; // gets overridden if we bounce
		if(!done && hitEnd) {
			// Alignment finsihed extending in one direction
			l2r = !l2r;
		}
		size_t dfsz = df.size();
		size_t pfsz = pf.size();
		TIndexOffU topf, botf, topb, botb;
		size_t d = posid_ + doff;
		if(e.e.isRefGap()) {
			d--; // might underflow
			if(doff == 0) {
				topf = topf_;
				botf = botf_;
				topb = topb_;
				botb = botb_;
				d = std::numeric_limits<size_t>::max();
				assert_eq(botf-topf, botb-topb);
			} else {
				assert_gt(al5pf_new, 0);
				assert_gt(d, 0);
				chr = pf[d].c;
				assert(pf[d].inited());
				assert_range(0, 3, chr);
				topf = pf[d].topf[chr];
				botf = pf[d].botf[chr];
				topb = pf[d].topb[chr];
				botb = pf[d].botb[chr];
				assert_eq(botf-topf, botb-topb);
			}
		} else {
			// A read gap or a mismatch
			assert(pf[d].inited());
			topf = pf[d].topf[chr];
			botf = pf[d].botf[chr];
			topb = pf[d].topb[chr];
			botb = pf[d].botb[chr];
			assert_eq(botf-topf, botb-topb);
		}
		assert_eq(d, e.d);
		assert_eq(topf, e.topf);
		assert_eq(botf, e.botf);
		assert_eq(topb, e.topb);
		assert_eq(botb, e.botb);
		if(done) {
			// Aligned the entire read end-to-end.  Presumably there's no need to
			// create a new Descent object.  We just report the alignment.
			alsink.reportAlignment(
								   q,        // query
								   gfmFw,    // forward index
								   gfmBw,    // backward index
								   topf,     // top of SA range in forward index
								   botf,     // bottom of SA range in forward index
								   topb,     // top of SA range in backward index
								   botb,     // bottom of SA range in backward index
								   descid_,  // Descent at the leaf
								   rid_,     // root id
								   e.e,      // extra edit, if necessary
								   best.pen, // penalty
								   df,       // factory with Descent
								   pf,       // factory with DescentPoss
								   rs,       // roots
								   cs);      // configs
			assert(alsink.repOk());
			return;
		}
		assert(al5pi_new != 0 || al5pf_new != q.length() - 1);
		TDescentId id = df.alloc();
		bool succ = df[id].init(
								q,         // query
								rid_,      // root id
								sc,        // scoring scheme
								minsc,     // minimum score
								maxpen,    // maximum penalty
								al5pi_new, // new near-5' extreme
								al5pf_new, // new far-5' extreme
								topf,      // SA range top in FW index
								botf,      // SA range bottom in FW index
								topb,      // SA range top in BW index
								botb,      // SA range bottom in BW index
								l2r,       // direction this descent will go in
								id,        // my ID
								descid_,   // parent ID
								best.pen,  // total penalties so far
								e.e,       // edit for incoming edge; uninitialized if bounced
								gfmFw,     // forward index
								gfmBw,     // mirror index
								re,        // redundancy checker
								df,        // Descent factory
								pf,        // DescentPos factory
								rs,        // DescentRoot list
								cs,        // DescentConfig list
								heap,      // heap
								alsink,    // alignment sink
								met,       // metrics
								prm);      // per-read metrics
		if(!succ) {
			// Reclaim memory we had used for this descent and its DescentPos info
			df.resize(dfsz);
			pf.resize(pfsz);
		}
		break;
	}
	if(!empty()) {
		// Re-insert this Descent with its new priority
		heap.insert(make_pair(out_.bestPri(), descid_));
	}
}

/**
 * Given the forward and backward indexes, and given topf/botf/topb/botb, get
 * tloc, bloc ready for the next step.
 */
template <typename index_t>
void Descent<index_t>::nextLocsBi(
								  const GFM<index_t>& gfmFw, // forward index
								  const GFM<index_t>& gfmBw, // mirror index
								  SideLocus<index_t>& tloc,    // top locus
								  SideLocus<index_t>& bloc,    // bot locus
								  index_t topf,     // top in BWT
								  index_t botf,     // bot in BWT
								  index_t topb,     // top in BWT'
								  index_t botb)     // bot in BWT'
{
	assert_gt(botf, 0);
	// Which direction are we going in next?
	if(l2r_) {
		// Left to right; use BWT'
		if(botb - topb == 1) {
			// Already down to 1 row; just init top locus
			tloc.initFromRow(topb, gfmBw.gh(), gfmBw.gfm());
			bloc.invalidate();
		} else {
			SideLocus<index_t>::initFromTopBot(
											   topb, botb, gfmBw.gh(), gfmBw.gfm(), tloc, bloc);
			assert(bloc.valid());
		}
	} else {
		// Right to left; use BWT
		if(botf - topf == 1) {
			// Already down to 1 row; just init top locus
			tloc.initFromRow(topf, gfmFw.gh(), gfmFw.gfm());
			bloc.invalidate();
		} else {
			SideLocus<index_t>::initFromTopBot(
											   topf, botf, gfmFw.gh(), gfmFw.gfm(), tloc, bloc);
			assert(bloc.valid());
		}
	}
	// Check if we should update the tracker with this refinement
	assert(botf - topf == 1 ||  bloc.valid());
	assert(botf - topf > 1  || !bloc.valid());
}

/**
 * Advance this descent by following read matches as far as possible.
 *
 * This routine doesn't have to consider the whole gamut of constraints on
 * which outgoing edges can be followed.  If it is a root descent, it does have
 * to know how deep the no-edit constraint goes, though, so we can decide
 * whether using the ftab would potentially jump over relevant branch points.
 * Apart from that, though, we simply proceed as far as it can go by matching
 * characters in the query, irrespective of the constraints.
 * recalcOutgoing(...) and followBestOutgoing(...) do have to consider these
 * constraints, though.
 *
 * Conceptually, as we make descending steps, we have:
 * 1. Before each step, a single range indicating how we departed the previous
 *    step
 * 2. As part of each step, a quad of ranges indicating what range would result
 *    if we proceeded on an a, c, g ot t
 *
 * Return true iff it is possible to branch from this descent.  If we haven't
 * exceeded the no-branch depth.
 */
template <typename index_t>
bool Descent<index_t>::followMatches(
									 const Read& q,     // query string
									 const Scoring& sc,         // scoring scheme
									 const GFM<index_t>& gfmFw,        // forward index
									 const GFM<index_t>& gfmBw,        // mirror index
									 DescentRedundancyChecker& re, // redundancy checker
									 EFactory<Descent<index_t> >& df,     // Descent factory
									 EFactory<DescentPos>& pf,  // DescentPos factory
									 const EList<DescentRoot>& rs,   // roots
									 const EList<DescentConfig>& cs, // configs
									 EHeap<TDescentPair>& heap, // heap
									 DescentAlignmentSink<index_t>& alsink, // alignment sink
									 DescentMetrics& met,       // metrics
									 PerReadMetrics& prm,       // per-read metrics
									 bool& branches,            // out: true -> there are > 0 ways to branch
									 bool& hitEnd,              // out: true -> hit read end with non-empty range
									 bool& done,                // out: true -> we made a full alignment
									 TReadOff& off5p_i,         // out: initial 5' offset
									 TIndexOffU& topf_bounce,    // out: top of SA range for fw idx for bounce
									 TIndexOffU& botf_bounce,    // out: bot of SA range for fw idx for bounce
									 TIndexOffU& topb_bounce,    // out: top of SA range for bw idx for bounce
									 TIndexOffU& botb_bounce)    // out: bot of SA range for bw idx for bounce
{
	// TODO: make these full-fledged parameters
	size_t nobranchDepth = 20;
	bool stopOnN = true;
	assert(q.repOk());
	assert(repOk(&q));
	assert_eq(gfmFw.eh().ftabChars(), gfmBw.gh().ftabChars());
#ifndef NDEBUG
	for(int i = 0; i < 4; i++) {
		assert_eq(gfmFw.fchr()[i], gfmBw.fchr()[i]);
	}
#endif
	SideLocus<index_t> tloc, bloc;
	TIndexOffU topf = topf_, botf = botf_, topb = topb_, botb = botb_;
    bool fw = rs[rid_].fw;
	bool toward3p;
	size_t off5p;
	assert_lt(al5pi_, q.length());
	assert_lt(al5pf_, q.length());
	while(true) {
		toward3p = (l2r_ == fw);
		assert_geq(al5pf_, al5pi_);
		assert(al5pi_ != 0 || al5pf_ != q.length() - 1);
		if(toward3p) {
			if(al5pf_ == q.length()-1) {
				l2r_ = !l2r_;
				continue;
			}
			if(al5pi_ == al5pf_ && root()) {
				off5p = off5p_i = al5pi_;
			} else {
				off5p = off5p_i = (al5pf_ + 1);
			}
		} else {
			if(al5pi_ == 0) {
				l2r_ = !l2r_;
				continue;
			}
			assert_gt(al5pi_, 0);
			if(al5pi_ == al5pf_ && root()) {
				off5p = off5p_i = al5pi_;
			} else {
				off5p = off5p_i = (al5pi_ - 1);
			}
		}
		break;
	}
	size_t off3p = q.length() - off5p - 1;
	assert_lt(off5p, q.length());
	assert_lt(off3p, q.length());
	bool firstPos = true;
	assert_eq(0, len_);
    
	// Number of times pf.alloc() is called.  So we can sanity check it.
	size_t nalloc = 0;
    // Set to true as soon as we encounter a branch point along this descent.
    branches = false;
    // hitEnd is set to true iff this edit pushes us to the extreme 5' or 3'
    // end of the alignment
    hitEnd = false;
    // done is set to true iff this edit aligns the only remaining character of
    // the read
    done = false;
	if(root()) {
        assert_eq(al5pi_, al5pf_);
		// Check whether/how far we can jump using ftab
		int ftabLen = gfmFw.gh().ftabChars();
		bool ftabFits = true;
		if(toward3p && ftabLen + off5p > q.length()) {
			ftabFits = false;
		} else if(!toward3p && off5p < (size_t)ftabLen) {
			ftabFits = false;
		}
		bool useFtab = ftabLen > 1 && (size_t)ftabLen <= nobranchDepth && ftabFits;
		bool ftabFailed = false;
		if(useFtab) {
			prm.nFtabs++;
			// Forward index: right-to-left
			size_t off_r2l = fw ? off5p : q.length() - off5p - 1;
			if(l2r_) {
				//
			} else {
				assert_geq((int)off_r2l, ftabLen - 1);
				off_r2l -= (ftabLen - 1);
			}
			bool ret = gfmFw.ftabLoHi(fw ? q.patFw : q.patRc, off_r2l,
                                      false, // reverse
                                      topf, botf);
			if(!ret) {
				// Encountered an N or something else that made it impossible
				// to use the ftab
				ftabFailed = true;
			} else {
				if(botf - topf == 0) {
					return false;
				}
				int c_r2l = fw ? q.patFw[off_r2l] : q.patRc[off_r2l];
				// Backward index: left-to-right
				size_t off_l2r = fw ? off5p : q.length() - off5p - 1;
				if(l2r_) {
					//
				} else {
					assert_geq((int)off_l2r, ftabLen - 1);
					off_l2r -= (ftabLen - 1);
				}
				ASSERT_ONLY(bool ret2 = )
				gfmBw.ftabLoHi(fw ? q.patFw : q.patRc, off_l2r,
                               false, // don't reverse
                               topb, botb);
				assert(ret == ret2);
				int c_l2r = fw ? q.patFw[off_l2r + ftabLen - 1] :
				q.patRc[off_l2r + ftabLen - 1];
				assert_eq(botf - topf, botb - topb);
				if(toward3p) {
					assert_geq((int)off3p, ftabLen - 1);
					off5p += ftabLen; off3p -= ftabLen;
				} else {
					assert_geq((int)off5p, ftabLen - 1);
					off5p -= ftabLen; off3p += ftabLen;
				}
				len_ += ftabLen;
				if(toward3p) {
					// By convention, al5pf_ and al5pi_ start out equal, so we only
					// advance al5pf_ by ftabLen - 1 (not ftabLen)
					al5pf_ += (ftabLen - 1); // -1 accounts for inclusive al5pf_
					if(al5pf_ == q.length() - 1) {
						hitEnd = true;
						done = (al5pi_ == 0);
					}
				} else {
					// By convention, al5pf_ and al5pi_ start out equal, so we only
					// advance al5pi_ by ftabLen - 1 (not ftabLen)
					al5pi_ -= (ftabLen - 1);
					if(al5pi_ == 0) {
						hitEnd = true;
						done = (al5pf_ == q.length()-1);
					}
				}
				// Allocate DescentPos data structures and leave them empty.  We
				// jumped over them by doing our lookup in the ftab, so we have no
				// info about outgoing edges from them, besides the matching
				// outgoing edge from the last pos which is in topf/botf and
				// topb/botb.
				size_t id = 0;
				if(firstPos) {
					posid_ = pf.alloc();
					pf[posid_].reset();
					firstPos = false;
					for(int i = 1; i < ftabLen; i++) {
						id = pf.alloc();
						pf[id].reset();
					}
				} else {
					for(int i = 0; i < ftabLen; i++) {
						id = pf.alloc();
						pf[id].reset();
					}
				}
				assert_eq(botf-topf, botb-topb);
				pf[id].c = l2r_ ? c_l2r : c_r2l;
				pf[id].topf[l2r_ ? c_l2r : c_r2l] = topf;
				pf[id].botf[l2r_ ? c_l2r : c_r2l] = botf;
				pf[id].topb[l2r_ ? c_l2r : c_r2l] = topb;
				pf[id].botb[l2r_ ? c_l2r : c_r2l] = botb;
				assert(pf[id].inited());
				nalloc += ftabLen;
			}
		}
		if(!useFtab || ftabFailed) {
			// Can't use ftab, use fchr instead
			int rdc = q.getc(off5p, fw);
			// If rdc is N, that's pretty bad!  That means we placed a root
			// right on an N.  The only thing we can reasonably do is to pick a
			// nucleotide at random and proceed.
			if(rdc > 3) {
				return false;
			}
			assert_range(0, 3, rdc);
			topf = topb = gfmFw.fchr()[rdc];
			botf = botb = gfmFw.fchr()[rdc+1];
			if(botf - topf == 0) {
				return false;
			}
			if(toward3p) {
				off5p++; off3p--;
			} else {
				off5p--; off3p++;
			}
			len_++;
            if(toward3p) {
                if(al5pf_ == q.length()-1) {
                    hitEnd = true;
                    done = (al5pi_ == 0);
                }
            } else {
                if(al5pi_ == 0) {
                    hitEnd = true;
                    done = (al5pf_ == q.length()-1);
                }
            }
			// Allocate DescentPos data structure.  We could fill it with the
			// four ranges from fchr if we wanted to, but that will never be
			// relevant.
			size_t id = 0;
			if(firstPos) {
				posid_ = id = pf.alloc();
                firstPos = false;
			} else {
				id = pf.alloc();
			}
			assert_eq(botf-topf, botb-topb);
			pf[id].c = rdc;
			pf[id].topf[rdc] = topf;
			pf[id].botf[rdc] = botf;
			pf[id].topb[rdc] = topb;
			pf[id].botb[rdc] = botb;
			assert(pf[id].inited());
			nalloc++;
		}
		assert_gt(botf, topf);
		assert_eq(botf - topf, botb - topb);
		// Check if this is redundant with an already-explored path
		if(!re.check(fw, l2r_, al5pi_, al5pf_, al5pf_ - al5pi_ + 1 + gapadd_,
		             topf, botf, pen_))
		{
			prm.nRedSkip++;
			return false;
		}
		prm.nRedFail++; // not pruned by redundancy list
		prm.nRedIns++;  // inserted into redundancy list
	}
    if(done) {
        Edit eempty;
        alsink.reportAlignment(
							   q,        // query
							   gfmFw,    // forward index
							   gfmBw,    // backward index
							   topf,     // top of SA range in forward index
							   botf,     // bottom of SA range in forward index
							   topb,     // top of SA range in backward index
							   botb,     // bottom of SA range in backward index
							   descid_,  // Descent at the leaf
							   rid_,     // root id
							   eempty,   // extra edit, if necessary
							   pen_,     // penalty
							   df,       // factory with Descent
							   pf,       // factory with DescentPoss
							   rs,       // roots
							   cs);      // configs
		assert(alsink.repOk());
        return true;
    } else if(hitEnd) {
		assert(botf > 0 || topf > 0);
        assert_gt(botf, topf);
        topf_bounce = topf;
        botf_bounce = botf;
        topb_bounce = topb;
        botb_bounce = botb;
        return true; // Bounced
    }
    // We just advanced either ftabLen characters, or 1 character,
    // depending on whether we used ftab or fchr.
    nextLocsBi(gfmFw, gfmBw, tloc, bloc, topf, botf, topb, botb);
    assert(tloc.valid());
	assert(botf - topf == 1 ||  bloc.valid());
	assert(botf - topf > 1  || !bloc.valid());
	TIndexOffU t[4], b[4];   // dest BW ranges
	TIndexOffU tp[4], bp[4]; // dest BW ranges for "prime" index
	ASSERT_ONLY(TIndexOff lasttot = botf - topf);
	bool fail = false;
	while(!fail && !hitEnd) {
        assert(!done);
		int rdc = q.getc(off5p, fw);
		int rdq = q.getq(off5p);
		assert_range(0, 4, rdc);
		assert_gt(botf, topf);
		assert(botf - topf == 1 ||  bloc.valid());
		assert(botf - topf > 1  || !bloc.valid());
		assert(tloc.valid());
        TIndexOffU width = botf - topf;
		bool ltr = l2r_;
		const GFM<index_t>& gfm = ltr ? gfmBw : gfmFw;
		t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
		int only = -1; // if we only get 1 non-empty range, this is the char
		size_t nopts = 1;
		if(bloc.valid()) {
			// Set up initial values for the primes
			if(ltr) {
				tp[0] = tp[1] = tp[2] = tp[3] = topf;
				bp[0] = bp[1] = bp[2] = bp[3] = botf;
			} else {
				tp[0] = tp[1] = tp[2] = tp[3] = topb;
				bp[0] = bp[1] = bp[2] = bp[3] = botb;
			}
			// Range delimited by tloc/bloc has size >1.  If size == 1,
			// we use a simpler query (see if(!bloc.valid()) blocks below)
			met.bwops++;
			met.bwops_bi++;
			prm.nSdFmops++;
			if(prm.doFmString) {
				prm.fmString.add(false, pen_, 1);
			}
			gfm.mapBiLFEx(tloc, bloc, t, b, tp, bp);
			// t, b, tp and bp now filled
			ASSERT_ONLY(TIndexOffU tot = (b[0]-t[0])+(b[1]-t[1])+(b[2]-t[2])+(b[3]-t[3]));
			ASSERT_ONLY(TIndexOffU totp = (bp[0]-tp[0])+(bp[1]-tp[1])+(bp[2]-tp[2])+(bp[3]-tp[3]));
			assert_eq(tot, totp);
			assert_leq(tot, lasttot);
			ASSERT_ONLY(lasttot = tot);
			fail = (rdc > 3 || b[rdc] <= t[rdc]);
			size_t nopts = 0;
			if(b[0] > t[0]) { nopts++; only = 0; }
			if(b[1] > t[1]) { nopts++; only = 1; }
			if(b[2] > t[2]) { nopts++; only = 2; }
			if(b[3] > t[3]) { nopts++; only = 3; }
            if(!fail && b[rdc] - t[rdc] < width) {
                branches = true;
            }
		} else {
			tp[0] = tp[1] = tp[2] = tp[3] = bp[0] = bp[1] = bp[2] = bp[3] = 0;
			// Range delimited by tloc/bloc has size 1
			TIndexOffU ntop = ltr ? topb : topf;
			met.bwops++;
			met.bwops_1++;
			prm.nSdFmops++;
			if(prm.doFmString) {
				prm.fmString.add(false, pen_, 1);
			}
			int cc = gfm.mapLF1(ntop, tloc);
			assert_range(-1, 3, cc);
            fail = (cc != rdc);
            if(fail) {
                branches = true;
            }
			if(cc >= 0) {
				only = cc;
				t[cc] = ntop; b[cc] = ntop+1;
				tp[cc] = ltr ? topf : topb;
				bp[cc] = ltr ? botf : botb;
			}
		}
		// Now figure out what to do with our N.
		int origRdc = rdc;
		if(rdc == 4) {
			fail = true;
		} else {
			topf = ltr ? tp[rdc] : t[rdc];
			botf = ltr ? bp[rdc] : b[rdc];
			topb = ltr ? t[rdc] : tp[rdc];
			botb = ltr ? b[rdc] : bp[rdc];
			assert_eq(botf - topf, botb - topb);
		}
		// The trouble with !stopOnN is that we don't have a way to store the N
		// edits.  There could be several per Descent.
		if(rdc == 4 && !stopOnN && nopts == 1) {
			fail = false;
			rdc = only;
			int pen = sc.n(rdq);
			assert_gt(pen, 0);
			pen_ += pen;
		}
		assert_range(0, 4, origRdc);
		assert_range(0, 4, rdc);
        // If 'fail' is true, we failed to align this read character.  We still
        // install the SA ranges into the DescentPos and increment len_ in this
        // case.
        
		// Convert t, tp, b, bp info tf, bf, tb, bb
		TIndexOffU *tf = ltr ? tp : t;
		TIndexOffU *bf = ltr ? bp : b;
		TIndexOffU *tb = ltr ? t : tp;
		TIndexOffU *bb = ltr ? b : bp;
		// Allocate DescentPos data structure.
		if(firstPos) {
			posid_ = pf.alloc();
            firstPos = false;
		} else {
			pf.alloc();
		}
		nalloc++;
		pf[posid_ + len_].reset();
        pf[posid_ + len_].c = origRdc;
		for(size_t i = 0; i < 4; i++) {
			pf[posid_ + len_].topf[i] = tf[i];
			pf[posid_ + len_].botf[i] = bf[i];
			pf[posid_ + len_].topb[i] = tb[i];
			pf[posid_ + len_].botb[i] = bb[i];
			assert_eq(pf[posid_ + len_].botf[i] - pf[posid_ + len_].topf[i],
			          pf[posid_ + len_].botb[i] - pf[posid_ + len_].topb[i]);
		}
		if(!fail) {
			// Check if this is redundant with an already-explored path
			size_t al5pf = al5pf_, al5pi = al5pi_;
			if(toward3p) {
				al5pf++;
			} else {
				al5pi--;
			}
			fail = !re.check(fw, l2r_, al5pi, al5pf,
			                 al5pf - al5pi + 1 + gapadd_, topf, botf, pen_);
			if(fail) {
				prm.nRedSkip++;
			} else {
				prm.nRedFail++; // not pruned by redundancy list
				prm.nRedIns++;  // inserted into redundancy list
			}
		}
		if(!fail) {
			len_++;
			if(toward3p) {
				al5pf_++;
				off5p++;
				off3p--;
				if(al5pf_ == q.length() - 1) {
					hitEnd = true;
					done = (al5pi_ == 0);
				}
			} else {
				assert_gt(al5pi_, 0);
				al5pi_--;
				off5p--;
				off3p++;
				if(al5pi_ == 0) {
					hitEnd = true;
					done = (al5pf_ == q.length() - 1);
				}
			}
		}
        if(!fail && !hitEnd) {
            nextLocsBi(gfmFw, gfmBw, tloc, bloc, tf[rdc], bf[rdc], tb[rdc], bb[rdc]);
        }
	}
	assert_geq(al5pf_, al5pi_);
	assert(!root() || al5pf_ - al5pi_ + 1 == nalloc || al5pf_ - al5pi_ + 2 == nalloc);
	assert_geq(pf.size(), nalloc);
    if(done) {
        Edit eempty;
        alsink.reportAlignment(
							   q,        // query
							   gfmFw,    // forward index
							   gfmBw,    // backward index
							   topf,     // top of SA range in forward index
							   botf,     // bottom of SA range in forward index
							   topb,     // top of SA range in backward index
							   botb,     // bottom of SA range in backward index
							   descid_,  // Descent at the leaf
							   rid_,     // root id
							   eempty,   // extra edit, if necessary
							   pen_,     // penalty
							   df,       // factory with Descent
							   pf,       // factory with DescentPoss
							   rs,       // roots
							   cs);      // configs
		assert(alsink.repOk());
        return true;
    } else if(hitEnd) {
        assert(botf > 0 || topf > 0);
        assert_gt(botf, topf);
        topf_bounce = topf;
        botf_bounce = botf;
        topb_bounce = topb;
        botb_bounce = botb;
        return true; // Bounced
    }
    assert(repOk(&q));
	assert(!hitEnd || topf_bounce > 0 || botf_bounce > 0);
	return true;
}

#endif
