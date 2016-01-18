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
 * group_walk.h
 *
 * Classes and routines for walking a set of BW ranges backwards from the edge
 * of a seed hit with the goal of resolving the offset of each row in each
 * range.  Here "offset" means offset into the concatenated string of all
 * references.  The main class is 'GroupWalk' and an important helper is
 * 'GWState'.
 *
 * For each combination of seed offset and orientation, there is an associated
 * QVal.  Each QVal describes a (possibly empty) set of suffix array ranges.
 * Call these "seed range sets."  Each range in the set is "backed" by a range
 * of the salist, represented as a PListSlice. Such a range is the origin of a
 * walk.
 *
 * When an offset is resolved, it is entered into the salist via the
 * PListSlice.  Note that other routines in this same thread might also be
 * setting elements of the salist, so routines here should expect that elements
 * can go from unresolved to resolved at any time.
 *
 * What bookkeeping do we have to do as we walk?  Before the first step, we
 * convert the initial QVal into a list of SATuples; the SATuples are our link
 * to the correpsonding ranges in the suffix array.  The list of SATuples is
 * then converted to a list of GWState objects; these keep track of where we
 * are in our walk (e.g. what 'top' and 'bot' are, how many steps have we gone,
 * etc) as well as how the elements in the current range correspond to elements
 * from the original range.
 *
 * The user asks the GroupWalk to resolve another offset by calling advance().
 * advance() can be called in various ways:
 *
 * (a) The user can request that the GroupWalk proceed until a
 *     *particular* element is resolved, then return that resolved
 *     element.  Other elements may be resolved along the way, but
 *     those results are buffered and may be dispensed in future calls
 *     to advance().
 *
 * (b) The user can request that the GroupWalk select an as-yet-
 *     unreported element at random and and proceed until that element
 *     is resolved and report it.  Again, other elements may be
 *     resolved along the way but they are buffered.
 *
 * (c) The user can request that the GroupWalk resolve elements in a
 *     particular BW range (with a particular offset and orientation)
 *     in an order of its choosing.  The GroupWalk in this case
 *     attempts to resolve as many offsets as possible as quickly as
 *     possible, and returns them as soon as they're found.  The res_
 *     buffer is used in this case.
 *
 * (d) Like (c) but resolving elements at a paritcular offset and
 *     orientation instead of at a specific BW range.  The res_ buffer
 *     is used in this case, since there's a chance that the 
 *
 * There are simple ways to heuristically reduce the problem size while
 * maintaining randomness.  For instance, the user put a ceiling on the
 * number of elements that we walk from any given seed offset or range.
 * We can then trim away random subranges to reduce the size of the
 * problem.  There is no need for the caller to do this for us.
 */

#ifndef GROUP_WALK_H_
#define GROUP_WALK_H_

#include <stdint.h>
#include <limits>
#include "ds.h"
#include "gfm.h"
#include "read.h"
#include "reference.h"
#include "mem_ids.h"

/**
 * Encapsulate an SA range and an associated list of slots where the resolved
 * offsets can be placed.
 */
template<typename T, typename index_t>
class SARangeWithOffs {

public:

	SARangeWithOffs() { reset(); };

	SARangeWithOffs(
                    index_t tf,
                    index_t bf,
                    index_t ntf,
                    index_t nbf,
                    const EList<pair<index_t, index_t> >& n_iedge_count,
                    size_t len,
                    const T& o) {
		init(tf, bf, ntf, nbf, n_iedge_count, len, o);
	}
	
	void init(
              index_t tf,
              index_t bf,
              index_t ntf,
              index_t nbf,
              const EList<pair<index_t, index_t> >& n_iedge_count,
              size_t len_,
              const T& o) {
        topf = tf;
        botf = bf;
        assert_lt(topf, botf);
        node_top = ntf;
        node_bot = nbf;
        assert_leq(node_bot - node_top, botf - topf);
        node_iedge_count = n_iedge_count;
        len = len_,
        offs = o;
	}

	/**
	 * Reset to uninitialized state.
	 */
	void reset() { topf = (index_t)INDEX_MAX; }
	
	/**
	 * Return true if this is initialized.
	 */
	bool inited() const {
		return topf != (index_t)INDEX_MAX;
	}
	
	/**
	 * Return the number of times this reference substring occurs in the
	 * reference, which is also the size of the 'offs' TSlice.
	 */
	size_t size() const { return offs.size(); }

	index_t topf;      // top in GBWT index
    index_t botf;
    index_t node_top;  // top node
    index_t node_bot;
    EList<pair<index_t, index_t> > node_iedge_count;
	size_t     len;        // length of the reference sequence involved
	T          offs;       // offsets
};

/**
 * A group of per-thread state that can be shared between all the GroupWalks
 * used in that thread.
 */
template <typename index_t>
struct GroupWalkState {

	GroupWalkState(int cat) : map(cat) {
		masks[0].setCat(cat);
		masks[1].setCat(cat);
		masks[2].setCat(cat);
		masks[3].setCat(cat);
	}

	EList<bool> masks[4];      // temporary list for masks; used in GWState
	EList<index_t, 16> map;   // temporary list of GWState maps
};

/**
 * Encapsulates counters that encode how much work the walk-left logic
 * has done.
 */
struct WalkMetrics {

	WalkMetrics() {
	    reset();
	}

	/**
	 * Sum each across this object and 'm'.  This is the only safe way
	 * to update a WalkMetrics shared by many threads.
	 */
	void merge(const WalkMetrics& m, bool getLock = false) {
		ThreadSafe ts(&mutex_m, getLock);
		bwops += m.bwops;
		branches += m.branches;
		resolves += m.resolves;
		refresolves += m.refresolves;
		reports += m.reports;
	}
	
	/**
	 * Set all to 0.
	 */
	void reset() {
		bwops = branches = resolves = refresolves = reports = 0;
	}

	uint64_t bwops;       // Burrows-Wheeler operations
	uint64_t branches;    // BW range branch-offs
	uint64_t resolves;    // # offs resolved with BW walk-left
	uint64_t refresolves; // # resolutions caused by reference scanning
	uint64_t reports;     // # offs reported (1 can be reported many times)
	MUTEX_T mutex_m;
};

/**
 * Coordinates for a BW element that the GroupWalk might resolve.
 */
template <typename index_t>
struct GWElt {

	GWElt() { reset(); }
	
	/**
	 * Reset GWElt to uninitialized state.
	 */
	void reset() {
		offidx = range = elt = len = (index_t)OFF_MASK;
		fw = false;
	}

	/**
	 * Initialize this WalkResult.
	 */
	void init(
		index_t oi,
		bool f,
		index_t r,
		index_t e,
		index_t l)
	{
		offidx = oi;
		fw = f;
		range = r;
		elt = e;
		len = l;
	}

	/**
	 * Return true iff this GWElt and the given GWElt refer to the same
	 * element.
	 */
	bool operator==(const GWElt& o) const {
		return offidx == o.offidx &&
		       fw == o.fw &&
		       range == o.range &&
		       elt == o.elt &&
		       len == o.len;
	}
	
	/**
	 * Return true iff this GWElt and the given GWElt refer to
	 * different elements.
	 */
	bool operator!=(const GWElt& o) const {
		return !(*this == o);
	}

	index_t offidx; // seed offset index
	bool    fw;     // strand
	index_t range;  // range
	index_t elt;    // element
	index_t len;    // length
};

/**
 * A record encapsulating the result of looking up one BW element in
 * the Bowtie index.
 */
template <typename index_t>
struct WalkResult {

	WalkResult() { reset(); }
	
	/**
	 * Reset GWElt to uninitialized state.
	 */
	void reset() {
		elt.reset();
		bwrow = toff = (index_t)OFF_MASK;
	}

	/**
	 * Initialize this WalkResult.
	 */
	void init(
		index_t oi,  // seed offset index
		bool f,       // strand
		index_t r,   // range
		index_t e,   // element
		index_t bwr, // BW row
		index_t len, // length
		index_t to)  // text offset
	{
		elt.init(oi, f, r, e, len);
		bwrow = bwr;
		toff = to;
	}

	GWElt<index_t> elt;   // element resolved
	index_t        bwrow; // SA row resolved
	index_t        toff;  // resolved offset from SA sample
};

/**
 * A GW hit encapsulates an SATuple describing a reference substring
 * in the cache, along with a bool indicating whether each element of
 * the hit has been reported yet.
 */
template<typename index_t, typename T>
class GWHit {

public:
	GWHit() :
		fmap(0, GW_CAT),
		offidx((index_t)OFF_MASK),
		fw(false),
		range((index_t)OFF_MASK),
		len((index_t)OFF_MASK),
		reported_(0, GW_CAT),
		nrep_(0)
	{
		assert(repOkBasic());
	}

	/**
	 * Initialize with a new SA range.  Resolve the done vector so that
	 * there's one bool per suffix array element.
	 */
	void init(
		SARangeWithOffs<T, index_t>& sa,
		index_t oi,
		bool f,
		index_t r)
	{
		nrep_ = 0;
		offidx = oi;
		fw = f;
		range = r;
		len = (index_t)sa.len;
		reported_.resize(sa.offs.size());
		reported_.fill(false);
		fmap.resize(sa.offs.size());
		fmap.fill(make_pair((index_t)OFF_MASK, (index_t)OFF_MASK));
	}
	
	/**
	 * Clear contents of sat and done.
	 */
	void reset() {
		reported_.clear();
		fmap.clear();
		nrep_ = 0;
		offidx = (index_t)OFF_MASK;
		fw = false;
		range = (index_t)OFF_MASK;
		len = (index_t)OFF_MASK;
	}
	
#ifndef NDEBUG
	/**
	 * Check that GWHit is internally consistent.  If a pointer to an
	 * EList of GWStates is given, we assume that it is the EList
	 * corresponding to this GWHit and check whether the forward and
	 * reverse mappings match up for the as-yet-unresolved elements.
	 */
	bool repOk(const SARangeWithOffs<T, index_t>& sa) const {
		assert_eq(reported_.size(), sa.offs.size());
		assert_eq(fmap.size(), sa.offs.size());
		// Shouldn't be any repeats among as-yet-unresolveds
		size_t nrep = 0;
		for(size_t i = 0; i < fmap.size(); i++) {
			if(reported_[i]) nrep++;
			if(sa.offs[i] != (index_t)OFF_MASK) {
				continue;
			}
			for(size_t j = i+1; j < fmap.size(); j++) {
				if(sa.offs[j] != (index_t)OFF_MASK) {
					continue;
				}
				assert(fmap[i] != fmap[j]);
			}
		}
		assert_eq(nrep_, nrep);
		return true;
	}

	/**
	 * Return true iff this GWHit is not obviously corrupt.
	 */
	bool repOkBasic() {
		return true;
	}
#endif
	
	/**
	 * Set the ith element to be reported.
	 */
	void setReported(index_t i) {
		assert(!reported_[i]);
		assert_lt(i, reported_.size());
		reported_[i] = true;
		nrep_++;
	}
	
	/**
	 * Return true iff element i has been reported.
	 */
	bool reported(index_t i) const {
		assert_lt(i, reported_.size());
		return reported_[i];
	}
	
	/**
	 * Return true iff all elements have been reported.
	 */
	bool done() const {
		assert_leq(nrep_, reported_.size());
		return nrep_ == reported_.size();
	}

	EList<std::pair<index_t, index_t>, 16> fmap; // forward map; to GWState & elt
	index_t offidx; // offset idx
	bool fw;         // orientation
	index_t range;  // original range index
	index_t len;    // length of hit

protected:

	EList<bool, 16> reported_; // per-elt bool indicating whether it's been reported
	index_t nrep_;
};

/**
 * Encapsulates the progress made along a particular path from the original
 * range.
 */
template<typename index_t, typename T>
class GWState {
	
public:

	GWState() : map_(0, GW_CAT) {
		reset(); assert(repOkBasic());
	}
	
	/**
	 * Initialize this GWState with new gfm, top, bot, step, and sat.
	 *
	 * We assume map is already set up.
	 *
	 * Returns true iff at least one elt was resolved.
	 */
	template<int S>
	pair<int, int> init(
		const GFM<index_t>& gfm,      // index to walk left in
		const BitPairReference& ref,  // bitpair-encoded reference
		SARangeWithOffs<T, index_t>& sa,       // SA range with offsets
		EList<GWState, S>& sts,       // EList of GWStates for range being advanced
		GWHit<index_t, T>& hit,       // Corresponding hit structure
		index_t range,                // which range is this?
		bool reportList,              // if true, "report" resolved offsets immediately by adding them to 'res' list
		EList<WalkResult<index_t>, 16>* res,   // EList where resolved offsets should be appended
		index_t tp,                   // top of range at this step
		index_t bt,                   // bot of range at this step
        index_t n_tp,                 // node at top
        index_t n_bt,                 // node at bot
        const EList<pair<index_t, index_t> >& n_iedge_count,
		index_t st,                   // # steps taken to get to this step
		WalkMetrics& met)
	{
		assert_gt(bt, tp);
		assert_lt(range, sts.size());
		top = tp;
		bot = bt;
        node_top = n_tp;
        node_bot = n_bt;
        node_iedge_count = n_iedge_count;
		step = st;
		assert(!inited_);
		ASSERT_ONLY(inited_ = true);
		ASSERT_ONLY(lastStep_ = step-1);
		return init(gfm, ref, sa, sts, hit, range, reportList, res, met);
	}

	/**
	 * Initialize this GWState.
	 *
	 * We assume map is already set up, and that 'step' is equal to the
	 * number of steps taken to get to the new top/bot pair *currently*
	 * in the top and bot fields.
	 *
	 * Returns a pair of numbers, the first being the number of
	 * resolved but unreported offsets found during this advance, the
	 * second being the number of as-yet-unresolved offsets.
	 */
	template<int S>
	pair<int, int> init(
		const GFM<index_t>& gfm,      // forward Bowtie index
		const BitPairReference& ref,  // bitpair-encoded reference
		SARangeWithOffs<T, index_t>& sa,       // SA range with offsets
		EList<GWState, S>& st,        // EList of GWStates for advancing range
		GWHit<index_t, T>& hit,       // Corresponding hit structure
		index_t range,                // range being inited
		bool reportList,              // report resolutions, adding to 'res' list?
		EList<WalkResult<index_t>, 16>* res,   // EList to append resolutions
		WalkMetrics& met)             // update these metrics
	{
		assert(inited_);
		assert_eq(step, lastStep_+1);
		ASSERT_ONLY(lastStep_++);
		assert_leq((index_t)step, gfm.gh().len());
		assert_lt(range, st.size());
        pair<int, int> ret = make_pair(0, 0);
		index_t trimBegin = 0, trimEnd = 0;
		bool empty = true; // assume all resolved until proven otherwise
		// Commit new information, if any, to the PListSlide.  Also,
		// trim and check if we're done.
        assert_eq(node_bot - node_top, map_.size());
        ASSERT_ONLY(index_t num_orig_iedges = 0, orig_e = 0);
        index_t num_iedges = 0, e = 0;
		for(size_t i = mapi_; i < map_.size(); i++) {
			bool resolved = (off((index_t)i, sa) != (index_t)OFF_MASK);
			if(!resolved) {
#ifndef NDEBUG
                while(orig_e < sa.node_iedge_count.size()) {
                    if(map((index_t)i) <= sa.node_iedge_count[orig_e].first) {
                        break;
                    }
                    num_orig_iedges += sa.node_iedge_count[orig_e].second;
                    orig_e++;
                }
#endif
                while(e < node_iedge_count.size()) {
                    if(i <= node_iedge_count[e].first) {
                        break;
                    }
                    num_iedges += node_iedge_count[e].second;
                    e++;
                }
				// Elt not resolved yet; try to resolve it now
				index_t bwrow = (index_t)(top + i + num_iedges);
                index_t node = (index_t)(node_top + i);
				index_t toff = gfm.tryOffset(bwrow, node);
                ASSERT_ONLY(index_t origBwRow = sa.topf + map((index_t)i) + num_orig_iedges);
                ASSERT_ONLY(index_t origNode = sa.node_top + map((index_t)i));
				assert_eq(bwrow, gfm.walkLeft(origBwRow, step));
				if(toff != (index_t)OFF_MASK) {
					// Yes, toff was resolvable
					assert_eq(toff, gfm.getOffset(bwrow, node));
					met.resolves++;
					toff += step;
                    assert_eq(toff, gfm.getOffset(origBwRow, origNode));
					setOff((index_t)i, toff, sa, met);
					if(!reportList) ret.first++;
#if 0
// used to be #ifndef NDEBUG, but since we no longer require that the reference
// string info be included, this is no longer relevant.

					// Sanity check that the reference characters under this
					// hit match the seed characters in hit.satup->key.seq.
					// This is NOT a check that we associated the exact right
					// text offset with the BW row.  This is an important
					// distinction because when resolved offsets are filled in
					// via refernce scanning, they are not necessarily the
					// exact right text offsets to associate with the
					// respective BW rows but they WILL all be correct w/r/t
					// the reference sequence underneath, which is what really
					// matters here.
					index_t tidx = (index_t)OFF_MASK, tof, tlen;
					bool straddled = false;
					gfm.joinedToTextOff(
						hit.len, // length of seed
						toff,    // offset in joined reference string
						tidx,    // reference sequence id
						tof,     // offset in reference coordinates
						tlen,    // length of reference sequence
						true,    // don't reject straddlers
						straddled);
					if(tidx != (index_t)OFF_MASK &&
					   hit.satup->key.seq != std::numeric_limits<uint64_t>::max())
					{
						// key: 2-bit characters packed into a 64-bit word with
						// the least significant bitpair corresponding to the
						// rightmost character on the Watson reference strand.
						uint64_t key = hit.satup->key.seq;
						for(int64_t j = tof + hit.len-1; j >= tof; j--) {
							// Get next reference base to the left
							int c = ref.getBase(tidx, j);
							assert_range(0, 3, c);
							// Must equal least significant bitpair of key
							if(c != (int)(key & 3)) {
								// Oops; when we jump to the piece of the
								// reference where the seed hit is, it doesn't
								// match the seed hit.  Before dying, check
								// whether we have the right spot in the joined
								// reference string
								SString<char> jref;
								gfm.restore(jref);
								uint64_t key2 = hit.satup->key.seq;
								for(int64_t k = toff + hit.len-1; k >= toff; k--) {
									int c = jref[k];
									assert_range(0, 3, c);
									assert_eq(c, (int)(key2 & 3));
									key2 >>= 2;
								}
								assert(false);
							}
							key >>= 2;
						}
					}
#endif
				}
			}
			// Is the element resolved?  We ask this regardless of how it was
			// resolved (whether this function did it just now, whether it did
			// it a while ago, or whether some other function outside GroupWalk
			// did it).
			if(off((index_t)i, sa) != (index_t)OFF_MASK) {
				if(reportList && !hit.reported(map((index_t)i))) {
					// Report it
					index_t toff = off((index_t)i, sa);
					assert(res != NULL);
					res->expand();
					index_t origBwRow = sa.topf + map((index_t)i);
					res->back().init(
						hit.offidx, // offset idx
						hit.fw,     // orientation
						hit.range,  // original range index
						map((index_t)i),     // original element offset
						origBwRow,  // BW row resolved
						hit.len,    // hit length
						toff);      // text offset
					hit.setReported(map((index_t)i));
					met.reports++;
				}
				// Offset resolved
				if(empty) {
					// Haven't seen a non-empty entry yet, so we
					// can trim this from the beginning.
					trimBegin++;
				} else {
					trimEnd++;
				}
			} else {
				// Offset not yet resolved
				ret.second++;
				trimEnd = 0;
				empty = false;
				// Set the forward map in the corresponding GWHit
				// object to point to the appropriate element of our
				// range
				assert_geq(i, mapi_);
				index_t bmap = map((index_t)i);
				hit.fmap[bmap].first = range;
				hit.fmap[bmap].second = (index_t)i;
#ifndef NDEBUG
				for(size_t j = 0; j < bmap; j++) {
					if(sa.offs[j] == (index_t)OFF_MASK &&
					   hit.fmap[j].first == range)
					{
						assert_neq(i, hit.fmap[j].second);
					}
				}
#endif
			}
		}
		// Trim from beginning
		assert_geq(trimBegin, 0);
		mapi_ += trimBegin;
        if(trimBegin > 0) {
            top += trimBegin;
            index_t e = 0;
            for(; e < node_iedge_count.size(); e++) {
                if(node_iedge_count[e].first >= trimBegin) break;
                assert_geq(top, node_iedge_count[e].second);
                top += node_iedge_count[e].second;
            }
            if(e > 0) node_iedge_count.erase(0, e);
            for(e = 0; e < node_iedge_count.size(); e++) {
                assert_geq(node_iedge_count[e].first, trimBegin);
                node_iedge_count[e].first -= trimBegin;
            }
        }
        node_top += trimBegin;
		if(trimEnd > 0) {
			// Trim from end
			map_.resize(map_.size() - trimEnd);
            bot -= trimEnd;
            index_t node_range = node_bot - node_top;
            while(node_iedge_count.size() > 0) {
                if(node_iedge_count.back().first < (node_range - trimEnd)) break;
                assert_geq(bot, node_iedge_count.back().second);
                bot -= node_iedge_count.back().second;
                node_iedge_count.pop_back();
            }
		}
        node_bot -= trimEnd;
#ifndef NDEBUG
        assert_leq(node_top, node_bot);
        index_t num_nodes = node_bot - node_top;
        index_t add = 0;
        for(index_t e = 0; e < node_iedge_count.size(); e++) {
            assert_lt(node_iedge_count[e].first, num_nodes);
            add += node_iedge_count[e].second;
        }
        assert_eq(bot - top, num_nodes + add);
        
#endif
		if(empty) {
			assert(done());
#ifndef NDEBUG
			// If range is done, all elements from map should be
			// resolved
			for(size_t i = mapi_; i < map_.size(); i++) {
				assert_neq((index_t)OFF_MASK, off((index_t)i, sa));
			}
			// If this range is done, then it should be the case that
			// all elements in the corresponding GWHit that point to
			// this range are resolved.
			for(size_t i = 0; i < hit.fmap.size(); i++) {
				if(sa.offs[i] == (index_t)OFF_MASK) {
					assert_neq(range, hit.fmap[i].first);
				}
			}
#endif
			return ret;
		} else {
			assert(!done());
		}
		// Is there a dollar sign in the middle of the range?
        tmp_zOffs.clear();
        for(index_t i = 0; i < gfm._zOffs.size(); i++) {
#ifndef NDEBUG
            if(i > 0) {
                assert_lt(gfm._zOffs[i-1], gfm._zOffs[i]);
            }
#endif
            assert_neq(top, gfm._zOffs[i]);
            // assert_neq(bot-1, gfm._zOffs[i]);
            if(gfm._zOffs[i] > top && gfm._zOffs[i] < bot) {
                tmp_zOffs.push_back(gfm._zOffs[i]);
            }
        }
        
        // Yes, the dollar sign is in the middle of this range.  We
        // must split it into the two ranges on either side of the
        // dollar.  Let 'bot' and 'top' delimit the portion of the
        // range prior to the dollar.
        if(tmp_zOffs.size() > 0) {
            tmp_gbwt_to_node.clear();
            index_t n = 0, e = 0;
            for(index_t r = 0; r < (bot - top); r++) {
                tmp_gbwt_to_node.push_back(n);
                if(e < node_iedge_count.size()) {
                    assert_leq(n, node_iedge_count[e].first);
                    if(n == node_iedge_count[e].first) {
                        for(index_t a = 0; a < node_iedge_count[e].second; a++) {
                            tmp_gbwt_to_node.push_back(n);
                            r++;
                        }
                        e++;
                    }
                }
                n++;
            }
            assert_eq(bot - top, tmp_gbwt_to_node.size());
            for(index_t i = 0; i < tmp_zOffs.size(); i++) {
                assert_lt(top, tmp_zOffs[i]);
                index_t diff = tmp_zOffs[i] - top;
                assert_lt(diff, tmp_gbwt_to_node.size());
                for(index_t j = diff + 1; j < tmp_gbwt_to_node.size(); j++) {
                    if(tmp_gbwt_to_node[i] == tmp_gbwt_to_node[j]) {
                        tmp_gbwt_to_node[j] = (index_t)INDEX_MAX;
                    } else {
                        break;
                    }
                }
                tmp_gbwt_to_node[diff] = (index_t)INDEX_MAX;
            }
            for(index_t i = 0; i < tmp_zOffs.size(); i++) {
                // Note: might be able to do additional trimming off the end.
                // Create a new range for the portion after the dollar.
                index_t new_top = tmp_zOffs[i] + 1;
                while(new_top - top < tmp_gbwt_to_node.size()) {
                    if(tmp_gbwt_to_node[new_top - top] != (index_t)INDEX_MAX) {
                        break;
                    }
                    new_top++;
                }
                assert_leq(new_top - top, tmp_gbwt_to_node.size());
                if(new_top - top == tmp_gbwt_to_node.size()) {
#if 0
                    if(node_iedge_count.size() > 0 &&
                       node_iedge_count.back().first + 1 == node_bot - node_top) {
                        assert_gt(node_iedge_count.back().second, 0);
                        node_iedge_count.back().second -= 1;
                        if(node_iedge_count.back().second == 0) {
                            node_iedge_count.resize(node_iedge_count.size()- 1);
                        }
                    }
#endif
                    break;
                }
                index_t new_node_top = tmp_gbwt_to_node[new_top - top] + node_top;
                assert_lt(new_node_top, node_bot);
                index_t new_bot;
                if(i + 1 < tmp_zOffs.size()) {
                    new_bot = tmp_zOffs[i+1];
                } else {
                    new_bot = bot;
                }
                index_t new_bot2 = new_bot;
                while(new_bot2 - top < tmp_gbwt_to_node.size()) {
                    if(tmp_gbwt_to_node[new_bot2 - top] != (index_t)INDEX_MAX) {
                        break;
                    }
                    new_bot2++;
                }
                index_t new_node_bot = node_bot;
                if(new_bot2 - top < tmp_gbwt_to_node.size()) {
                    new_node_bot = node_top + tmp_gbwt_to_node[new_bot2 - top];
                }
                tmp_node_iedge_count.clear();
                if(new_top >= new_bot) continue;
                for(index_t j = new_top - top; j + 1 < new_bot - top;) {
                    index_t n = tmp_gbwt_to_node[j];
                    index_t j2 = j + 1;
                    while(j2 < new_bot - top) {
                        if(n != tmp_gbwt_to_node[j2]) {
                            break;
                        }
                        j2++;
                    }
                    if(j + 1 < j2) {
                        tmp_node_iedge_count.expand();
                        assert_lt(node_top, new_node_top);
                        tmp_node_iedge_count.back().first = n - (new_node_top - node_top);
                        tmp_node_iedge_count.back().second = j2 - j - 1;
                    }
                    j = j2;
                }
                st.expand();
                st.back().reset();
                st.back().initMap(new_node_bot - new_node_top);
                for(index_t j = new_node_top; j < new_node_bot; j++) {
                    st.back().map_[j - new_node_top] = map(j - node_top + mapi_);
                }
                st.back().init(
                               gfm,
                               ref,
                               sa,
                               st,
                               hit,
                               (index_t)st.size()-1,
                               reportList,
                               res,
                               new_top,
                               new_bot,
                               new_node_top,
                               new_node_bot,
                               tmp_node_iedge_count,
                               step,
                               met);
            }
            assert_eq((index_t)map_.size(), node_bot - node_top + mapi_);
            bot = tmp_zOffs[0];
            assert_lt(bot - top, tmp_gbwt_to_node.size());
            node_bot = tmp_gbwt_to_node[bot - top - 1] + node_top + 1;
            map_.resize(node_bot - node_top + mapi_);
            index_t width = node_bot - node_top;
            for(index_t e = 0; e < node_iedge_count.size(); e++) {
                if(node_iedge_count[e].first >= node_bot - node_top) {
                    node_iedge_count.resize(e);
                    break;
                }
                width += node_iedge_count[e].second;
            }
            if(width != bot - top) {
                assert_eq(width, bot - top + 1);
                assert_gt(node_iedge_count.size(), 0);
                assert_gt(node_iedge_count.back().second, 0);
                node_iedge_count.back().second -= 1;
                if(node_iedge_count.back().second == 0) {
                    node_iedge_count.resize(node_iedge_count.size()- 1);
                }
            }
        }
		assert_gt(bot, top);
		// Prepare SideLocus's for next step
		if(bot-top > 1) {
			SideLocus<index_t>::initFromTopBot(top, bot, gfm.gh(), gfm.gfm(), tloc, bloc);
			assert(tloc.valid()); assert(tloc.repOk(gfm.gh()));
			assert(bloc.valid()); assert(bloc.repOk(gfm.gh()));
		} else {
			tloc.initFromRow(top, gfm.gh(), gfm.gfm());
			assert(tloc.valid()); assert(tloc.repOk(gfm.gh()));
			bloc.invalidate();
		}
		return ret;
	}
	
#ifndef NDEBUG
	/**
	 * Check if this GWP is internally consistent.
	 */
	bool repOk(
		const GFM<index_t>& gfm,
		GWHit<index_t, T>& hit,
		index_t range) const
	{
		assert(done() || bot > top);
		assert(doneResolving(hit) || (tloc.valid() && tloc.repOk(gfm.gh())));
		assert(doneResolving(hit) || bot == top+1 || (bloc.valid() && bloc.repOk(gfm.gh())));
		assert_eq(map_.size()-mapi_, bot-top);
		// Make sure that 'done' is compatible with whether we have >=
		// 1 elements left to resolve.
		int left = 0;
		for(size_t i = mapi_; i < map_.size(); i++) {
			ASSERT_ONLY(index_t row = (index_t)(top + i - mapi_));
			ASSERT_ONLY(index_t origRow = hit.satup->topf + map(i));
			assert(step == 0 || row != origRow);
			assert_eq(row, gfm.walkLeft(origRow, step));
			assert_lt(map_[i], hit.satup->offs.size());
			if(off(i, hit) == (index_t)OFF_MASK) left++;
		}
		assert(repOkMapRepeats());
		assert(repOkMapInclusive(hit, range));
		return true;
	}
	
	/**
	 * Return true iff this GWState is not obviously corrupt.
	 */
	bool repOkBasic() {
		assert_geq(bot, top);
		return true;
	}

	/**
	 * Check that the fmap elements pointed to by our map_ include all
	 * of the fmap elements that point to this range.
	 */
	bool repOkMapInclusive(GWHit<index_t, T>& hit, index_t range) const {
		for(size_t i = 0; i < hit.fmap.size(); i++) {
			if(hit.satup->offs[i] == (index_t)OFF_MASK) {
				if(range == hit.fmap[i].first) {
					ASSERT_ONLY(bool found = false);
					for(size_t j = mapi_; j < map_.size(); j++) {
						if(map(j) == i) {
							ASSERT_ONLY(found = true);
							break;
						}
					}
					assert(found);
				}
			}
		}
		return true;
	}
	
	/**
	 * Check that no two elements in map_ are the same.
	 */
	bool repOkMapRepeats() const {
		for(size_t i = mapi_; i < map_.size(); i++) {
			for(size_t j = i+1; j < map_.size(); j++) {
				assert_neq(map_[i], map_[j]);
			}
		}
		return true;
	}
#endif
	
	/**
	 * Return the offset currently assigned to the ith element.  If it
	 * has not yet been resolved, return 0xffffffff.
	 */
	index_t off(
				index_t i,
				const SARangeWithOffs<T, index_t>& sa)
	{
		assert_geq(i, mapi_);
		assert_lt(i, map_.size());
		assert_lt(map_[i], sa.offs.size());
		return sa.offs.get(map_[i]);
	}

	/**
	 * Return the offset of the element within the original range's
	 * PListSlice that the ith element of this range corresponds to.
	 */
	index_t map(index_t i) const {
		assert_geq(i, mapi_);
		assert_lt(i, map_.size());
		return map_[i];
	}

	/**
	 * Return the offset of the first untrimmed offset in the map.
	 */
	index_t mapi() const {
		return mapi_;
	}

	/**
	 * Return number of active elements in the range being tracked by
	 * this GWState.
	 */
	index_t size() const {
		return (index_t)(map_.size() - mapi_);
	}
	
	/**
	 * Return true iff all elements in this leaf range have been
	 * resolved.
	 */
	bool done() const {
		return size() == 0;
	}

	/**
	 * Set the PListSlice element that corresponds to the ith element
	 * of 'map' to the specified offset.
	 */
	void setOff(
		index_t i,
		index_t off,
		SARangeWithOffs<T, index_t>& sa,
		WalkMetrics& met)
	{
		assert_lt(i + mapi_, map_.size());
		assert_lt(map_[i + mapi_], sa.offs.size());
		size_t saoff = map_[i + mapi_];
		sa.offs[saoff] = off;
		assert_eq(off, sa.offs[saoff]);
	}

	/**
	 * Advance this GWState by one step (i.e. one BW operation).  In
	 * the event of a "split", more elements are added to the EList
	 * 'st', which must have room for at least 3 more elements without
	 * needing another expansion.  If an expansion of 'st' is
	 * triggered, this GWState object becomes invalid.
	 *
	 * Returns a pair of numbers, the first being the number of
	 * resolved but unreported offsets found during this advance, the
	 * second being the number of as-yet-unresolved offsets.
	 */
	template <int S>
	pair<int, int> advance(
		const GFM<index_t>& gfm,     // the forward Bowtie index, for stepping left
		const BitPairReference& ref, // bitpair-encoded reference
		SARangeWithOffs<T, index_t>& sa,      // SA range with offsets
		GWHit<index_t, T>& hit,      // the associated GWHit object
		index_t range,               // which range is this?
		bool reportList,             // if true, "report" resolved offsets immediately by adding them to 'res' list
		EList<WalkResult<index_t>, 16>* res,  // EList where resolved offsets should be appended
		EList<GWState, S>& st,       // EList of GWStates for range being advanced
		GroupWalkState<index_t>& gws,         // temporary storage for masks
		WalkMetrics& met,
		PerReadMetrics& prm)
	{
		ASSERT_ONLY(index_t origTop = top);
		ASSERT_ONLY(index_t origBot = bot);
		assert_geq(step, 0);
		assert_eq(step, lastStep_);
		// assert_geq(st.capacity(), st.size() + 4);
		assert(tloc.valid()); assert(tloc.repOk(gfm.gh()));
		assert_eq(node_bot-node_top, (index_t)(map_.size()-mapi_));
		pair<int, int> ret = make_pair(0, 0);
		assert_eq(top, tloc.toBWRow(gfm.gh()));
		if(bot - top > 1) {
            bool first = true;
            ASSERT_ONLY(index_t sum = 0);
            index_t newtop = 0, newbot = 0;
            index_t new_node_top = 0, new_node_bot = 0;
            gws.map.clear();
			// Still multiple elements being tracked
            index_t curtop = top, curbot = bot;
            index_t cur_node_top = node_top, cur_node_bot = node_bot;
            for(index_t e = 0; e < node_iedge_count.size() + 1; e++) {
                if(e >= node_iedge_count.size()) {
                    if(e > 0) {
                        curtop = curbot + node_iedge_count[e-1].second;
                        curbot = bot;
                        if(curtop >= curbot) {
                            assert_eq(curtop, curbot);
                            break;
                        }
                        cur_node_top = cur_node_bot;
                        cur_node_bot = node_bot;
                    }
                } else {
                    if(e > 0) {
                        curtop = curbot + node_iedge_count[e-1].second;
                        assert_lt(node_iedge_count[e-1].first, node_iedge_count[e].first);
                        curbot = curtop + (node_iedge_count[e].first - node_iedge_count[e-1].first);
                        cur_node_top = cur_node_bot;
                    } else {
                        curbot = curtop + node_iedge_count[e].first + 1;
                    }
                    cur_node_bot = node_top + node_iedge_count[e].first + 1;
                }
                assert_lt(curtop, curbot);
                index_t upto[4], in[4];
                upto[0] = in[0] = upto[1] = in[1] =
                upto[2] = in[2] = upto[3] = in[3] = 0;
                // assert_eq(bot, bloc.toBWRow(gfm.gh()));
                met.bwops++;
                prm.nExFmops++;
                // Assert that there's not a dollar sign in the middle of
                // this range
#ifndef NDEBUG
                for(index_t i = 0; i < gfm._zOffs.size(); i++) {
                    assert(curbot <= gfm._zOffs[i] || curtop > gfm._zOffs[i]);
                }
#endif
                SideLocus<index_t> curtloc, curbloc;
                SideLocus<index_t>::initFromTopBot(curtop, curbot, gfm.gh(), gfm.gfm(), curtloc, curbloc);
                gfm.mapLFRange(curtloc, curbloc, curbot-curtop, upto, in, gws.masks);
#ifndef NDEBUG
                for(int i = 0; i < 4; i++) {
                    assert_eq(curbot-curtop, (index_t)(gws.masks[i].size()));
                }
#endif
                
                for(int i = 0; i < 4; i++) {
                    if(in[i] > 0) {
                        // Non-empty range resulted
                        if(first) {
                            // For the first one,
                            first = false;
                            pair<index_t, index_t> range, node_range;
                            backup_node_iedge_count.clear();
                            SideLocus<index_t>::initFromTopBot(curtop, curbot, gfm.gh(), gfm.gfm(), curtloc, curbloc);
                            range = gfm.mapGLF(curtloc, curbloc, i, &node_range, &backup_node_iedge_count, cur_node_bot - cur_node_top);
                            newtop = range.first;
                            newbot = range.second;
                            new_node_top = node_range.first;
                            new_node_bot = node_range.second;
                            // Range narrowed so we have to look at the masks
                            for(size_t j = 0; j < gws.masks[i].size(); j++) {
                                assert_lt(j+mapi_+(cur_node_top - node_top), map_.size());
                                if(gws.masks[i][j]) {
                                    gws.map.push_back(map_[j+mapi_+(cur_node_top - node_top)]);
                                    assert(gws.map.size() <= 1 || gws.map.back() != gws.map[gws.map.size()-2]);
#if 0
                                    // If this element is not yet resolved,
                                    // then check that it really is the
                                    // expected number of steps to the left
                                    // of the corresponding element in the
                                    // root range
                                    assert_lt(gws.map.back(), sa.size());
                                    if(sa.offs[gws.map.back()] == (index_t)OFF_MASK) {
                                        assert_eq(newtop + gws.map.size() - 1,
                                                  gfm.walkLeft(sa.topf + gws.map.back(), step+1));
                                    }
#endif
                                }
                            }
                            assert_lt(new_node_top, new_node_bot);
                            if(new_node_bot - new_node_top < gws.map.size()) {
                                assert_eq(curbot - curtop, cur_node_bot - cur_node_top);
                                SideLocus<index_t> tmptloc, tmpbloc;
                                pair<index_t, index_t> tmp_node_range;
                                index_t j1 = 0, j2 = 0;
                                for(index_t c = 0; c < gws.masks[i].size(); c++) {
                                    if(gws.masks[i][c]) {
                                        j1 = c;
                                        break;
                                    }
                                }
                                for(index_t j = 0; j + 1 < gws.map.size(); j++) {
                                    for(index_t c = j1 + 1; c < gws.masks[i].size(); c++) {
                                        if(gws.masks[i][c]) {
                                            j2 = c;
                                            break;
                                        }
                                    }
                                    assert_lt(j1, j2);
                                    SideLocus<index_t>::initFromTopBot(curtop + j1, curtop + j2 + 1, gfm.gh(), gfm.gfm(), tmptloc, tmpbloc);
                                    gfm.mapGLF(tmptloc, tmpbloc, i, &tmp_node_range);
                                    assert_gt(tmp_node_range.second - tmp_node_range.first, 0);
                                    if(tmp_node_range.second - tmp_node_range.first == 1) {
                                        index_t jmap = gws.map[j];
                                        assert_lt(jmap, sa.offs.size());
                                        sa.offs[jmap] = gws.map[j];
                                        gws.map[j] = (index_t)OFF_MASK;
                                    }
                                    j1 = j2;
                                    j2 = 0;
                                }
                                for(index_t j = 0; j < gws.map.size();) {
                                    if(gws.map[j] == (index_t)OFF_MASK) {
                                        gws.map.erase(j);
                                    } else j++;
                                }
#ifndef NDEBUG
                                for(index_t j = 0; j < gws.map.size(); j++) {
                                    assert_neq(gws.map[j], (index_t)OFF_MASK);
                                }
#endif
                            }
                            assert_eq(new_node_bot - new_node_top, (index_t)(gws.map.size()));
                        } else {
                            // For each beyond the first, create a new
                            // GWState and add it to the GWState list.
                            // NOTE: this can cause the underlying list to
                            // be expanded which in turn might leave 'st'
                            // pointing to bad memory.
                            st.expand();
                            st.back().reset();
                            tmp_node_iedge_count.clear();
                            pair<index_t, index_t> range, node_range;
                            SideLocus<index_t>::initFromTopBot(curtop, curbot, gfm.gh(), gfm.gfm(), curtloc, curbloc);
                            range = gfm.mapGLF(curtloc, curbloc, i, &node_range, &tmp_node_iedge_count, cur_node_bot - cur_node_top);
                            assert_geq(range.second - range.first, node_range.second - node_range.first);
                            index_t ntop = range.first;
                            index_t nbot = range.second;
                            st.back().mapi_ = 0;
                            st.back().map_.clear();
                            met.branches++;
                            // Range narrowed so we have to look at the masks
                            for(size_t j = 0; j < gws.masks[i].size(); j++) {
                                if(gws.masks[i][j]) st.back().map_.push_back(map_[j+mapi_+(cur_node_top - node_top)]);
                            }
                            assert_lt(node_range.first, node_range.second);
                            if(node_range.second - node_range.first < st.back().map_.size()) {
                                assert_eq(curbot - curtop, cur_node_bot - cur_node_top);
                                SideLocus<index_t> tmptloc, tmpbloc;
                                pair<index_t, index_t> tmp_node_range;
                                index_t j1 = 0, j2 = 0;
                                for(index_t c = 0; c < gws.masks[i].size(); c++) {
                                    if(gws.masks[i][c]) {
                                        j1 = c;
                                        break;
                                    }
                                }
                                for(index_t j = 0; j + 1 < st.back().map_.size(); j++) {
                                    for(index_t c = j1 + 1; c < gws.masks[i].size(); c++) {
                                        if(gws.masks[i][c]) {
                                            j2 = c;
                                            break;
                                        }
                                    }
                                    assert_lt(j1, j2);
                                    SideLocus<index_t>::initFromTopBot(curtop + j1, curtop + j2 + 1, gfm.gh(), gfm.gfm(), tmptloc, tmpbloc);
                                    gfm.mapGLF(tmptloc, tmpbloc, i, &tmp_node_range);
                                    assert_gt(tmp_node_range.second - tmp_node_range.first, 0);
                                    if(tmp_node_range.second - tmp_node_range.first == 1) {
                                        index_t jmap = st.back().map_[j];
                                        assert_lt(jmap, sa.offs.size());
                                        sa.offs[jmap] = st.back().map_[j];
                                        st.back().map_[j] = (index_t)OFF_MASK;
                                    }
                                    j1 = j2;
                                    j2 = 0;
                                }
                                for(index_t j = 0; j < st.back().map_.size();) {
                                    if(st.back().map_[j] == (index_t)OFF_MASK) {
                                        st.back().map_.erase(j);
                                    } else j++;
                                }
#ifndef NDEBUG
                                for(index_t j = 0; j < st.back().map_.size(); j++) {
                                    assert_neq(st.back().map_[j], (index_t)OFF_MASK);
                                }
#endif
                            }
                            assert_eq(node_range.second - node_range.first, st.back().map_.size());
                            pair<int, int> rret =
                            st.back().init(
                                           gfm,         // forward Bowtie index
                                           ref,         // bitpair-encodede reference
                                           sa,          // SA range with offsets
                                           st,          // EList of all GWStates associated with original range
                                           hit,         // associated GWHit object
                                           (index_t)st.size()-1, // range offset
                                           reportList,  // if true, report hits to 'res' list
                                           res,         // report hits here if reportList is true
                                           ntop,        // BW top of new range
                                           nbot,        // BW bot of new range
                                           node_range.first,
                                           node_range.second,
                                           tmp_node_iedge_count,
                                           step+1,      // # steps taken to get to this new range
                                           met);        // update these metrics
                            ret.first += rret.first;
                            ret.second += rret.second;
                        }
                        ASSERT_ONLY(sum += in[i]);
                    }
                }
            }
            mapi_ = 0;
            // assert_eq(new_node_bot-new_node_top, sum);
            assert_gt(newbot, newtop);
            assert(top != newtop || bot != newbot);
            //assert(!(newtop < top && newbot > top));
            top = newtop;
            bot = newbot;
            node_top = new_node_top;
            node_bot = new_node_bot;
            node_iedge_count = backup_node_iedge_count;
            backup_node_iedge_count.clear();
            if(!gws.map.empty()) {
                map_ = gws.map;
            }
            //assert(repOkMapRepeats());
            //assert(repOkMapInclusive(hit, range));
            assert_eq(node_bot-node_top, (index_t)map_.size());
        } else {
            // Down to one element
			assert_eq(bot, top+1);
			assert_eq(1, map_.size()-mapi_);
			// Sets top, returns char walked through (which we ignore)
			ASSERT_ONLY(index_t oldtop = top);
			met.bwops++;
			prm.nExFmops++;
            pair<index_t, index_t> node_range(0, 0);
			pair<index_t, index_t> range = gfm.mapGLF1(top, tloc, &node_range);
            top = range.first;
			assert_neq(top, oldtop);
			bot = top+1;
            node_top = node_range.first;
            node_bot = node_range.second;
			if(mapi_ > 0) {
				map_[0] = map_[mapi_];
				mapi_ = 0;
			}
			map_.resize(1);
		}
		assert(top != origTop || bot != origBot);
		step++;
		assert_gt(step, 0);
		assert_leq((index_t)step, gfm.gh().len());
		pair<int, int> rret =
		init<S>(
			gfm,        // forward GFM index
			ref,        // bitpair-encodede reference
			sa,         // SA range with offsets
			st,         // EList of all GWStates associated with original range
			hit,        // associated GWHit object
			range,      // range offset
			reportList, // if true, report hits to 'res' list
			res,        // report hits here if reportList is true
			met);       // update these metrics
		ret.first += rret.first;
		ret.second += rret.second;
		return ret;
	}

	/**
	 * Clear all state in preparation for the next walk.
	 */
	void reset() {
		top = bot = node_top = node_bot = step = mapi_ = 0;
		ASSERT_ONLY(lastStep_ = -1);
		ASSERT_ONLY(inited_ = false);
		tloc.invalidate();
		bloc.invalidate();
		map_.clear();
        node_iedge_count.clear();
        backup_node_iedge_count.clear();
        tmp_node_iedge_count.clear();
	}
	
	/**
	 * Resize the map_ field to the given size.
	 */
	void initMap(size_t newsz) {
		mapi_ = 0;
		map_.resize(newsz);
		for(size_t i = 0; i < newsz; i++) {
			map_[i] = (index_t)i;
		}
	}

	/**
	 * Return true iff all rows corresponding to this GWState have been
	 * resolved and reported.
	 */
	bool doneReporting(const GWHit<index_t, T>& hit) const {
		for(size_t i = mapi_; i < map_.size(); i++) {
			if(!hit.reported(map(i))) return false;
		}
		return true;
	}

	/**
	 * Return true iff all rows corresponding to this GWState have been
	 * resolved (but not necessarily reported).
	 */
	bool doneResolving(const SARangeWithOffs<T, index_t>& sa) const {
		for(size_t i = mapi_; i < map_.size(); i++) {
			if(sa.offs[map((index_t)i)] == (index_t)OFF_MASK) return false;
		}
		return true;
	}

	SideLocus<index_t> tloc;      // SideLocus for top
	SideLocus<index_t> bloc;      // SideLocus for bottom
	index_t            top;       // top elt of range in GBWT
	index_t            bot;       // bot elt of range in GBWT
    index_t            node_top;
    index_t            node_bot;
    EList<pair<index_t, index_t> > node_iedge_count;
	int                step;      // how many steps have we walked to the left so far
    
    // temporary
    EList<pair<index_t, index_t> > backup_node_iedge_count;
    EList<pair<index_t, index_t> > tmp_node_iedge_count;
    
    EList<index_t>                 tmp_zOffs;
    EList<index_t>                 tmp_gbwt_to_node;

protected:
	
	ASSERT_ONLY(bool inited_);
	ASSERT_ONLY(int lastStep_);
	EList<index_t, 16> map_; // which elts in range 'range' we're tracking
	index_t mapi_;           // first untrimmed element of map
};

template<typename index_t, typename T, int S>
class GroupWalk2S {
public:
	typedef EList<GWState<index_t, T>, S> TStateV;

	GroupWalk2S() : st_(8, GW_CAT) {
		reset();
	}
	
	/**
	 * Reset the GroupWalk in preparation for the next SeedResults.
	 */
	void reset() {
		elt_ = rep_ = 0;
		ASSERT_ONLY(inited_ = false);
	}

	/**
	 * Initialize a new group walk w/r/t a QVal object.
	 */
	void init(
		const GFM<index_t>& gfmFw,   // forward Bowtie index for walking left
		const BitPairReference& ref, // bitpair-encoded reference
		SARangeWithOffs<T, index_t>& sa,      // SA range with offsets
		RandomSource& rnd,           // pseudo-random generator for sampling rows
		WalkMetrics& met)            // update metrics here
	{
		reset();
#ifndef NDEBUG
		inited_ = true;
#endif
		// Init GWHit
		hit_.init(sa, 0, false, 0);
		// Init corresponding GWState
		st_.resize(1);
		st_.back().reset();
		assert(st_.back().repOkBasic());
		index_t top = sa.topf;
		index_t bot = sa.botf;
        index_t node_top = sa.node_top;
        index_t node_bot = (index_t)(node_top + sa.size());
		st_.back().initMap(sa.size());
		st_.ensure(4);
		st_.back().init(
			gfmFw,              // Bowtie index
			ref,                // bitpair-encoded reference
			sa,                 // SA range with offsets
			st_,                // EList<GWState>
			hit_,               // GWHit
			0,                  // range 0
			false,              // put resolved elements into res_?
			NULL,               // put resolved elements here
			top,                // GBW row at top
			bot,                // GBW row at bot
            node_top,           // node at top
            node_bot,           // node at bot
            sa.node_iedge_count,
			0,                  // # steps taken
			met);               // update metrics here
		elt_ += sa.size();
		assert(hit_.repOk(sa));
	}

	//
	// ELEMENT-BASED
	//

	/**
	 * Advance the GroupWalk until all elements have been resolved.
	 */
	void resolveAll(WalkMetrics& met, PerReadMetrics& prm) {
		WalkResult<index_t> res; // ignore results for now
		for(size_t i = 0; i < elt_; i++) {
			advanceElement((index_t)i, res, met, prm);
		}
	}

	/**
	 * Advance the GroupWalk until the specified element has been
	 * resolved.
	 */
	bool advanceElement(
		index_t elt,                  // element within the range
		const GFM<index_t>& gfmFw,    // forward Bowtie index for walking left
		const BitPairReference& ref,  // bitpair-encoded reference
		SARangeWithOffs<T, index_t>& sa,       // SA range with offsets
		GroupWalkState<index_t>& gws, // GroupWalk state; scratch space
		WalkResult<index_t>& res,     // put the result here
		WalkMetrics& met,             // metrics
		PerReadMetrics& prm)          // per-read metrics
	{
		assert(inited_);
		assert(!done());
		assert(hit_.repOk(sa));
		assert_lt(elt, sa.size()); // elt must fall within range
		// Until we've resolved our element of interest...
		while(sa.offs[elt] == (index_t)OFF_MASK) {
			// Get the GWState that contains our element of interest
			size_t range = hit_.fmap[elt].first;
            assert_lt(range, st_.size());
			st_.ensure(st_[range].node_bot - st_[range].node_top);
            // st_.ensure(4);
			GWState<index_t, T>& st = st_[range];
			assert(!st.doneResolving(sa));
			// Returns a pair of numbers, the first being the number of
			// resolved but unreported offsets found during this advance, the
			// second being the number of as-yet-unresolved offsets.
			st.advance(
				gfmFw,
				ref,
				sa,
				hit_,
				(index_t)range,
				false,
				NULL,
				st_,
				gws,
				met,
				prm);
			assert(sa.offs[elt] != (index_t)OFF_MASK ||
			       !st_[hit_.fmap[elt].first].doneResolving(sa));
		}
		assert_neq((index_t)OFF_MASK, sa.offs[elt]);
		// Report it!
		if(!hit_.reported(elt)) {
			hit_.setReported(elt);
		}
		met.reports++;
		res.init(
			0,              // seed offset
			false,          // orientation
			0,              // range
			elt,            // element
			sa.topf + elt,  // bw row
			(index_t)sa.len, // length of hit
			sa.offs[elt]);  // resolved text offset
		rep_++;
		return true;
	}

	/**
	 * Return true iff all elements have been resolved and reported.
	 */
	bool done() const { return rep_ == elt_; }
	
#ifndef NDEBUG
	/**
	 * Check that GroupWalk is internally consistent.
	 */
	bool repOk(const SARangeWithOffs<T, index_t>& sa) const {
		assert(hit_.repOk(sa));
		assert_leq(rep_, elt_);
		// This is a lot of work
		size_t resolved = 0, reported = 0;
		// For each element
		const size_t sz = sa.size();
		for(size_t m = 0; m < sz; m++) {
			// Is it resolved?
			if(sa.offs[m] != (index_t)OFF_MASK) {
				resolved++;
			} else {
				assert(!hit_.reported(m));
			}
			// Is it reported?
			if(hit_.reported(m)) {
				reported++;
			}
			assert_geq(resolved, reported);
		}
		assert_geq(resolved, reported);
		assert_eq(rep_, reported);
		assert_eq(elt_, sz);
		return true;
	}
#endif

	/**
	 * Return the number of BW elements that we can resolve.
	 */
	index_t numElts() const { return elt_; }
	
	/**
	 * Return the size occupied by this GroupWalk and all its constituent
	 * objects.
	 */
	size_t totalSizeBytes() const {
		return 2 * sizeof(size_t) + st_.totalSizeBytes() + sizeof(GWHit<index_t, T>);
	}
	/**
	 * Return the capacity of this GroupWalk and all its constituent objects.
	 */
	size_t totalCapacityBytes() const {
		return 2 * sizeof(size_t) + st_.totalCapacityBytes() + sizeof(GWHit<index_t, T>);
	}
	
#ifndef NDEBUG
	bool initialized() const { return inited_; }
#endif
	
protected:

	ASSERT_ONLY(bool inited_);    // initialized?
	
	index_t elt_;    // # BW elements under the control of the GropuWalk
	index_t rep_;    // # BW elements reported

	// For each orientation and seed offset, keep a GWState object that
	// holds the state of the walk so far.
	TStateV st_;

	// For each orientation and seed offset, keep an EList of GWHit.
	GWHit<index_t, T> hit_;
};

#endif /*GROUP_WALK_H_*/
