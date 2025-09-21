/*
 * Copyright 2015, Daehwan Kim <infphilo@gmail.com>
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

#ifndef HI_ALIGNER_H_
#define HI_ALIGNER_H_

#include <iostream>
#include <utility>
#include <limits>
#include "qual.h"
#include "ds.h"
#include "sstring.h"
#include "alphabet.h"
#include "edit.h"
#include "read.h"
// Threading is necessary to synchronize the classes that dump
// intermediate alignment results to files.  Otherwise, all data herein
// is constant and shared, or per-thread.
#include "threading.h"
#include "aligner_result.h"
#include "aligner_cache.h"
#include "scoring.h"
#include "mem_ids.h"
#include "simple_func.h"
#include "aligner_driver.h"
#include "aligner_sw_driver.h"
#include "group_walk.h"
#include "tp.h"
#include "gp.h"

#if defined(_MSC_VER)
#pragma warning(disable:4146)
#endif

// Allow longer introns for long anchored reads involving canonical splice sites
inline uint32_t MaxIntronLen(uint32_t anchor, uint32_t minAnchorLen) {
    uint32_t intronLen = 0;
    if(anchor >= minAnchorLen) {
        if(anchor < 2) anchor = 2;
        uint32_t shift = (anchor << 1) - 4;
        shift = min<uint32_t>(max<uint32_t>(shift, 13), 30);
        intronLen = 1 << shift;
    }
    return intronLen;
}

inline float intronLen_prob(uint32_t anchor, uint32_t intronLen, uint32_t maxIntronLen) {
    uint32_t expected_intron_len = maxIntronLen;
    if(anchor < 14) expected_intron_len = 1 << ((anchor << 1) + 4);
    if(expected_intron_len > maxIntronLen) expected_intron_len = maxIntronLen;
    assert_gt(expected_intron_len, 0);
    float result = ((float)intronLen) / ((float)expected_intron_len);
    if(result > 1.0f) result = 1.0f;
    return result;
}

// Allow longer introns for long anchored reads involving non-canonical splice sites
inline uint32_t MaxIntronLen_noncan(uint32_t anchor, uint32_t minAnchorLen_noncan) {
    uint32_t intronLen = 0;
    if(anchor >= minAnchorLen_noncan) {
        if(anchor < 5) anchor = 5;
        uint32_t shift = (anchor << 1) - 10;
        shift = min<uint32_t>(shift, 30);
        intronLen = 1 << shift;
    }
    return intronLen;
}

inline float intronLen_prob_noncan(uint32_t anchor, uint32_t intronLen, uint32_t maxIntronLen) {
    uint32_t expected_intron_len = maxIntronLen;
    if(anchor < 16) expected_intron_len = 1 << (anchor << 1);
    if(expected_intron_len > maxIntronLen) expected_intron_len = maxIntronLen;
    assert_gt(expected_intron_len, 0);
    float result = ((float)intronLen) / ((float)expected_intron_len);
    if(result > 1.0f) result = 1.0f;
    return result;
}

/**
 * Hit types for BWTHit class below
 * Three hit types to anchor a read on the genome
 *
 */
enum {
    CANDIDATE_HIT = 1,
    PSEUDOGENE_HIT,
    ANCHOR_HIT,
};

/**
 * Simple struct for holding a partial alignment for the read
 * The alignment locations are represented by FM offsets [top, bot),
 * and later genomic offsets are calculated when necessary
 */
template <typename index_t>
struct BWTHit {
	
	BWTHit() { reset(); }
	
	void reset() {
		_top = _bot = 0;
        _node_top = _node_bot = 0;
        _node_iedge_count.clear();
		_fw = true;
		_bwoff = (index_t)INDEX_MAX;
		_len = 0;
		_coords.clear();
        _anchor_examined = false;
        _hit_type = CANDIDATE_HIT;
	}
	
	void init(
			  index_t top,
			  index_t bot,
              index_t node_top,
              index_t node_bot,
              const EList<pair<index_t, index_t> >& node_iedge_count,
  			  bool fw,
			  uint32_t bwoff,
			  uint32_t len,
              index_t hit_type = CANDIDATE_HIT)
	{
        assert_leq(node_bot - node_top, bot - top);
#ifndef NDEBUG
        if(node_bot - node_top < bot - top) {
            assert_gt(node_iedge_count.size(), 0);
        }
#endif
		_top = top;
        _bot = bot;
        _node_top = node_top;
        _node_bot = node_bot;
        _node_iedge_count = node_iedge_count;
		_fw = fw;
		_bwoff = bwoff;
		_len = len;
        _coords.clear();
        _anchor_examined = false;
        _hit_type = hit_type;
	}
    
    bool hasGenomeCoords() const { return !_coords.empty(); }
	
	/**
	 * Return true iff there is no hit.
	 */
	bool empty() const {
		return _bot <= _top;
	}
	
	/**
	 * Higher score = higher priority.
	 */
	bool operator<(const BWTHit& o) const {
		return _len > o._len;
	}
	
	/**
	 * Return the size of the alignments SA ranges.
	 */
	index_t size() const {
        assert_leq(_top, _bot);
        return _bot - _top;
    }
    
    index_t len() const {
        assert_gt(_len, 0);
        return _len;
    }
	
#ifndef NDEBUG
	/**
	 * Check that hit is sane w/r/t read.
	 */
	bool repOk(const Read& rd) const {
		assert_gt(_bot, _top);
		assert_neq(_bwoff, (index_t)INDEX_MAX);
		assert_gt(_len, 0);
		return true;
	}
#endif
	
	index_t         _top;               // start of the range in the FM index
	index_t         _bot;               // end of the range in the FM index
    index_t         _node_top;
    index_t         _node_bot;
    EList<pair<index_t, index_t> >  _node_iedge_count;
	bool            _fw;                // whether read is forward or reverse complemented
	index_t         _bwoff;             // current base of a read to search from the right end
	index_t         _len;               // read length
	
    EList<Coord>    _coords;            // genomic offsets corresponding to [_top, _bot)
    
    bool            _anchor_examined;   // whether or not this hit is examined
    index_t         _hit_type;          // hit type (anchor hit, pseudogene hit, or candidate hit)
};


/**
 * Simple struct for holding alignments for the read
 * The alignments are represented by chains of BWTHits
 */
template <typename index_t>
struct ReadBWTHit {
	
	ReadBWTHit() { reset(); }
	
	void reset() {
        _fw = true;
		_len = 0;
        _cur = 0;
        _done = false;
        _numPartialSearch = 0;
        _numUniqueSearch = 0;
        _repeat = false;
        _partialHits.clear();
	}

	void init(
			  bool fw,
              index_t len)
	{
        _fw = fw;
        assert_gt(len, 0);
        _len = len;
        _cur = 0;
        _done = false;
        _numPartialSearch = 0;
        _numUniqueSearch = 0;
        _repeat = false;
        _partialHits.clear();
	}
    
    bool done() {
#ifndef NDEBUG
        assert_gt(_len, 0);
        if(_cur >= _len) {
            assert(_done);
        }
#endif
        return _done;
    }
    
    void done(bool done) {
        assert(!_done);
        assert(done);
        _done = done;
    }
    
    index_t len() const { return _len; }
    index_t cur() const { return _cur; }
    bool    repeat() const { return _repeat; }
    
    index_t offsetSize()             { return (index_t)_partialHits.size(); }
    size_t  numPartialSearch()       { return _numPartialSearch; }
    index_t numActualPartialSearch()
    {
        assert_leq(_numUniqueSearch, _numPartialSearch);
        return (index_t)(_numPartialSearch - _numUniqueSearch);
    }
    
    bool width(index_t offset_) {
        assert_lt(offset_, _partialHits.size());
        return _partialHits[offset_].size();
    }
    
    bool hasGenomeCoords(index_t offset_) {
        assert_lt(offset_, _partialHits.size());
        index_t width_ = width(offset_);
        if(width_ == 0) {
            return true;
        } else {
            return _partialHits[offset_].hasGenomeCoords();
        }
    }
    
    bool hasAllGenomeCoords() {
        if(_cur < _len) return false;
        if(_partialHits.size() <= 0) return false;
        for(size_t oi = 0; oi < _partialHits.size(); oi++) {
            if(!_partialHits[oi].hasGenomeCoords())
                return false;
        }
        return true;
    }
    
    /**
     *
     */
    index_t minWidth(index_t& offset) const {
        index_t minWidth_ = (index_t)INDEX_MAX;
        index_t minWidthLen_ = 0;
        for(size_t oi = 0; oi < _partialHits.size(); oi++) {
            const BWTHit<index_t>& hit = _partialHits[oi];
            if(hit.empty()) continue;
            // if(!hit.hasGenomeCoords()) continue;
            assert_gt(hit.size(), 0);
            if((minWidth_ > hit.size()) ||
               (minWidth_ == hit.size() && minWidthLen_ < hit.len())) {
                minWidth_ = hit.size();
                minWidthLen_ = hit.len();
                offset = (index_t)oi;
            }
        }
        return minWidth_;
    }
    
    // add policy for calculating a search score
    int64_t searchScore(index_t minK) {
        int64_t score = 0;
        const int64_t penaltyPerOffset = minK * minK;
        for(size_t i = 0; i < _partialHits.size(); i++) {
            index_t len = _partialHits[i]._len;
            score += (len * len);
        }
        
        assert_geq(_numPartialSearch, _partialHits.size());
        index_t actualPartialSearch = numActualPartialSearch();
        score -= (actualPartialSearch * penaltyPerOffset);
        score -= (1 << (actualPartialSearch << 1));
        return score;
    }
    
    BWTHit<index_t>& getPartialHit(index_t offset_) {
        assert_lt(offset_, _partialHits.size());
        return _partialHits[offset_];
    }
    
    bool adjustOffset(index_t minK) {
        assert_gt(_partialHits.size(), 0);
        const BWTHit<index_t>& hit = _partialHits.back();
        if(hit.len() >= minK + 3) {
            return false;
        }
        assert_geq(_cur, hit.len());
        index_t origCur = _cur - hit.len();
        _cur = origCur + max(hit.len(), minK + 1) - minK;
        _partialHits.pop_back();
        return true;
    }
    
    void setOffset(index_t offset) {
        assert_lt(offset, _len);
        _cur = offset;
    }
    
#ifndef NDEBUG
	/**
	 */
	bool repOk() const {
        for(size_t i = 0; i < _partialHits.size(); i++) {
            if(i == 0) {
                assert_geq(_partialHits[i]._bwoff, 0);
            }
            
            if(i + 1 < _partialHits.size()) {
                assert_leq(_partialHits[i]._bwoff + _partialHits[i]._len, _partialHits[i+1]._bwoff);
            } else {
                assert_eq(i+1, _partialHits.size());
                assert_eq(_partialHits[i]._bwoff + _partialHits[i]._len, _cur);
            }
        }
		return true;
	}
#endif
	
	bool     _fw;
	index_t  _len;
    index_t  _cur;
    bool     _done;
    index_t  _numPartialSearch;
    index_t  _numUniqueSearch;
    index_t  _cur_local;
    bool     _repeat;
    
    EList<BWTHit<index_t> >       _partialHits;
};


/**
 * this is per-thread data, which are shared by GenomeHit classes
 * the main purpose of this struct is to avoid extensive use of memory related functions
 * such as new and delete - those are really slow and lock based
 */
template <typename index_t>
struct SharedTempVars {
    SStringExpandable<char> raw_refbuf;
    SStringExpandable<char> raw_refbuf2;
    EList<int64_t> temp_scores;
    EList<int64_t> temp_scores2;
    
    // Align with alternatives
    EList<pair<index_t, int> >      ssOffs;
    EList<pair<index_t, int> >      offDiffs;
    EList<SStringExpandable<char> > raw_refbufs;
    EList<Edit>                     alt_edits;
    ELList<Edit, 128, 4>            candidate_edits;
    ELList<pair<index_t, index_t> > ht_llist;
    Haplotype<index_t>              cmp_ht;
    
    ASSERT_ONLY(SStringExpandable<uint32_t> destU32);
    
    ASSERT_ONLY(BTDnaString editstr);
    ASSERT_ONLY(BTDnaString partialseq);
    ASSERT_ONLY(BTDnaString refstr);
    ASSERT_ONLY(EList<index_t> reflens);
    ASSERT_ONLY(EList<index_t> refoffs);
    
    LinkedEList<EList<Edit> > raw_edits;
    LinkedEList<EList<pair<index_t, index_t> > > raw_ht_lists;
};

/**
 * GenomeHit represents read alignment or alignment of a part of a read
 * Two GenomeHits that represents alignments of different parts of a read
 * can be combined together.  Also, GenomeHit can be extended in both directions.
 */
template <typename index_t>
struct GenomeHit {
	GenomeHit() :
    _fw(false),
    _rdoff((index_t)INDEX_MAX),
    _len((index_t)INDEX_MAX),
    _trim5(0),
    _trim3(0),
    _tidx((index_t)INDEX_MAX),
    _toff((index_t)INDEX_MAX),
    _joinedOff((index_t)INDEX_MAX),
    _repeat(false),
    _edits(NULL),
    _ht_list(NULL),
    _score(MIN_I64),
    _localscore(MIN_I64),
    _hitcount(1),
    _edits_node(NULL),
    _ht_list_node(NULL),
    _sharedVars(NULL)
    {
    }
    
    GenomeHit(const GenomeHit& otherHit) :
    _edits(NULL),
    _ht_list(NULL),
    _hitcount(1),
    _edits_node(NULL),
    _ht_list_node(NULL),
    _sharedVars(NULL)
    {
        init(otherHit._fw,
             otherHit._rdoff,
             otherHit._len,
             otherHit._trim5,
             otherHit._trim3,
             otherHit._tidx,
             otherHit._toff,
             otherHit._joinedOff,
             *(otherHit._sharedVars),
             otherHit._repeat,
             otherHit._edits,
             otherHit._ht_list,
             otherHit._score,
             otherHit._localscore,
             otherHit._splicescore);
    }
    
    GenomeHit<index_t>& operator=(const GenomeHit<index_t>& otherHit) {
        if(this == &otherHit) return *this;
        init(otherHit._fw,
             otherHit._rdoff,
             otherHit._len,
             otherHit._trim5,
             otherHit._trim3,
             otherHit._tidx,
             otherHit._toff,
             otherHit._joinedOff,
             *(otherHit._sharedVars),
             otherHit._repeat,
             otherHit._edits,
             otherHit._ht_list,
             otherHit._score,
             otherHit._localscore,
             otherHit._splicescore);
        
        return *this;
    }
    
    ~GenomeHit() {
        if(_edits_node != NULL) {
            assert(_edits != NULL);
            assert(_sharedVars != NULL);
            _sharedVars->raw_edits.delete_node(_edits_node);
            _edits = NULL;
            _edits_node = NULL;
        }
        if(_ht_list_node != NULL) {
            assert(_ht_list != NULL);
            assert(_sharedVars != NULL);
            _sharedVars->raw_ht_lists.delete_node(_ht_list_node);
            _ht_list = NULL;
            _ht_list_node = NULL;
        }
        _sharedVars = NULL;
    }
	
	void init(
              bool                      fw,
			  index_t                   rdoff,
			  index_t                   len,
              index_t                   trim5,
              index_t                   trim3,
              index_t                   tidx,
              index_t                   toff,
              index_t                   joinedOff,
              SharedTempVars<index_t>&  sharedVars,
              bool                      repeat = false,
              EList<Edit>*              edits = NULL,
              EList<pair<index_t, index_t> >* ht_list = NULL,
              int64_t                   score = 0,
              int64_t                   localscore = 0,
              double                    splicescore = 0.0)
	{
		_fw = fw;
		_rdoff = rdoff;
		_len = len;
        _trim5 = trim5;
        _trim3 = trim3;
        _tidx = tidx;
        _toff = toff;
        _joinedOff = joinedOff;
        _repeat = repeat;
		_score = score;
        _localscore = localscore;
        _splicescore = splicescore;
        
        assert(_sharedVars == NULL || _sharedVars == &sharedVars);
        _sharedVars = &sharedVars;
        if(_edits == NULL) {
            assert(_edits_node == NULL);
            _edits_node = _sharedVars->raw_edits.new_node();
            assert(_edits_node != NULL);
            _edits = &(_edits_node->payload);
        }
        assert(_edits != NULL);
        _edits->clear();
        if(edits != NULL) *_edits = *edits;
        
        if(_ht_list == NULL) {
            assert(_ht_list_node == NULL);
            _ht_list_node = _sharedVars->raw_ht_lists.new_node();
            assert(_ht_list_node != NULL);
            _ht_list = &(_ht_list_node->payload);
        }
        assert(_ht_list != NULL);
        _ht_list->clear();
        if(ht_list != NULL) *_ht_list = *ht_list;
        
        _hitcount = 1;
	}
    
    bool inited() const {
        return _len >= 0 && _len < (index_t)INDEX_MAX;
    }
    
    /**
     * Check if it is compatible with another GenomeHit with respect to indels or introns
     */
    bool compatibleWith(
                        const GenomeHit<index_t>& otherHit,
                        index_t minIntronLen,
                        index_t maxIntronLen,
                        bool no_spliced_alignment = false) const;
    
    /**
     * Combine itself with another GenomeHit
     */
    bool combineWith(
                     const GenomeHit&           otherHit,
                     const Read&                rd,
                     const GFM<index_t>&        gfm,
                     const BitPairReference&    ref,
                     const ALTDB<index_t>&      altdb,
                     const RepeatDB<index_t>&   repeatdb,
                     SpliceSiteDB&              ssdb,
                     SwAligner&                 swa,
                     SwMetrics&                 swm,
                     const Scoring&             sc,
                     TAlScore                   minsc,
                     RandomSource&              rnd,                  // pseudo-random source
                     index_t                    minK_local,
                     index_t                    minIntronLen,
                     index_t                    maxIntronLen,
                     index_t                    minAnchorLen,         // minimum anchor length for canonical splice site
                     index_t                    minAnchorLen_noncan,  // minimum anchor length for non-canonical splice site
                     const index_t              maxAltsTried,
                     const SpliceSite*          spliceSite = NULL,    // penalty for splice site
                     bool                       no_spliced_alignment = false);
    
    /**
     * Extend the partial alignment (GenomeHit) bidirectionally
     */
    bool extend(
                const Read&                rd,
                const GFM<index_t>&        gfm,
                const BitPairReference&    ref,
                const ALTDB<index_t>&      altdb,
                const RepeatDB<index_t>&   repeatdb,
                SpliceSiteDB&              ssdb,
                SwAligner&                 swa,
                SwMetrics&                 swm,
                PerReadMetrics&            prm,
                const Scoring&             sc,
                TAlScore                   minsc,
                RandomSource&              rnd,           // pseudo-random source
                index_t                    minK_local,
                const TranscriptomePolicy& tpol,
                const GraphPolicy&         gpol,
                index_t&                   leftext,
                index_t&                   rightext,
                index_t                    mm = 0);
    
    /**
     * Adjust alignment with respect to SNPs, usually updating Edits
     *
     */
    static bool adjustWithALT(
                              index_t                     rdoff,
                              index_t                     len,
                              const Coord&                coord,
                              SharedTempVars<index_t>&    sharedVars,
                              EList<GenomeHit<index_t> >& genomeHits,
                              const Read&                 rd,
                              const GFM<index_t>&         gfm,
                              const ALTDB<index_t>&       altdb,
                              const BitPairReference&     ref,
                              const GraphPolicy&          gpol);
    
    /**
     * Adjust alignment with respect to SNPs, usually updating Edits
     *
     */
    bool adjustWithALT(
                       const Read&             rd,
                       const GFM<index_t>&     gfm,
                       const ALTDB<index_t>&   altdb,
                       const BitPairReference& ref,
                       const GraphPolicy&      gpol);
   
    /*
     *
     */
    static void findSSOffs(
                           const GFM<index_t>&         gfm,
                           const ALTDB<index_t>&       altdb,
                           index_t                     start,
                           index_t                     end,
                           EList<pair<index_t, int> >& ssOffs);
    
    /*
     * Find offset differences due to deletions
     */
    static index_t findOffDiffs(
                                const GFM<index_t>&         gfm,
                                const ALTDB<index_t>&       altdb,
                                index_t                     start,
                                index_t                     end,
                                EList<pair<index_t, int> >& offDiffs);
    
    /*
     *
     */
    static index_t alignWithALTs(
                                 const EList<ALT<index_t> >&       alts,
                                 const EList<Haplotype<index_t> >& haplotypes,
                                 const EList<index_t>&             haplotype_maxrights,
                                 index_t                           joinedOff,
                                 const BTDnaString&                rdseq,
                                 index_t                           base_rdoff,
                                 index_t                           rdoff,
                                 index_t                           rdlen,
                                 const BitPairReference&           ref,
                                 SharedTempVars<index_t>&          sharedVar,
                                 index_t                           tidx,
                                 int                               rfoff,
                                 index_t                           rflen,
                                 bool                              left,
                                 const GraphPolicy&                gpol,
                                 EList<Edit>&                      edits,
                                 ELList<pair<index_t, index_t> >&  ht_llist,
                                 EList<pair<index_t, index_t> >&   ht_list,
                                 Haplotype<index_t>&               cmp_ht,
                                 ELList<Edit, 128, 4>*             candidate_edits = NULL,
                                 index_t                           mm = 0,
                                 index_t*                          numNs = NULL)
    {
        int best_rdoff = (int)rdoff;
        if(numNs != NULL) *numNs = 0;
        index_t numALTsTried = 0;
        EList<Edit>& alt_edits = sharedVar.alt_edits;
        alt_edits = edits;
        index_t nedits = (index_t)edits.size();
        if(candidate_edits != NULL) candidate_edits->clear();
        ht_llist.clear();
        // ht_llist.expand();
        // ht_llist[0] = ht_list;
        alignWithALTs_recur(
                            alts,
                            haplotypes,
                            haplotype_maxrights,
                            joinedOff,
                            rdseq,
                            rdoff - base_rdoff,
                            rdoff,
                            rdlen,
                            ref,
                            sharedVar.raw_refbufs,
                            ASSERT_ONLY(sharedVar.destU32,)
                            alt_edits,
                            best_rdoff,
                            NULL, /* rfseq */
                            tidx,
                            rfoff,
                            rflen,
                            left,
                            edits,
                            mm,
                            ht_llist,
                            cmp_ht,
                            candidate_edits,
                            0, /* tmp_numNs */
                            numNs,
                            0,    /* dep */
                            gpol,
                            numALTsTried);
        index_t extlen = 0;
        if(left) {
            assert_geq(best_rdoff, -1);
            assert_leq(best_rdoff, (int)rdoff);
            extlen = rdoff - best_rdoff;
        } else {
            assert_leq(best_rdoff, (int)(rdoff + rdlen));
            assert_geq(best_rdoff, (int)rdoff);
            extlen = best_rdoff - rdoff;
        }
        if(extlen > 0 && edits.size() > 0) {
            const Edit& f = edits.front();
            if(f.pos + extlen == base_rdoff + 1) {
                if(f.type == EDIT_TYPE_READ_GAP ||
                   f.type == EDIT_TYPE_REF_GAP ||
                   f.type == EDIT_TYPE_SPL) {
                    extlen = 0;
                }
                if(f.type == EDIT_TYPE_MM && f.chr == 'N') {
                    extlen = 0;
                }
            }
            const Edit& b = edits.back();
            if(extlen > 0 && b.pos == rdoff - base_rdoff + extlen - 1) {
                if(b.type == EDIT_TYPE_READ_GAP ||
                   b.type == EDIT_TYPE_REF_GAP) {
                    extlen = 0;
                }
            }
            if(extlen == 0 && edits.size() > nedits) {
                if(left) {
                    edits.erase(0, edits.size() - nedits);
                } else {
                    edits.resize(nedits);
                }
            }
        }
        return extlen;
    }
    
    /*
     *
     */
    static index_t alignWithALTs_recur(
                                       const EList<ALT<index_t> >&       alts,
                                       const EList<Haplotype<index_t> >& haplotypes,
                                       const EList<index_t>&             haplotype_maxrights,
                                       index_t                           joinedOff,
                                       const BTDnaString&                rdseq,
                                       index_t                           rdoff_add,
                                       index_t                           rdoff,
                                       index_t                           rdlen,
                                       const BitPairReference&           ref,
                                       EList<SStringExpandable<char> >&  raw_refbufs,
                                       ASSERT_ONLY(SStringExpandable<uint32_t> destU32,)
                                       EList<Edit>&                      tmp_edits,
                                       int&                              best_rdoff,
                                       const char*                       rfseq,
                                       index_t                           tidx,
                                       int                               rfoff,
                                       index_t                           rflen,
                                       bool                              left,
                                       EList<Edit>&                      edits,
                                       index_t                           mm,
                                       ELList<pair<index_t, index_t> >&  ht_llist,
                                       Haplotype<index_t>&               cmp_ht,
                                       ELList<Edit, 128, 4>*             candidate_edits,
                                       index_t                           tmp_numNs,
                                       index_t*                          numNs,
                                       index_t                           dep,
                                       const GraphPolicy&                gpol,
                                       index_t&                          numALTsTried,
                                       ALT_TYPE                          prev_alt_type = ALT_NONE);
    
    /**
     * For alignment involving indel, move the indels
     * to the left most possible position
     */
    void leftAlign(const Read& rd);
    
    index_t rdoff() const { return _rdoff; }
    index_t len()   const { return _len; }
    index_t trim5() const { return _trim5; }
    index_t trim3() const { return _trim3; }
    
    void trim5(index_t                 trim5,
               const Read&             rd,
               SpliceSiteDB&           ssdb,
               const Scoring&          sc,
               index_t                 minK_local,
               index_t                 minIntronLen,
               index_t                 maxIntronLen,
               index_t                 minAnchorLen,
               index_t                 minAnchorLen_noncan,
               const BitPairReference& ref)
    {
        assert_eq(_rdoff, trim5);
        assert_eq(_trim5, 0);
        _trim5 = trim5;
        calculateScore(rd,
                       ssdb,
                       sc,
                       minK_local,
                       minIntronLen,
                       maxIntronLen,
                       minAnchorLen,
                       minAnchorLen_noncan,
                       ref);
    }
    void trim3(index_t                 trim3,
               const Read&             rd,
               SpliceSiteDB&           ssdb,
               const Scoring&          sc,
               index_t                 minK_local,
               index_t                 minIntronLen,
               index_t                 maxIntronLen,
               index_t                 minAnchorLen,
               index_t                 minAnchorLen_noncan,
               const BitPairReference& ref)
    {
        _trim3 = trim3;
        calculateScore(rd,
                       ssdb,
                       sc,
                       minK_local,
                       minIntronLen,
                       maxIntronLen,
                       minAnchorLen,
                       minAnchorLen_noncan,
                       ref);
    }
    
    // for repeat alignments
    //  reverse fw
    void reverse(const Read& rd)
    {
        _fw = !_fw;
        index_t end = _trim5 + _rdoff + _len;
        assert_leq(end, rd.length());
        _rdoff = rd.length() - end;
        index_t tmp_trim = _trim5;
        _trim5 = _trim3;
        _trim3 = tmp_trim;
        Edit::invertPoss(*_edits, rd.length());
        // complements
        Edit::complement(*_edits);
    }
    
    index_t ref()    const { return _tidx; }
    index_t refoff() const { return _toff; }
    index_t fw()     const { return _fw; }
    
    bool    repeat() const { return _repeat; }
    void    repeat(bool repeat) { _repeat = repeat;}
    
    index_t hitcount() const { return _hitcount; }
    
    /**
     * Leftmost coordinate
     */
    Coord coord() const {
        return Coord(_tidx, _toff, _fw);
    }
    
    int64_t score()       const { return _score; }
    int64_t localscore()  const { return _localscore; }
    double  splicescore() const { return _splicescore; }
    
    const EList<Edit>& edits() const { return *_edits; }
    
    /**
     * Retrieve the partial alignment from the left until indel or intron
     */
    void getLeft(index_t&       rdoff,
                 index_t&       len,
                 index_t&       toff,
                 int64_t*       score = NULL,
                 const Read*    rd = NULL,
                 const Scoring* sc = NULL) const
    {
        assert(inited());
        toff = _toff, rdoff = _rdoff, len = _len;
        const BTString* qual = NULL;
        if(score != NULL) {
            assert(rd != NULL);
            assert(sc != NULL);
            *score = 0;
            qual = &(_fw ? rd->qual : rd->qualRev);
        }
        for(index_t i = 0; i < _edits->size(); i++) {
            const Edit& edit = (*_edits)[i];
            if(edit.type == EDIT_TYPE_SPL ||
               edit.type == EDIT_TYPE_READ_GAP ||
               edit.type == EDIT_TYPE_REF_GAP ||
               (edit.type == EDIT_TYPE_MM && edit.snpID != (index_t)INDEX_MAX)) {
                len = edit.pos;
                break;
            }
            if(score != NULL) {
                if(edit.type == EDIT_TYPE_MM) {
                    assert(qual != NULL);
                    if(edit.snpID == (index_t)INDEX_MAX) {
                        *score += sc->score(
                                            dna2col[edit.qchr] - '0',
                                            asc2dnamask[edit.chr],
                                            (*qual)[this->_rdoff + edit.pos] - 33);
                    }
                }
            }
        }
        assert_geq(len, 0);
    }
    
    /**
     * Retrieve the partial alignment from the right until indel or intron
     */
    void getRight(index_t&       rdoff,
                  index_t&       len,
                  index_t&       toff,
                  int64_t*       score = NULL,
                  const Read*    rd = NULL,
                  const Scoring* sc = NULL) const
    {
        assert(inited());
        toff = _toff, rdoff = _rdoff, len = _len;
        const BTString* qual = NULL;
        if(score != NULL) {
            assert(rd != NULL);
            assert(sc != NULL);
            *score = 0;
            qual = &(_fw ? rd->qual  : rd->qualRev);
        }
        if(_edits->size() == 0) return;
        for(int i = (int)_edits->size() - 1; i >= 0; i--) {
            const Edit& edit = (*_edits)[i];
            if(edit.type == EDIT_TYPE_SPL ||
               edit.type == EDIT_TYPE_READ_GAP ||
               edit.type == EDIT_TYPE_REF_GAP ||
               (edit.type == EDIT_TYPE_MM && edit.snpID != (index_t)INDEX_MAX)) {
                rdoff = _rdoff + edit.pos;
                assert_lt(edit.pos, _len);
                len = _len - edit.pos;
                if(edit.type == EDIT_TYPE_REF_GAP) {
                    assert_lt(edit.pos + 1, _len);
                    assert_gt(len, 1);
                    rdoff++;
                    len--;
                } else if(edit.type == EDIT_TYPE_MM) {
                    assert_leq(edit.pos + 1, _len);
                    assert_geq(len, 1);
                    rdoff++;
                    len--;
                }
                toff = getRightOff() - len;
                break;
            }
            if(score != NULL) {
                if(edit.type == EDIT_TYPE_MM) {
                    assert(qual != NULL);
                    if(edit.snpID == (index_t)INDEX_MAX) {
                        *score += sc->score(
                                            dna2col[edit.qchr] - '0',
                                            asc2dnamask[edit.chr],
                                            (*qual)[this->_rdoff + edit.pos] - 33);
                    }
                }
            }
        }
        assert_geq(len, 0);
    }
    
    /**
     * Retrieve the genomic offset of the right end
     */
    index_t getRightOff() const {
        assert(inited());
        index_t toff = _toff + _len;
        for(index_t i = 0; i < _edits->size(); i++) {
            const Edit& ed = (*_edits)[i];
            if(ed.type == EDIT_TYPE_SPL) {
                toff += ed.splLen;
            } else if(ed.type == EDIT_TYPE_READ_GAP) {
                toff++;
            } else if(ed.type == EDIT_TYPE_REF_GAP) {
                assert_gt(toff, 0);
                toff--;
            }
        }
        return toff;
    }
    
    /**
     * Retrieve left anchor length and number of edits in the anchor
     */
    void getLeftAnchor(index_t& leftanchor,
                       index_t& nedits) const
    {
        assert(inited());
        leftanchor = _len;
        nedits = 0;
        for(index_t i = 0; i < _edits->size(); i++) {
            const Edit& edit = (*_edits)[i];
            if(edit.type == EDIT_TYPE_SPL) {
                leftanchor = edit.pos;
                break;
            } else if(edit.type == EDIT_TYPE_MM ||
                      edit.type == EDIT_TYPE_READ_GAP ||
                      edit.type == EDIT_TYPE_REF_GAP) {
                nedits++;
            }
        }
    }
    
    /**
     * Retrieve right anchor length and number of edits in the anchor
     */
    void getRightAnchor(index_t& rightanchor,
                        index_t& nedits) const
    {
        rightanchor = _len;
        nedits = 0;
        if(_edits->size() == 0) return;
        for(int i = (int)_edits->size() - 1; i >= 0; i--) {
            const Edit& edit = (*_edits)[i];
            if(edit.type == EDIT_TYPE_SPL) {
                rightanchor = _len - edit.pos - 1;
                break;
            } else if(edit.type == EDIT_TYPE_MM ||
                      edit.type == EDIT_TYPE_READ_GAP ||
                      edit.type == EDIT_TYPE_REF_GAP) {
                nedits++;
            }
        }
    }
    
    
    /**
     * Is it spliced alignment?
     * Return: first is spliced-alignment, second is spliced-alignment to known transcripts
     */
    pair<bool, bool> spliced() const {
        pair<bool, bool> result(false, true);
        for(index_t i = 0; i < _edits->size(); i++) {
            const Edit& e = (*_edits)[i];
            if(e.type == EDIT_TYPE_SPL) {
                result.first = true;
                result.second &= e.knownSpl;
            }
        }
        result.second &= result.first;
        return result;
    }
    
    /**
     *
     */
    bool spliced_consistently() const {
        int splDir = SPL_UNKNOWN;
        for(index_t i = 0; i < _edits->size(); i++) {
            const Edit& edit = (*_edits)[i];
            if(edit.type != EDIT_TYPE_SPL) continue;
            if(splDir != SPL_UNKNOWN) {
                if(edit.splDir != SPL_UNKNOWN) {
                    if(splDir == SPL_FW || splDir == SPL_SEMI_FW) {
                        if(edit.splDir != SPL_FW && edit.splDir != SPL_SEMI_FW)
                            return false;
                    }
                    if(splDir == SPL_RC || splDir == SPL_SEMI_RC) {
                        if(edit.splDir != SPL_RC && edit.splDir != SPL_SEMI_RC)
                            return false;
                    }
                }
            } else {
                splDir = edit.splDir;
            }
        }
        return true;
    }
    
    /**
     * return one of EDIT_SPL_FW, EDIT_SPL_RC, EDIT_SPL_UNKNOWN
     */
    int splicing_dir() const {
        int splDir = SPL_UNKNOWN;
        for(index_t i = 0; i < _edits->size(); i++) {
            const Edit& edit = (*_edits)[i];
            if(edit.type != EDIT_TYPE_SPL) continue;
            if(splDir != SPL_UNKNOWN) {
                if(edit.splDir != SPL_UNKNOWN) {
                    if(splDir == SPL_FW || splDir == SPL_SEMI_FW) {
                        if(edit.splDir != SPL_FW && edit.splDir != SPL_SEMI_FW)
                            return SPL_UNKNOWN;
                    }
                    if(splDir == SPL_RC || splDir == SPL_SEMI_RC) {
                        if(edit.splDir != SPL_RC && edit.splDir != SPL_SEMI_RC)
                            return SPL_UNKNOWN;
                    }
                }
            } else {
                splDir = edit.splDir;
            }
        }
        if(splDir == SPL_FW || splDir == SPL_SEMI_FW)
            return SPL_FW;
        else if(splDir == SPL_RC || splDir == SPL_SEMI_RC)
            return SPL_RC;
        else
            return SPL_UNKNOWN;
    }
    
    bool operator== (const GenomeHit<index_t>& other) const {
        if(_fw != other._fw ||
           _rdoff != other._rdoff ||
           _len != other._len ||
           _tidx != other._tidx ||
           _toff != other._toff ||
           _trim5 != other._trim5 ||
           _trim3 != other._trim3) {
            return false;
        }
        if(_edits->size() != other._edits->size()) return false;
        for(index_t i = 0; i < _edits->size(); i++) {
            const Edit& e = (*_edits)[i];
            const Edit& oe = (*other._edits)[i];
            if(e.isReadGap()) {
                if(!oe.isReadGap()) return false;
            } else if(e.isRefGap()) {
                if(!oe.isRefGap()) return false;
            } else {
                if(!(e == oe)) {
                    return false;
                }
            }
        }
        // daehwan - this may not be true when some splice sites are provided from outside
        // assert_eq(_score, other._score);
        return true;
    }
    
    bool contains(const GenomeHit<index_t>& other) const {
        return (*this) == other;
    }

	/**
	 * Return number of mismatches in the alignment.
	 */
	int mms() const {
#if 0
		if     (_e2.inited()) return 2;
		else if(_e1.inited()) return 1;
		else                  return 0;
#endif
        return 0;
	}
	
	/**
	 * Return the number of Ns involved in the alignment.
	 */
	int ns() const {
#if 0
		int ns = 0;
		if(_e1.inited() && _e1.hasN()) {
			ns++;
			if(_e2.inited() && _e2.hasN()) {
				ns++;
			}
		}
		return ns;
#endif
        return 0;
	}
    
    int ngaps() const {
        return 0;
    }
	
#ifndef NDEBUG
	/**
	 * Check that hit is sane w/r/t read.
	 */
	bool repOk(const Read& rd, const BitPairReference& ref);
#endif
    
    void replace_edits_with_alts(const Read&      rd,
                                 const EList<ALT<index_t> >& alts,
                                 SpliceSiteDB&    ssdb,
                                 const Scoring&   sc,
                                 index_t          minK_local,
                                 index_t          minIntronLen,
                                 index_t          maxIntronLen,
                                 index_t          minAnchorLen,
                                 index_t          minAnchorLen_noncan,
                                 const            BitPairReference& ref) {
        assert(inited());
        if(alts.size() <= 0)
            return;
        if(_edits->size() <= 0)
            return;
        
        index_t joinedOff = _joinedOff;
        int offset = 0;
        size_t i = 0, next_i;
        while(i < _edits->size()) {
            next_i = i + 1;
            Edit& ed = (*_edits)[i];
            if(ed.type == EDIT_TYPE_SPL) {
                assert(false);
            } else if(ed.type == EDIT_TYPE_READ_GAP || ed.type == EDIT_TYPE_REF_GAP) {
                for(; next_i < _edits->size(); next_i++) {
                    Edit& next_ed = (*_edits)[next_i];
                    if(ed.type != next_ed.type) break;
                }
            }

            if(ed.snpID == (index_t)INDEX_MAX) {
                ALT<index_t> cmp_alt;
                cmp_alt.pos = joinedOff + ed.pos + offset;
                index_t alt_i = (index_t)alts.bsearchLoBound(cmp_alt);
                for(; alt_i < alts.size(); alt_i++) {
                    const ALT<index_t>& alt = alts[alt_i];
                    if(alt.left > cmp_alt.pos) break;
                    if(ed.type == EDIT_TYPE_MM) {
                        if(alt.type != ALT_SNP_SGL) continue;
                        if("ACGT"[alt.seq] == ed.qchr) {
                            ed.snpID = alt_i;
                            break;
                        }
                    } else {
                        size_t gap = next_i - i;
                        if(ed.type == EDIT_TYPE_READ_GAP) {
                            if(alt.type != ALT_SNP_DEL) continue;
                            if(alt.len == gap) {
                                for(size_t ii = i; ii < next_i; ii++) {
                                    Edit& ii_ed = (*_edits)[ii];
                                    ii_ed.snpID = alt_i;
                                }
                                break;
                            }
                        } else {
                            assert_eq(ed.type, EDIT_TYPE_REF_GAP);
                            if(alt.type != ALT_SNP_INS) continue;
                            if(alt.len == gap) {
                                uint64_t seq = 0;
                                for(size_t ii = i; ii < next_i; ii++) {
                                    Edit& ii_ed = (*_edits)[ii];
                                    seq = (seq << 2) | asc2dna[ii_ed.qchr];
                                }
                                if(alt.seq == seq) {
                                    for(size_t ii = i; ii < next_i; ii++) {
                                        Edit& ii_ed = (*_edits)[ii];
                                        ii_ed.snpID = alt_i;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            
            if(ed.type == EDIT_TYPE_SPL) {
                offset += ed.splLen;
            } else if(ed.type == EDIT_TYPE_READ_GAP || ed.type == EDIT_TYPE_REF_GAP) {
                size_t gap = next_i - i;
                if(ed.type == EDIT_TYPE_READ_GAP) {
                    assert_gt(joinedOff, gap);
                    offset += gap;
                } else {
                    offset -= gap;
                }
            }
            
            i = next_i;
        }
        
        calculateScore(rd,
                       ssdb,
                       sc,
                       minK_local,
                       minIntronLen,
                       maxIntronLen,
                       minAnchorLen,
                       minAnchorLen_noncan,
                       ref);
    }
    
private:
    /**
	 * Calculate alignment score
	 */
    int64_t calculateScore(
                           const Read&      rd,
                           SpliceSiteDB&    ssdb,
                           const Scoring&   sc,
                           index_t          minK_local,
                           index_t          minIntronLen,
                           index_t          maxIntronLen,
                           index_t          minAnchorLen,
                           index_t          minAnchorLen_noncan,
                           const            BitPairReference& ref);
    
public:
	bool            _fw;
	index_t         _rdoff;
	index_t         _len;
    index_t         _trim5;
    index_t         _trim3;
    
    index_t         _tidx;
    index_t         _toff;
    index_t         _joinedOff;
    bool            _repeat;
	EList<Edit>*    _edits;
    EList<pair<index_t, index_t> >* _ht_list;
    int64_t         _score;
    int64_t         _localscore;
    double          _splicescore;
    
    index_t         _hitcount;  // for selection purposes
    
    LinkedEListNode<EList<Edit> >*  _edits_node;
    LinkedEListNode<EList<pair<index_t, index_t> > >*  _ht_list_node;
    SharedTempVars<index_t>* _sharedVars;
};

/**
 * Check if it is compatible with another GenomeHit with respect to indels or introns
 */
template <typename index_t>
bool GenomeHit<index_t>::compatibleWith(
                                        const GenomeHit<index_t>& otherHit,
                                        index_t minIntronLen,
                                        index_t maxIntronLen,
                                        bool no_spliced_alignment) const
{
    if(this == &otherHit) return false;
    // check if they are on the same strand and on the same contig
    if(_fw != otherHit._fw || _tidx != otherHit._tidx) return false;
    // make sure itself is closer to the left end of read than otherHit
    if(_rdoff > otherHit._rdoff) return false;
    // do not consider a case itself (read portion) includes otherHit
    if(_rdoff + _len > otherHit._rdoff + otherHit._len) return false;
    // make sure itself comes before otherHit wrt. genomic positions
    if(_toff > otherHit._toff) return false;
    
    index_t this_rdoff, this_len, this_toff;
    this->getRight(this_rdoff, this_len, this_toff);
    assert_geq(this_len, 0);
    index_t other_rdoff, other_len, other_toff;
    otherHit.getLeft(other_rdoff, other_len, other_toff);
    assert_geq(other_len, 0);
    
    if(this_rdoff > other_rdoff) return false;
    if(this_rdoff + this_len > other_rdoff + other_len) return false;
    if(this_toff > other_toff) return false;
    
    index_t refdif = other_toff - this_toff;
    index_t rddif = other_rdoff - this_rdoff;
    
    // check if there is a deletion, an insertion, or a potential intron
    // between the two partial alignments
    if(!no_spliced_alignment) {
        if(refdif > rddif + maxIntronLen) {
            return false;
        }
    }
    return true;
}

/**
 * Combine itself with another GenomeHit
 * while allowing mismatches, an insertion, a deletion, or an intron
 */
template <typename index_t>
bool GenomeHit<index_t>::combineWith(
                                     const GenomeHit&           otherHit,
                                     const Read&                rd,
                                     const GFM<index_t>&        gfm,
                                     const BitPairReference&    ref,
                                     const ALTDB<index_t>&      altdb,
                                     const RepeatDB<index_t>&   repeatdb,
                                     SpliceSiteDB&              ssdb,
                                     SwAligner&                 swa,
                                     SwMetrics&                 swm,
                                     const Scoring&             sc,
                                     TAlScore                   minsc,
                                     RandomSource&              rnd,                    // pseudo-random source
                                     index_t                    minK_local,
                                     index_t                    minIntronLen,
                                     index_t                    maxIntronLen,
                                     index_t                    minAnchorLen,           // minimum anchor length for canonical splice site
                                     index_t                    minAnchorLen_noncan,    // minimum anchor length for non-canonical splice site
                                     const index_t              maxAltsTried,
                                     const SpliceSite*          spliceSite,             // penalty for splice site
                                     bool                       no_spliced_alignment)
{
    if(this == &otherHit) return false;
    assert(compatibleWith(otherHit, minIntronLen, maxIntronLen, no_spliced_alignment));
    assert_eq(this->_tidx, otherHit._tidx);
    assert_lt(this->_tidx, ref.numRefs());
    
    // get the partial part of the alignment from the right
    // until an indel or splice sites
    index_t this_rdoff, this_len, this_toff;
    int64_t this_score;
    this->getRight(this_rdoff, this_len, this_toff, &this_score, &rd, &sc);
    assert_geq(this_len, 0);
    assert_leq(this_score, 0);
    assert_geq(this_score, this->_score);
    
    // get the partial part of the other alignment from the left
    // until an indel or splice sites
    index_t other_rdoff, other_len, other_toff;
    int64_t other_score;
    otherHit.getLeft(other_rdoff, other_len, other_toff, &other_score, &rd, &sc);
    assert_geq(other_len, 0);
    assert_leq(other_score, 0);
    assert_geq(other_score, otherHit._score);
    
    assert_leq(this_rdoff, other_rdoff);
    if(this_len != 0 &&
       other_len != 0 &&
       this_rdoff + this_len > other_rdoff + other_len) return false;
    assert_leq(this_rdoff + this_len, other_rdoff + other_len);
    index_t len = other_rdoff - this_rdoff + other_len;
    const index_t reflen = ref.approxLen(_tidx);
    if(this_toff + len > reflen) return false;
    assert_leq(this_toff + len, reflen);
    
    // check if an indel or an intron is necessary
    index_t refdif = other_toff - this_toff;
    index_t rddif = other_rdoff - this_rdoff;
    bool spliced = false, ins = false, del = false;
    if(refdif != rddif) {
        if(refdif > rddif) {
            if(!no_spliced_alignment && refdif - rddif >= minIntronLen) {
                assert_leq(refdif - rddif, maxIntronLen);
                spliced = true;
            } else {
                del = true;
            }
        } else {
            ins = true;
        }
    }
#ifndef NDEBUG
    if(ins) {
        assert(!spliced && !del);
    } else {
        if(spliced) assert(!del);
        else        assert(!spliced);
    }
#endif
    
    if(no_spliced_alignment) {
        if(spliced) return false;
    }
    
    // if the combination of the two alignments does not involve an indel or an intron,
    // then simply combine them and return
    if(!spliced && !ins && !del && this_rdoff + this_len == other_rdoff) {
        index_t addoff = otherHit._rdoff - this->_rdoff;
        for(index_t i = 0; i < otherHit._edits->size(); i++) {
            _edits->push_back((*otherHit._edits)[i]);
            _edits->back().pos += addoff;
        }
        _len += otherHit._len;
        calculateScore(
                       rd,
                       ssdb,
                       sc,
                       minK_local,
                       minIntronLen,
                       maxIntronLen,
                       minAnchorLen,
                       minAnchorLen_noncan,
                       ref);
        assert(repOk(rd, ref));
        return true;
    }
    
    // calculate the maximum gap lengths based on the current score and the mimumimu alignment score to be reported
    const BTDnaString& seq = this->_fw ? rd.patFw : rd.patRc;
    const BTString& qual = this->_fw ? rd.qual : rd.qualRev;
    index_t rdlen = (index_t)seq.length();
    int64_t remainsc = minsc - (_score - this_score) - (otherHit._score - other_score);
    if(remainsc > 0) remainsc = 0;
    int read_gaps = 0, ref_gaps = 0;
    if(!spliced) {
        read_gaps = sc.maxReadGaps(remainsc + sc.canSpl(), rdlen);
        ref_gaps = sc.maxRefGaps(remainsc + sc.canSpl(), rdlen);
    }
    if(ins) {
        if(refdif + ref_gaps < rddif) {
            return false;
        }
    } else if(del) {
        if(rddif + read_gaps < refdif) {
            return false;
        }
    }
    int this_ref_ext = read_gaps;
    if(spliced) this_ref_ext += (int)intronic_len;
    if(this_toff + len > reflen) return false;
    if(this_toff + len + this_ref_ext > reflen) this_ref_ext = reflen - (this_toff + len);
    assert(_sharedVars != NULL);
    SStringExpandable<char>& raw_refbuf = _sharedVars->raw_refbuf;
    EList<int64_t>& temp_scores = _sharedVars->temp_scores;
    EList<int64_t>& temp_scores2 = _sharedVars->temp_scores2;
    ASSERT_ONLY(SStringExpandable<uint32_t>& destU32 = _sharedVars->destU32);
    raw_refbuf.resize(len + this_ref_ext + 16);
    int off = ref.getStretch(
                             reinterpret_cast<uint32_t*>(raw_refbuf.wbuf()),
                             (size_t)this->_tidx,
                             (size_t)this_toff,
                             len + this_ref_ext
                             ASSERT_ONLY(, destU32));
    assert_lt(off, 16);
    char *refbuf = raw_refbuf.wbuf() + off, *refbuf2 = NULL;
    
    // discover a splice site, an insertion, or a deletion
    index_t maxscorei = (index_t)INDEX_MAX;
    int64_t maxscore = MIN_I64;
    uint32_t maxspldir = SPL_UNKNOWN;
    float maxsplscore = 0.0f;
    // allow an indel near a splice site
    index_t splice_gap_maxscorei = (index_t)INDEX_MAX;
    int64_t donor_seq = 0, acceptor_seq = 0;
    int splice_gap_off = 0;
    if(spliced || ins || del) {
        int other_ref_ext = min<int>(read_gaps + (int)intronic_len, other_toff + other_len - len);
        SStringExpandable<char>& raw_refbuf2 = _sharedVars->raw_refbuf2;
        raw_refbuf2.resize(len + other_ref_ext + 16);
        int off2 = ref.getStretch(
                                  reinterpret_cast<uint32_t*>(raw_refbuf2.wbuf()),
                                  (size_t)otherHit._tidx,
                                  (size_t)(other_toff + other_len - len - other_ref_ext),
                                  len + other_ref_ext
                                  ASSERT_ONLY(, destU32));
        refbuf2 = raw_refbuf2.wbuf() + off2 + other_ref_ext;
        temp_scores.resize(len);
        temp_scores2.resize(len);
        if(spliced) {
            static const char GT   = 0x23, AG   = 0x02;
            static const char GTrc = 0x01, AGrc = 0x13;
            static const char GC   = 0x21, GCrc = 0x21;
            static const char AT   = 0x03, AC   = 0x01;
            static const char ATrc = 0x03, ACrc = 0x20;
            int i;
            for(i = 0; i < (int)len; i++) {
                int rdc = seq[this_rdoff + i], rfc = refbuf[i];
                if(i > 0) {
                    temp_scores[i] = temp_scores[i-1];
                } else {
                    temp_scores[i] = 0;
                }
                if(rdc != rfc) {
                    temp_scores[i] += sc.score(rdc, 1 << rfc, qual[this_rdoff + i] - 33);
                }
                if(temp_scores[i] < remainsc) {
                    break;
                }
            }
            int i_limit = min<int>(i, len);
            int i2;
            for(i2 = len - 1; i2 >= 0; i2--) {
                int rdc = seq[this_rdoff + i2], rfc = refbuf2[i2];
                if((index_t)(i2 + 1) < len) {
                    temp_scores2[i2] = temp_scores2[i2+1];
                } else {
                    temp_scores2[i2] = 0;
                }
                if(rdc != rfc) {
                    temp_scores2[i2] += sc.score(rdc, 1 << rfc, qual[this_rdoff + i2] - 33);
                }
                if(temp_scores2[i2] < remainsc) {
                    break;
                }
            }
            int i2_limit = max<int>(i2, 0);
            if(spliceSite != NULL){
                assert_leq(this_toff, (int)spliceSite->left());
                if(i2_limit <= (int)(spliceSite->left() - this_toff)) {
                    i2_limit = (int)(spliceSite->left() - this_toff);
                    i_limit = i2_limit + 1;
                } else {
                    i_limit = i2_limit;
                }
            }
            for(i = i2_limit, i2 = i2_limit + 1;
                i < i_limit && i2 < (int)len;
                i++, i2++) {
                int64_t tempscore = temp_scores[i] + temp_scores2[i2];
                char donor = 0xff, acceptor = 0xff;
                if((index_t)(i + 2) < len + this_ref_ext) {
                    donor = refbuf[i + 1];
                    donor = (donor << 4) | refbuf[i + 2];
                }
                if(i2 - 2 >= -other_ref_ext) {
                    acceptor = refbuf2[i2 - 2];
                    acceptor = (acceptor << 4) | refbuf2[i2 - 1];
                }
                bool canonical = false, semi_canonical = false;
                uint32_t spldir = SPL_UNKNOWN;
                if((donor == GT && acceptor == AG) /* || (donor == AT && acceptor == AC) */) {
                    spldir = SPL_FW;
                    canonical = true;
                } else if((donor == AGrc && acceptor == GTrc) /* || (donor == ACrc && acceptor == ATrc) */) {
                    spldir = SPL_RC;
                    canonical = true;
                } else if((donor == GC && acceptor == AG) || (donor == AT && acceptor == AC)) {
                    spldir = SPL_SEMI_FW;
                    semi_canonical = true;
                } else if((donor == AGrc && acceptor == GCrc) || (donor == ACrc && acceptor == ATrc)) {
                    spldir = SPL_SEMI_RC;
                    semi_canonical = true;
                }
                tempscore -= (canonical ? sc.canSpl() : sc.noncanSpl());
                int64_t temp_donor_seq = 0, temp_acceptor_seq = 0;
                float splscore = 0.0f;
                if(canonical) {
                    // in case of canonical splice site, extract donor side sequence and acceptor side sequence
                    //    to calculate a score of the splicing event.
                    if(spldir == SPL_FW) {
                        if(i + 1 >= (int)donor_exonic_len &&
                           (int)(len + this_ref_ext) > i + (int)donor_intronic_len &&
                           i2 + (int)other_ref_ext >= (int)acceptor_intronic_len &&
                           (int)len > i2 + (int)acceptor_exonic_len - 1) {
                            int from = i + 1 - (int)donor_exonic_len;
                            int to = i + (int)donor_intronic_len;
                            for(int j = from; j <= to; j++) {
                                assert_geq(j, 0);
                                assert_lt(j, (int)(len + this_ref_ext));
                                int base = refbuf[j];
                                if(base > 3) base = 0;
                                temp_donor_seq = temp_donor_seq << 2 | base;
                            }
                            from = i2 - acceptor_intronic_len;
                            to = i2 + acceptor_exonic_len - 1;
                            for(int j = from; j <= to; j++) {
                                assert_geq(j, -(int)other_ref_ext);
                                assert_lt(j, (int)len);
                                int base = refbuf2[j];
                                if(base > 3) base = 0;
                                temp_acceptor_seq = temp_acceptor_seq << 2 | base;
                            }
                        }
                    } else if(spldir == SPL_RC) {
                        if(i + 1 >= (int)acceptor_exonic_len &&
                           (int)(len + this_ref_ext) > i + (int)acceptor_intronic_len &&
                           i2 + (int)other_ref_ext >= (int)donor_intronic_len &&
                           (int)len > i2 + (int)donor_exonic_len - 1) {
                            int from = i + 1 - (int)acceptor_exonic_len;
                            int to = i + (int)acceptor_intronic_len;
                            for(int j = to; j >= from; j--) {
                                assert_geq(j, 0);
                                assert_lt(j, (int)(len + this_ref_ext));
                                int base = refbuf[j];
                                if(base > 3) base = 0;
                                temp_acceptor_seq = temp_acceptor_seq << 2 | (base ^ 0x3);
                            }
                            from = i2 - donor_intronic_len;
                            to = i2 + donor_exonic_len - 1;
                            for(int j = to; j >= from; j--) {
                                assert_geq(j, -(int)other_ref_ext);
                                assert_lt(j, (int)len);
                                int base = refbuf2[j];
                                if(base > 3) base = 0;
                                temp_donor_seq = temp_donor_seq << 2 | (base ^ 0x3);
                            }
                        }
                    }
                    
                    splscore = SpliceSiteDB::probscore(temp_donor_seq, temp_acceptor_seq);
                }
                // daehwan - for debugging purposes
                // choose a splice site with the better score
                if((maxspldir == SPL_UNKNOWN && spldir == SPL_UNKNOWN && maxscore < tempscore) ||
                   (maxspldir == SPL_UNKNOWN && spldir == SPL_UNKNOWN && maxscore == tempscore && semi_canonical) ||
                   (maxspldir != SPL_UNKNOWN && spldir != SPL_UNKNOWN && (maxscore < tempscore || (maxscore == tempscore && maxsplscore < splscore))) ||
                   (maxspldir == SPL_UNKNOWN && spldir != SPL_UNKNOWN)) {
                    maxscore = tempscore;
                    maxscorei = i;
                    maxspldir = spldir;
                    maxsplscore = splscore;
                    if(maxspldir != SPL_UNKNOWN) {
                        donor_seq = temp_donor_seq;
                        acceptor_seq = temp_acceptor_seq;
                    } else {
                        donor_seq = 0;
                        acceptor_seq = 0;
                    }
                }
            }
        } else {
            // discover an insertion or a deletion
            assert(ins || del);
            int inslen = (ins ? rddif - refdif : 0);
            int dellen = (del ? refdif - rddif : 0);
            int64_t gap_penalty;
            if(ins) {
                gap_penalty = -(sc.refGapOpen() + sc.refGapExtend() * (inslen - 1));
            } else {
                assert(del);
                gap_penalty = -(sc.readGapOpen() + sc.readGapExtend() * (dellen - 1));
            }
            if(gap_penalty < remainsc) return false;
            int i;
            for(i = 0; i < (int)len; i++) {
                int rdc = seq[this_rdoff + i], rfc = refbuf[i];
                if(i > 0) {
                    temp_scores[i] = temp_scores[i-1];
                } else {
                    temp_scores[i] = 0;
                }
                if(rdc != rfc) {
                    temp_scores[i] += sc.score(rdc, 1 << rfc, qual[this_rdoff + i] - 33);
                }
                if(temp_scores[i] + gap_penalty < remainsc) {
                    break;
                }
            }
            int i_limit = min<int>(i, len);
            int i2;
            for(i2 = len - 1; i2 >= 0; i2--) {
                int rdc = seq[this_rdoff + i2], rfc = refbuf2[i2];
                if((index_t)(i2 + 1) < len) {
                    temp_scores2[i2] = temp_scores2[i2+1];
                } else {
                    temp_scores2[i2] = 0;
                }
                if(rdc != rfc) {
                    temp_scores2[i2] += sc.score(rdc, 1 << rfc, qual[this_rdoff + i2] - 33);
                }
                if(temp_scores2[i2] + gap_penalty < remainsc) {
                    break;
                }
            }
            int i2_limit = (i2 < inslen ? 0 : i2 - inslen);
            for(i = i2_limit, i2 = i2_limit + 1 + inslen;
                i < i_limit && i2 < (int)len;
                i++, i2++) {
                int64_t tempscore = temp_scores[i] + temp_scores2[i2] + gap_penalty;
                if(maxscore < tempscore) {
                    maxscore = tempscore;
                    maxscorei = i;
                }
            }
        }
        if(maxscore == MIN_I64) return false;
        assert_lt(maxscorei, len);
        if(spliced && spliceSite == NULL) {
            uint32_t shorter_anchor_len = min<uint32_t>(maxscorei + 1, len - maxscorei - 1);
            assert_leq(this_toff, other_toff);
            if(maxspldir == SPL_SEMI_FW || maxspldir == SPL_SEMI_RC || maxspldir == SPL_UNKNOWN) {
                if(shorter_anchor_len < minAnchorLen_noncan) {
                    float intronLenProb = intronLen_prob_noncan(shorter_anchor_len, other_toff - this_toff, maxIntronLen);
                    if(intronLenProb > 0.01f)
                        return false;
                }
            } else {
                if(shorter_anchor_len < minAnchorLen) {
                    float intronLenProb = intronLen_prob(shorter_anchor_len, other_toff - this_toff, maxIntronLen);
                    if(intronLenProb > 0.01f)
                        return false;
                }
            }
        }
        if(maxscore < remainsc)
            return false;
    }
    
    bool clear = true;
    for(int i = (int)_edits->size() - 1; i >= 0; i--) {
        const Edit& edit = (*_edits)[i];
        if(edit.type == EDIT_TYPE_SPL ||
           edit.type == EDIT_TYPE_READ_GAP ||
           edit.type == EDIT_TYPE_REF_GAP ||
           (edit.type == EDIT_TYPE_MM && edit.snpID != (index_t)INDEX_MAX)) {
            _edits->resize(i+1);
            clear = false;
            break;
        }
    }
    if(clear) this->_edits->clear();
    // combine two alignments while updating edits
    if(spliced) {
        assert_geq(this_rdoff, this->_rdoff);
        index_t addoff = this_rdoff - this->_rdoff;
        int rd_gap_off = -min<int>(splice_gap_off, 0);
        int ref_gap_off = max<int>(splice_gap_off, 0);
        for(int i = 0; i < (int)len; i++) {
            assert_lt(this_rdoff + i, rdlen);
            int rdc = seq[this_rdoff + i];
            assert_range(0, 4, rdc);
            int rfc;
            if(splice_gap_maxscorei <= maxscorei) {
	      if(i <= (int)splice_gap_maxscorei) {
                    rfc = refbuf[i];
	      } else if(i <= (int)maxscorei) {
                    rfc = refbuf[i - ref_gap_off + rd_gap_off];
                } else {
                    rfc = refbuf2[i];
                }
            } else {
	      if(i <= (int)maxscorei) {
                    rfc = refbuf[i];
	      } else if(i <= (int)splice_gap_maxscorei) {
                    rfc = refbuf2[i + ref_gap_off - rd_gap_off];
                } else {
                    rfc = refbuf2[i];
                }
            }
            assert_range(0, 4, rfc);
            if(rdc != rfc) {
                Edit e((uint32_t)(i + addoff), rfc, rdc, EDIT_TYPE_MM, false);
                _edits->push_back(e);
            }
            if(i == (int)maxscorei) {
                index_t left = this_toff + i + 1;
                if(splice_gap_maxscorei <= maxscorei) {
                    left = left - ref_gap_off + rd_gap_off;
                }
                index_t right = other_toff + other_len - (len - i - 1);
                if(splice_gap_maxscorei > maxscorei) {
                    right = right + ref_gap_off - rd_gap_off;
                }
                index_t skipLen = 0;
                assert_lt(left, right);
                skipLen = right - left;
                Edit e((uint32_t)(i + 1 + addoff), 0, 0, EDIT_TYPE_SPL, skipLen, maxspldir, spliceSite != NULL, false);
                e.donor_seq = donor_seq;
                e.acceptor_seq = acceptor_seq;
                _edits->push_back(e);
            }
            if(i == (int)splice_gap_maxscorei && splice_gap_off != 0) {
                if(rd_gap_off > 0) {
                    assert_lt(left, right);
                    for(index_t j = 0; j < (index_t)rd_gap_off; j++) {
                        int temp_rfc_off = i + 1 + j;
                        int temp_rfc;
                        if(i < (int)maxscorei) {
                            temp_rfc = refbuf[temp_rfc_off];
                        } else {
                            temp_rfc = refbuf2[temp_rfc_off - rd_gap_off];
                        }
                        assert_range(0, 4, temp_rfc);
                        Edit e((uint32_t)(i + 1 + addoff), "ACGTN"[temp_rfc], '-', EDIT_TYPE_READ_GAP);
                        _edits->push_back(e);
                    }
                } else {
                    assert_gt(ref_gap_off, 0);
                    for(index_t j = 0; j < (index_t)ref_gap_off; j++) {
                        assert_lt(this_rdoff + i + 1 + j, rdlen);
                        int temp_rdc = seq[this_rdoff + i + 1 + j];
                        assert_range(0, 4, temp_rdc);
                        Edit e((uint32_t)(i + 1 + j + addoff), '-', "ACGTN"[temp_rdc], EDIT_TYPE_REF_GAP);
                        _edits->push_back(e);
                    }
                    i += ref_gap_off;
                }
            }
        }
    } else {
        index_t ins_len = 0;
        for(index_t i = 0; i < len; i++) {
            char rdc = seq[this_rdoff + i];
            char rfc = (i <= maxscorei ? refbuf[i] : refbuf2[i]);
            assert_geq(this_rdoff, this->_rdoff);
            index_t addoff = this_rdoff - this->_rdoff;
            if(rdc != rfc) {
                ALT<index_t> cmp_alt;
                assert_geq(this_toff, this->_toff);
                cmp_alt.pos = this->_joinedOff + i + (this_toff - this->_toff) - ins_len;
                index_t alt_i = (index_t)altdb.alts().bsearchLoBound(cmp_alt);
                index_t add_alt_i = std::numeric_limits<index_t>::max();
                for(; alt_i < altdb.alts().size(); alt_i++) {
                    const ALT<index_t>& alt = altdb.alts()[alt_i];
                    if(alt.left > cmp_alt.pos) break;
                    if(alt.type != ALT_SNP_SGL) continue;
                    if(alt.seq == rdc) {
                        add_alt_i = alt_i;
                        break;
                    }
                }
                
                Edit e((uint32_t)(i + addoff), rfc, rdc, EDIT_TYPE_MM, false, add_alt_i);
                _edits->push_back(e);
            }
            if(i == maxscorei) {
                index_t left = this_toff + i + 1;
                if(other_toff + other_len < len - i - 1)
                    return false;
                index_t right = other_toff + other_len - (len - i - 1);
                index_t skipLen = 0;
                if(del) {
                    assert_lt(left, right);
                    skipLen = right - left;
                    for(index_t j = 0; j < skipLen; j++) {
                        int temp_rfc;
                        if(i + 1 + j < len) temp_rfc = refbuf[i + 1 + j];
                        else                temp_rfc = ref.getBase(this->_tidx, this_toff + i + 1 + j);
                        assert_range(0, 4, temp_rfc);
                        Edit e((uint32_t)(i + 1 + addoff), "ACGTN"[temp_rfc], '-', EDIT_TYPE_READ_GAP);
                        _edits->push_back(e);
                    }
                } else {
                    assert(ins);
                    assert_lt(right, left);
                    skipLen = left - right;
                    for(index_t j = 0; j < skipLen; j++) {
                        assert_lt(this_rdoff + i + 1 + j, seq.length());
                        int temp_rdc = seq[this_rdoff + i + 1 + j];
                        assert_range(0, 4, temp_rdc);
                        Edit e((uint32_t)(i + 1 + j + addoff), '-', "ACGTN"[temp_rdc], EDIT_TYPE_REF_GAP);
                        _edits->push_back(e);
                    }
                    i += skipLen;
                    ins_len += skipLen;
                }
            }
        }
    }
    index_t fsi = (index_t)otherHit._edits->size();
    for(index_t i = 0; i < otherHit._edits->size(); i++) {
        const Edit& edit = (*otherHit._edits)[i];
        if(edit.type == EDIT_TYPE_SPL ||
           edit.type == EDIT_TYPE_READ_GAP ||
           edit.type == EDIT_TYPE_REF_GAP ||
           (edit.type == EDIT_TYPE_MM && edit.snpID != (index_t)INDEX_MAX)) {
            fsi = i;
            break;
        }
    }
    assert_leq(this->_rdoff, otherHit._rdoff);
    index_t addoff = otherHit._rdoff - this->_rdoff;
    for(index_t i = fsi; i < otherHit._edits->size(); i++) {
        _edits->push_back((*otherHit._edits)[i]);
        _edits->back().pos += addoff;
    }
    // for alignment involving indel, left align so that
    // indels go to the left most of the combined alignment
    if(ins || del || (spliced && splice_gap_off != 0)) {
        leftAlign(rd);
    }
    
    // update alignment score, trims
    assert_leq(this->_rdoff + this->_len, otherHit._rdoff + otherHit._len);
    _len = otherHit._rdoff + otherHit._len - this->_rdoff;
    assert_eq(_trim3, 0);
    _trim3 += otherHit._trim3;
    calculateScore(
                   rd,
                   ssdb,
                   sc,
                   minK_local,
                   minIntronLen,
                   maxIntronLen,
                   minAnchorLen,
                   minAnchorLen_noncan,
                   ref);
#ifndef NDEBUG
    if(_joinedOff != (index_t)INDEX_MAX) {
        ASSERT_ONLY(bool straddled = false);
        ASSERT_ONLY(index_t tmp_tidx = 0, tmp_toff = 0, tmp_tlen = 0);
        gfm.joinedToTextOff(
                            0,
                            _joinedOff,
                            tmp_tidx,
                            tmp_toff,
                            tmp_tlen,
                            true,        // reject straddlers?
                            straddled);  // straddled?
        assert_eq(tmp_tidx, _tidx);
        assert_eq(tmp_toff, _toff);
    }
#endif
    assert(repOk(rd, ref));
    return true;
}

/**
 * Extend the partial alignment (GenomeHit) bidirectionally
 */
template <typename index_t>
bool GenomeHit<index_t>::extend(
                                const Read&                rd,
                                const GFM<index_t>&        gfm,
                                const BitPairReference&    ref,
                                const ALTDB<index_t>&      altdb,
                                const RepeatDB<index_t>&   repeatdb,
                                SpliceSiteDB&              ssdb,
                                SwAligner&                 swa,
                                SwMetrics&                 swm,
                                PerReadMetrics&            prm,
                                const Scoring&             sc,
                                TAlScore                   minsc,
                                RandomSource&              rnd,           // pseudo-random source
                                index_t                    minK_local,
                                const TranscriptomePolicy& tpol,
                                const GraphPolicy&         gpol,
                                index_t&                   leftext,
                                index_t&                   rightext,
                                index_t                    mm)
{
    assert_lt(this->_tidx, ref.numRefs());
    index_t max_leftext = leftext, max_rightext = rightext;
    assert(max_leftext > 0 || max_rightext > 0);
    leftext = 0, rightext = 0;
    index_t rdlen = (index_t)rd.length();
    bool doLeftAlign = false;
    assert(_sharedVars != NULL);
    
    const index_t minIntronLen = tpol.minIntronLen();
    const index_t maxIntronLen = tpol.maxIntronLen();
    const index_t minAnchorLen = tpol.minAnchorLen();
    const index_t minAnchorLen_noncan = tpol.minAnchorLen_noncan();
    
    // extend the alignment further in the left direction
    // with 'mm' mismatches allowed
    const BTDnaString& seq = _fw ? rd.patFw : rd.patRc;
    if(max_leftext > 0 && _rdoff > 0) {
        assert_gt(_rdoff, 0);
        index_t left_rdoff, left_len, left_toff;
        this->getLeft(left_rdoff, left_len, left_toff);
        assert_eq(left_rdoff, _rdoff);
        assert_eq(left_toff, _toff);
        if(_toff <= 0) return false;
        int rl = (int)_toff - (int)_rdoff;
        assert_geq(_score, minsc);
        index_t reflen = _rdoff + 10;
        rl -= (reflen - _rdoff);
        if(rl < 0) {
            reflen += rl;
            rl = 0;
        }
        index_t numNs = 0;
        index_t num_prev_edits = (index_t)_edits->size();
        index_t best_ext = alignWithALTs(
                                         altdb.alts(),
                                         altdb.haplotypes(),
                                         altdb.haplotype_maxrights(),
                                         this->_joinedOff,
                                         seq,
                                         this->_rdoff - 1,
                                         this->_rdoff - 1,
                                         this->_rdoff,
                                         ref,
                                         *_sharedVars,
                                         _tidx,
                                         rl,
                                         reflen,
                                         true, /* left? */
                                         gpol,
                                         *this->_edits,
                                         _sharedVars->ht_llist,
                                         *this->_ht_list,
                                         _sharedVars->cmp_ht,
                                         NULL,
                                         mm,
                                         &numNs);
        // Do not allow for any edits including known snps and splice sites when extending zero-length hit
        if(_len == 0 && mm == 0 && _edits->size() > 0) {
            _edits->clear();
            return false;
        }
        if(best_ext > 0) {
            leftext = best_ext;
            assert_leq(num_prev_edits, _edits->size());
            index_t added_edits = (index_t)_edits->size() - num_prev_edits;
            int ref_ext = (int)best_ext;
            for(index_t i = 0; i < added_edits; i++) {
                const Edit& edit = (*_edits)[i];
                if(edit.type == EDIT_TYPE_REF_GAP)       ref_ext--;
                else if(edit.type == EDIT_TYPE_READ_GAP) ref_ext++;
                else if(edit.type == EDIT_TYPE_SPL)      ref_ext += edit.splLen;
            }
            assert_leq(best_ext, _rdoff);
            _rdoff -= best_ext;
            assert_leq(ref_ext, _toff);
            _toff -= ref_ext;
            _len += best_ext;
            assert_leq(_len, rdlen);
            assert_leq((int)numNs, ref_ext);
            assert_leq(ref_ext - (int)numNs, _joinedOff);
            _joinedOff -= (ref_ext - (int)numNs);
            for(index_t i = 0; i < _edits->size(); i++) {
                if(i < added_edits) {
                    assert_geq((*_edits)[i].pos, _rdoff);
                    (*_edits)[i].pos -= _rdoff;
                } else {
                    (*_edits)[i].pos += best_ext;
                }
            }
        }
    }
    
    // extend the alignment further in the right direction
    // with 'mm' mismatches allowed
    if(max_rightext > 0 && _rdoff + _len < rdlen) {
        index_t right_rdoff, right_len, right_toff;
        this->getRight(right_rdoff, right_len, right_toff);
        index_t rl = right_toff + right_len;
        assert_eq(_rdoff + _len, right_rdoff + right_len);
        index_t rr = rdlen - (right_rdoff + right_len);
        index_t tlen = ref.approxLen(_tidx);
        if(rl < tlen) {
            index_t reflen = rr + 10;
            if(rl + reflen > tlen) {
                reflen = tlen - rl;
            }
            int ref_ext = (int)_len;
            for(index_t ei = 0; ei < _edits->size(); ei++) {
                const Edit& e = (*_edits)[ei];
                if(e.type == EDIT_TYPE_REF_GAP)       ref_ext--;
                else if(e.type == EDIT_TYPE_READ_GAP) ref_ext++;
                else if(e.type == EDIT_TYPE_SPL)      ref_ext += e.splLen;
                else if(e.type == EDIT_TYPE_MM && e.chr == 'N') ref_ext--;
            }
            index_t best_ext = alignWithALTs(
                                             altdb.alts(),
                                             altdb.haplotypes(),
                                             altdb.haplotype_maxrights(),
                                             this->_joinedOff + ref_ext,
                                             seq,
                                             this->_rdoff,
                                             this->_rdoff + this->_len,
                                             rdlen - (this->_rdoff + this->_len),
                                             ref,
                                             *_sharedVars,
                                             _tidx,
                                             (int)rl,
                                             reflen,
                                             false,
                                             gpol,
                                             *this->_edits,
                                             _sharedVars->ht_llist,
                                             *this->_ht_list,
                                             _sharedVars->cmp_ht,
                                             NULL,
                                             mm);
            // Do not allow for any edits including known snps and splice sites when extending zero-length hit
            if(_len == 0 && mm == 0 && _edits->size() > 0) {
                _edits->clear();
                return false;
            }
            if(best_ext > 0) {
                rightext = best_ext;
                _len += best_ext;
            }
        }
    }
    
#ifndef NDEBUG
    if(_joinedOff != (index_t)INDEX_MAX && seq[_rdoff] < 4) {
        ASSERT_ONLY(bool straddled = false);
        ASSERT_ONLY(index_t tmp_tidx = 0, tmp_toff = 0, tmp_tlen = 0);
        gfm.joinedToTextOff(
                            0,
                            _joinedOff,
                            tmp_tidx,
                            tmp_toff,
                            tmp_tlen,
                            true,        // reject straddlers?
                            straddled);  // straddled?
        if(!gfm.repeat()) {
            assert_eq(tmp_tidx, _tidx);
        }
        assert_eq(tmp_toff, _toff);
    }
#endif
    
    if(doLeftAlign) leftAlign(rd);
    assert_leq(_rdoff + _len, rdlen);
    calculateScore(
                   rd,
                   ssdb,
                   sc,
                   minK_local,
                   minIntronLen,
                   maxIntronLen,
                   minAnchorLen,
                   minAnchorLen_noncan,
                   ref);
    assert(repOk(rd, ref));
    return leftext > 0 || rightext > 0;
}

/**
 * Adjust alignment with respect to SNPs, usually updating Edits
 *
 */
template <typename index_t>
bool GenomeHit<index_t>::adjustWithALT(
                                       index_t                     rdoff,
                                       index_t                     len,
                                       const Coord&                coord,
                                       SharedTempVars<index_t>&    sharedVars,
                                       EList<GenomeHit<index_t> >& genomeHits,
                                       const Read&                 rd,
                                       const GFM<index_t>&         gfm,
                                       const ALTDB<index_t>&       altdb,
                                       const BitPairReference&     ref,
                                       const GraphPolicy&          gpol)
{
    if(gfm.gh().linearFM()) {
        genomeHits.expand();
        genomeHits.back().init(
                               coord.orient(),
                               rdoff,
                               len,
                               0, // trim5
                               0, // trim3
                               (index_t)coord.ref(),
                               (index_t)coord.off(),
                               (index_t)coord.joinedOff(),
                               sharedVars);
        return true;
    }
    index_t width = 1 << (gfm.gh()._offRate + 2);
    EList<pair<index_t, int> >& ssOffs = sharedVars.ssOffs;
    findSSOffs(gfm, altdb, (coord.joinedOff() >= width ? (index_t)(coord.joinedOff() - width) : 0), (index_t)(coord.joinedOff() + width), ssOffs);
    assert_gt(ssOffs.size(), 0);
    bool found = false;
    for(index_t s = 0; s < ssOffs.size(); s++) {
        index_t off = (index_t)coord.off();
        index_t joinedOff = (index_t)coord.joinedOff();
        pair<index_t, int>& ssOff = ssOffs[s];
        if(ssOff.first > 0) {
            assert_neq(ssOff.second, 0);
            if(ssOff.second > 0) {
                off += ssOff.first;
                joinedOff += ssOff.first;
            } else {
                off -= ssOff.first;
                joinedOff -= ssOff.first;
            }
        }
        size_t numGenomeHits = genomeHits.size();
        genomeHits.expand();
        genomeHits.back().init(
                               coord.orient(),
                               rdoff,
                               len,
                               0, // trim5
                               0, // trim3
                               (index_t)coord.ref(),
                               off,
                               joinedOff,
                               sharedVars);
        GenomeHit<index_t>& genomeHit = genomeHits.back();
        EList<pair<index_t, int> >& offDiffs = sharedVars.offDiffs;
        const index_t single_offDiffs_size = findOffDiffs(gfm,
                                                          altdb,
                                                          (genomeHit._joinedOff >= width ? genomeHit._joinedOff - width : 0),
                                                          genomeHit._joinedOff + width,
                                                          offDiffs);
        assert_leq(single_offDiffs_size, offDiffs.size());
        
        const BTDnaString& seq = genomeHit._fw ? rd.patFw : rd.patRc;
        const EList<ALT<index_t> >& alts = altdb.alts();
        
        index_t orig_joinedOff = genomeHit._joinedOff;
        index_t orig_toff = genomeHit._toff;
        bool found2 = false;
        // maxAltsTried is not directly related to the size of offDiffs,
        // but let's make the size of offDiffs is determined by maxAltsTried
        const index_t max_offDiffs_size = max<index_t>(4, gpol.maxAltsTried() / 4);
        if(offDiffs.size() - single_offDiffs_size > max_offDiffs_size) offDiffs.resize(single_offDiffs_size + max_offDiffs_size);
        for(index_t o = 0; o < offDiffs.size() && !found2; o++) {
            const pair<index_t, int>& offDiff = offDiffs[o];
#ifndef NDEBUG
            if(o == 0) {
                assert_eq(offDiff.first, 0);
                assert_eq(offDiff.second, 0);
            }
#endif
            if(offDiff.second >= 0) {
                genomeHit._joinedOff = orig_joinedOff + offDiff.first;
                genomeHit._toff = orig_toff + offDiff.first;
            } else {
                if(orig_toff < offDiff.first) continue;
                assert_geq(orig_joinedOff, offDiff.first);
                genomeHit._joinedOff = orig_joinedOff - offDiff.first;
                genomeHit._toff = orig_toff - offDiff.first;
            }
            
            genomeHit._edits->clear();
            ELList<Edit, 128, 4>& candidate_edits = sharedVars.candidate_edits;
            candidate_edits.clear();
            index_t reflen = genomeHit._len + 10;
            index_t alignedLen = alignWithALTs(
                                               alts,
                                               altdb.haplotypes(),
                                               altdb.haplotype_maxrights(),
                                               genomeHit._joinedOff,
                                               seq,
                                               genomeHit._rdoff,
                                               genomeHit._rdoff,
                                               genomeHit._len,
                                               ref,
                                               sharedVars,
                                               genomeHit._tidx,
                                               (int)genomeHit._toff,
                                               reflen,
                                               false, /* left? */
                                               gpol,
                                               *genomeHit._edits,
                                               sharedVars.ht_llist,
                                               *genomeHit._ht_list,
                                               sharedVars.cmp_ht,
                                               &candidate_edits);
            if(alignedLen == genomeHit._len) {
                found2 = true;
                assert(genomeHit.repOk(rd, ref));
                for(index_t i = 0; i < genomeHits.size() - 1; i++) {
                    if(genomeHits[i] == genomeHits.back()) {
                        found2 = false;
                    }
                }
                if(found2) {
                    for(index_t e = 0; e < candidate_edits.size(); e++) {
                        genomeHits.expand();
                        genomeHits.back() = genomeHits[genomeHits.size() - 2];
                        *(genomeHits.back()._edits) = candidate_edits[e];
                        assert(genomeHits.back().repOk(rd, ref));
                        for(size_t i = 0; i < genomeHits.size() - 1; i++) {
                            if(genomeHits[i] == genomeHits.back()) {
                                genomeHits.pop_back();
                                break;
                            }
                        }
                    }
                }
            } else {
                genomeHit._edits->clear();
            }
        }
        if(!found2) genomeHits.pop_back();
        found = genomeHits.size() > numGenomeHits;
    }
    return found;
}

/**
 * Adjust alignment with respect to SNPs, usually updating Edits
 *
 */
template <typename index_t>
bool GenomeHit<index_t>::adjustWithALT(
                                       const Read&             rd,
                                       const GFM<index_t>&     gfm,
                                       const ALTDB<index_t>&   altdb,
                                       const BitPairReference& ref,
                                       const GraphPolicy&      gpol)
{
    if(gfm.gh().linearFM()) return true;
    assert_lt(this->_tidx, ref.numRefs());
    
    assert(_sharedVars != NULL);
    EList<pair<index_t, int> >& offDiffs = _sharedVars->offDiffs;
    index_t width = 1 << (gfm.gh()._offRate + 2);
    const index_t single_offDiffs_size = findOffDiffs(gfm,
                                                      altdb,
                                                      (this->_joinedOff >= width ? this->_joinedOff - width : 0),
                                                      this->_joinedOff + width,
                                                      offDiffs);
    assert_leq(single_offDiffs_size, offDiffs.size());
    
    const BTDnaString& seq = _fw ? rd.patFw : rd.patRc;
    const EList<ALT<index_t> >& alts = altdb.alts();
    
    index_t orig_joinedOff = this->_joinedOff;
    index_t orig_toff = this->_toff;
    bool found = false;
    // maxAltsTried is not directly related to the size of offDiffs,
    // but let's make the size of offDiffs is determined by maxAltsTried
    const index_t max_offDiffs_size = max<index_t>(4, gpol.maxAltsTried() / 4);
    if(offDiffs.size() - single_offDiffs_size > max_offDiffs_size) offDiffs.resize(single_offDiffs_size + max_offDiffs_size);
    for(index_t o = 0; o < offDiffs.size() && !found; o++) {
        const pair<index_t, int>& offDiff = offDiffs[o];
#ifndef NDEBUG
        if(o == 0) {
            assert_eq(offDiff.first, 0);
            assert_eq(offDiff.second, 0);
        }
#endif
        if(offDiff.second >= 0) {
            this->_joinedOff = orig_joinedOff + offDiff.first;
            this->_toff = orig_toff + offDiff.first;
        } else {
            if(orig_toff < offDiff.first) continue;
            assert_geq(orig_joinedOff, offDiff.first);
            this->_joinedOff = orig_joinedOff - offDiff.first;
            this->_toff = orig_toff - offDiff.first;
        }
        index_t reflen = this->_len + 10;
        index_t alignedLen = alignWithALTs(
                                           alts,
                                           altdb.haplotypes(),
                                           altdb.haplotype_maxrights(),
                                           this->_joinedOff,
                                           seq,
                                           this->_rdoff,
                                           this->_rdoff,
                                           this->_len,
                                           ref,
                                           *_sharedVars,
                                           this->_tidx,
                                           (int)this->_toff,
                                           reflen,
                                           false, /* left? */
                                           gpol,
                                           *this->_edits,
                                           _sharedVars->ht_llist,
                                           *this->_ht_list,
                                           _sharedVars->cmp_ht,
                                           &_sharedVars->candidate_edits);
        if(alignedLen == this->_len) {
            found = true;
        } else {
            this->_edits->clear();
        }
    }
#ifndef NDEBUG
    if(found) {
        assert(repOk(rd, ref));
    }
#endif
    return found;
}

/*
 * Find offset differences due to splice sites
 */
template <typename index_t>
void GenomeHit<index_t>::findSSOffs(
                                    const GFM<index_t>&         gfm,
                                    const ALTDB<index_t>&       altdb,
                                    index_t                     start,
                                    index_t                     end,
                                    EList<pair<index_t, int> >& ssOffs)
{
    ssOffs.clear();
    ssOffs.expand();
    ssOffs.back().first = ssOffs.back().second = 0;
    if(gfm.gh().linearFM() || !altdb.hasSpliceSites()) return;
    const EList<ALT<index_t> >& alts = altdb.alts();
    
    // Find splice sites included in this region
    ALT<index_t> alt_search;
    alt_search.left = start;
    for(index_t i = (index_t)alts.bsearchLoBound(alt_search); i < alts.size(); i++) {
        const ALT<index_t>& alt = alts[i];
        if(alt.left >= end) break;
        if(!alt.splicesite()) continue;
        //
        if(alt.left < alt.right) {
            ssOffs.expand();
            ssOffs.back().first = alt.right - alt.left + 1;
            ssOffs.back().second = 1;
            
            const index_t relax = 5;
            if(alt.right > relax) alt_search.left = alt.right - relax;
            else                  alt_search.left = 0;
            for(index_t j = (index_t)alts.bsearchLoBound(alt_search); j < alts.size(); j++) {
                const ALT<index_t>& alt2 = alts[j];
                if(!alt2.splicesite()) continue;
                if(alt2.left < alt2.right) continue;
                if(alt2.left + alt2.right == alt.left + alt.right) continue;
                if(alt2.left > alt.right + relax) break;
                ssOffs.expand();
                if(alt2.right < alt.left) {
                    ssOffs.back().first = alt.left - alt2.right;
                    ssOffs.back().second = -1;
                } else {
                    ssOffs.back().first = alt2.right - alt.left;
                    ssOffs.back().second = 1;
                }
            }
        } else {
            ssOffs.expand();
            ssOffs.back().first = alt.left - alt.right + 1;
            ssOffs.back().second = -1;
        }
    }

    if(ssOffs.size() > 1) {
        ssOffs.sort();
        index_t new_size = (index_t)(unique(ssOffs.begin(), ssOffs.end()) - ssOffs.begin());
        ssOffs.resize(new_size);
    }
}


/*
 * Find offset differences due to indels
 */
template <typename index_t>
index_t GenomeHit<index_t>::findOffDiffs(
                                         const GFM<index_t>&         gfm,
                                         const ALTDB<index_t>&       altdb,
                                         index_t                     start,
                                         index_t                     end,
                                         EList<pair<index_t, int> >& offDiffs)
{
    offDiffs.clear();
    offDiffs.expand();
    offDiffs.back().first = offDiffs.back().second = 0;
    if(gfm.gh().linearFM()) return offDiffs.size();
    const EList<ALT<index_t> >& alts = altdb.alts();
    pair<index_t, index_t> alt_range;
    
    // Find SNPs included in this region
    {
        ALT<index_t> alt_search;
        alt_search.pos = start;
        alt_range.first = alt_range.second = (index_t)alts.bsearchLoBound(alt_search);
        for(alt_range.second = alt_range.first; alt_range.second < alts.size(); alt_range.second++) {
            const ALT<index_t>& alt = alts[alt_range.second];
            if(alt.splicesite() && alt.left > alt.right) continue;
            if(alt.deletion() && alt.reversed) continue;
            if(alt.pos >= end) break;
        }
    }
    if(alt_range.first >= alt_range.second) return offDiffs.size();
    
    for(index_t second = alt_range.second; second > alt_range.first; second--) {
        assert_leq(second, alts.size());
        const ALT<index_t>& alt = alts[second - 1];
        if(!alt.gap() || alt.splicesite() || (alt.deletion() && alt.reversed))
            continue;
        int off = 0;
        if(alt.type == ALT_SNP_DEL) {
            off = alt.len;
        } else {
            assert_eq(alt.type, ALT_SNP_INS);
            off = -alt.len;
        }
        assert_neq(off, 0);
        offDiffs.expand();
        offDiffs.back().first = abs(off);
        offDiffs.back().second = (off > 0 ? 1 : -1);
    }

    if(offDiffs.size() > 1) {
        offDiffs.sort();
        index_t new_size = (index_t)(unique(offDiffs.begin(), offDiffs.end()) - offDiffs.begin());
        offDiffs.resize(new_size);
    }
    
    const index_t single_offDiffs_size = offDiffs.size();
    for(index_t second = alt_range.second; second > alt_range.first; second--) {
        assert_leq(alt_range.second, alts.size());
        const ALT<index_t>& alt = alts[second - 1];
        if(!alt.gap() || alt.splicesite() || (alt.deletion() && alt.reversed))
            continue;
        int off = 0;
        if(alt.type == ALT_SNP_DEL) {
            off = alt.len;
        } else {
            assert_eq(alt.type, ALT_SNP_INS);
            off = -alt.len;
        }
        for(index_t second2 = second - 1; second2 > alt_range.first; second2--) {
            const ALT<index_t>& alt2 = alts[second2 - 1];
            if(!alt2.gap() || alt2.splicesite() || (alt2.deletion() && alt2.reversed))
                continue;
            if(alt2.type == ALT_SNP_DEL) {
                if(alt2.pos + alt2.len >= alt.pos)
                    continue;
                off += alt2.len;
            } else {
                assert_eq(alt2.type, ALT_SNP_INS);
                if(alt2.pos >= alt.pos)
                    continue;
                off -= alt2.len;
            }
            bool found = false;
            for(index_t i = 0; i < offDiffs.size(); i++) {
                int off_cmp = offDiffs[i].first * offDiffs[i].second;
                if(off == off_cmp) {
                    found = true;
                    break;
                }
            }
            if(!found) {
                offDiffs.expand();
                offDiffs.back().first = abs(off);
                offDiffs.back().second = (off > 0 ? 1 : -1);
            }
        }
    }
    
    return single_offDiffs_size;
}


/*
 *
 */
template <typename index_t>
void add_haplotypes(
                    const EList<ALT<index_t> >&       alts,
                    const EList<Haplotype<index_t> >& haplotypes,
                    const EList<index_t>&             haplotype_maxrights,
                    Haplotype<index_t>&               cmp_ht,
                    EList<pair<index_t, index_t> >&   ht_list,
                    index_t                           rdlen,
                    bool                              left_ext = true,
                    bool                              initial = false)
{
    pair<int, int> ht_range;
    ht_range.first = ht_range.second = (int)haplotypes.bsearchLoBound(cmp_ht);
    if(ht_range.first >= haplotypes.size())
        return;
    
    if(left_ext) {
        for(; ht_range.first >= 0; ht_range.first--) {
            const Haplotype<index_t>& ht = haplotypes[ht_range.first];
            if(!initial) {
                if(ht.right >= cmp_ht.left) continue;
            }
            index_t ht_maxright = haplotype_maxrights[ht_range.first];
            assert_geq(ht_maxright, ht.right);
            if(ht_maxright + rdlen - 1 < cmp_ht.left) break;
            if(ht.alts.size() <= 0) continue;
            bool added = false;
            for(index_t h = 0; h < ht_list.size(); h++) {
                if(ht_list[h].first == ht_range.first) {
                    added = true;
                    break;
                }
            }
            if(added) continue;
            ht_list.expand();
            ht_list.back().first = ht_range.first;
            assert_gt(ht.alts.size(), 0);
            if(ht.right < cmp_ht.left) {
                ht_list.back().second = ht.alts.size() - 1;
            } else {
                assert(initial);
                ht_list.back().second = ht.alts.size();
                for(int a = (int)ht.alts.size() - 1; a >= 0; a--) {
                    index_t alti = ht.alts[a];
                    assert_lt(alti, alts.size());
                    const ALT<index_t>& alt = alts[alti];
                    assert(alt.snp());
                    ht_list.back().second = (index_t)a;
                    if(cmp_ht.left > alt.pos) break;
                }
                if(ht_list.back().second == ht.alts.size()) {
                    ht_list.pop_back();
                }
            }
        }
    } else {
        if(initial) {
            for(; ht_range.first >= 0; ht_range.first--) {
                const Haplotype<index_t>& ht = haplotypes[ht_range.first];
                index_t ht_maxright = haplotype_maxrights[ht_range.first];
                assert_geq(ht_maxright, ht.right);
                if(ht_maxright < cmp_ht.left) break;
                if(ht.right < cmp_ht.left || ht.left > cmp_ht.left) continue;
                if(ht.alts.size() <= 0) continue;
                bool added = false;
                for(index_t h = 0; h < ht_list.size(); h++) {
                    if(ht_list[h].first == ht_range.first) {
                        added = true;
                        break;
                    }
                }
                if(added) continue;
                ht_list.expand();
                ht_list.back().first = ht_range.first;
                assert_gt(ht.alts.size(), 0);
                ht_list.back().second = ht.alts.size();
                for(index_t a = 0; a < ht.alts.size(); a++) {
                    index_t alti = ht.alts[a];
                    assert_lt(alti, alts.size());
                    const ALT<index_t>& alt = alts[alti];
                    assert(alt.snp());
                    ht_list.back().second = a;
                    if(cmp_ht.left <= alt.pos) break;
                }
                if(ht_list.back().second == ht.alts.size()) {
                    ht_list.pop_back();
                }
            }
        }

        for(; ht_range.second < haplotypes.size(); ht_range.second++) {
                const Haplotype<index_t>& ht = haplotypes[ht_range.second];
                if(ht.left < cmp_ht.right) continue;
                if(ht.left >= cmp_ht.right + rdlen) break;
                if(ht.alts.size() <= 0) continue;
                bool added = false;
                for(index_t h = 0; h < ht_list.size(); h++) {
                    if(ht_list[h].first == ht_range.second) {
                        added = true;
                        break;
                    }
                }
                if(added) continue;
                ht_list.expand();
                ht_list.back().first = ht_range.second;
                assert_gt(ht.alts.size(), 0);
                ht_list.back().second = 0;
        }
    }
}


/*
 *
 */
template <typename index_t>
index_t GenomeHit<index_t>::alignWithALTs_recur(
                                                const EList<ALT<index_t> >&       alts,
                                                const EList<Haplotype<index_t> >& haplotypes,
                                                const EList<index_t>&             haplotype_maxrights,
                                                index_t                           joinedOff,
                                                const BTDnaString&                rdseq,
                                                index_t                           rdoff_add,
                                                index_t                           rdoff,
                                                index_t                           rdlen,
                                                const BitPairReference&           ref,
                                                EList<SStringExpandable<char> >&  raw_refbufs,
                                                ASSERT_ONLY(SStringExpandable<uint32_t> destU32,)
                                                EList<Edit>&                      tmp_edits,
                                                int&                              best_rdoff,
                                                const char*                       rfseq,
                                                index_t                           tidx,
                                                int                               rfoff,
                                                index_t                           rflen,
                                                bool                              left,
                                                EList<Edit>&                      edits,
                                                index_t                           mm,
                                                ELList<pair<index_t, index_t> >&  ht_llist,
                                                Haplotype<index_t>&               cmp_ht,
                                                ELList<Edit, 128, 4>*             candidate_edits,
                                                index_t                           tmp_numNs,
                                                index_t*                          numNs,
                                                index_t                           dep,
                                                const GraphPolicy&                gpol,
                                                index_t&                          numALTsTried,
                                                ALT_TYPE                          prev_alt_type)
{
    if(numALTsTried > gpol.maxAltsTried() + dep) return 0;
    assert_gt(rdlen, 0);
    assert_gt(rflen, 0);
    if(ht_llist.size() <= dep) ht_llist.expand();
    if(raw_refbufs.size() <= dep) raw_refbufs.expand();
    if(rfoff < -16) return 0;
    size_t contig_len = ref.approxLen(tidx);
    if(rfoff >= contig_len) return 0;
    if(rfoff >= 0 && rfoff + rflen > contig_len) {
        rflen = contig_len - rfoff;
    } else if(rfoff < 0 && rflen > contig_len) {
        rflen = contig_len;
    }
    if(rflen == 0) return 0;
    if(rfseq == NULL) {
        SStringExpandable<char>& raw_refbuf = raw_refbufs[dep];
        raw_refbuf.resize(rflen + 16 + 16);
        raw_refbuf.fill(0x4);
        int off = ref.getStretch(
                                 reinterpret_cast<uint32_t*>(raw_refbuf.wbuf() + 16),
                                 tidx,
                                 max<int>(rfoff, 0),
                                 rfoff > 0 ? rflen : rflen + rfoff
                                 ASSERT_ONLY(, destU32));
        assert_lt(off, 16);
        rfseq = raw_refbuf.wbuf() + 16 + off + min<int>(rfoff, 0);
    }
    
    if(left) {
        index_t tmp_mm = 0;
        int min_rd_i = (int)rdoff;
        int mm_min_rd_i = (int)rdoff;
        index_t mm_tmp_numNs = 0;
        for(int rf_i = (int)rflen - 1; rf_i >= 0 && mm_min_rd_i >= 0; rf_i--, mm_min_rd_i--) {
            int rf_bp = rfseq[rf_i];
            int rd_bp = rdseq[mm_min_rd_i];
            if(rf_bp != rd_bp || rd_bp == 4) {
                if(tmp_mm == 0) {
                    min_rd_i = mm_min_rd_i;
                }
                if(tmp_mm >= mm) break;
                tmp_mm++;
                Edit e(
                       mm_min_rd_i,
                       "ACGTN"[rf_bp],
                       "ACGTN"[rd_bp],
                       EDIT_TYPE_MM);
                tmp_edits.insert(e, 0);
            }
            if(rf_bp == 4) {
                if(tmp_mm == 0) tmp_numNs++;
                mm_tmp_numNs++;
            }
        }
        if(tmp_mm == 0) {
            min_rd_i = mm_min_rd_i;
        }
        if(mm_min_rd_i < best_rdoff) {
            best_rdoff = mm_min_rd_i;
            edits = tmp_edits;
            if(numNs != NULL) *numNs = mm_tmp_numNs;
        }
        if(mm_min_rd_i < 0) return rdlen;
        if(tmp_mm > 0) {
            tmp_edits.erase(0, tmp_mm);
            tmp_mm = 0;
        }
        
        // Find SNPs included in this region
        pair<int, int> alt_range(0, 0);
        if(alts.size() > 0) {
            ALT<index_t> cmp_alt;
            const index_t minK = 16;
            assert_leq(mm_min_rd_i, rdoff);
            index_t rd_diff = rdoff - mm_min_rd_i;
            rd_diff = (rd_diff > minK ? rd_diff - minK : 0);
            if(gpol.enableCODIS()) {
                rd_diff = 0;
            }
            if(rd_diff >= joinedOff) {
                cmp_alt.pos = joinedOff;
            } else {
                cmp_alt.pos = joinedOff - rd_diff;
            }
            alt_range.first = alt_range.second = (int)alts.bsearchLoBound(cmp_alt);
            if(alt_range.first >= alts.size()) {
                assert_gt(alts.size(), 0);
                alt_range.first = alt_range.second = alt_range.second - 1;
            }
            for(; alt_range.first >= 0; alt_range.first--) {
                const ALT<index_t>& alt = alts[alt_range.first];
                if(alt.snp()) {
                    if(alt.deletion() && !alt.reversed) continue;
                    if(alt.pos + rdlen < joinedOff) break;
                } else if(alt.splicesite()) {
                    if(alt.left < alt.right) continue;
                    if(alt.left + rdlen - 1 < joinedOff) break;
                } else {
                    assert(alt.exon());
                    continue;
                }
            }
        }
        
        // Update and find Haplotypes
        EList<pair<index_t, index_t> >& ht_list = ht_llist[dep];
        ht_list.clear();
        if(gpol.useHaplotype() && haplotypes.size() > 0) {
            if(dep > 0) {
                EList<pair<index_t, index_t> >& ht_prev_list = ht_llist[dep-1];
                for(index_t p = 0; p < ht_prev_list.size(); p++) {
                    const pair<index_t, index_t>& ht_ref = ht_prev_list[p];
                    const Haplotype<index_t>& ht = haplotypes[ht_ref.first];
                    assert_lt(ht_ref.second, ht.alts.size());
                    index_t alt_id = ht.alts[ht_ref.second];
                    assert_gt(tmp_edits.size(), 0);
                    const ALT<index_t>& alt = alts[tmp_edits[0].snpID];
                    const ALT<index_t>& ht_alt = alts[alt_id];
                    if(!alt.isSame(ht_alt)) continue;
                    if(ht_ref.second == 0) {
                        cmp_ht.left = cmp_ht.right = joinedOff;
                        add_haplotypes(alts,
                                       haplotypes,
                                       haplotype_maxrights,
                                       cmp_ht,
                                       ht_list,
                                       rdlen);
                    } else {
                        ht_list.push_back(ht_ref);
                        ht_list.back().second--;
                    }
                }
            }
            if(ht_list.size() <= 0) {
                cmp_ht.left = cmp_ht.right = joinedOff;
                add_haplotypes(alts,
                               haplotypes,
                               haplotype_maxrights,
                               cmp_ht,
                               ht_list,
                               rdlen,
                               true, // left_ext?
                               dep == 0); // initial?
            }
        }
        
        assert_geq(rdoff, 0);
        const index_t orig_nedits = (index_t)tmp_edits.size();
        for(; alt_range.second > alt_range.first; alt_range.second--) {
            ALT<index_t> alt = alts[alt_range.second];
            if(alt.pos >= joinedOff) continue;
            if(alt.splicesite()) {
                if(alt.left < alt.right) continue;
                index_t tmp = alt.left;
                alt.left = alt.right;
                alt.right = tmp;
            }
            if(alt.deletion()) {
                if(!alt.reversed) continue;
                alt.pos = alt.pos - alt.len + 1;
            }
            if(alt.exon()) continue;
            bool alt_compatible = false;
            int rf_i = (int)rflen - 1, rd_i = (int)rdoff;
            int diff = 0;
            if(alt.type == ALT_SNP_SGL) {
                diff = joinedOff - alt.pos - 1;
            } else if(alt.type == ALT_SNP_DEL) {
                if(alt.pos + alt.len >= joinedOff) continue;
                diff = joinedOff - (alt.pos + alt.len);
            } else if(alt.type == ALT_SNP_INS) {
                diff = joinedOff - alt.pos;
            } else {
                assert(alt.splicesite());
                diff = joinedOff - (alt.right + 1);
            }
            if(rf_i < diff || rd_i < diff) continue;
            rf_i -= diff;
            rd_i -= diff;
            int rd_bp = rdseq[rd_i];
            if(rd_i < min_rd_i) {
                if(alt.type == ALT_SNP_INS) {
                    if(rd_i + 1 >= min_rd_i) continue;
                }
                break;
            }
            
            // Check to see if there is a haplotype that supports this alt
            if(ht_list.size() > 0 && alt.snp()) {
                bool ht_found = false;
                for(index_t h = 0; h < ht_list.size(); h++) {
                    const pair<index_t, index_t>& ht_ref = ht_list[h];
                    const Haplotype<index_t>& ht = haplotypes[ht_ref.first];
                    assert_lt(ht_ref.second, ht.alts.size());
                    index_t ht_alti = ht.alts[ht_ref.second];
                    const ALT<index_t>& ht_alt = alts[ht_alti];
                    if(alts[alt_range.second].isSame(ht_alt)) {
                        ht_found = true;
                        break;
                    }
                }
                if(!ht_found) continue;
            }

            if(alt.type == ALT_SNP_SGL) {
                if(rd_bp == (int)alt.seq) {
                    int rf_bp = rfseq[rf_i];
                    Edit e(
                           rd_i,
                           "ACGTN"[rf_bp],
                           "ACGTN"[rd_bp],
                           EDIT_TYPE_MM,
                           true, /* chars? */
                           alt_range.second);
                    tmp_edits.insert(e, 0);
                    rd_i--;
                    rf_i--;
                    alt_compatible = true;
                }
            } else if(alt.type == ALT_SNP_DEL) {
                if(rfoff + rf_i > (int)alt.len) {
                    if(rf_i > (int)alt.len) {
                        for(index_t i = 0; i < alt.len; i++) {
                            int rf_bp = rfseq[rf_i - i];
                            Edit e(
                                   rd_i + 1,
                                   "ACGTN"[rf_bp],
                                   '-',
                                   EDIT_TYPE_READ_GAP,
                                   true, /* chars? */
                                   alt_range.second);
                            tmp_edits.insert(e, 0);
                        }
                        
                    } else {
                        // long deletions
                        int new_rfoff = rfoff - alt.len;
                        index_t new_rflen = rf_i + alt.len + 10;
                        if(raw_refbufs.size() <= dep + 1) raw_refbufs.expand();
                        SStringExpandable<char>& raw_refbuf = raw_refbufs[dep + 1];
                        raw_refbuf.resize(new_rflen + 16 + 16);
                        raw_refbuf.fill(0x4);
                        int off = ref.getStretch(
                                                 reinterpret_cast<uint32_t*>(raw_refbuf.wbuf() + 16),
                                                 tidx,
                                                 max<int>(new_rfoff, 0),
                                                 new_rfoff > 0 ? new_rflen : new_rflen + new_rfoff
                                                 ASSERT_ONLY(, destU32));
                        assert_lt(off, 16);
                        const char* new_rfseq = raw_refbuf.wbuf() + 16 + off + min<int>(new_rfoff, 0);
                        for(int i = 0; i < alt.len; i++) {
                            int rf_bp = new_rfseq[rf_i - i + alt.len];
                            Edit e(
                                   rd_i + 1,
                                   "ACGTN"[rf_bp],
                                   '-',
                                   EDIT_TYPE_READ_GAP,
                                   true, /* chars? */
                                   alt_range.second);
                            tmp_edits.insert(e, 0);
                        }
                    }
                    rf_i -= (int)alt.len;
                    alt_compatible = true;
                }
            } else if(alt.type == ALT_SNP_INS) {
                if(rd_i > (int)alt.len) {
                    bool same_seq = true;
                    for(index_t i = 0; i < alt.len; i++) {
                        rd_bp = rdseq[rd_i - i];
                        int snp_bp = (alt.seq >> (i << 1)) & 0x3;
                        if(rd_bp != snp_bp) {
                            same_seq = false;
                            break;
                        }
                        Edit e(
                               rd_i - i,
                               '-',
                               "ACGTN"[rd_bp],
                               EDIT_TYPE_REF_GAP,
                               true, /* chars? */
                               alt_range.second);
                        tmp_edits.insert(e, 0);
                    }
                    if(same_seq) {
                        rd_i -= (int)alt.len;
                        alt_compatible = true;
                    }
                }
            } else if(alt.type == ALT_SPLICESITE) {
                bool add_splicesite = true;
                if(rd_i == rdoff && prev_alt_type == ALT_SPLICESITE) {
                    add_splicesite = false;
                }
                if(add_splicesite) {
                    assert_lt(rd_i, rflen);
                    assert_lt(alt.left, alt.right);
                    index_t intronLen = alt.right - alt.left + 1;
                    Edit e(rd_i + 1,
                           0,
                           0,
                           EDIT_TYPE_SPL,
                           intronLen,
                           alt.fw ? SPL_FW : SPL_RC,
                           true,   /* known splice site? */
                           false); /* chrs? */
                    tmp_edits.insert(e, 0);
                    alt_compatible = true;
                }
            }
            if(alt_compatible) {
                numALTsTried++;
                assert_leq(rd_i, (int)rdoff);
                if(rd_i < 0) {
                    best_rdoff = rd_i;
                    edits = tmp_edits;
                    return rdlen;
                }
                index_t next_joinedOff = alt.pos;
                int next_rfoff = rfoff, next_rdoff = rd_i;
                const char* next_rfseq = rfseq;
                int next_rflen = rf_i + 1, next_rdlen = rd_i + 1;
                if(alt.splicesite()) {
                    assert_lt(alt.left, alt.right);
                    next_joinedOff = alt.left;
                    index_t intronLen = alt.right - alt.left + 1;
                    assert_geq(next_rfoff, intronLen);
                    next_rfoff -= intronLen;
                    next_rfseq = NULL;
                }
                if(next_rflen < next_rdlen) {
                    int add_len = next_rdlen + 10 - next_rflen;
                    if(next_rfoff < add_len) add_len = next_rfoff;
                    next_rfoff -= add_len;
                    next_rflen += add_len;
                    next_rfseq = NULL;
                }
                index_t alignedLen = alignWithALTs_recur(
                                                         alts,
                                                         haplotypes,
                                                         haplotype_maxrights,
                                                         next_joinedOff,
                                                         rdseq,
                                                         rdoff_add,
                                                         next_rdoff,
                                                         next_rdlen,
                                                         ref,
                                                         raw_refbufs,
                                                         ASSERT_ONLY(destU32,)
                                                         tmp_edits,
                                                         best_rdoff,
                                                         next_rfseq,
                                                         tidx,
                                                         next_rfoff,
                                                         next_rflen,
                                                         left,
                                                         edits,
                                                         mm,
                                                         ht_llist,
                                                         cmp_ht,
                                                         candidate_edits,
                                                         tmp_numNs,
                                                         numNs,
                                                         dep + 1,
                                                         gpol,
                                                         numALTsTried,
                                                         alt.type);
                if(alignedLen == next_rdlen) return rdlen;
            }
            // Restore to the earlier state
            assert_leq(orig_nedits, tmp_edits.size());
            if(orig_nedits < tmp_edits.size()) tmp_edits.erase(0, tmp_edits.size() - orig_nedits);
        }
        return 0;
    } else {
        index_t tmp_mm = 0;
        index_t max_rd_i = 0;
        index_t mm_max_rd_i = 0;
        index_t mm_tmp_numNs = 0;
        for(index_t rf_i = 0; rf_i < rflen && mm_max_rd_i < rdlen; rf_i++, mm_max_rd_i++) {
            int rf_bp = rfseq[rf_i];
            int rd_bp = rdseq[rdoff + mm_max_rd_i];
            if(rf_bp != rd_bp || rd_bp == 4) {
                if(tmp_mm == 0) {
                    max_rd_i = mm_max_rd_i;
                }
                if(tmp_mm >= mm) break;
                tmp_mm++;
                Edit e(
                       mm_max_rd_i + rdoff_add,
                       "ACGTN"[rf_bp],
                       "ACGTN"[rd_bp],
                       EDIT_TYPE_MM);
                tmp_edits.push_back(e);
            }
            if(rf_bp == 4) {
                if(tmp_mm == 0) tmp_numNs++;
                mm_tmp_numNs++;
            }
        }
        if(tmp_mm == 0) {
            max_rd_i = mm_max_rd_i;
        }
        if(mm_max_rd_i + rdoff > best_rdoff) {
            best_rdoff = mm_max_rd_i + rdoff;
            edits = tmp_edits;
            if(numNs != NULL) *numNs = mm_tmp_numNs;
            if(candidate_edits != NULL) candidate_edits->clear();
        } else if(mm_max_rd_i + rdoff == best_rdoff) {
            if(candidate_edits != NULL) {
                candidate_edits->expand();
                candidate_edits->back() = tmp_edits;
            }
        }
        if(mm_max_rd_i == rflen) {
            return mm_max_rd_i;
        }
        
        // Find SNPs included in this region
        pair<index_t, index_t> alt_range;
        {
            ALT<index_t> cmp_alt;
            const index_t minK = 16;
            index_t rd_diff = (max_rd_i > minK ? max_rd_i - minK : 0);
            if(gpol.enableCODIS()) {
                rd_diff = 0;
            }
            cmp_alt.pos = joinedOff + rd_diff;
            alt_range.first = alt_range.second = (index_t)alts.bsearchLoBound(cmp_alt);
            if(alt_range.first >= alts.size()) return 0;
            for(; alt_range.second < alts.size(); alt_range.second++) {
                const ALT<index_t>& alt = alts[alt_range.second];
                if(alt.splicesite()) {
                    if(alt.left > alt.right) continue;
                }
                if(alt.deletion()) {
                    if(alt.reversed) continue;
                }
                if(alt.left > joinedOff + max_rd_i) break;
            }
        }
        if(mm_max_rd_i == rdlen) {
            bool further_search = false;
            for(index_t s = alt_range.first; s < alt_range.second; s++) {
                const ALT<index_t>& alt = alts[s];
                if(alt.splicesite() && alt.left < alt.right) {
                    further_search = true;
                    break;
                }
            }
            if(!further_search) return mm_max_rd_i;
        }
        if(tmp_mm > 0) {
            tmp_edits.resize(tmp_edits.size() - tmp_mm);
            tmp_mm = 0;
        }
        
        // Update and find Haplotypes
        EList<pair<index_t, index_t> >& ht_list = ht_llist[dep];
        ht_list.clear();
        if(gpol.useHaplotype() && haplotypes.size() > 0) {
            if(dep > 0) {
                EList<pair<index_t, index_t> >& ht_prev_list = ht_llist[dep-1];
                for(index_t p = 0; p < ht_prev_list.size(); p++) {
                    const pair<index_t, index_t>& ht_ref = ht_prev_list[p];
                    const Haplotype<index_t>& ht = haplotypes[ht_ref.first];
                    if(ht_ref.second < ht.alts.size()) {
                        index_t alt_id = ht.alts[ht_ref.second];
                        assert_gt(tmp_edits.size(), 0);
                        const ALT<index_t>& alt = alts[tmp_edits.back().snpID];
                        const ALT<index_t>& ht_alt = alts[alt_id];
                        if(!alt.isSame(ht_alt)) continue;
                    }
                    if(ht_ref.second + 1 >= ht.alts.size() && joinedOff > ht.right) {
                        cmp_ht.left = cmp_ht.right = joinedOff;
                        add_haplotypes(alts,
                                       haplotypes,
                                       haplotype_maxrights,
                                       cmp_ht,
                                       ht_list,
                                       rdlen,
                                       false); // left_ext?
                    } else {
                        ht_list.push_back(ht_ref);
                        ht_list.back().second++;
                    }
                }
            }
            if(ht_list.size() <= 0) {
                cmp_ht.left = cmp_ht.right = joinedOff;
                add_haplotypes(alts,
                               haplotypes,
                               haplotype_maxrights,
                               cmp_ht,
                               ht_list,
                               rdlen,
                               false, // left_ext?
                               dep == 0 && rdoff_add == 0); // initial?
            }
        }
        
        const index_t orig_nedits = (index_t)tmp_edits.size();
        for(; alt_range.first < alt_range.second; alt_range.first++) {
            const ALT<index_t>& alt = alts[alt_range.first];
            if(alt.splicesite()) {
                if(alt.left > alt.right) continue;
            }
            if(alt.exon()) continue;
            if(alt.deletion()) {
                if(alt.reversed) continue;
            }
            bool alt_compatible = false;
            assert_leq(joinedOff, alt.pos);
            index_t rf_i, rd_i;
            rf_i = rd_i = alt.pos - joinedOff;
            if(rd_i >= rdlen) continue;
            assert_leq(rd_i, max_rd_i);
            int rf_bp = rfseq[rf_i];
            int rd_bp = rdseq[rdoff + rd_i];
            
            // Check to see if there is a haplotype that supports this alt
            if(ht_list.size() > 0 && alt.snp()) {
                bool ht_found = false;
                for(index_t h = 0; h < ht_list.size(); h++) {
                    const pair<index_t, index_t>& ht_ref = ht_list[h];
                    const Haplotype<index_t>& ht = haplotypes[ht_ref.first];
                    if(ht_ref.second >= ht.alts.size())
                        continue;
                    index_t ht_alti = ht.alts[ht_ref.second];
                    const ALT<index_t>& ht_alt = alts[ht_alti];
                    if(alts[alt_range.first].isSame(ht_alt)) {
                        ht_found = true;
                        break;
                    }
                }
                if(!ht_found) continue;
            }
            
            if(alt.type == ALT_SNP_SGL) {
                if(rd_bp == (int)alt.seq) {
                    Edit e(
                           rd_i + rdoff_add,
                           "ACGTN"[rf_bp],
                           "ACGTN"[rd_bp],
                           EDIT_TYPE_MM,
                           true, /* chars? */
                           alt_range.first);
                    tmp_edits.push_back(e);
                    rd_i++;
                    rf_i++;
                    alt_compatible = true;
                }
            } else if(alt.type == ALT_SNP_DEL) {
                bool try_del = rd_i > 0;
                if(rd_i == 0 && dep > 0) {
                    // Avoid consecutive deletions
                    assert_gt(tmp_edits.size(), 0);
                    const Edit& e = tmp_edits.back();
                    if(e.type != EDIT_TYPE_READ_GAP) {
                        try_del = true;
                    }
                }
                if(try_del) {
                    if(rf_i + alt.len <= rflen) {
                        for(index_t i = 0; i < alt.len; i++) {
                            rf_bp = rfseq[rf_i + i];
                            Edit e(
                                   rd_i + rdoff_add,
                                   "ACGTN"[rf_bp],
                                   '-',
                                   EDIT_TYPE_READ_GAP,
                                   true, /* chars? */
                                   alt_range.first);
                            tmp_edits.push_back(e);
                        }
                        rf_i += alt.len;
                        alt_compatible = true;
                    } else {
                        // long deletions
                        index_t new_rflen = rf_i + alt.len + 10;
                        if(raw_refbufs.size() <= dep + 1) raw_refbufs.expand();
                        SStringExpandable<char>& raw_refbuf = raw_refbufs[dep + 1];
                        raw_refbuf.resize(new_rflen + 16 + 16);
                        raw_refbuf.fill(0x4);
                        int off = ref.getStretch(
                                                 reinterpret_cast<uint32_t*>(raw_refbuf.wbuf() + 16),
                                                 tidx,
                                                 max<int>(rfoff, 0),
                                                 rfoff > 0 ? new_rflen : new_rflen + rfoff
                                                 ASSERT_ONLY(, destU32));
                        assert_lt(off, 16);
                        const char* new_rfseq = raw_refbuf.wbuf() + 16 + off + min<int>(rfoff, 0);
                        for(index_t i = 0; i < alt.len; i++) {
                            rf_bp = new_rfseq[rf_i + i];
                            Edit e(
                                   rd_i + rdoff_add,
                                   "ACGTN"[rf_bp],
                                   '-',
                                   EDIT_TYPE_READ_GAP,
                                   true, /* chars? */
                                   alt_range.first);
                            tmp_edits.push_back(e);
                        }
                        rf_i += alt.len;
                        alt_compatible = true;
                    }
                }
            } else if(alt.type == ALT_SNP_INS) {
                if(rd_i + alt.len <= rdlen && rf_i > 0) {
                    bool same_seq = true;
                    for(index_t i = 0; i < alt.len; i++) {
                        rd_bp = rdseq[rdoff + rd_i + i];
                        int snp_bp = (alt.seq >> ((alt.len - i - 1) << 1)) & 0x3;
                        if(rd_bp != snp_bp) {
                            same_seq = false;
                            break;
                        }
                        Edit e(
                               rd_i + i + rdoff_add,
                               '-',
                               "ACGTN"[rd_bp],
                               EDIT_TYPE_REF_GAP,
                               true, /* chars? */
                               alt_range.first);
                        tmp_edits.push_back(e);
                    }
                    if(same_seq) {
                        rd_i += alt.len;
                        alt_compatible = true;
                    }
                }
            } else if(alt.type == ALT_SPLICESITE) {
                bool try_splice = rd_i > 0;
                if(rd_i == 0 && dep > 0) {
                    // Avoid consecutive introns
                    assert_gt(tmp_edits.size(), 0);
                    const Edit& e = tmp_edits.back();
                    if(e.type != EDIT_TYPE_SPL) {
                        try_splice = true;
                    }
                }
                if(try_splice) {
                    assert_lt(rd_i, rflen);
                    index_t intronLen = alt.right - alt.left + 1;
                    Edit e(rd_i + rdoff_add,
                           0,
                           0,
                           EDIT_TYPE_SPL,
                           intronLen,
                           alt.fw ? SPL_FW : SPL_RC,
                           true,   /* known splice site? */
                           false); /* chrs? */
                    tmp_edits.push_back(e);
                    alt_compatible = true;
                }
            }
            if(alt_compatible) {
                numALTsTried++;
                if(rd_i == rdlen) {
                    assert_leq(best_rdoff, rdoff + rd_i);
                    if(best_rdoff < rdoff + rd_i) {
                        if(candidate_edits != NULL) candidate_edits->clear();
                    }
                    if(candidate_edits != NULL) {
                        candidate_edits->expand();
                        candidate_edits->back() = tmp_edits;
                    }
                    best_rdoff = rdoff + rd_i;
                    edits = tmp_edits;
                    return rd_i;
                }
                index_t next_joinedOff = 0;
                int next_rfoff = rfoff + rf_i, next_rdoff = rdoff + rd_i;
                const char* next_rfseq = rfseq + rf_i;
                index_t next_rflen = rflen - rf_i, next_rdlen = rdlen - rd_i;
                if(alt.type == ALT_SNP_SGL) {
                    next_joinedOff = alt.pos + 1;
                } else if(alt.type == ALT_SNP_DEL) {
                    next_joinedOff = alt.pos + alt.len;
                    if(rflen <= rf_i) {
                        next_rflen = 0; // Will reset next_rfseq and next_rflen below
                    }
                } else if(alt.type == ALT_SNP_INS) {
                    next_joinedOff = alt.pos;
                } else if(alt.type == ALT_SPLICESITE) {
                    next_joinedOff = alt.right + 1;
                    index_t intronLen = alt.right - alt.left + 1;
                    next_rfoff += intronLen;
                    next_rfseq = NULL;
                } else {
                    assert(false);
                }
                if(next_rflen < next_rdlen) {
                    next_rflen = next_rdlen + 10;
                    next_rfseq = NULL;
                }
                index_t alignedLen = alignWithALTs_recur(
                                                         alts,
                                                         haplotypes,
                                                         haplotype_maxrights,
                                                         next_joinedOff,
                                                         rdseq,
                                                         rdoff_add + rd_i,
                                                         next_rdoff,
                                                         next_rdlen,
                                                         ref,
                                                         raw_refbufs,
                                                         ASSERT_ONLY(destU32,)
                                                         tmp_edits,
                                                         best_rdoff,
                                                         next_rfseq,
                                                         tidx,
                                                         next_rfoff,
                                                         next_rflen,
                                                         left,
                                                         edits,
                                                         mm,
                                                         ht_llist,
                                                         cmp_ht,
                                                         candidate_edits,
                                                         tmp_numNs,
                                                         numNs,
                                                         dep + 1,
                                                         gpol,
                                                         numALTsTried,
                                                         alt.type);
                if(alignedLen > 0) {
                    assert_leq(rdoff + rd_i + alignedLen, best_rdoff);
                    bool search_further = false;
                    if(alt.splicesite()) {
                        for(index_t sf = alt_range.first + 1; sf < alt_range.second; sf++) {
                            const ALT<index_t>& alt2 = alts[sf];
                            if(alt2.splicesite() && alt2.left < alt2.right) {
                                search_further = true;
                                break;
                            }
                        }
                    }
                    if(!search_further) {
                        if(rd_i + alignedLen == rdlen) {
                            return rd_i + alignedLen;
                        }
                    }
                }
            }
            
            // Restore to the earlier state
            assert_leq(orig_nedits, tmp_edits.size());
            if(orig_nedits < tmp_edits.size()) tmp_edits.resize(orig_nedits);
        }
        return 0;
    }
}


/**
 * For alignment involving indel, move the indels
 * to the left most possible position
 */
template <typename index_t>
void GenomeHit<index_t>::leftAlign(const Read& rd)
{
    ASSERT_ONLY(const index_t rdlen = (index_t)rd.length());
    const BTDnaString& seq = _fw ? rd.patFw : rd.patRc;
    for(index_t ei = 0; ei < _edits->size(); ei++) {
        Edit& edit = (*_edits)[ei];
        if(edit.type != EDIT_TYPE_READ_GAP && edit.type != EDIT_TYPE_REF_GAP)
            continue;
        if(edit.snpID != (index_t)INDEX_MAX)
            continue;
        index_t ei2 = ei + 1;
        for(; ei2 < _edits->size(); ei2++) {
            const Edit& edit2 = (*_edits)[ei2];
            if(edit2.type != edit.type) break;
            if(edit.type == EDIT_TYPE_READ_GAP) {
                if(edit.pos != edit2.pos) break;
            } else {
                assert_eq(edit.type, EDIT_TYPE_REF_GAP);
                if(edit.pos + ei2 - ei != edit2.pos) break;
            }
        }
        assert_gt(ei2, 0);
        ei2 -= 1;
        Edit& edit2 = (*_edits)[ei2];
        int b = 0;
        if(ei > 0) {
            const Edit& prev_edit = (*_edits)[ei - 1];
            b = prev_edit.pos;
        }
        int l = edit.pos - 1;
        while(l > b) {
            assert_lt(l, (int)rdlen);
            int rdc = seq[_rdoff + l];
            assert_range(0, 4, rdc);
            char rfc = (edit.type == EDIT_TYPE_READ_GAP ? edit2.chr : edit2.qchr);
            if(rfc != "ACGTN"[rdc]) break;
            for(int ei3 = ei2; ei3 > (int)ei; ei3--) {
                if(edit.type == EDIT_TYPE_READ_GAP) {
                    (*_edits)[ei3].chr = (*_edits)[ei3 - 1].chr;
                } else {
                    (*_edits)[ei3].qchr = (*_edits)[ei3 - 1].qchr;
                }
                (*_edits)[ei3].pos -= 1;
            }
            rdc = seq[_rdoff + l];
            assert_range(0, 4, rdc);
            if(edit.type == EDIT_TYPE_READ_GAP) {
                edit.chr = "ACGTN"[rdc];
            } else {
                edit.qchr = "ACGTN"[rdc];
            }
            edit.pos -= 1;
            l--;
        }
        ei = ei2;
    }
}

#ifndef NDEBUG
/**
 * Check that hit is sane w/r/t read.
 */
template <typename index_t>
bool GenomeHit<index_t>::repOk(const Read& rd, const BitPairReference& ref)
{
    if(_len <= 0) return true;
    assert(_sharedVars != NULL);
    SStringExpandable<char>& raw_refbuf = _sharedVars->raw_refbuf;
    SStringExpandable<uint32_t>& destU32 = _sharedVars->destU32;
    
    BTDnaString& editstr = _sharedVars->editstr;
    BTDnaString& partialseq = _sharedVars->partialseq;
    BTDnaString& refstr = _sharedVars->refstr;
    EList<index_t>& reflens = _sharedVars->reflens;
    EList<index_t>& refoffs = _sharedVars->refoffs;
    
    editstr.clear(); partialseq.clear(); refstr.clear();
    reflens.clear(); refoffs.clear();
    
    const BTDnaString& seq = _fw ? rd.patFw : rd.patRc;
    partialseq.install(seq.buf() + this->_rdoff, (size_t)this->_len);
    Edit::toRef(partialseq, *_edits, editstr);
    
    index_t refallen = 0;
    int64_t reflen = 0;
    int64_t refoff = this->_toff;
    refoffs.push_back((index_t)refoff);
    size_t eidx = 0;
    for(size_t i = 0; i < _len; i++, reflen++, refoff++) {
        while(eidx < _edits->size() && (*_edits)[eidx].pos == i) {
            const Edit& edit = (*_edits)[eidx];
            if(edit.isReadGap()) {
                reflen++;
                refoff++;
            } else if(edit.isRefGap()) {
                reflen--;
                refoff--;
            }
            if(edit.isSpliced()) {
                assert_gt(reflen, 0);
                refallen += reflen;
                reflens.push_back((index_t)reflen);
                reflen = 0;
                refoff += edit.splLen;
                assert_gt(refoff, 0);
                refoffs.push_back((index_t)refoff);
            }
            eidx++;
        }
    }
    assert_gt(reflen, 0);
    refallen += (index_t)reflen;
    reflens.push_back((index_t)reflen);
    assert_gt(reflens.size(), 0);
    assert_gt(refoffs.size(), 0);
    assert_eq(reflens.size(), refoffs.size());
    refstr.clear();
    for(index_t i = 0; i < reflens.size(); i++) {
        assert_gt(reflens[i], 0);
        if(i > 0) {
            assert_gt(refoffs[i], refoffs[i-1]);
        }
        raw_refbuf.resize(reflens[i] + 16);
        raw_refbuf.clear();
        int off = ref.getStretch(
                                 reinterpret_cast<uint32_t*>(raw_refbuf.wbuf()),
                                 (size_t)this->_tidx,
                                 (size_t)max<TRefOff>(refoffs[i], 0),
                                 reflens[i],
                                 destU32);
        assert_leq(off, 16);
        for(index_t j = 0; j < reflens[i]; j++) {
            char rfc = *(raw_refbuf.buf()+off+j);
            refstr.append(rfc);
        }
    }
    if(refstr != editstr) {
        cerr << "Decoded nucleotides and edits don't match reference:" << endl;
        //cerr << "           score: " << score.score()
        //<< " (" << gaps << " gaps)" << endl;
        cerr << "           edits: ";
        Edit::print(cerr, *_edits);
        cerr << endl;
        cerr << "    decoded nucs: " << partialseq << endl;
        cerr << "     edited nucs: " << editstr << endl;
        cerr << "  reference nucs: " << refstr << endl;
        assert(0);
    }

    return true;
}
#endif

/**
 * Calculate alignment score
 */
template <typename index_t>
int64_t GenomeHit<index_t>::calculateScore(
                                           const Read&             rd,
                                           SpliceSiteDB&           ssdb,
                                           const Scoring&          sc,
                                           index_t                 minK_local,
                                           index_t                 minIntronLen,
                                           index_t                 maxIntronLen,
                                           index_t                 minAnchorLen,
                                           index_t                 minAnchorLen_noncan,
                                           const BitPairReference& ref)
{
    int64_t score = 0;
    double splicescore = 0;
    int64_t localscore = 0;
    index_t numsplices = 0;
    index_t mm = 0;
    const BTDnaString& seq = _fw ? rd.patFw : rd.patRc;
    const BTString& qual = _fw ? rd.qual : rd.qualRev;
    index_t rdlen = (index_t)seq.length();
    int64_t toff_base = _toff;
    bool conflict_splicesites = false;
    uint8_t whichsense = SPL_UNKNOWN;
    for(index_t i = 0; i < _edits->size(); i++) {
        const Edit& edit = (*_edits)[i];
        assert_lt(edit.pos, _len);
        if(edit.type == EDIT_TYPE_MM) {
            if(edit.snpID == std::numeric_limits<uint32_t>::max()) {
                int pen = sc.score(
                                   dna2col[edit.qchr] - '0',
                                   asc2dnamask[edit.chr],
                                   qual[this->_rdoff + edit.pos] - 33);
                score += pen;
                mm++;
            }
        } else if(edit.type == EDIT_TYPE_SPL) {
            // int left = toff_base + edit.pos - 1;
            // assert_geq(left, 0);
            // int right = left + edit.splLen + 1;
            // assert_geq(right, 0);
            if(!edit.knownSpl) {
                int left_anchor_len = _rdoff + edit.pos;
                assert_gt(left_anchor_len, 0);
                assert_lt(left_anchor_len, (int)rdlen);
                int right_anchor_len = rdlen - left_anchor_len;
                index_t mm2 = 0;
                for(index_t j = i + 1; j < _edits->size(); j++) {
                    const Edit& edit2 = (*_edits)[j];
                    if(edit2.type == EDIT_TYPE_MM ||
                       edit2.type == EDIT_TYPE_READ_GAP ||
                       edit2.type == EDIT_TYPE_REF_GAP) mm2++;
                }
                left_anchor_len -= (mm * 2);
                right_anchor_len -= (mm2 * 2);
                int shorter_anchor_len = min<int>(left_anchor_len, right_anchor_len);
                if(shorter_anchor_len <= 0) shorter_anchor_len = 1;
                assert_gt(shorter_anchor_len, 0);
                uint32_t intronLen_thresh = ((edit.splDir == SPL_FW || edit.splDir == SPL_RC) ?
                                             MaxIntronLen(shorter_anchor_len, minAnchorLen) :
                                             MaxIntronLen_noncan(shorter_anchor_len, minAnchorLen_noncan));
                if(intronLen_thresh < maxIntronLen) {
                    if(edit.splLen > intronLen_thresh) {
                        score += MIN_I32;
                    }
                    
                    if(edit.splDir == SPL_FW || edit.splDir == SPL_RC) {
                        float probscore = ssdb.probscore(edit.donor_seq, edit.acceptor_seq);
                        
                        float probscore_thresh = 0.8f;
                        if(edit.splLen >> 16) probscore_thresh = 0.99f;
                        else if(edit.splLen >> 15) probscore_thresh = 0.97f;
                        else if(edit.splLen >> 14) probscore_thresh = 0.94f;
                        else if(edit.splLen >> 13) probscore_thresh = 0.91f;
                        else if(edit.splLen >> 12) probscore_thresh = 0.88f;
                        if(probscore < probscore_thresh) score += MIN_I32;
                    }
                    if(shorter_anchor_len == left_anchor_len) {
                        if(_trim5 > 0) score += MIN_I32;
                        for(int j = (int)i - 1; j >= 0; j--) {
                            if((*_edits)[j].type == EDIT_TYPE_MM ||
                               (*_edits)[j].type == EDIT_TYPE_READ_GAP ||
                               (*_edits)[j].type == EDIT_TYPE_REF_GAP)
                                score += MIN_I32;
                        }
                    } else {
                        if(_trim3 > 0) score += MIN_I32;
                        for(index_t j = i + 1; j < _edits->size(); j++) {
                            if((*_edits)[j].type == EDIT_TYPE_MM ||
                               (*_edits)[j].type == EDIT_TYPE_READ_GAP ||
                               (*_edits)[j].type == EDIT_TYPE_REF_GAP)
                                score += MIN_I32;
                        }
                    }
                }

                if(edit.snpID == std::numeric_limits<uint32_t>::max()) {
                    if(edit.splDir == SPL_FW || edit.splDir == SPL_RC) {
                        score -= sc.canSpl((int)edit.splLen);
                    } else {
                        score -= sc.noncanSpl((int)edit.splLen);
                    }
                }
                
                // daehwan - for debugging purposes
                if(shorter_anchor_len <= 15) {
                    numsplices += 1;
                    splicescore += (double)edit.splLen;
                }
            }
            
            if(!conflict_splicesites) {
                if(whichsense == SPL_UNKNOWN) {
                    whichsense = edit.splDir;
                } else if(edit.splDir != SPL_UNKNOWN) {
                    assert_neq(whichsense, SPL_UNKNOWN);
                    if(edit.splDir == SPL_FW || edit.splDir == SPL_SEMI_FW) {
                        if(whichsense != SPL_FW && whichsense != SPL_SEMI_FW) {
                            conflict_splicesites = true;
                        }
                    }
                    if(edit.splDir == SPL_RC || edit.splDir == SPL_SEMI_RC) {
                        if(whichsense != SPL_RC && whichsense != SPL_SEMI_RC) {
                            conflict_splicesites = true;
                        }
                    }
                }
            }
            
            toff_base += edit.splLen;
        } else if(edit.type == EDIT_TYPE_READ_GAP) {
            bool open = true;
            if(i > 0 &&
               (*_edits)[i-1].type == EDIT_TYPE_READ_GAP &&
               (*_edits)[i-1].pos == edit.pos) {
                open = false;
            }
            if(edit.snpID == std::numeric_limits<uint32_t>::max()) {
                if(open)    score -= sc.readGapOpen();
                else        score -= sc.readGapExtend();
            }
            toff_base++;
        } else if(edit.type == EDIT_TYPE_REF_GAP) {
            bool open = true;
            if(i > 0 &&
               (*_edits)[i-1].type == EDIT_TYPE_REF_GAP &&
               (*_edits)[i-1].pos + 1 == edit.pos) {
                open = false;
            }
            if(edit.snpID == std::numeric_limits<uint32_t>::max()) {
                if(open)    score -= sc.refGapOpen();
                else        score -= sc.refGapExtend();
            }
            toff_base--;
        }
#ifndef NDEBUG
        else {
            assert(false);
        }
#endif
    }
    
    // Penalty for soft-clipping
    for(index_t i = 0; i < _trim5; i++) {
        score -= sc.sc(qual[i]);
    }
    
    for(index_t i = 0; i < _trim3; i++) {
        score -= sc.sc(qual[i]);
    }
    
    if(conflict_splicesites) {
        score -= sc.conflictSpl();
    }
    
    if (numsplices > 1) splicescore /= (double)numsplices;
    score += (_len - mm) * sc.match();
    _score = score;
    _splicescore = splicescore;
    _localscore = localscore;
    
    return score;
}

/**
 * Encapsulates counters that measure how much work has been done by
 * hierarchical indexing
 */
struct HIMetrics {
    
	HIMetrics() : mutex_m() {
	    reset();
	}
    
	void reset() {
		anchoratts = 0;
        localatts = 0;
        localindexatts = 0;
        localextatts = 0;
        localsearchrecur = 0;
        globalgenomecoords = 0;
        localgenomecoords = 0;
	}
	
	void init(
              uint64_t localatts_,
              uint64_t anchoratts_,
              uint64_t localindexatts_,
              uint64_t localextatts_,
              uint64_t localsearchrecur_,
              uint64_t globalgenomecoords_,
              uint64_t localgenomecoords_)
	{
        localatts = localatts_;
        anchoratts = anchoratts_;
        localindexatts = localindexatts_;
        localextatts = localextatts_;
        localsearchrecur = localsearchrecur_;
        globalgenomecoords = globalgenomecoords_;
        localgenomecoords = localgenomecoords_;
    }
	
	/**
	 * Merge (add) the counters in the given HIMetrics object into this
	 * object.  This is the only safe way to update a HIMetrics shared
	 * by multiple threads.
	 */
	void merge(const HIMetrics& r, bool getLock = false) {
        ThreadSafe ts(&mutex_m, getLock);
        localatts += r.localatts;
        anchoratts += r.anchoratts;
        localindexatts += r.localindexatts;
        localextatts += r.localextatts;
        localsearchrecur += r.localsearchrecur;
        globalgenomecoords += r.globalgenomecoords;
        localgenomecoords += r.localgenomecoords;
    }
	   
    uint64_t localatts;      // # attempts of local search
    uint64_t anchoratts;     // # attempts of anchor search
    uint64_t localindexatts; // # attempts of local index search
    uint64_t localextatts;   // # attempts of extension search
    uint64_t localsearchrecur;
    uint64_t globalgenomecoords;
    uint64_t localgenomecoords;
	
	MUTEX_T mutex_m;
};

/**
 * With a hierarchical indexing, SplicedAligner provides several alignment strategies
 * , which enable effective alignment of RNA-seq reads
 */
template <typename index_t, typename local_index_t>
class HI_Aligner {

public:
	
	/**
	 * Initialize with index.
	 */
	HI_Aligner(
               const GFM<index_t>& gfm,
               bool anchorStop = true,
               uint64_t threads_rids_mindist = 0) :
    _anchorStop(anchorStop),
    _gwstate(GW_CAT),
    _gwstate_local(GW_CAT),
    _thread_rids_mindist(threads_rids_mindist)
    {
        index_t genomeLen = gfm.gh().len();
        _minK = 0;
        while(genomeLen > 0) {
            genomeLen >>= 2;
            _minK++;
        }
        _minK_local = 8;
    }
    
    HI_Aligner() {
    }
    
    /**
     */
    void initRead(Read *rd, bool nofw, bool norc, TAlScore minsc, TAlScore maxpen, bool rightendonly = false) {
        assert(rd != NULL);
        _rds[0] = rd;
        _rds[1] = NULL;
		_paired = false;
        _rightendonly = rightendonly;
        _nofw[0] = nofw;
        _nofw[1] = true;
        _norc[0] = norc;
        _norc[1] = true;
        _minsc[0] = minsc;
        _minsc[1] = INDEX_MAX;
        _maxpen[0] = maxpen;
        _maxpen[1] = INDEX_MAX;
        for(size_t fwi = 0; fwi < 2; fwi++) {
            bool fw = (fwi == 0);
            _hits[0][fwi].init(fw, (index_t)_rds[0]->length());
        }
        _genomeHits.clear();
        _genomeHits_rep[0].clear();
        _hits_searched[0].clear();
        assert(!_paired);
    }
    
    /**
     */
    void initReads(Read *rds[2], bool nofw[2], bool norc[2], TAlScore minsc[2], TAlScore maxpen[2]) {
        assert(rds[0] != NULL && rds[1] != NULL);
		_paired = true;
        _rightendonly = false;
        for(size_t rdi = 0; rdi < 2; rdi++) {
            _rds[rdi] = rds[rdi];
            _nofw[rdi] = nofw[rdi];
            _norc[rdi] = norc[rdi];
            _minsc[rdi] = minsc[rdi];
            _maxpen[rdi] = maxpen[rdi];
            for(size_t fwi = 0; fwi < 2; fwi++) {
                bool fw = (fwi == 0);
		        _hits[rdi][fwi].init(fw, (index_t)_rds[rdi]->length());
            }
            _hits_searched[rdi].clear();
        }
        _genomeHits.clear();
        _genomeHits_rep[0].clear();
        _genomeHits_rep[1].clear();
        _concordantIdxInspected.first = _concordantIdxInspected.second = 0;
        assert(_paired);
        assert(!_rightendonly);
    }
    
    /**
     * Aligns a read or a pair
     * This funcion is called per read or pair
     */
    virtual
    int go(
           const Scoring&             sc,
           const PairedEndPolicy&     pepol, // paired-end policy
           const TranscriptomePolicy& tpol,
           const GraphPolicy&         gpol,
           const GFM<index_t>&        gfm,
           const GFM<index_t>*        rgfm,
           const ALTDB<index_t>&      altdb,
           const RepeatDB<index_t>&   repeatdb,
           const ALTDB<index_t>&      raltdb,
           const BitPairReference&    ref,
           const BitPairReference*    rref,
           SwAligner&                 swa,
           SpliceSiteDB&              ssdb,
           WalkMetrics&               wlm,
           PerReadMetrics&            prm,
           SwMetrics&                 swm,
           HIMetrics&                 him,
           RandomSource&              rnd,
           AlnSinkWrap<index_t>&      sink)
    {
        const ReportingParams& rp = sink.reportingParams();
        
        index_t rdi;
        bool fw;
        bool found[2][2] = {{true, true}, {this->_paired, this->_paired}};
        // given read and its reverse complement
        //  (and mate and the reverse complement of mate in case of pair alignment),
        // pick up one with best partial alignment
        while(nextBWT(sc, pepol, tpol, gpol, gfm, altdb, ref, rdi, fw, wlm, prm, him, rnd, sink)) {
            // given the partial alignment, try to extend it to full alignments
            index_t fwi = (fw == true ? 0 : 1);
            found[rdi][fwi] = align(sc, pepol, tpol, gpol, gfm, altdb, repeatdb, ref, swa, ssdb, rdi, fw, wlm, prm, swm, him, rnd, sink);
            if(!found[0][0] && !found[0][1] && !found[1][0] && !found[1][1]) {
                break;
            }
            
            // try to combine this alignment with some of mate alignments
            // to produce pair alignment
            if(this->_paired) {
                pairReads(sc, pepol, tpol, gpol, gfm, altdb, repeatdb, ref, wlm, prm, him, rnd, sink);
                // if(sink.bestPair() >= _minsc[0] + _minsc[1]) break;
            }
        }
        
        // if no concordant pair is found, try to use alignment of one-end
        // as an anchor to align the other-end
        if(this->_paired) {
            if(sink.numPair() == 0 &&
               (sink.bestUnp1() >= _minsc[0] || sink.bestUnp2() >= _minsc[1])) {
                bool mate_found = false;
                const EList<AlnRes> *rs[2] = {NULL, NULL};
                sink.getUnp1(rs[0]); assert(rs[0] != NULL);
                sink.getUnp2(rs[1]); assert(rs[1] != NULL);
                index_t rs_size[2] = {(index_t)rs[0]->size(), (index_t)rs[1]->size()};
                for(index_t i = 0; i < 2; i++) {
                    for(index_t j = 0; j < rs_size[i]; j++) {
                        const AlnRes& res = (*rs[i])[j];
                        bool fw = (res.orient() == 1);
                        mate_found |= alignMate(
                                                sc,
                                                pepol,
                                                tpol,
                                                gpol,
                                                gfm,
                                                altdb,
                                                repeatdb,
                                                ref,
                                                swa,
                                                ssdb,
                                                i,
                                                fw,
                                                wlm,
                                                prm,
                                                swm,
                                                him,
                                                rnd,
                                                sink,
                                                (index_t)res.refid(),
                                                (index_t)res.refoff());
                    }
                }
                
                if(mate_found) {
                    pairReads(
                              sc,
                              pepol,
                              tpol,
                              gpol,
                              gfm,
                              altdb,
                              repeatdb,
                              ref,
                              wlm,
                              prm,
                              him,
                              rnd,
                              sink);
                }
            }
        }
        
        
        // Determine whether reads map to repetitive sequences
        bool repeat[2][2] = {{false, false}, {false, false}};
        bool perform_repeat_alignment = false;
        
        index_t indexIdx[2] = {0, 0};
#if 1
        if(rgfm != NULL && !((RFM<index_t>*)rgfm)->empty()) {
            // use repeat index to decide whether a read or a pair is from repetitive sequences
            indexIdx[0] = ((RFM<index_t>*)rgfm)->getLocalRFM_idx((*_rds)[0].length());
            if(_paired) {
                indexIdx[1] = ((RFM<index_t>*)rgfm)->getLocalRFM_idx((*_rds)[1].length());
            }
            LocalRFM<index_t>& rfm = ((RFM<index_t>*)rgfm)->getLocalRFM(indexIdx[0]);
            bool skip_repeat[2][2] = {{false, false}, {false, false}};
            if(_paired) {
                const EList<AlnRes> *rs[2] = {NULL, NULL};
                sink.getPair(rs[0], rs[1]);
                assert_eq(rs[0]->size(), rs[1]->size());
                TAlScore bestScore[2][2] = {{_minsc[rdi], _minsc[rdi]}, {_minsc[rdi], _minsc[rdi]}};
                for(size_t r = 0; r < rs[0]->size(); r++) {
                    const AlnRes& rs1 = (*rs[0])[r];
                    const AlnRes& rs2 = (*rs[1])[r];
                    TAlScore score = rs1.score().score() + rs2.score().score();
                    int fwi[2] = {rs1.fw() ? 0 : 1, rs2.fw() ? 0 : 1};
                    if(score > bestScore[fwi[0]][fwi[1]]) {
                        bestScore[fwi[0]][fwi[1]] = score;
                    }
                }
                for(index_t fwi = 0; fwi < 2; fwi++) {
                    for(index_t fwi2 = 0; fwi2 < 2; fwi2++) {
                        if(bestScore[fwi][fwi2] < 0)
                            continue;
                        
                        ReadBWTHit<index_t>& hit = _hits[0][fwi];
                        bool unique = false;
                        for(size_t hi = 0; hi < hit.offsetSize(); hi++) {
                            BWTHit<index_t>& partialHit = hit.getPartialHit(hi);
                            if(partialHit.len() >= _minK + 8 && partialHit.size() == 1) {
                                unique = true;
                                break;
                            }
                        }
                        if(!unique)
                            continue;
                        
                        bool unique2 = false;
                        ReadBWTHit<index_t>& hit2 = _hits[1][fwi2];
                        for(size_t hi = 0; hi < hit2.offsetSize(); hi++) {
                            BWTHit<index_t>& partialHit = hit2.getPartialHit(hi);
                            if(partialHit.len() >= _minK + 8 && partialHit.size() == 1) {
                                unique2 = true;
                                break;
                            }
                        }
                        if(!unique2)
                            continue;
                        
                        skip_repeat[0][fwi] = skip_repeat[1][fwi2] = true;
                    }
                }
            } else {
                const EList<AlnRes> *rs = NULL;
                if(rdi == 0) sink.getUnp1(rs);
                else sink.getUnp2(rs);
                
                TAlScore bestScore[2] = {_minsc[rdi], _minsc[rdi]};
                for(index_t r = 0; r < rs->size(); r++) {
                    TAlScore score = (*rs)[r].score().score();
                    if((*rs)[r].fw()) {
                        if(score > bestScore[0]) {
                            bestScore[0] = score;
                        }
                    } else {
                        if(score > bestScore[1]) {
                            bestScore[1] = score;
                        }
                    }
                }
                for(index_t fwi = 0; fwi < 2; fwi++) {
                    if(bestScore[fwi] < 0)
                        continue;
                    
                    ReadBWTHit<index_t>& hit = _hits[rdi][fwi];
                    index_t offsetSize = hit.offsetSize();
                    for(size_t hi = 0; hi < offsetSize; hi++) {
                        BWTHit<index_t>& partialHit = hit.getPartialHit(hi);
                        if(partialHit.len() >= _minK + 8 && partialHit.size() == 1) {
                            skip_repeat[rdi][fwi] = true;
                            break;
                        }
                    }
                    if(skip_repeat[rdi][fwi]) break;
                }
            }
            
            for(size_t rdi = 0; rdi < (_paired ? 2 : 1); rdi++) {
                for(size_t fwi = 0; fwi < 2; fwi++) {
                    if(skip_repeat[rdi][fwi]) continue;
                    bool fw = (fwi == 0);
                    _hits[rdi][fwi].init(fw, (index_t)_rds[rdi]->length());
                }
            }
            
            while(nextBWT(sc, pepol, tpol, gpol, rfm, altdb, ref, rdi, fw, wlm, prm, him, rnd, sink));
            for(size_t rdi = 0; rdi < (_paired ? 2 : 1); rdi++) {
                for(size_t fwi = 0; fwi < 2; fwi++) {
                    if(skip_repeat[rdi][fwi]) continue;
                    ReadBWTHit<index_t>& hit = _hits[rdi][fwi];
                    index_t offsetSize = hit.offsetSize();
                    assert_gt(offsetSize, 0);
                    for(size_t hi = 0; hi < offsetSize; hi++) {
                        BWTHit<index_t>& partialHit = hit.getPartialHit(hi);
                        if(partialHit.len() >= (_minK << 1)) {
                            repeat[rdi][fwi] = true;
                            perform_repeat_alignment = true;
                            break;
                        }
                    }
                }
            }
        }
#else
        // use minimizer to decide whether a read or a pair is from repetitive sequences
        perform_repeat_alignment = false;
        for(size_t rdi = 0; rdi < (_paired ? 2 : 1); rdi++) {
            Read& read = *_rds[rdi];
            for(size_t fwi = 0; fwi < 2; fwi++) {
                const BTDnaString& seq = (fwi == 0 ? read.patFw : read.patRc);
                repeat[rdi][fwi] = repeat_kmertable.isRepeat(seq, _tmp_minimizers);
                perform_repeat_alignment |= repeat[rdi][fwi];
            }
        }
#endif
        
        // Handle alignment to repetitive regions
        if(rgfm != NULL &&
           perform_repeat_alignment) {
            LocalRFM<index_t>& rfm = ((RFM<index_t>*)rgfm)->getLocalRFM(indexIdx[0]);
            RB_KmerTable& repeatKmertable = ((RFM<index_t>*)rgfm)->getKmertable(indexIdx[0]);
            
            _repeatConcordant.clear();
            index_t prev_align_size[2] = {0, 0};
            for(size_t rdi = 0; rdi < (_paired ? 2 : 1); rdi++) {
                const EList<AlnRes> *rs = NULL;
                if(rdi == 0) sink.getUnp1(rs);
                else         sink.getUnp2(rs);
                prev_align_size[rdi] = rs->size();
            }
            
            for(size_t rdi = 0; rdi < (_paired ? 2 : 1); rdi++) {
                for(size_t fwi = 0; fwi < 2; fwi++) {
                    if(!repeat[rdi][fwi]) continue;
                    
                    // choose candidate partial alignments for further alignment
                    index_t maxsize = max<index_t>(rp.khits, rp.kseeds);
      
#if 0
                    ReadBWTHit<index_t>& hit = _hits[rdi][fwi];
                    if(!hit.done()) continue;
                    getAnchorHits(rfm,
                                  pepol,
                                  tpol,
                                  gpol,
                                  altdb,
                                  repeatdb,
                                  *rref,
                                  rnd,
                                  rdi,
                                  fwi == 0, // fw
                                  _genomeHits_rep[rdi],
                                  _genomeHits_rep[rdi].size() + maxsize,
                                  _sharedVars,
                                  wlm,
                                  prm,
                                  him,
                                  true); // repeat?
#else
                    
                    getRepeatHits(rfm,
                                  pepol,
                                  tpol,
                                  gpol,
                                  raltdb,
                                  repeatdb,
                                  repeatKmertable,
                                  *rref,
                                  rnd,
                                  rdi,
                                  fwi == 0, // fw
                                  indexIdx[0],
                                  _genomeHits_rep[rdi],
                                  _genomeHits_rep[rdi].size() + maxsize,
                                  _sharedVars,
                                  sc,
                                  swa,
                                  ssdb,
                                  swm,
                                  wlm,
                                  prm,
                                  him,
                                  sink);
#endif
                }
            }
            
            EList<pair<RepeatCoord<index_t>, RepeatCoord<index_t> > >& positions = _positions;
            for(size_t rdi = 0; rdi < (_paired ? 2 : 1); rdi++) {
                for(size_t i = 0; i < _genomeHits_rep[rdi].size(); i++) {
                    if(_genomeHits_rep[rdi][i].len() < (_minK << 1)) continue;
                    
                    // DK - debugging purposes
#if 0
                    positions.clear();
                    repeatdb.getCoords(_genomeHits_rep[rdi][i]._tidx,
                                       _genomeHits_rep[rdi][i]._joinedOff,
                                       _genomeHits_rep[rdi][i]._joinedOff + _genomeHits_rep[rdi][i].len(),
                                       _snpIDs,
                                       raltdb,
                                       positions,
                                       rp.khits * 1000);
#endif
                    
                    const EList<AlnRes> *rs = NULL;
                    if(rdi == 0) sink.getUnp2(rs);
                    else         sink.getUnp1(rs);
                    assert(rs != NULL);
                    
                    bool candidate_found = false;
                    for(size_t j = 0; j < prev_align_size[1-rdi]; j++) {
                        const AlnRes& res = (*rs)[j];
                        if(res.repeat())
                            continue;
                        
                        TAlScore estScore = res.score().score() + _genomeHits_rep[rdi][i].score();
                        if(sink.bestPair() >= estScore && sink.numBestPair().first > rp.khits)
                            break;
                        
                        positions.clear();
                        index_t joinedOff = 0;
                        gfm.textOffToJoined(res.refid(), res.refoff(), joinedOff);
                        repeatdb.findCoords(joinedOff,
                                            joinedOff + res.refExtent(),
                                            _genomeHits_rep[rdi][i]._tidx,
                                            _genomeHits_rep[rdi][i]._joinedOff,
                                            _genomeHits_rep[rdi][i]._joinedOff + _genomeHits_rep[rdi][i].len(),
                                            _snpIDs,
                                            raltdb,
                                            positions,
                                            rp.khits * 10);
                        if(positions.size() <= 0)
                            continue;
                        
                        for(size_t p = 0; p < positions.size(); p++) {
                            if(positions[p].first.tid != res.refid()) continue;
                            if(positions[p].first.toff + 1000 < res.refoff() ||
                               res.refoff() + 1000 < positions[p].first.toff) continue;
                            if(sink.bestPair() >= estScore && sink.numBestPair().first > rp.khits)
                                break;
                            
                            candidate_found = true;
                            
                            _genomeHits.clear();
                            _genomeHits.expand();
                            _genomeHits.back() = _genomeHits_rep[rdi][i];
                            _genomeHits.back()._tidx = positions[p].first.tid;
                            _genomeHits.back()._toff = positions[p].first.toff;
                            _genomeHits.back()._joinedOff = positions[p].first.joinedOff;
                            if(!positions[p].first.fw) {
                                _genomeHits.back().reverse(*_rds[rdi]);
                            }
                            
                            // extend the partial alignments bidirectionally using
                            // local search, extension, and (less often) global search
                            hybridSearch(sc,
                                         pepol,
                                         tpol,
                                         gpol,
                                         gfm,
                                         altdb,
                                         repeatdb,
                                         ref,
                                         swa,
                                         ssdb,
                                         rdi,
                                         _genomeHits.back()._fw,
                                         wlm,
                                         prm,
                                         swm,
                                         him,
                                         rnd,
                                         sink);
                        }
                        
                        if(candidate_found) {
                            pairReads(
                                      sc,
                                      pepol,
                                      tpol,
                                      gpol,
                                      gfm,
                                      altdb,
                                      repeatdb,
                                      ref,
                                      wlm,
                                      prm,
                                      him,
                                      rnd,
                                      sink);
                        }
                    }
                    
                    if(rdi == 0 && _paired) {
                        for(size_t j = 0; j < _genomeHits_rep[1].size(); j++) {
                            if(_genomeHits_rep[1][j].len() < (_minK << 1)) continue;
                            
                            TAlScore estScore = _genomeHits_rep[0][i].score() + _genomeHits_rep[1][j].score();
                            // if(sink.bestPair() >= estScore && sink.numBestPair().first > rp.khits)
                            //    break;

                            positions.clear();
                            repeatdb.findCommonCoords(_genomeHits_rep[0][i]._tidx,
                                                      _genomeHits_rep[0][i]._joinedOff,
                                                      _genomeHits_rep[0][i]._joinedOff + _genomeHits_rep[0][i].len(),
                                                      _snpIDs,
                                                      _genomeHits_rep[1][j]._tidx,
                                                      _genomeHits_rep[1][j]._joinedOff,
                                                      _genomeHits_rep[1][j]._joinedOff + _genomeHits_rep[1][j].len(),
                                                      _snpIDs2,
                                                      raltdb,
                                                      positions,
                                                      rp.khits * 10);
                            if(positions.size() <= 0) continue;
                            
                            _repeatConcordant.expand();
                            _repeatConcordant.back().first = _genomeHits_rep[0][i]._joinedOff;
                            _repeatConcordant.back().second = _genomeHits_rep[1][j]._joinedOff;
                            
                            for(size_t p = 0; p < positions.size(); p++) {
                                if(sink.bestPair() >= estScore && sink.numBestPair().first > rp.khits)
                                    break;
                                
                                _genomeHits.clear();
                                _genomeHits.expand();
                                _genomeHits.back() = _genomeHits_rep[0][i];
                                _genomeHits.back()._tidx = positions[p].first.tid;
                                _genomeHits.back()._toff = positions[p].first.toff;
                                _genomeHits.back()._joinedOff = positions[p].first.joinedOff;
                                if(!positions[p].first.fw) {
                                    _genomeHits.back().reverse(*_rds[0]);
                                }
                                
                                // extend the partial alignments bidirectionally using
                                // local search, extension, and (less often) global search
                                hybridSearch(sc,
                                             pepol,
                                             tpol,
                                             gpol,
                                             gfm,
                                             altdb,
                                             repeatdb,
                                             ref,
                                             swa,
                                             ssdb,
                                             0,
                                             _genomeHits.back()._fw,
                                             wlm,
                                             prm,
                                             swm,
                                             him,
                                             rnd,
                                             sink);
                                
                                _genomeHits.clear();
                                _genomeHits.expand();
                                _genomeHits.back() = _genomeHits_rep[1][j];
                                _genomeHits.back()._tidx = positions[p].second.tid;
                                _genomeHits.back()._toff = positions[p].second.toff;
                                _genomeHits.back()._joinedOff = positions[p].second.joinedOff;
                                if(!positions[p].second.fw) {
                                    _genomeHits.back().reverse(*_rds[1]);
                                }
                                
                                // extend the partial alignments bidirectionally using
                                // local search, extension, and (less often) global search
                                hybridSearch(sc,
                                             pepol,
                                             tpol,
                                             gpol,
                                             gfm,
                                             altdb,
                                             repeatdb,
                                             ref,
                                             swa,
                                             ssdb,
                                             1,
                                             _genomeHits.back()._fw,
                                             wlm,
                                             prm,
                                             swm,
                                             him,
                                             rnd,
                                             sink);
                            }
                            
                            if(positions.size() > 0) {
                                pairReads(
                                          sc,
                                          pepol,
                                          tpol,
                                          gpol,
                                          gfm,
                                          altdb,
                                          repeatdb,
                                          ref,
                                          wlm,
                                          prm,
                                          him,
                                          rnd,
                                          sink);
                            }
                        }
                    } // if(rdi == 0)
                } // for(size_t i = 0; i < _genomeHits_rep[rdi].size()
                
                bool align2repeat = false;
                if(_paired) {
                    index_t numBestPair = sink.numBestPair().first;
                    align2repeat = (numBestPair == 0 || numBestPair > rp.khits);
                } else {
                    const EList<AlnRes> *rs = NULL;
                    if(rdi == 0) sink.getUnp1(rs);
                    else         sink.getUnp2(rs);
                    assert(rs != NULL);
                    align2repeat = (rs->size() == 0 || sink.numBestUnp(rdi).first > rp.khits);
                }
                
                if(align2repeat) {
                    for(size_t i = 0; i < _genomeHits_rep[rdi].size(); i++) {
                        _genomeHits.clear();
                        _genomeHits.expand();
                        _genomeHits.back() = _genomeHits_rep[rdi][i];
                        _genomeHits.back()._repeat = true;
                        hybridSearch(sc,
                                     pepol,
                                     tpol,
                                     gpol,
                                     rfm,
                                     altdb,
                                     repeatdb,
                                     *rref,
                                     swa,
                                     ssdb,
                                     rdi,
                                     _genomeHits.back()._fw,
                                     wlm,
                                     prm,
                                     swm,
                                     him,
                                     rnd,
                                     sink);
                    }
                    
                    if(_paired && rdi == 1) {
                        if(sink.numBestUnp(rdi).first > rp.khits) {
                            pairReads(
                                      sc,
                                      pepol,
                                      tpol,
                                      gpol,
                                      gfm,
                                      altdb,
                                      repeatdb,
                                      ref,
                                      wlm,
                                      prm,
                                      him,
                                      rnd,
                                      sink);
                        }
                    }
                }
                
            } // for(size_t rdi = 0
        } // repeat
        
        return EXTEND_POLICY_FULFILLED;
    }
    
    /**
     * Given a read or its reverse complement (or mate),
     * align the unmapped portion using the global FM index
     */
    virtual
    bool nextBWT(
                 const Scoring&             sc,
                 const PairedEndPolicy&     pepol, // paired-end policy
                 const TranscriptomePolicy& tpol,
                 const GraphPolicy&         gpol,
                 const GFM<index_t>&        gfm,
                 const ALTDB<index_t>&      altdb,
                 const BitPairReference&    ref,
                 index_t&                   rdi,
                 bool&                      fw,
                 WalkMetrics&               wlm,
                 PerReadMetrics&            prm,
                 HIMetrics&                 him,
                 RandomSource&              rnd,
                 AlnSinkWrap<index_t>&      sink)
    {
        const ReportingParams& rp = sink.reportingParams();
        
        // Pick up a candidate from a read or its reverse complement
        // (for pair, also consider mate and its reverse complement)
        while(pickNextReadToSearch(rdi, fw)) {
            size_t mineFw = 0, mineRc = 0;
            index_t fwi = (fw ? 0 : 1);
            ReadBWTHit<index_t>& hit = _hits[rdi][fwi];
            assert(!hit.done());
            bool pseudogeneStop = gfm.gh().linearFM() && !tpol.no_spliced_alignment();
            bool anchorStop = _anchorStop && !gfm.repeat();
            if(!rp.secondary) {
                index_t numSearched = hit.numActualPartialSearch();
                int64_t bestScore = 0;
                if(rdi == 0) {
                    bestScore = sink.bestUnp1();
                    if(bestScore >= _minsc[rdi]) {
                        // do not further align this candidate
                        // unless it may be at least as good as the alignment of its reverse complement
                        index_t maxmm = (index_t)((-bestScore + sc.mmpMax - 1) / sc.mmpMax);
                        if(numSearched > maxmm + sink.bestSplicedUnp1() + 1) {
                            hit.done(true);
                            if(_paired) {
                                if(sink.bestUnp2() >= _minsc[1-rdi] &&
                                   sink.numPair() > 0) return false;
                                else continue;
                            } else {
                                return false;
                            }
                        }
                    }
                } else {
                    assert(_paired);
                    assert_eq(rdi, 1);
                    bestScore = sink.bestUnp2();
                    if(bestScore >= _minsc[rdi]) {
                        // Do not further extend this alignment
                        // unless it may be at least as good as the previous alignemnt
                        index_t maxmm = (index_t)((-bestScore + sc.mmpMax - 1) / sc.mmpMax);
                        if(numSearched > maxmm + sink.bestSplicedUnp2() + 1) {
                            hit.done(true);
                            if(_paired) {
                                if(sink.bestUnp1() >= _minsc[1-rdi] &&
                                   sink.numPair() > 0) return false;
                                else continue;
                            } else {
                                return false;
                            }
                        }
                    }
                }
                
                ReadBWTHit<index_t>& rchit = _hits[rdi][1-fwi];
                if(rchit.done() && bestScore < _minsc[rdi]) {
                    if(numSearched > rchit.numActualPartialSearch() + (anchorStop ? 1 : 0)) {
                        hit.done(true);
                        return false;
                    }
                }
            }
            
            // Align this read beginning from previously stopped base
            // stops when it is uniquelly mapped with at least 28bp or
            // it may involve processed pseudogene
            partialSearch(
                          gfm,
                          *_rds[rdi],
                          sc,
                          sink.reportingParams(),
                          fw,
                          0,
                          mineFw,
                          mineRc,
                          hit,
                          rnd,
                          pseudogeneStop,
                          anchorStop);
            
            assert(hit.repOk());
            if(hit.done()) return true;
            // Advance hit._cur by 1
            if(!pseudogeneStop) {
                if(hit._cur + 1 < hit._len) hit._cur++;
            }
            if(anchorStop) {
                hit.done(true);
                return true;
            }
            // hit.adjustOffset(_minK);
        }
        
        return false;
    }
    
    /**
     * Given partial alignments of a read, try to further extend
     * the alignment bidirectionally
     */
    virtual
    bool align(
               const Scoring&                   sc,
               const PairedEndPolicy&           pepol, // paired-end policy
               const TranscriptomePolicy&       tpol,
               const GraphPolicy&               gpol,
               const GFM<index_t>&              gfm,
               const ALTDB<index_t>&            altdb,
               const RepeatDB<index_t>&         repeatdb,
               const BitPairReference&          ref,
               SwAligner&                       swa,
               SpliceSiteDB&                    ssdb,
               index_t                          rdi,
               bool                             fw,
               WalkMetrics&                     wlm,
               PerReadMetrics&                  prm,
               SwMetrics&                       swm,
               HIMetrics&                       him,
               RandomSource&                    rnd,
               AlnSinkWrap<index_t>&            sink);
    
    /**
     * Given the alignment of its mate as an anchor,
     * align the read
     */
    virtual
    bool alignMate(
                   const Scoring&                   sc,
                   const PairedEndPolicy&           pepol, // paired-end policy
                   const TranscriptomePolicy&       tpol,
                   const GraphPolicy&               gpol,
                   const GFM<index_t>&              gfm,
                   const ALTDB<index_t>&            altdb,
                   const RepeatDB<index_t>&         repeatdb,
                   const BitPairReference&          ref,
                   SwAligner&                       swa,
                   SpliceSiteDB&                    ssdb,
                   index_t                          rdi,
                   bool                             fw,
                   WalkMetrics&                     wlm,
                   PerReadMetrics&                  prm,
                   SwMetrics&                       swm,
                   HIMetrics&                       him,
                   RandomSource&                    rnd,
                   AlnSinkWrap<index_t>&            sink,
                   index_t                          tidx,
                   index_t                          toff);
    
    /**
     * Given a partial alignment of a read, try to further extend
     * the alignment bidirectionally using a combination of
     * local search, extension, and global search
     */
    virtual
    void hybridSearch(
                      const Scoring&                     sc,
                      const PairedEndPolicy&             pepol, // paired-end policy
                      const TranscriptomePolicy&         tpol,
                      const GraphPolicy&                 gpol,
                      const GFM<index_t>&                gfm,
                      const ALTDB<index_t>&              altdb,
                      const RepeatDB<index_t>&           repeatdb,
                      const BitPairReference&            ref,
                      SwAligner&                         swa,
                      SpliceSiteDB&                      ssdb,
                      index_t                            rdi,
                      bool                               fw,
                      WalkMetrics&                       wlm,
                      PerReadMetrics&                    prm,
                      SwMetrics&                         swm,
                      HIMetrics&                         him,
                      RandomSource&                      rnd,
                      AlnSinkWrap<index_t>&              sink)
    {}
    
    /**
     * Given a partial alignment of a read, try to further extend
     * the alignment bidirectionally using a combination of
     * local search, extension, and global search
     */
    virtual
    int64_t hybridSearch_recur(
                               const Scoring&                   sc,
                               const PairedEndPolicy&           pepol, // paired-end policy
                               const TranscriptomePolicy&       tpol,
                               const GraphPolicy&               gpol,
                               const GFM<index_t>&              gfm,
                               const ALTDB<index_t>&            altdb,
                               const RepeatDB<index_t>&         repeatdb,
                               const BitPairReference&          ref,
                               SwAligner&                       swa,
                               SpliceSiteDB&                    ssdb,
                               index_t                          rdi,
                               const GenomeHit<index_t>&        hit,
                               index_t                          hitoff,
                               index_t                          hitlen,
                               WalkMetrics&                     wlm,
                               PerReadMetrics&                  prm,
                               SwMetrics&                       swm,
                               HIMetrics&                       him,
                               RandomSource&                    rnd,
                               AlnSinkWrap<index_t>&            sink,
                               bool                             alignMate = false,
                               index_t                          dep = 0)
    { return numeric_limits<int64_t>::min(); }
    
    /**
     * Choose a candidate for alignment from a read or its reverse complement
     * (also from a mate or its reverse complement for pair)
     */
    bool pickNextReadToSearch(index_t& rdi, bool& fw) {
        rdi = 0; fw = true;
        bool picked = false;
        int64_t maxScore = std::numeric_limits<int64_t>::min();
        for(index_t rdi2 = 0; rdi2 < (_paired ? 2 : 1); rdi2++) {
            assert(_rds[rdi2] != NULL);
            for(index_t fwi = 0; fwi < 2; fwi++) {
                if     (fwi == 0 && _nofw[rdi2]) continue;
                else if(fwi == 1 && _norc[rdi2]) continue;

                if(_hits[rdi2][fwi].done()) continue;
                int64_t curScore = _hits[rdi2][fwi].searchScore((index_t)_minK);
                if(_hits[rdi2][fwi].cur() == 0) {
                    curScore = std::numeric_limits<int64_t>::max();
                }
                assert_gt(curScore, std::numeric_limits<int64_t>::min());
                if(curScore > maxScore) {
                    maxScore = curScore;
                    rdi = rdi2;
                    fw = (fwi == 0);
                    picked = true;
                }
            }
        }
        
        return picked;
    }

	/**
     * Align a part of a read without any edits
	 */
    index_t partialSearch(
                          const GFM<index_t>&     gfm,     // GFM index
                          const Read&             read,    // read to align
                          const Scoring&          sc,      // scoring scheme
                          const ReportingParams&  rp,
                          bool                    fw,      // don't align forward read
                          size_t                  mineMax, // don't care about edit bounds > this
                          size_t&                 mineFw,  // minimum # edits for forward read
                          size_t&                 mineRc,  // minimum # edits for revcomp read
                          ReadBWTHit<index_t>&    hit,     // holds all the seed hits (and exact hit)
                          RandomSource&           rnd,
                          bool&                   pseudogeneStop,  // stop if mapped to multiple locations due to processed pseudogenes
                          bool&                   anchorStop,
                          index_t                 maxHitLen = (index_t)INDEX_MAX);
    
    /**
     * Global FM index search
     */
    index_t globalGFMSearch(
                            const GFM<index_t>&    gfm,   // GFM index
                            const Read&            read,  // read to align
                            const Scoring&         sc,    // scoring scheme
                            const ReportingParams& rp,
                            bool                   fw,
                            index_t                hitoff,
                            index_t&               hitlen,
                            index_t&               top,
                            index_t&               bot,
                            index_t&               node_top,
                            index_t&               node_bot,
                            EList<pair<index_t, index_t> >& node_iedge_count,
                            RandomSource&          rnd,
                            bool&                  uniqueStop,
                            index_t                maxHitLen = (index_t)INDEX_MAX);
    
    /**
     * Local FM index search
     */
    index_t localGFMSearch(
                           const LocalGFM<local_index_t, index_t>&  gfm,  // GFM index
                           const Read&                      read,    // read to align
                           const Scoring&                   sc,      // scoring scheme
                           const ReportingParams&           rp,
                           bool                             fw,
                           index_t                          rdoff,
                           index_t&                         hitlen,
                           local_index_t&                   top,
                           local_index_t&                   bot,
                           local_index_t&                   node_top,
                           local_index_t&                   node_bot,
                           EList<pair<local_index_t, local_index_t> >& local_node_iedge_count,
                           RandomSource&                    rnd,
                           bool&                            uniqueStop,
                           local_index_t                    minUniqueLen,
                           local_index_t                    maxHitLen = (local_index_t)INDEX_MAX,
                           local_index_t                    maxHits = 0);
    
    /**
     * Convert FM offsets to the corresponding genomic offset (chromosome id, offset)
     **/
    bool getGenomeCoords(
                         const GFM<index_t>&        gfm,
                         const ALTDB<index_t>&      altdb,
                         const BitPairReference&    ref,
                         RandomSource&              rnd,
                         index_t                    top,
                         index_t                    bot,
                         index_t                    node_top,
                         index_t                    node_bot,
                         const EList<pair<index_t, index_t> >& node_iedge_count,
                         bool                       fw,
                         index_t                    maxelt,
                         index_t                    rdoff,
                         index_t                    rdlen,
                         EList<Coord>&              coords,
                         WalkMetrics&               met,
                         PerReadMetrics&            prm,
                         HIMetrics&                 him,
                         bool                       rejectStraddle,
                         bool&                      straddled);
    
    /**
     * Convert FM offsets to the corresponding genomic offset (chromosome id, offset)
     **/
    bool getGenomeCoords_local(
                               const GFM<local_index_t>&    gfm,
                               const ALTDB<index_t>&        altdb,
                               const BitPairReference&      ref,
                               RandomSource&                rnd,
                               local_index_t                top,
                               local_index_t                bot,
                               local_index_t                node_top,
                               local_index_t                node_bot,
                               const EList<pair<local_index_t, local_index_t> >& node_iedge_count,
                               bool                         fw,
                               index_t                      rdoff,
                               index_t                      rdlen,
                               EList<Coord>&                coords,
                               WalkMetrics&                 met,
                               PerReadMetrics&              prm,
                               HIMetrics&                   him,
                               bool                         rejectStraddle,
                               bool&                        straddled);
    
    /**
     * Given a set of partial alignments for a read,
     * choose some that are longer and mapped to fewer places
     */
    index_t getAnchorHits(
                          const GFM<index_t>&               gfm,
                          const PairedEndPolicy&            pepol, // paired-end policy
                          const TranscriptomePolicy&        tpol,
                          const GraphPolicy&                gpol,
                          const ALTDB<index_t>&             altdb,
                          const RepeatDB<index_t>&          repeatdb,
                          const BitPairReference&           ref,
                          RandomSource&                     rnd,
                          index_t                           rdi,
                          bool                              fw,
                          EList<GenomeHit<index_t> >&       genomeHits,
                          index_t                           maxGenomeHitSize,
                          SharedTempVars<index_t>&          sharedVars,
                          WalkMetrics&                      wlm,
                          PerReadMetrics&                   prm,
                          HIMetrics&                        him,
                          bool                              repeat = false)
    {
        index_t fwi = (fw ? 0 : 1);
        assert_lt(rdi, 2);
        assert(_rds[rdi] != NULL);
        ReadBWTHit<index_t>& hit = _hits[rdi][fwi];
        assert(hit.done());
        index_t offsetSize = hit.offsetSize();
        assert_gt(offsetSize, 0);
        for(size_t hi = 0; hi < offsetSize; hi++) {
            index_t hj = 0;
            for(; hj < offsetSize; hj++) {
                BWTHit<index_t>& partialHit_j = hit.getPartialHit(hj);
                if(partialHit_j.empty() ||
                   partialHit_j.hasGenomeCoords() ||
                   partialHit_j.len() <= _minK + 2) continue;
                else break;
            }
            if(hj >= offsetSize) break;
            for(index_t hk = hj + 1; hk < offsetSize; hk++) {
                BWTHit<index_t>& partialHit_j = hit.getPartialHit(hj);
                BWTHit<index_t>& partialHit_k = hit.getPartialHit(hk);
                if(partialHit_k.empty() ||
                   partialHit_k.hasGenomeCoords() ||
                   partialHit_k.len() <= _minK + 2) continue;
                
                if(partialHit_j._hit_type == partialHit_k._hit_type) {
                    if((partialHit_j.size() > partialHit_k.size()) ||
                       (partialHit_j.size() == partialHit_k.size() && partialHit_j.len() < partialHit_k.len())) {
                        hj = hk;
                    }
                } else {
                    if(partialHit_k._hit_type > partialHit_j._hit_type) {
                        hj = hk;
                    }
                }
            }
            BWTHit<index_t>& partialHit = hit.getPartialHit(hj);
            assert(!partialHit.hasGenomeCoords());
            
            // Retrieve genomic coordinates
            //  If there are too many genomic coordinates to get,
            //  then we randomly choose and retrieve a small set of them
            assert_leq(genomeHits.size(), maxGenomeHitSize);
            index_t remainedGenomeHitSize = maxGenomeHitSize - genomeHits.size();
            if(remainedGenomeHitSize <= 0)
                break;
            index_t expectedNumCoords = partialHit._node_bot - partialHit._node_top;
            bool straddled = false;
            if(expectedNumCoords <= remainedGenomeHitSize) {
                getGenomeCoords(
                                gfm,
                                altdb,
                                ref,
                                rnd,
                                partialHit._top,
                                partialHit._bot,
                                partialHit._node_top,
                                partialHit._node_bot,
                                partialHit._node_iedge_count,
                                fw,
                                partialHit._bot - partialHit._top,
                                hit._len - partialHit._bwoff - partialHit._len,
                                partialHit._len,
                                partialHit._coords,
                                wlm,
                                prm,
                                him,
                                false, // reject straddled
                                straddled);
            } else {
                index_t edgeIdx = 0;
                index_t top = partialHit._top;
                index_t added = 0;
                for(index_t node = partialHit._node_top; node < partialHit._node_bot; node++, expectedNumCoords--) {
                    index_t bot = top + 1;
                    _tmp_node_iedge_count.clear();
                    if(edgeIdx < partialHit._node_iedge_count.size()) {
                        assert_leq(node - partialHit._node_top, partialHit._node_iedge_count[edgeIdx].first);
                        if(node - partialHit._node_top == partialHit._node_iedge_count[edgeIdx].first) {
                            bot += partialHit._node_iedge_count[edgeIdx].second;
                            _tmp_node_iedge_count.expand();
                            _tmp_node_iedge_count.back().first = 0;
                            _tmp_node_iedge_count.back().second = partialHit._node_iedge_count[edgeIdx].second;
                            edgeIdx++;
                        }
                    }
                    assert_lt(added, remainedGenomeHitSize);
                    uint32_t rndi = rnd.nextU32() % expectedNumCoords;
                    if(rndi < remainedGenomeHitSize - added) {
                        getGenomeCoords(
                                        gfm,
                                        altdb,
                                        ref,
                                        rnd,
                                        top,
                                        bot,
                                        node,
                                        node + 1,
                                        _tmp_node_iedge_count,
                                        fw,
                                        partialHit._bot - partialHit._top,
                                        hit._len - partialHit._bwoff - partialHit._len,
                                        partialHit._len,
                                        partialHit._coords,
                                        wlm,
                                        prm,
                                        him,
                                        false, // reject straddled
                                        straddled);
                        added++;
                        if(added >= remainedGenomeHitSize) break;
                    }
                    top = bot;
                }
            }
            
            if(!partialHit.hasGenomeCoords()) continue;
            EList<Coord>& coords = partialHit._coords;
            assert_gt(coords.size(), 0);
            const index_t genomeHit_size = (index_t)genomeHits.size();
            if(genomeHit_size + coords.size() > maxGenomeHitSize) {
                coords.shufflePortion(0, coords.size(), rnd);
            }
            for(index_t k = 0; k < coords.size(); k++) {
                const Coord& coord = coords[k];
                if(coord.ref() == numeric_limits<index_t>::max())
                    continue;
                index_t len = partialHit._len;
                index_t rdoff = hit._len - partialHit._bwoff - len;
                bool overlapped = false;
                for(index_t l = 0; l < genomeHit_size; l++) {
                    GenomeHit<index_t>& genomeHit = genomeHits[l];
                    if(genomeHit.ref() != (index_t)coord.ref() || genomeHit.fw() != coord.fw()) continue;
                    assert_lt(genomeHit.rdoff(), hit._len);
                    assert_lt(rdoff, hit._len);
                    index_t hitoff = genomeHit.refoff() + hit._len - genomeHit.rdoff();
                    index_t hitoff2 = (index_t)coord.off() + hit._len - rdoff;
                    int64_t hitoff_diff = (tpol.no_spliced_alignment() ? 0 : tpol.maxIntronLen());
                    if(abs((int64_t)hitoff - (int64_t)hitoff2) <= hitoff_diff) {
                        overlapped = true;
                        genomeHit._hitcount++;
                        break;
                    }
                }
                if(repeat) {
                    if(!repeatdb.repeatExist(coord.ref(), coord.off(), coord.off() + len)) {
                        continue;
                    }
                }
                if(!overlapped) {
                    GenomeHit<index_t>::adjustWithALT(
                                                      rdoff,
                                                      len,
                                                      coord,
                                                      _sharedVars,
                                                      genomeHits,
                                                      *_rds[rdi],
                                                      gfm,
                                                      altdb,
                                                      ref,
                                                      gpol);
                }
                if(partialHit._hit_type == CANDIDATE_HIT && genomeHits.size() >= maxGenomeHitSize) break;
            }
            if(partialHit._hit_type == CANDIDATE_HIT && genomeHits.size() >= maxGenomeHitSize) break;
        }

        return (index_t)genomeHits.size();
    }
    
    /**
     *
     */
    index_t getRepeatHits(
                          const GFM<index_t>&               gfm,
                          const PairedEndPolicy&            pepol, // paired-end policy
                          const TranscriptomePolicy&        tpol,
                          const GraphPolicy&                gpol,
                          const ALTDB<index_t>&             altdb,
                          const RepeatDB<index_t>&          repeatdb,
                          const RB_KmerTable&               repeat_kmertable,
                          const BitPairReference&           ref,
                          RandomSource&                     rnd,
                          index_t                           rdi,
                          bool                              fw,
                          index_t                           repID,
                          EList<GenomeHit<index_t> >&       genomeHits,
                          index_t                           maxGenomeHitSize,
                          SharedTempVars<index_t>&          sharedVars,
                          const Scoring&                    sc,
                          SwAligner&                        swa,
                          SpliceSiteDB&                     ssdb,
                          SwMetrics&                        swm,
                          WalkMetrics&                      wlm,
                          PerReadMetrics&                   prm,
                          HIMetrics&                        him,
                          AlnSinkWrap<index_t>&             sink)
    {
        assert_lt(rdi, 2);
        assert(_rds[rdi] != NULL);
        Read& rd = *_rds[rdi];
        const BTDnaString& seq = fw ? rd.patFw : rd.patRc;
        repeat_kmertable.findAlignments(seq,
                                        _tmp_minimizers,
                                        _tmp_position2D,
                                        _tmp_alignments);
        
        const TAlScore cushion = sc.mmpMax;
     
        TAlScore bestScore = _minsc[rdi];
        size_t prev_numHits = genomeHits.size();
        for(index_t i = 0; i < _tmp_alignments.size(); i++) {
            const RB_Alignment& coord = _tmp_alignments[i];
            index_t len = seq.length();
            index_t rdoff = 0;
            if(!repeatdb.repeatExist(repID, coord.pos, coord.pos + len)) {
                continue;
            }
            
            genomeHits.expand();
            GenomeHit<index_t>& genomeHit = genomeHits.back();
            genomeHit.init(fw,
                           rdoff,
                           0,
                           0, // trim5
                           0, // trim3
                           repID, // ref,
                           coord.pos,
                           coord.pos,
                           this->_sharedVars);
            
            index_t maxmm = (index_t)(-bestScore / sc.mmpMax);
            index_t leftext = 0, rightext = len;
            genomeHit.extend(rd,
                             gfm,
                             ref,
                             altdb,
                             repeatdb,
                             ssdb,
                             swa,
                             swm,
                             prm,
                             sc,
                             this->_minsc[rdi],
                             rnd,
                             (index_t)this->_minK_local,
                             tpol,
                             gpol,
                             leftext,
                             rightext,
                             maxmm);
            
            if(genomeHit.len() < len) {
                genomeHits.pop_back();
                continue;
            }

            if(genomeHit.score() > bestScore) {
                bestScore = genomeHit.score();
                size_t remove_count = 0;
                size_t k = prev_numHits;
                for(size_t j = prev_numHits; j < genomeHits.size(); j++) {
                    if(genomeHits[j].score() >= max(_minsc[rdi], bestScore - cushion)) {
                        assert_leq(k, j);
                        if(k < j) {
                            genomeHits[k] = genomeHits[j];
                        }
                        k++;
                    } else {
                        remove_count++;
                    }
                }
                assert_eq(k + remove_count, genomeHits.size());
                assert_leq(prev_numHits + remove_count, genomeHits.size());
                if(remove_count > 0) {
                    genomeHits.resize(genomeHits.size() - remove_count);
                }
            } else if(genomeHit.score() < max(_minsc[rdi], bestScore - cushion)) {
                genomeHits.pop_back();
            }
        }
        
        return (index_t)genomeHits.size();
    }
    
    bool pairReads(
                   const Scoring&             sc,
                   const PairedEndPolicy&     pepol, // paired-end policy
                   const TranscriptomePolicy& tpol,
                   const GraphPolicy&         gpol,
                   const GFM<index_t>&        gfm,
                   const ALTDB<index_t>&      altdb,
                   const RepeatDB<index_t>&   repeatdb,
                   const BitPairReference&    ref,
                   WalkMetrics&               wlm,
                   PerReadMetrics&            prm,
                   HIMetrics&                 him,
                   RandomSource&              rnd,
                   AlnSinkWrap<index_t>&      sink);

    /**
     *
     **/
    bool reportHit(
                   const Scoring&                   sc,
                   const PairedEndPolicy&           pepol, // paired-end policy
                   const TranscriptomePolicy&       tpol,
                   const GraphPolicy&               gpol,
                   const GFM<index_t>&              gfm,
                   const ALTDB<index_t>&            altdb,
                   const RepeatDB<index_t>&         repeatdb,
                   const BitPairReference&          ref,
                   const SpliceSiteDB&              ssdb,
                   AlnSinkWrap<index_t>&            sink,
                   index_t                          rdi,
                   const GenomeHit<index_t>&        hit,
                   bool                             alignMate = false,
                   const GenomeHit<index_t>*        ohit = NULL);
    
    /**
     * check this alignment is already examined
     **/
    bool redundant(
                   AlnSinkWrap<index_t>&    sink,
                   index_t                  rdi,
                   index_t                  tidx,
                   index_t                  toff);
    
    /**
     * check this alignment is already examined
     **/
    bool redundant(
                   AlnSinkWrap<index_t>&            sink,
                   index_t                          rdi,
                   const GenomeHit<index_t>&        hit);
    
    
    /**
     *
     **/
    bool isSearched(
                    const GenomeHit<index_t>&       hit,
                    index_t                         rdi);
    
    /**
     *
     **/
    void addSearched(const GenomeHit<index_t>&       hit,
                     index_t                         rdi);
    
    
protected:
  
    Read *   _rds[2];
    bool     _paired;
    bool     _rightendonly;
    bool     _nofw[2];
    bool     _norc[2];
    TAlScore _minsc[2];
    TAlScore _maxpen[2];
    
    bool     _anchorStop;
    
    ReadBWTHit<index_t> _hits[2][2];
    
    EList<index_t, 16>                                 _offs;
    SARangeWithOffs<EListSlice<index_t, 16>, index_t>  _sas;
    GroupWalk2S<index_t, EListSlice<index_t, 16>, 16>  _gws;
    GroupWalkState<index_t>                            _gwstate;
    
    EList<local_index_t, 16>                                       _offs_local;
    SARangeWithOffs<EListSlice<local_index_t, 16>, local_index_t>  _sas_local;
    GroupWalk2S<local_index_t, EListSlice<local_index_t, 16>, 16>  _gws_local;
    GroupWalkState<local_index_t>                                  _gwstate_local;
            
    // temporary and shared variables used for GenomeHit
    // this should be defined before _genomeHits and _hits_searched
    SharedTempVars<index_t> _sharedVars;
    
    // temporary and shared variables for AlnRes
    LinkedEList<EList<Edit> > _rawEdits;
    
    // temporary
    EList<GenomeHit<index_t> >     _genomeHits;
    EList<GenomeHit<index_t> >     _genomeHits_rep[2];
    EList<index_t>                 _snpIDs;
    EList<index_t>                 _snpIDs2;
    EList<bool>                    _genomeHits_done;
    ELList<Coord>                  _coords;
    EList<pair<RepeatCoord<index_t>, RepeatCoord<index_t> > >     _positions;
    ELList<SpliceSite>             _spliceSites;
    
    pair<index_t, index_t>         _concordantIdxInspected;
    EList<pair<index_t, index_t> > _repeatConcordant;
    
    size_t _minK; // log4 of the size of a genome
    size_t _minK_local; // log4 of the size of a local index (8)

    ELList<GenomeHit<index_t> >     _local_genomeHits;
    EList<uint8_t>                  _anchors_added;
    uint64_t max_localindexatts;
    
	uint64_t bwops_;                    // Burrows-Wheeler operations
	uint64_t bwedits_;                  // Burrows-Wheeler edits
    
    //
    EList<GenomeHit<index_t> >     _hits_searched[2];

    uint64_t   _thread_rids_mindist;
    
    //
    EList<pair<index_t, index_t> > _node_iedge_count;
    EList<pair<index_t, index_t> > _tmp_node_iedge_count;
    
    EList<pair<local_index_t, local_index_t> > _local_node_iedge_count;
    EList<pair<local_index_t, local_index_t> > _tmp_local_node_iedge_count;
    
    EList<pair<uint64_t, size_t> > _tmp_minimizers;
    ELList<RB_Alignment>           _tmp_position2D;
    EList<RB_Alignment>            _tmp_alignments;

    // For AlnRes::matchesRef
	ASSERT_ONLY(EList<bool> raw_matches_);
	ASSERT_ONLY(BTDnaString tmp_rf_);
	ASSERT_ONLY(BTDnaString tmp_rdseq_);
	ASSERT_ONLY(BTString tmp_qseq_);
};

#define HIER_INIT_LOCS(top, bot, tloc, bloc, e) { \
	if(bot - top == 1) { \
		tloc.initFromRow(top, (e).gh(), (e).gfm()); \
		bloc.invalidate(); \
	} else { \
		SideLocus<index_t>::initFromTopBot(top, bot, (e).gh(), (e).gfm(), tloc, bloc); \
		assert(bloc.valid()); \
	} \
}

#define HIER_SANITY_CHECK_4TUP(t, b, tp, bp) { \
	ASSERT_ONLY(cur_index_t tot = (b[0]-t[0])+(b[1]-t[1])+(b[2]-t[2])+(b[3]-t[3])); \
	ASSERT_ONLY(cur_index_t totp = (bp[0]-tp[0])+(bp[1]-tp[1])+(bp[2]-tp[2])+(bp[3]-tp[3])); \
	assert_eq(tot, totp); \
}

#define LOCAL_INIT_LOCS(top, bot, tloc, bloc, e) { \
    if(bot - top == 1) { \
        tloc.initFromRow(top, (e).gh(), (e).gfm()); \
        bloc.invalidate(); \
    } else { \
        SideLocus<local_index_t>::initFromTopBot(top, bot, (e).gh(), (e).gfm(), tloc, bloc); \
        assert(bloc.valid()); \
    } \
}

/**
 * Given partial alignments of a read, try to further extend
 * the alignment bidirectionally
 */
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::align(
                                               const Scoring&                   sc,
                                               const PairedEndPolicy&           pepol, // paired-end policy
                                               const TranscriptomePolicy&       tpol,
                                               const GraphPolicy&               gpol,
                                               const GFM<index_t>&              gfm,
                                               const ALTDB<index_t>&            altdb,
                                               const RepeatDB<index_t>&         repeatdb,
                                               const BitPairReference&          ref,
                                               SwAligner&                       swa,
                                               SpliceSiteDB&                    ssdb,
                                               index_t                          rdi,
                                               bool                             fw,
                                               WalkMetrics&                     wlm,
                                               PerReadMetrics&                  prm,
                                               SwMetrics&                       swm,
                                               HIMetrics&                       him,
                                               RandomSource&                    rnd,
                                               AlnSinkWrap<index_t>&            sink)
{
    const ReportingParams& rp = sink.reportingParams();
    index_t fwi = (fw ? 0 : 1);
    assert_lt(rdi, 2);
    assert(_rds[rdi] != NULL);
    ReadBWTHit<index_t>& hit = _hits[rdi][fwi];
    assert(hit.done());
    index_t minOff = 0;
    if(hit.minWidth(minOff) == std::numeric_limits<index_t>::max()) return false;
    
    // Don't try to align if the potential alignment for this read might be
    // worse than the best alignment of its reverse complement
    int64_t bestScore = (rdi == 0 ? sink.bestUnp1() : sink.bestUnp2());
    index_t num_spliced = (rdi == 0 ? sink.bestSplicedUnp1() : sink.bestSplicedUnp2());
    if(bestScore < _minsc[rdi]) bestScore = _minsc[rdi];
    index_t maxmm = (index_t)((-bestScore + sc.mmpMax - 1) / sc.mmpMax);
    index_t numActualPartialSearch = hit.numActualPartialSearch();
    if(!rp.secondary && numActualPartialSearch > maxmm + num_spliced + 1) return true;
    
    // choose candidate partial alignments for further alignment
    const index_t maxsize = max<index_t>(rp.khits, rp.kseeds);
    _genomeHits.clear();
    index_t numHits = getAnchorHits(gfm,
                                    pepol,
                                    tpol,
                                    gpol,
                                    altdb,
                                    repeatdb,
                                    ref,
                                    rnd,
                                    rdi,
                                    fw,
                                    _genomeHits,
                                    maxsize,
                                    _sharedVars,
                                    wlm,
                                    prm,
                                    him,
                                    gfm.repeat());
    if(numHits <= 0) return false;
   
    // limit the number of local index searches used for alignment of the read
    uint64_t add = 0;
    if(rp.secondary) add = (-_minsc[rdi] / sc.mmpMax) * numHits * 2;
    else           add = (-_minsc[rdi] / sc.mmpMax) * numHits;
    max_localindexatts = him.localindexatts + max<uint64_t>(10, add);
    // extend the partial alignments bidirectionally using
    // local search, extension, and (less often) global search
    hybridSearch(sc,
                 pepol,
                 tpol,
                 gpol,
                 gfm,
                 altdb,
                 repeatdb,
                 ref,
                 swa,
                 ssdb,
                 rdi,
                 fw,
                 wlm,
                 prm,
                 swm,
                 him,
                 rnd,
                 sink);
    
    return true;
}


/**
 * Given the alignment of its mate as an anchor,
 * align the read
 */
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::alignMate(
                                                   const Scoring&                   sc,
                                                   const PairedEndPolicy&           pepol, // paired-end policy
                                                   const TranscriptomePolicy&       tpol,
                                                   const GraphPolicy&               gpol,
                                                   const GFM<index_t>&              gfm,
                                                   const ALTDB<index_t>&            altdb,
                                                   const RepeatDB<index_t>&         repeatdb,
                                                   const BitPairReference&          ref,
                                                   SwAligner&                       swa,
                                                   SpliceSiteDB&                    ssdb,
                                                   index_t                          rdi,
                                                   bool                             fw,
                                                   WalkMetrics&                     wlm,
                                                   PerReadMetrics&                  prm,
                                                   SwMetrics&                       swm,
                                                   HIMetrics&                       him,
                                                   RandomSource&                    rnd,
                                                   AlnSinkWrap<index_t>&            sink,
                                                   index_t                          tidx,
                                                   index_t                          toff)
{
    const ReportingParams& rp = sink.reportingParams();
    
    assert_lt(rdi, 2);
    index_t ordi = 1 - rdi;
    bool ofw = (fw == gMate2fw ? gMate1fw : gMate2fw);
    assert(_rds[ordi] != NULL);
    const Read& ord = *_rds[ordi];
    index_t rdlen = (index_t)ord.length();
    assert_gt(rdlen, 0);
    
    _genomeHits.clear();
    if(_coords.size() == 0) {
        _coords.expand();
    }
    EList<Coord>& coords = _coords.front();
    
    // local search to find anchors
    const HGFM<index_t, local_index_t>* hGFM = (const HGFM<index_t, local_index_t>*)(&gfm);
    const LocalGFM<local_index_t, index_t>* lGFM = hGFM->getLocalGFM(tidx, toff);
    bool first = true;
    index_t count = 0;
    index_t max_hitlen = 0;
    while(count++ < 2) {
        if(first) {
            first = false;
        } else {
            if(_genomeHits.size() > 0) break;
            if(fw) {
                lGFM = hGFM->nextLocalGFM(lGFM);
            } else {
                lGFM = hGFM->prevLocalGFM(lGFM);
            }
            if(lGFM == NULL || lGFM->empty()) break;
        }
        index_t hitoff = rdlen - 1;
        while(hitoff >= _minK_local - 1) {
            index_t hitlen = 0;
            local_index_t top = (local_index_t)INDEX_MAX, bot = (local_index_t)INDEX_MAX;
            local_index_t node_top = (local_index_t)INDEX_MAX, node_bot = (local_index_t)INDEX_MAX;
            _local_node_iedge_count.clear();
            bool uniqueStop = false;
            index_t nelt = localGFMSearch(
                                          *lGFM,   // GFM index
                                          ord,     // read to align
                                          sc,      // scoring scheme
                                          sink.reportingParams(),
                                          ofw,
                                          hitoff,
                                          hitlen,
                                          top,
                                          bot,
                                          node_top,
                                          node_bot,
                                          _local_node_iedge_count,
                                          rnd,
                                          uniqueStop,
                                          _minK_local);
            assert_leq(top, bot);
            assert_eq(nelt, (index_t)(node_bot - node_top));
            assert_leq(hitlen, hitoff + 1);
            if(nelt > 0 && nelt <= rp.kseeds && hitlen > max_hitlen) {
                coords.clear();
                bool straddled = false;
                getGenomeCoords_local(
                                      *lGFM,
                                      altdb,
                                      ref,
                                      rnd,
                                      top,
                                      bot,
                                      node_top,
                                      node_bot,
                                      _local_node_iedge_count,
                                      ofw,
                                      hitoff - hitlen + 1,
                                      hitlen,
                                      coords,
                                      wlm,
                                      prm,
                                      him,
                                      true, // reject straddled?
                                      straddled);
                assert_leq(coords.size(), nelt);
                _genomeHits.clear();
                for(index_t ri = 0; ri < coords.size(); ri++) {
                    const Coord& coord = coords[ri];
                    if(tpol.no_spliced_alignment()) {
                        if(coord.off() + pepol.maxFragLen() * 2 < toff || toff + pepol.maxFragLen() * 2 < coord.off())
                            continue;
                    }
                    
                    GenomeHit<index_t>::adjustWithALT(
                                                      hitoff - hitlen + 1,
                                                      hitlen,
                                                      coord,
                                                      _sharedVars,
                                                      _genomeHits,
                                                      *this->_rds[ordi],
                                                      gfm,
                                                      altdb,
                                                      ref,
                                                      gpol);
                }
                 max_hitlen = hitlen;
            }
            assert_leq(hitlen, hitoff + 1);
            if(hitlen > 0) hitoff -= (hitlen - 1);
            if(hitoff > 0) hitoff -= 1;
        } // while(hitoff >= minHitLen - 1)
    } // while(count++ < 2)
    
    // randomly select
    const index_t maxsize = rp.kseeds;
    if(_genomeHits.size() > maxsize) {
        _genomeHits.shufflePortion(0, _genomeHits.size(), rnd);
        _genomeHits.resize(maxsize);
    }
    
    // local search using the anchor
    for(index_t hi = 0; hi < _genomeHits.size(); hi++) {
        him.anchoratts++;
        GenomeHit<index_t>& genomeHit = _genomeHits[hi];
        index_t leftext = (index_t)INDEX_MAX, rightext = (index_t)INDEX_MAX;
        genomeHit.extend(
                         ord,
                         gfm,
                         ref,
                         altdb,
                         repeatdb,
                         ssdb,
                         swa,
                         swm,
                         prm,
                         sc,
                         _minsc[ordi],
                         rnd,
                         (index_t)_minK_local,
                         tpol,
                         gpol,
                         leftext,
                         rightext);
        hybridSearch_recur(
                           sc,
                           pepol,
                           tpol,
                           gpol,
                           gfm,
                           altdb,
                           repeatdb,
                           ref,
                           swa,
                           ssdb,
                           ordi,
                           genomeHit,
                           genomeHit.rdoff(),
                           genomeHit.len(),
                           wlm,
                           prm,
                           swm,
                           him,
                           rnd,
                           sink,
                           true); // alignMate?
    }
    
    return true;
}


/**
 * convert FM offsets to the corresponding genomic offset (chromosome id, offset)
 **/
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::getGenomeCoords(
                                                         const GFM<index_t>&        gfm,
                                                         const ALTDB<index_t>&      altdb,
                                                         const BitPairReference&    ref,
                                                         RandomSource&              rnd,
                                                         index_t                    top,
                                                         index_t                    bot,
                                                         index_t                    node_top,
                                                         index_t                    node_bot,
                                                         const EList<pair<index_t, index_t> >& node_iedge_count,
                                                         bool                       fw,
                                                         index_t                    maxelt,
                                                         index_t                    rdoff,
                                                         index_t                    rdlen,
                                                         EList<Coord>&              coords,
                                                         WalkMetrics&               met,
                                                         PerReadMetrics&            prm,
                                                         HIMetrics&                 him,
                                                         bool                       rejectStraddle,
                                                         bool&                      straddled)
{
    straddled = false;
    assert_gt(bot, top);
    assert_leq(node_bot - node_top, bot - top);
    index_t nelt = node_bot - node_top;
    nelt = min<index_t>(nelt, maxelt);
    him.globalgenomecoords += nelt;
    _offs.resize(nelt);
    _offs.fill((index_t)INDEX_MAX);
    _sas.init(
              top,
              bot,
              node_top,
              node_bot,
              node_iedge_count,
              rdlen,
              EListSlice<index_t, 16>(_offs, 0, nelt));
    _gws.init(gfm, ref, _sas, rnd, met);
    
    for(index_t off = 0; off < nelt; off++) {
        WalkResult<index_t> wr;
        index_t tidx = 0, toff = 0, tlen = 0;
        _gws.advanceElement(
                            off,
                            gfm,          // forward Bowtie index for walking left
                            ref,          // bitpair-encoded reference
                            _sas,         // SA range with offsets
                            _gwstate,     // GroupWalk state; scratch space
                            wr,           // put the result here
                            met,          // metrics
                            prm);         // per-read metrics
        assert_neq(wr.toff, (index_t)INDEX_MAX);
        bool straddled2 = false;
        gfm.joinedToTextOff(
                            wr.elt.len,
                            wr.toff,
                            tidx,
                            toff,
                            tlen,
                            rejectStraddle,        // reject straddlers?
                            straddled2);  // straddled?
        
        straddled |= straddled2;
        
        if(tidx == (index_t)INDEX_MAX) {
            // The seed hit straddled a reference boundary so the seed
            // hit isn't valid
            return false;
        }
        index_t global_toff = toff, global_tidx = tidx;
        
        // Coordinate of the seed hit w/r/t the pasted reference string
        coords.expand();
        if(!straddled2) {
            coords.back().init(global_tidx, (int64_t)global_toff, fw, wr.toff);
        } else {
            coords.back().init(numeric_limits<index_t>::max(), (int64_t)global_toff, fw, wr.toff);
        }
    }
    
    return true;
}

/**
 * convert FM offsets to the corresponding genomic offset (chromosome id, offset)
 **/
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::getGenomeCoords_local(
                                                               const GFM<local_index_t>&    gfm,
                                                               const ALTDB<index_t>&        altdb,
                                                               const BitPairReference&      ref,
                                                               RandomSource&                rnd,
                                                               local_index_t                top,
                                                               local_index_t                bot,
                                                               local_index_t                node_top,
                                                               local_index_t                node_bot,
                                                               const EList<pair<local_index_t, local_index_t> >& node_iedge_count,
                                                               bool                         fw,
                                                               index_t                      rdoff,
                                                               index_t                      rdlen,
                                                               EList<Coord>&                coords,
                                                               WalkMetrics&                 met,
                                                               PerReadMetrics&              prm,
                                                               HIMetrics&                   him,
                                                               bool                         rejectStraddle,
                                                               bool&                        straddled)
{
    straddled = false;
    assert_gt(bot, top);
    assert_leq(node_bot - node_top, bot - top);
    index_t nelt = node_bot - node_top;
    him.localgenomecoords += nelt;
    _offs_local.resize(nelt);
    _offs_local.fill((local_index_t)INDEX_MAX);
    _sas_local.init(
                    top,
                    bot,
                    node_top,
                    node_bot,
                    node_iedge_count,
                    rdlen,
                    EListSlice<local_index_t, 16>(_offs_local, 0, nelt));
    _gws_local.init(gfm, ref, _sas_local, rnd, met);
    
    for(local_index_t off = 0; off < nelt; off++) {
        WalkResult<local_index_t> wr;
        local_index_t tidx = 0, toff = 0, tlen = 0;
        _gws_local.advanceElement(
                                  off,
                                  gfm,          // forward Bowtie index for walking left
                                  ref,          // bitpair-encoded reference
                                  _sas_local,   // SA range with offsets
                                  _gwstate_local, // GroupWalk state; scratch space
                                  wr,           // put the result here
                                  met,          // metrics
                                  prm);         // per-read metrics
        assert_neq(wr.toff, (local_index_t)INDEX_MAX);
        bool straddled2 = false;
        bool result = gfm.joinedToTextOff(
                                          wr.elt.len,
                                          wr.toff,
                                          tidx,
                                          toff,
                                          tlen,
                                          rejectStraddle,        // reject straddlers?
                                          straddled2);  // straddled?
        if(!result) continue;
        
        straddled |= straddled2;
        
        if(tidx == (local_index_t)INDEX_MAX) {
            // The seed hit straddled a reference boundary so the seed
            // hit isn't valid
            return false;
        }
        LocalGFM<local_index_t, index_t>* localGFM = (LocalGFM<local_index_t, index_t>*)&gfm;
        index_t global_tidx = localGFM->_tidx;
        index_t global_toff = toff + localGFM->_localOffset;
        index_t joinedOff = wr.toff + localGFM->_joinedOffset;
        if(global_toff < rdoff) continue;
        
        // Coordinate of the seed hit w/r/t the pasted reference string
        coords.expand();
        coords.back().init(global_tidx, (int64_t)global_toff, fw, joinedOff);
    }
    
    return true;
}


/**
 * examine alignments of left and right reads to produce concordant pair alignment
 **/
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::pairReads(
                                                   const Scoring&             sc,
                                                   const PairedEndPolicy&     pepol, // paired-end policy
                                                   const TranscriptomePolicy& tpol,
                                                   const GraphPolicy&         gpol,
                                                   const GFM<index_t>&        gfm,
                                                   const ALTDB<index_t>&      altdb,
                                                   const RepeatDB<index_t>&   repeatdb,
                                                   const BitPairReference&    ref,
                                                   WalkMetrics&               wlm,
                                                   PerReadMetrics&            prm,
                                                   HIMetrics&                 him,
                                                   RandomSource&              rnd,
                                                   AlnSinkWrap<index_t>&      sink)
{
    const ReportingParams& rp = sink.reportingParams();
    assert(_paired);
    const EList<AlnRes> *rs1 = NULL, *rs2 = NULL;
    sink.getUnp1(rs1); assert(rs1 != NULL);
    sink.getUnp2(rs2); assert(rs2 != NULL);
    index_t start_i = _concordantIdxInspected.first, start_j = _concordantIdxInspected.second;
    _concordantIdxInspected.first = rs1->size();
    _concordantIdxInspected.second = rs2->size();
    for(index_t i = 0; i < rs1->size(); i++) {
        for(index_t j = (i >= start_i ? 0 : start_j); j < rs2->size(); j++) {
            if(sink.state().doneConcordant()) {
                return true;
            }
            const AlnRes& r1 = (*rs1)[i];
            Coord left = r1.refcoord(), right = r1.refcoord_right();
            assert_eq(left.ref(), right.ref());
            const AlnRes& r2 = (*rs2)[j];
            Coord left2 = r2.refcoord(), right2 = r2.refcoord_right();
            assert_eq(left2.ref(), right2.ref());
            
            if(r1.repeat() != r2.repeat())
                continue;
            
            bool dna_frag_pass = true;
            if(r1.repeat() && r2.repeat()) {
                bool found = false;
                for(size_t r = 0; r < _repeatConcordant.size(); r++) {
                    if(_repeatConcordant[r].first == left.off() &&
                       _repeatConcordant[r].second == left2.off()) {
                        found = true;
                        break;
                    }
                }
                if(!found) {
                    dna_frag_pass = false;
                }
            } else{
                if(left.ref() != left2.ref()) continue;
                assert_eq(left.orient(), right.orient());
                assert_eq(left2.orient(), right2.orient());
                if(left.orient() == gMate1fw) {
                    if(left2.orient() != gMate2fw) continue;
                } else {
                    if(left2.orient() == gMate2fw) continue;
                    Coord temp = left; left = left2; left2 = temp;
                    temp = right; right = right2; right2 = temp;
                }
                if(left.off() > left2.off()) continue;
                if(right.off() > right2.off()) continue;
                if(right.off() + (int)tpol.maxIntronLen() < left2.off()) continue;
                assert_geq(r1.score().score(), _minsc[0]);
                assert_geq(r2.score().score(), _minsc[1]);
                if(tpol.no_spliced_alignment()){
                    int pairCl = PE_ALS_DISCORD;
                    assert_eq(r1.refid(), r2.refid());
                    index_t off1, off2, len1, len2;
                    bool fw1, fw2;
                    if(r1.refoff() < r2.refoff()) {
                        off1 = r1.refoff(); off2 = r2.refoff();
                        len1 = r1.refExtent(); len2 = r2.refExtent();
                        fw1 = r1.fw(); fw2 = r2.fw();
                    } else {
                        off1 = r2.refoff(); off2 = r1.refoff();
                        len1 = r2.refExtent(); len2 = r1.refExtent();
                        fw1 = r2.fw(); fw2 = r1.fw();
                    }
                    // Check that final mate alignments are consistent with
                    // paired-end fragment constraints
                    pairCl = pepol.peClassifyPair(
                                                  off1,
                                                  len1,
                                                  fw1,
                                                  off2,
                                                  len2,
                                                  fw2);
                    dna_frag_pass = (pairCl != PE_ALS_DISCORD);
                }
            }
        
            if(!tpol.no_spliced_alignment() || dna_frag_pass) {
                TAlScore threshold = sink.bestPair();
                if(sink.bestUnp1() >= _minsc[0] && sink.bestUnp2() >= _minsc[1]) {
                    TAlScore tmp = sink.bestUnp1() + sink.bestUnp2() - (r1.readLength() + r2.readLength()) * 0.03 * sc.mm(255);
                    if(tmp > threshold) {
                        threshold = tmp;
                    }
                }
                if(r1.score().score() + r2.score().score() >= threshold || rp.secondary) {
                    sink.report(0, &r1, &r2);
                }
            }
        }
    }
    return true;
}


/**
 * report read (or pair) alignment
 **/
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::reportHit(
                                                   const Scoring&                   sc,
                                                   const PairedEndPolicy&           pepol, // paired-end policy
                                                   const TranscriptomePolicy&       tpol,
                                                   const GraphPolicy&               gpol,
                                                   const GFM<index_t>&              gfm,
                                                   const ALTDB<index_t>&            altdb,
                                                   const RepeatDB<index_t>&         repeatdb,
                                                   const BitPairReference&          ref,
                                                   const SpliceSiteDB&              ssdb,
                                                   AlnSinkWrap<index_t>&            sink,
                                                   index_t                          rdi,
                                                   const GenomeHit<index_t>&        hit,
                                                   bool                             alignMate,
                                                   const GenomeHit<index_t>*        ohit)
{
    assert_lt(rdi, 2);
    assert(_rds[rdi] != NULL);
    const Read& rd = *_rds[rdi];
    index_t rdlen = (index_t)rd.length();
    if(hit.rdoff() - hit.trim5() > 0 || hit.len() + hit.trim5() + hit.trim3() < rdlen) return false;
    if(hit.score() < _minsc[rdi]) return false;
    if(!sink.reportingParams().repeat && hit.repeat()) return false;
    
    // Edits are represented from 5' end of read to 3' end, not an alignment of read
    EList<Edit>& edits = const_cast<EList<Edit>&>(hit.edits());
    if(hit.trim5() > 0) {
        for(size_t i = 0; i < edits.size(); i++) {
            edits[i].pos += hit.trim5();
        }
    }
    if(!hit.fw()) {
        Edit::invertPoss(edits, rdlen, false);
    }
    // in case of multiple exonic alignments, choose the ones near (known) splice sites
    // this helps eliminate cases of reads being mapped to pseudogenes
    pair<bool, bool> spliced = hit.spliced(); // pair<spliced, spliced_to_known>
    if(tpol.xs_only() && spliced.first) {
        if(hit.splicing_dir() == SPL_UNKNOWN)
            return false;
    }
    if(!tpol.no_spliced_alignment() && tpol.avoid_pseudogene()) {
        if(!spliced.first) {
            assert(!spliced.second);
            const index_t max_exon_size = 10000;
            index_t left = 0;
            if(hit.refoff() > max_exon_size) {
                left = hit.refoff() - max_exon_size;
            }
            index_t right = hit.refoff() + hit.len() + max_exon_size;
            spliced.first = ssdb.hasSpliceSites(
                                                hit.ref(),
                                                left,
                                                right,
                                                left,
                                                right,
                                                true); // include novel splice sites
            if(altdb.hasExons()) {
                spliced.second = ssdb.insideExon(hit.ref(), hit.refoff(), hit.refoff() + hit.len() - 1);
            }
        }
    }
    if(tpol.transcriptome_mapping_only() && !spliced.second)
        return false;
    
    AlnScore asc(
                 hit.score(),       // numeric score
                 hit.ns(),          // # Ns
                 hit.ngaps(),       // # gaps
                 hit.repeat(),
                 hit.splicescore(), // splice score
                 spliced.second,    // mapped to known transcripts?
                 spliced.first,     // spliced alignment or near splice sites (novel)?
                 hit.trim5(),       // left trim length
                 hit.trim3());      // right trim length
    bool softTrim = hit.trim5() > 0 || hit.trim3() > 0;
    AlnRes rs;
    rs.init(
            rdlen,                      // # chars after hard trimming
            rd.rdid,                    // read ID
            asc,                        // alignment score
            &hit.edits(),               // nucleotide edits array
            0,                          // nucleotide edits first pos
            hit.edits().size(),         // nucleotide edits last pos
            NULL,                       // ambig base array
            0,                          // ambig base first pos
            0,                          // ambig base last pos
            hit.coord(),                // coord of leftmost aligned char in ref
            hit.repeat() ? gfm.plen()[0] : gfm.plen()[hit.ref()], // length of reference aligned to
            &_rawEdits,
            -1,                         // # seed mms allowed
            -1,                         // seed length
            -1,                         // seed interval
            0,                          // minimum score for valid alignment (daehwan)
            -1,                         // nuc5p (for colorspace)
            -1,                         // nuc3p (for colorspace)
            false,                      // soft pre-trimming?
            0,                          // 5p pre-trimming
            0,                          // 3p pre-trimming
            softTrim,                   // soft trimming?
            hit.fw() ? hit.trim5() : hit.trim3(),  // 5p trimming
            hit.fw() ? hit.trim3() : hit.trim5(),  // 3p trimming
            hit.repeat());              // repeat?
    if(!hit.fw()) {
        Edit::invertPoss(edits, rdlen, false);
    }
    if(hit.trim5() > 0) {
        for(size_t i = 0; i < edits.size(); i++) {
            edits[i].pos -= hit.trim5();
        }
    }
    

    //rs.setRefNs(nrefn);
    assert(rs.matchesRef(
                         rd,
                         ref,
                         tmp_rf_,
                         tmp_rdseq_,
                         tmp_qseq_,
                         _sharedVars.raw_refbuf,
                         _sharedVars.destU32,
                         raw_matches_,
                         _sharedVars.raw_refbuf2,
                         _sharedVars.reflens,
                         _sharedVars.refoffs));
    if(ohit == NULL) {
        bool done;
        if(rdi == 0 && !_rightendonly) {
            done = sink.report(0, &rs, NULL, alignMate);
        } else {
            done = sink.report(0, NULL, &rs, alignMate);
        }
        return done;
    }
    
    assert(ohit != NULL);
    const Read& ord = *_rds[1-rdi];
    index_t ordlen = (index_t)ord.length();
    if(ohit->rdoff() - ohit->trim5() > 0 || ohit->len() + ohit->trim5() + ohit->trim3() < ordlen) return false;
    if(ohit->score() < _minsc[1-rdi]) return false;
    EList<Edit>& oedits = const_cast<EList<Edit>&>(ohit->edits());
    if(ohit->trim5() > 0) {
        for(size_t i = 0; i < oedits.size(); i++) {
            oedits[i].pos += ohit->trim5();
        }
    }
    if(!ohit->fw()) {
        Edit::invertPoss(oedits, ordlen, false);
    }
    AlnScore oasc(
                  ohit->score(),   // numeric score
                  ohit->ns(),      // # Ns
                  ohit->ngaps(),   // # gaps
                  ohit->repeat()); // repeat?
    bool osoftTrim = ohit->trim5() > 0 || ohit->trim3() > 0;
    AlnRes ors;
    ors.init(
             ordlen,                     // # chars after hard trimming
             ord.rdid,                   // read ID
             oasc,                       // alignment score
             &ohit->edits(),             // nucleotide edits array
             0,                          // nucleotide edits first pos
             ohit->edits().size(),       // nucleotide edits last pos
             NULL,                       // ambig base array
             0,                          // ambig base first pos
             0,                          // ambig base last pos
             ohit->coord(),              // coord of leftmost aligned char in ref
             gfm.plen()[ohit->ref()],    // length of reference aligned to
             &_rawEdits,
             -1,                         // # seed mms allowed
             -1,                         // seed length
             -1,                         // seed interval
             0,                          // minimum score for valid alignment (daehwan)
             -1,                         // nuc5p (for colorspace)
             -1,                         // nuc3p (for colorspace)
             false,                      // soft pre-trimming?
             0,                          // 5p pre-trimming
             0,                          // 3p pre-trimming
             osoftTrim,                  // soft trimming?
             ohit->fw() ? ohit->trim5() : ohit->trim3(),  // 5p trimming
             ohit->fw() ? ohit->trim3() : ohit->trim5(),  // 3p trimming
             ohit->repeat());            // repeat?
    if(!ohit->fw()) {
        Edit::invertPoss(oedits, ordlen, false);
    }
    if(ohit->trim5() > 0) {
        for(size_t i = 0; i < oedits.size(); i++) {
            oedits[i].pos -= ohit->trim5();
        }
    }
    //rs.setRefNs(nrefn);
    assert(ors.matchesRef(
                          ord,
                          ref,
                          tmp_rf_,
                          tmp_rdseq_,
                          tmp_qseq_,
                          _sharedVars.raw_refbuf,
                          _sharedVars.destU32,
                          raw_matches_,
                          _sharedVars.raw_refbuf2,
                          _sharedVars.reflens,
                          _sharedVars.refoffs));
    
    bool done;
    if(rdi == 0) {
        done = sink.report(0, &rs, &ors);
    } else {
        done = sink.report(0, &ors, &rs);
    }
    return done;
}

/**
 * check this alignment is already examined
 **/
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::redundant(
                                                   AlnSinkWrap<index_t>&    sink,
                                                   index_t                  rdi,
                                                   index_t                  tidx,
                                                   index_t                  toff)
{
    assert_lt(rdi, 2);
    const EList<AlnRes>* rs = NULL;
    if(rdi == 0) sink.getUnp1(rs);
    else         sink.getUnp2(rs);
    assert(rs != NULL);
    for(index_t i = 0; i < rs->size(); i++) {
        Coord coord_left = (*rs)[i].refcoord(), coord_right = (*rs)[i].refcoord_right();
        assert_eq(coord_left.ref(), coord_right.ref());
        assert_lt(coord_left.off(), coord_right.off());
        assert_eq(coord_left.orient(), coord_right.orient());
        
        if(tidx != coord_left.ref()) continue;
        if(toff >= coord_left.off() && toff <= coord_right.off()) return true;
    }
    
    return false;
}


/**
 * check this alignment is already examined
 **/
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::redundant(
                                                   AlnSinkWrap<index_t>&            sink,
                                                   index_t                          rdi,
                                                   const GenomeHit<index_t>&        hit)
{
    assert_lt(rdi, 2);
    assert(_rds[rdi] != NULL);
    index_t rdlen = (index_t)_rds[rdi]->length();
    const EList<AlnRes>* rs = NULL;
    if(rdi == 0) sink.getUnp1(rs);
    else         sink.getUnp2(rs);
    assert(rs != NULL);
    for(index_t i = 0; i < rs->size(); i++) {
        const AlnRes& rsi = (*rs)[i];
        if(rsi.refcoord() == hit.coord()) {
            const EList<Edit>& editsi = rsi.ned();
            const EList<Edit>& edits = hit.edits();
            if(editsi.size() == edits.size()) {
                size_t eidx = 0;
                if(!hit.fw()) {
                    Edit::invertPoss(const_cast<EList<Edit>&>(edits), rdlen, false);
                }
                // daehwan: this is a temporary solution to compare edits
                for(; eidx < editsi.size(); eidx++) {
                    if(!(editsi[eidx] == edits[eidx])) {
                        break;
                    }
                }
                if(!hit.fw()) {
                    Edit::invertPoss(const_cast<EList<Edit>&>(edits), rdlen, false);
                }
                if(eidx >= editsi.size()) {
                    assert_eq(eidx, editsi.size());
                    return true;
                }
            }
        }
    }
    
    return false;
}


/**
 * Sweep right-to-left and left-to-right using exact matching.  Remember all
 * the SA ranges encountered along the way.  Report exact matches if there are
 * any.  Calculate a lower bound on the number of edits in an end-to-end
 * alignment.
 */
template <typename index_t, typename local_index_t>
index_t HI_Aligner<index_t, local_index_t>::partialSearch(
                                                          const GFM<index_t>&       gfm,    // BWT index
                                                          const Read&               read,    // read to align
                                                          const Scoring&            sc,      // scoring scheme
                                                          const ReportingParams&    rp,
                                                          bool                      fw,
                                                          size_t                    mineMax, // don't care about edit bounds > this
                                                          size_t&                   mineFw,  // minimum # edits for forward read
                                                          size_t&                   mineRc,  // minimum # edits for revcomp read
                                                          ReadBWTHit<index_t>&      hit,     // holds all the seed hits (and exact hit)
                                                          RandomSource&             rnd,     // pseudo-random source
                                                          bool&                     pseudogeneStop,
                                                          bool&                     anchorStop,
                                                          index_t                   maxHitLen)
{
    bool pseudogeneStop_ = pseudogeneStop, anchorStop_ = anchorStop;
    pseudogeneStop = anchorStop = false;
	const index_t ftabLen = gfm.gh().ftabChars();
    const bool linearFM = gfm.gh().linearFM();
	SideLocus<index_t> tloc, bloc;
	const index_t len = (index_t)read.length();
    const BTDnaString& seq = fw ? read.patFw : read.patRc;
    assert(!seq.empty());
    
    size_t nelt = 0;
    EList<BWTHit<index_t> >& partialHits = hit._partialHits;
    index_t& cur = hit._cur;
    assert_lt(cur, hit._len);
    
    hit._numPartialSearch++;
    
    index_t offset = cur;
    index_t dep = offset;
    pair<index_t, index_t> range(0, 0);
    pair<index_t, index_t> rangeTemp(0, 0);
    pair<index_t, index_t> node_range(0, 0);
    pair<index_t, index_t> node_rangeTemp(0, 0);
    _node_iedge_count.clear();
    _tmp_node_iedge_count.clear();
    index_t left = len - dep;
    assert_gt(left, 0);
    if(left < ftabLen + 1) {
        cur = hit._len;
        partialHits.expand();
        partialHits.back().init((index_t)INDEX_MAX,
                                (index_t)INDEX_MAX,
                                (index_t)INDEX_MAX,
                                (index_t)INDEX_MAX,
                                _node_iedge_count,
                                fw,
                                (index_t)offset,
                                (index_t)(cur - offset));
        hit.done(true);
		return 0;
    }
    // Does N interfere with use of Ftab?
    for(index_t i = 0; i < ftabLen; i++) {
        int c = seq[len-dep-1-i];
        if(c > 3) {
            cur += (i+1);
            partialHits.expand();
            partialHits.back().init((index_t)INDEX_MAX,
                                    (index_t)INDEX_MAX,
                                    (index_t)INDEX_MAX,
                                    (index_t)INDEX_MAX,
                                    _node_iedge_count,
                                    fw,
                                    (index_t)offset,
                                    (index_t)(cur - offset));
            if(cur >= hit._len) {
                hit.done(true);
            }
			return 0;
        }
    }
    
    // Use ftab
    gfm.ftabLoHi(seq, len - dep - ftabLen, false, range.first, range.second);
    dep += ftabLen;
    if(range.first >= range.second) {
        cur = dep;
        partialHits.expand();
        partialHits.back().init((index_t)INDEX_MAX,
                                (index_t)INDEX_MAX,
                                (index_t)INDEX_MAX,
                                (index_t)INDEX_MAX,
                                _node_iedge_count,
                                fw,
                                (index_t)offset,
                                (index_t)(cur - offset));
        if(cur >= hit._len) {
            hit.done(true);
        }
        return 0;
    }
    index_t same_range = 0, similar_range = 0;
    HIER_INIT_LOCS(range.first, range.second, tloc, bloc, gfm);
    // Keep going
    while(dep < len && dep - offset < maxHitLen) {
        int c = seq[len-dep-1];
        if(c > 3) {
            rangeTemp.first = rangeTemp.second = 0;
            node_rangeTemp.first = node_rangeTemp.second = 0;
            _tmp_node_iedge_count.clear();
        } else {
            if(bloc.valid()) {
                bwops_ += 2;
                if(linearFM) {
                    rangeTemp = gfm.mapLF(tloc, bloc, c, &node_rangeTemp);
                } else {
                    rangeTemp = gfm.mapGLF(tloc, bloc, c, &node_rangeTemp, &_tmp_node_iedge_count, (index_t)rp.kseeds);
                }
            } else {
                bwops_++;
                rangeTemp = gfm.mapGLF1(range.first, tloc, c, &node_rangeTemp);
                if(rangeTemp.first + 1 < rangeTemp.second) {
                    assert_eq(node_rangeTemp.first + 1, node_rangeTemp.second);
                    _tmp_node_iedge_count.clear();
                    _tmp_node_iedge_count.expand();
                    _tmp_node_iedge_count.back().first = 0;
                    _tmp_node_iedge_count.back().second = rangeTemp.second - rangeTemp.first - 1;
                }
            }
        }
        if(rangeTemp.first >= rangeTemp.second) {
            break;
        }
        if(pseudogeneStop_) {
            if(node_rangeTemp.second - node_rangeTemp.first < node_range.second - node_range.first && node_range.second - node_range.first <= min<index_t>(5, (index_t)rp.khits)) {
                static const index_t minLenForPseudogene = (index_t)_minK + 6;
                if(dep - offset >= minLenForPseudogene && similar_range >= 5) {
                    hit._numUniqueSearch++;
                    pseudogeneStop = true;
                    break;
                }
            }
            if(node_rangeTemp.second - node_rangeTemp.first != 1) {
                if(node_rangeTemp.second - node_rangeTemp.first + 2 >= node_range.second - node_range.first) similar_range++;
                else if(node_rangeTemp.second - node_rangeTemp.first + 4 < node_range.second - node_range.first) similar_range = 0;
            } else {
                pseudogeneStop_ = false;
            }
        }
        
        if(anchorStop_) {
            if(node_rangeTemp.second - node_rangeTemp.first != 1 && node_range.second - node_range.first == node_rangeTemp.second - node_rangeTemp.first) {
                same_range++;
                if(same_range >= 5) {
                    anchorStop_ = false;
                }
            } else {
                same_range = 0;
            }
        
            if(dep - offset >= _minK + 8 && node_rangeTemp.second - node_rangeTemp.first >= 4) {
                anchorStop_ = false;
            }
        }
        
        range = rangeTemp;
        node_range = node_rangeTemp;
        if(_tmp_node_iedge_count.size() > 0) {
            _node_iedge_count = _tmp_node_iedge_count;
            _tmp_node_iedge_count.clear();
        } else {
            _node_iedge_count.clear();
        }
        dep++;
        
        if(anchorStop_) {
            if(dep - offset >= _minK + 12 && range.second - range.first == 1) {
                hit._numUniqueSearch++;
                anchorStop = true;
                break;
            }
        }
        
        HIER_INIT_LOCS(range.first, range.second, tloc, bloc, gfm);
    }
    
    // Done
    if(range.first < range.second) {
        assert_leq(node_range.second - node_range.first, range.second - range.first);
        assert_gt(dep, offset);
        assert_leq(dep, len);
        partialHits.expand();
        index_t hit_type = CANDIDATE_HIT;
        if(anchorStop) hit_type = ANCHOR_HIT;
        else if(pseudogeneStop) hit_type = PSEUDOGENE_HIT;
        bool report = node_range.first < node_range.second;
        if(node_range.second - node_range.first < range.second - range.first) {
            if(_node_iedge_count.size() == 0) report = false;
        }
        if(report) {
#ifndef NDEBUG
            if(node_range.second - node_range.first < range.second - range.first) {
                ASSERT_ONLY(index_t add = 0);
                for(index_t e = 0; e < _node_iedge_count.size(); e++) {
                    if(e > 0) {
                        assert_lt(_node_iedge_count[e-1].first, _node_iedge_count[e].first);
                    }
                    assert_gt(_node_iedge_count[e].second, 0);
                    add += _node_iedge_count[e].second;
                }
                assert_eq(node_range.second - node_range.first + add, range.second - range.first);
            } else {
                assert(_node_iedge_count.empty());
            }
#endif
            partialHits.back().init(range.first,
                                    range.second,
                                    node_range.first,
                                    node_range.second,
                                    _node_iedge_count,
                                    fw,
                                    (index_t)offset,
                                    (index_t)(dep - offset),
                                    hit_type);
        } else {
            _node_iedge_count.clear();
            partialHits.back().init(INDEX_MAX,
                                    INDEX_MAX,
                                    INDEX_MAX,
                                    INDEX_MAX,
                                    _node_iedge_count,
                                    fw,
                                    (index_t)offset,
                                    (index_t)(dep - offset),
                                    hit_type);
        }
        
        nelt += (node_range.second - node_range.first);
        cur = dep;
        if(cur >= hit._len) {
            if(hit_type == CANDIDATE_HIT) hit._numUniqueSearch++;
            hit.done(true);
        }
    }
    return (index_t)nelt;
}


/**
 */
template <typename index_t, typename local_index_t>
index_t HI_Aligner<index_t, local_index_t>::globalGFMSearch(
                                                            const GFM<index_t>&    gfm,  // BWT index
                                                            const Read&            read,  // read to align
                                                            const Scoring&         sc,    // scoring scheme
                                                            const ReportingParams& rp,
                                                            bool                   fw,
                                                            index_t                hitoff,
                                                            index_t&               hitlen,
                                                            index_t&               top,
                                                            index_t&               bot,
                                                            index_t&               node_top,
                                                            index_t&               node_bot,
                                                            EList<pair<index_t, index_t> >& node_iedge_count,
                                                            RandomSource&          rnd,
                                                            bool&                  uniqueStop,
                                                            index_t                maxHitLen)
{
    bool uniqueStop_ = uniqueStop;
    uniqueStop = false;
    const index_t ftabLen = gfm.gh().ftabChars();
    const bool linearFM = gfm.gh().linearFM();
	SideLocus<index_t> tloc, bloc;
	const index_t len = (index_t)read.length();
    
	size_t nelt = 0;
    const BTDnaString& seq = fw ? read.patFw : read.patRc;
    assert(!seq.empty());
    
    index_t offset = len - hitoff - 1;
    index_t dep = offset;
    pair<index_t, index_t> range(0, 0);
    pair<index_t, index_t> rangeTemp(0, 0);
    pair<index_t, index_t> node_range(0, 0);
    pair<index_t, index_t> node_rangeTemp(0, 0);
    node_iedge_count.clear();
    _tmp_node_iedge_count.clear();
    index_t left = len - dep;
    assert_gt(left, 0);
    if(left < ftabLen + 1) {
        hitlen = left;
        return 0;
    }
    
    // Does N interfere with use of Ftab?
    for(index_t i = 0; i < ftabLen; i++) {
        int c = seq[len-dep-1-i];
        if(c > 3) {
            hitlen = (i+1);
            return 0;
        }
    }
    
    // Use ftab
    gfm.ftabLoHi(seq, len - dep - ftabLen, false, range.first, range.second);
    dep += ftabLen;
    if(range.first >= range.second) {
        hitlen = ftabLen;
        return 0;
    }
    
    HIER_INIT_LOCS(range.first, range.second, tloc, bloc, gfm);
    // Keep going
    while(dep < len) {
        int c = seq[len-dep-1];
        if(c > 3) {
            rangeTemp.first = rangeTemp.second = 0;
            node_rangeTemp.first = node_rangeTemp.second = 0;
            _tmp_node_iedge_count.clear();
        } else {
            if(bloc.valid()) {
                bwops_ += 2;
                if(linearFM) {
                    rangeTemp = gfm.mapLF(tloc, bloc, c, &node_rangeTemp);
                } else {
                    rangeTemp = gfm.mapGLF(tloc, bloc, c, &node_rangeTemp, &_tmp_node_iedge_count, (index_t)rp.kseeds);
                }
            } else {
                bwops_++;
                rangeTemp = gfm.mapGLF1(range.first, tloc, c, &node_rangeTemp);
                if(rangeTemp.first + 1 < rangeTemp.second) {
                    assert_eq(node_rangeTemp.first + 1, node_rangeTemp.second);
                    _tmp_node_iedge_count.clear();
                    _tmp_node_iedge_count.expand();
                    _tmp_node_iedge_count.back().first = 0;
                    _tmp_node_iedge_count.back().second = rangeTemp.second - rangeTemp.first - 1;
                }
            }
        }
        if(rangeTemp.first >= rangeTemp.second) {
            break;
        }
        
        range = rangeTemp;
        node_range = node_rangeTemp;
        if(_tmp_node_iedge_count.size() > 0) {
            node_iedge_count = _tmp_node_iedge_count;
            _tmp_node_iedge_count.clear();
        } else {
            node_iedge_count.clear();
        }
        dep++;
        
        if(uniqueStop_) {
            if(range.second - range.first == 1 && dep - offset >= _minK) {
                uniqueStop = true;
                break;
            }
        }
        
        HIER_INIT_LOCS(range.first, range.second, tloc, bloc, gfm);
    }
    
    // Done
    if(node_range.first < node_range.second && node_range.second - node_range.first <= rp.kseeds) {
        assert_leq(node_range.second - node_range.first, range.second - range.first);
#ifndef NDEBUG
        if(node_range.second - node_range.first < range.second - range.first) {
            ASSERT_ONLY(index_t add = 0);
            for(index_t e = 0; e < node_iedge_count.size(); e++) {
                if(e > 0) {
                    assert_lt(node_iedge_count[e-1].first, node_iedge_count[e].first);
                }
                assert_gt(node_iedge_count[e].second, 0);
                add += node_iedge_count[e].second;
            }
            assert_eq(node_range.second - node_range.first + add, range.second - range.first);
        } else {
            assert(node_iedge_count.empty());
        }
#endif
        assert_gt(dep, offset);
        assert_leq(dep, len);
        top = range.first; bot = range.second;
        node_top = node_range.first; node_bot = node_range.second;
        nelt += (node_bot - node_top);
        hitlen = dep - offset;
    }
    return (index_t)nelt;
}


/**
 *
 **/
template <typename index_t, typename local_index_t>
index_t HI_Aligner<index_t, local_index_t>::localGFMSearch(
                                                           const LocalGFM<local_index_t, index_t>&  gfm,  // GFM index
                                                           const Read&                      read,    // read to align
                                                           const Scoring&                   sc,      // scoring scheme
                                                           const ReportingParams&           rp,
                                                           bool                             fw,
                                                           index_t                          rdoff,
                                                           index_t&                         hitlen,
                                                           local_index_t&                   top,
                                                           local_index_t&                   bot,
                                                           local_index_t&                   node_top,
                                                           local_index_t&                   node_bot,
                                                           EList<pair<local_index_t, local_index_t> >& local_node_iedge_count,
                                                           RandomSource&                    rnd,
                                                           bool&                            uniqueStop,
                                                           local_index_t                    minUniqueLen,
                                                           local_index_t                    maxHitLen,
                                                           local_index_t                    maxHits)
{
    maxHits = max<local_index_t>(maxHits, rp.kseeds);
    bool uniqueStop_ = uniqueStop;
    uniqueStop = false;
    const local_index_t ftabLen = (local_index_t)gfm.gh().ftabChars();
    const bool linearFM = gfm.gh().linearFM();
	SideLocus<local_index_t> tloc, bloc;
	const local_index_t len = (local_index_t)read.length();
	size_t nelt = 0;
    
    const BTDnaString& seq = fw ? read.patFw : read.patRc;
    assert(!seq.empty());
    
    local_index_t offset = len - rdoff - 1;
    local_index_t dep = offset;
    pair<local_index_t, local_index_t> range(0, 0);
    pair<local_index_t, local_index_t> rangeTemp(0, 0);
    pair<local_index_t, local_index_t> node_range(0, 0);
    pair<local_index_t, local_index_t> node_rangeTemp(0, 0);
    top = bot = node_top = node_bot = 0;
    local_node_iedge_count.clear();
    _tmp_local_node_iedge_count.clear();
    local_index_t left = len - dep;
    assert_gt(left, 0);
    if(left < ftabLen + 1) {
        hitlen = left;
		return 0;
    }
    // Does N interfere with use of Ftab?
    for(local_index_t i = 0; i < ftabLen; i++) {
        int c = seq[len-dep-1-i];
        if(c > 3) {
            hitlen = i + 1;
			return 0;
        }
    }
    
    gfm.ftabLoHi(seq, len - dep - ftabLen, false, range.first, range.second);
    dep += ftabLen;
    if(range.first >= range.second) {
        hitlen = ftabLen;
        return 0;
    }
    LOCAL_INIT_LOCS(range.first, range.second, tloc, bloc, gfm);
    // Keep going
    while(dep < len) {
        int c = seq[len-dep-1];
        if(c > 3) {
            rangeTemp.first = rangeTemp.second = 0;
            node_rangeTemp.first = node_rangeTemp.second = 0;
            _tmp_local_node_iedge_count.clear();
        } else {
            if(bloc.valid()) {
                bwops_ += 2;
                if(linearFM) {
                    rangeTemp = gfm.mapLF(tloc, bloc, c, &node_rangeTemp);
                } else {
                    rangeTemp = gfm.mapGLF(tloc, bloc, c, &node_rangeTemp, &_tmp_local_node_iedge_count, rp.kseeds);
                }
            } else {
                bwops_++;
                rangeTemp = gfm.mapGLF1(range.first, tloc, c, &node_rangeTemp);
                if(rangeTemp.first + 1 < rangeTemp.second) {
                    assert_eq(node_rangeTemp.first + 1, node_rangeTemp.second);
                    _tmp_local_node_iedge_count.clear();
                    _tmp_local_node_iedge_count.expand();
                    _tmp_local_node_iedge_count.back().first = 0;
                    _tmp_local_node_iedge_count.back().second = rangeTemp.second - rangeTemp.first - 1;
                }
            }
        }
        if(rangeTemp.first >= rangeTemp.second) {
            break;
        }
        
        range = rangeTemp;
        node_range = node_rangeTemp;
        if(_tmp_local_node_iedge_count.size() > 0) {
            local_node_iedge_count = _tmp_local_node_iedge_count;
            _tmp_local_node_iedge_count.clear();
        } else {
            local_node_iedge_count.clear();
        }
        dep++;

        if(uniqueStop_) {
            if(range.second - range.first == 1 && dep - offset >= minUniqueLen) {
                uniqueStop = true;
                break;
            }
        }
        
        if(dep - offset >= maxHitLen) break;
        LOCAL_INIT_LOCS(range.first, range.second, tloc, bloc, gfm);
    }
    
    // Done
    if(node_range.first < node_range.second && node_range.second - node_range.first <= maxHits) {
        assert_leq(node_range.second - node_range.first, range.second - range.first);
#ifndef NDEBUG
        if(node_range.second - node_range.first < range.second - range.first) {
            ASSERT_ONLY(index_t add = 0);
            for(index_t e = 0; e < local_node_iedge_count.size(); e++) {
                if(e > 0) {
                    assert_lt(local_node_iedge_count[e-1].first, local_node_iedge_count[e].first);
                }
                assert_gt(local_node_iedge_count[e].second, 0);
                add += local_node_iedge_count[e].second;
            }
            assert_eq(node_range.second - node_range.first + add, range.second - range.first);
        } else {
            assert(local_node_iedge_count.empty());
        }
#endif
        assert_gt(dep, offset);
        assert_leq(dep, len);
        top = range.first; bot = range.second;
        node_top = node_range.first; node_bot = node_range.second;
        nelt += (node_bot - node_top);
        hitlen = dep - offset;
    }

    return (index_t)nelt;
}

/**
 *
 **/
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::isSearched(
                                                    const GenomeHit<index_t>&   hit,
                                                    index_t                     rdi)
{
    assert_lt(rdi, 2);
    EList<GenomeHit<index_t> >& searchedHits = _hits_searched[rdi];
    for(index_t i = 0; i < searchedHits.size(); i++) {
        if(searchedHits[i].contains(hit)) return true;
    }
    return false;
}

/**
 *
 **/
template <typename index_t, typename local_index_t>
void HI_Aligner<index_t, local_index_t>::addSearched(
                                                     const GenomeHit<index_t>&   hit,
                                                     index_t                     rdi)
{
    assert_lt(rdi, 2);
    assert(!isSearched(hit, rdi));
    EList<GenomeHit<index_t> >& searchedHits = _hits_searched[rdi];
    searchedHits.push_back(hit);
}

#endif /*HI_ALIGNER_H_*/
