/*
 * Copyright 2013, Daehwan Kim <infphilo@gmail.com>
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

#ifndef SPLICE_SITE_H_
#define SPLICE_SITE_H_

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <limits>
#include "assert_helpers.h"
#include "mem_ids.h"
#include "ref_coord.h"
#include "ds.h"
#include "read.h"
#include "reference.h"
#include "hier_idx_common.h"

using namespace std;

// #define NEW_PROB_MODEL
// the following parameters are borrowed from Yeo and Burge 2004 in Journal of Computational Biology
const size_t donor_exonic_len = 3;
const size_t donor_intronic_len = 6;
const size_t donor_len = donor_exonic_len + donor_intronic_len;

#if defined(NEW_PROB_MODEL)
const size_t acceptor_intronic_len = 20;
const size_t acceptor_exonic_len = 3;
const size_t acceptor_len = acceptor_intronic_len + acceptor_exonic_len;
#else
// the following parameters are borrowed from Ch 3. in Bioinformatics - From Genomes to Drugs. Volume I: Basic Technologies (Victor Solovyev)
const size_t acceptor_intronic_len = 14;
const size_t acceptor_exonic_len = 1;
const size_t acceptor_len = acceptor_intronic_len + acceptor_exonic_len;
const size_t acceptor_len1 = acceptor_len / 2;
const size_t acceptor_len2 = acceptor_len - acceptor_len1;

// the following parameters are borrowed from Yeo and Burge 2004 in Journal of Computational Biology
const float background_prob[4] = {0.27f, 0.23f, 0.23f, 0.27f};

extern float donor_prob[4][donor_len];
extern float acceptor_prob[4][acceptor_len];

extern float donor_prob_sum[1 << (donor_len << 1)];
extern float acceptor_prob_sum1[1 << (acceptor_len1 << 1)];
extern float acceptor_prob_sum2[1 << (acceptor_len2 << 1)];

#endif
const size_t intronic_len = max<size_t>(donor_intronic_len, acceptor_intronic_len);

extern void init_junction_prob();

/**
 *
 */
class SpliceSitePos {
    
public:
    
	SpliceSitePos() { reset(); }
    
	SpliceSitePos(const SpliceSitePos& c) { init(c); }
	
	SpliceSitePos(uint32_t ref, uint32_t left, uint32_t right, bool fw, bool canonical) { init(ref, left, right, fw, canonical); }
    
	/**
	 * Copy given fields into this Coord.
	 */
	virtual void init(uint32_t ref, uint32_t left, uint32_t right, bool fw, bool canonical) {
		_ref = ref;
		_left = left;
        _right = right;
		_fw = fw;
        _canonical = canonical;
	}
    
	/**
	 * Copy contents of given Coord into this one.
	 */
	void init(const SpliceSitePos& c) {
		_ref = c._ref;
		_left = c._left;
		_right = c._right;
        _fw = c._fw;
        _canonical = c._canonical;
	}
	
	/**
	 * Return true iff this Coord is identical to the given Coord.
	 */
	bool operator==(const SpliceSitePos& o) const {
		assert(inited());
		assert(o.inited());
		return _ref == o._ref &&
        _left == o._left &&
        _right == o._right &&
        _fw == o._fw &&
        _canonical == o._canonical;
	}
    
	/**
	 * Return true iff this Coord is less than the given Coord.  One Coord is
	 * less than another if (a) its reference id is less, (b) its orientation is
	 * less, or (c) its offset is less.
	 */
	bool operator<(const SpliceSitePos& o) const {
		if(_ref < o._ref) return true;
		if(_ref > o._ref) return false;
		if(_left < o._left) return true;
		if(_left > o._left) return false;
        if(_right < o._right) return true;
		if(_right > o._right) return false;
        if(_fw != o._fw) return _fw;
        if(_canonical != o._canonical) return _canonical;
		return false;
	}
	
	/**
	 * Return the opposite result from operator<.
	 */
	bool operator>=(const SpliceSitePos& o) const {
		return !((*this) < o);
	}
	
	/**
	 * Return true iff this Coord is greater than the given Coord.  One Coord
	 * is greater than another if (a) its reference id is greater, (b) its
	 * orientation is greater, or (c) its offset is greater.
	 */
	bool operator>(const SpliceSitePos& o) const {
		if(_ref > o._ref) return true;
		if(_ref < o._ref) return false;
		if(_left > o._left) return true;
		if(_left < o._left) return false;
        if(_right > o._right) return true;
		if(_right < o._right) return false;
        if(_fw != o._fw) return !_fw;
        if(_canonical != o._canonical) return !_canonical;
		return false;
	}
	
	/**
	 * Return the opposite result from operator>.
	 */
	bool operator<=(const SpliceSitePos& o) const {
		return !((*this) > o);
	}
	
	/**
	 * Reset this coord to uninitialized state.
	 */
	virtual void reset() {
		_ref = std::numeric_limits<uint32_t>::max();
		_left = std::numeric_limits<uint32_t>::max();
        _right = std::numeric_limits<uint32_t>::max();
        _fw = true;
        _canonical = true;
	}
	
	/**
	 * Return true iff this Coord is initialized (i.e. ref and off have both
	 * been set since the last call to reset()).
	 */
	bool inited() const {
		if(_ref != std::numeric_limits<uint32_t>::max() &&
		   _left != std::numeric_limits<uint32_t>::max() &&
           _right != std::numeric_limits<uint32_t>::max())
		{
			return true;
		}
		return false;
	}
	
#ifndef NDEBUG
	/**
	 * Check that coord is internally consistent.
	 */
	bool repOk() const {
		if(_ref == std::numeric_limits<uint32_t>::max() ||
		   _left == std::numeric_limits<uint32_t>::max() ||
           _right == std::numeric_limits<uint32_t>::max()) {
            return false;
		}
		return true;
	}
#endif
	
	/**
	 * Check whether an interval defined by this coord and having
	 * length 'len' is contained within an interval defined by
	 * 'inbegin' and 'inend'.
	 */
#if 0
	bool within(int64_t len, int64_t inbegin, int64_t inend) const {
		return off_ >= inbegin && off_ + len <= inend;
	}
#endif
	
	uint32_t ref()          const { return _ref; }
	uint32_t left()         const { return _left; }
    uint32_t right()        const { return _right; }
	bool     fw()           const { return _fw; }
    bool     canonical()    const { return _canonical; }
    uint32_t intron_len()   const { return _right - _left - 1; }
    
protected:
    
	uint32_t  _ref;             // which reference?
	uint32_t  _left;            // 0-based offset of the right most base of the left flanking exon
    uint32_t  _right;           // 0-based offset of the left most base of the right flanking exon
	bool      _fw;              // true -> Watson strand
    bool      _canonical;
};

/**
 *
 */
class SpliceSite : public SpliceSitePos {

public:

	SpliceSite() { reset(); }

	SpliceSite(const SpliceSite& c) { init(c); }
	
	SpliceSite(uint32_t ref, uint32_t left, uint32_t right, bool fw, bool canonical) { init(ref, left, right, fw, canonical); }

	/**
	 * Copy given fields into this Coord.
	 */
	virtual void init(uint32_t ref, uint32_t left, uint32_t right, bool fw, bool canonical) {
        SpliceSitePos::init(ref, left, right, fw, canonical);
    
        // _donordint = 0;
        // _acceptordint = 0;
        // _leftseq = 0;
        // _rightseq = 0;
        _leftext = 0;
        _rightext = 0;
        _numreads = 0;
        _editdist = 0;
        // _probscore = 0.0f;
        _readid = 0;
        _fromfile = false;
	}

	/**
	 * Copy contents of given Coord into this one.
	 */
	void init(const SpliceSite& c) {
        SpliceSitePos::init(c);
        
        // _donordint = c._donordint;
        // _acceptordint = c._acceptordint;
        // _leftseq = c._leftseq;
        // _rightseq = c._rightseq;
        _leftext = c._leftext;
        _rightext = c._rightext;
        _numreads = c._numreads;
        _editdist = c._editdist;
        // _probscore = c._probscore;
        _readid = c._readid;
        _fromfile = c._fromfile;
	}
	
	/**
	 * Reset this coord to uninitialized state.
	 */
	virtual void reset() {
        SpliceSitePos::reset();
        
        // _donordint = 0;
        // _acceptordint = 0;
        // _leftseq = 0;
        // _rightseq = 0;
        _leftext = 0;
        _rightext = 0;
        _numreads = 0;
        _editdist = 0;
        // _probscore = 0.0;
        _readid = 0;
        _fromfile = false;
	}
	
    // uint8_t  donordint()    const { return _donordint; }
    // uint8_t  acceptordint() const { return _acceptordint; }
    // uint64_t leftseq()      const { return _leftseq; }
    // uint64_t rightseq()     const { return _rightseq; }
    
    uint32_t leftext()      const { return _leftext; }
    uint32_t rightext()     const { return _rightext; }
    
    uint64_t numreads()     const { return _numreads; }
    uint32_t editdist()     const { return _editdist; }
    // float    probscore()    const { return _probscore; }

public:

    // uint8_t   _donordint;       // 3' dinucleotide on Watson strand
    // uint8_t   _acceptordint;    // 5' dinucleotide on Watson strand
    // uint64_t  _leftseq;         // left 32bp flanking seq on Watson strand
    // uint64_t  _rightseq;        // right 32bp flanking seq on Watson strand
    
    uint32_t  _leftext;
    uint32_t  _rightext;
    
    uint64_t  _numreads;        // number of reads spanning this splice site
    
    uint32_t  _editdist;
    
    // float     _probscore;
    
    uint64_t _readid;
    bool     _fromfile;
};

std::ostream& operator<<(std::ostream& out, const SpliceSite& c);

class AlnRes;

// #define SPLIT_DB

#if defined(SPLIT_DB)

class SpliceSiteDB {
public:
    typedef RedBlackNode<SpliceSitePos, uint32_t> Node;
    
public:
    SpliceSiteDB(
                 const BitPairReference& refs,
                 const EList<string>& refnames,
                 bool threadSafe = true,
                 bool write = false,
                 bool read = false);
    ~SpliceSiteDB();
    
    bool addSpliceSite(
                       const Read& rd,
                       const AlnRes& rs,
                       uint32_t minAnchorLen = 15);
    
   static float probscore(
                          int64_t donor_seq,
                          int64_t acceptor_seq);

#if 0
    size_t size(uint64_t ref) const;
    bool empty(uint64_t ref) const;
#endif
    
    bool empty() { return _empty; }
    
    bool write() const { return _write; }
    bool read() const { return _read; }
    
    bool getSpliceSite(SpliceSite& ss) const;
    void getLeftSpliceSites(uint32_t ref, uint32_t left, uint32_t range, EList<SpliceSite>& spliceSites) const;
    void getRightSpliceSites(uint32_t ref, uint32_t right, uint32_t range, EList<SpliceSite>& spliceSites) const;

    void print(ofstream& out);
    void read(ifstream& in, bool novel = true);
    
private:
    void getSpliceSites_recur(
                              const RedBlackNode<SpliceSitePos, uint32_t> *node,
			      const EList<SpliceSite>& spliceSites_db,
                              uint32_t ref,
                              uint32_t left,
                              uint32_t right,
                              EList<SpliceSite>& spliceSites) const;
    
    const RedBlackNode<SpliceSitePos, uint32_t>* getSpliceSite_temp(const SpliceSitePos& ssp) const;
    
    void print_recur(
                     const RedBlackNode<SpliceSitePos, uint32_t> *node,
                     ofstream& out,
                     EList<SpliceSitePos>& ss_list);
    
    Pool& pool(uint64_t ref, uint64_t ref_seg);
    
    void print_impl(
                    ofstream& out,
                    EList<SpliceSitePos>& ss_list,
                    const SpliceSitePos* ss = NULL);
    
private:
    uint64_t                            _numRefs;
    EList<string>                       _refnames;
    
    ELList<RedBlack<SpliceSitePos, uint32_t>* >   _fwIndex;
    ELList<RedBlack<SpliceSitePos, uint32_t>* >   _bwIndex;
    
    ELLList<SpliceSite>                  _spliceSites;
    
    ELLList<Pool*>                       _pool;   // dispenses memory pages
    
    bool                                _write;
    bool                                _read;
    
    ELList<MUTEX_T>                     _mutex;
    bool                                _threadSafe;
    
    SStringExpandable<char>             raw_refbuf;
    ASSERT_ONLY(SStringExpandable<uint32_t> destU32);
    BTDnaString                         donorstr;
    BTDnaString                         acceptorstr;
    
    bool                                _empty;
};

#else

class SpliceSiteDB {
public:
    typedef RedBlackNode<SpliceSitePos, uint32_t> Node;
    
public:
    SpliceSiteDB(
                 const BitPairReference& refs,
                 const EList<string>& refnames,
                 bool threadSafe = true,
                 bool write = false,
                 bool read = false);
    ~SpliceSiteDB();
    
    bool addSpliceSite(
                       const Read& rd,
                       const AlnRes& rs,
                       uint32_t minAnchorLen = 15);
    
    static float probscore(
                           int64_t donor_seq,
                           int64_t acceptor_seq);
    
    size_t size(uint64_t ref) const;
    bool empty(uint64_t ref) const;
    
    bool empty() { return _empty; }
    
    bool write() const { return _write; }
    bool read() const { return _read; }
    
    bool getSpliceSite(SpliceSite& ss) const;
    void getLeftSpliceSites(uint32_t ref, uint32_t left, uint32_t range, EList<SpliceSite>& spliceSites) const;
    void getRightSpliceSites(uint32_t ref, uint32_t right, uint32_t range, EList<SpliceSite>& spliceSites) const;
    
    void print(ofstream& out);
    void read(ifstream& in, bool novel = true);
    
private:
    void getSpliceSites_recur(
                              const RedBlackNode<SpliceSitePos, uint32_t> *node,
                              uint32_t ref,
                              uint32_t left,
                              uint32_t right,
                              EList<SpliceSite>& spliceSites) const;
    
    const RedBlackNode<SpliceSitePos, uint32_t>* getSpliceSite_temp(const SpliceSitePos& ssp) const;
    
    void print_recur(
                     const RedBlackNode<SpliceSitePos, uint32_t> *node,
                     ofstream& out,
                     const uint32_t numreads_cutoff,
                     const uint32_t numreads_cutoff2,
                     EList<SpliceSite>& ss_list);
    
    Pool& pool(uint64_t ref);
    
    void print_impl(
                    ofstream& out,
                    EList<SpliceSite>& ss_list,
                    const SpliceSite* ss = NULL);
    
private:
    uint64_t                            _numRefs;
    EList<string>                       _refnames;
    
    EList<RedBlack<SpliceSitePos, uint32_t>* >   _fwIndex;
    EList<RedBlack<SpliceSitePos, uint32_t>* >   _bwIndex;
    
    ELList<SpliceSite>                   _spliceSites;
    
    ELList<Pool*>                        _pool;   // dispenses memory pages
    
    bool                                _write;
    bool                                _read;
    
    EList<MUTEX_T>                      _mutex;
    bool                                _threadSafe;
    
    SStringExpandable<char>             raw_refbuf;
    ASSERT_ONLY(SStringExpandable<uint32_t> destU32);
    BTDnaString                         donorstr;
    BTDnaString                         acceptorstr;
    
    bool                                _empty;
};

#endif

#endif /*ifndef SPLICE_SITE_H_*/
