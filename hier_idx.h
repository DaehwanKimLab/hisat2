/*
 * Copyright 2013, Daehwan Kim <infphilo@gmail.com>
 *
 * This file is part of Beast.  Beast is based on Bowtie 2.
 *
 * Beast is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Beast is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Beast.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HIEREBWT_H_
#define HIEREBWT_H_

#include "hier_idx_common.h"
#include "bt2_idx.h"
#include "bt2_io.h"
#include "bt2_util.h"

/**
 * Extended Burrows-Wheeler transform data.
 * LocalEbwt is a specialized Ebwt index that represents ~64K bps
 * and therefore uses two bytes as offsets within 64K bps.
 * This class has only two additional member variables to denote the genomic sequenuce it represents:
 * (1) the contig index and (2) the offset within the contig.
 *
 */
template <typename index_t = uint16_t, typename full_index_t = uint32_t>
class LocalEbwt : public Ebwt<index_t> {
	typedef Ebwt<index_t> PARENT_CLASS;
public:
	/// Construct an Ebwt from the given input file
	LocalEbwt(const string& in,
			  FILE *in5,
			  FILE *in6,
			  char *mmFile5,
			  char *mmFile6,
			  full_index_t& tidx,
			  full_index_t& localOffset,
			  bool switchEndian,
			  size_t& bytesRead,
			  int color,
			  int needEntireReverse,
			  bool fw,
			  int32_t overrideOffRate, // = -1,
			  int32_t offRatePlus, // = -1,
			  uint32_t lineRate,
			  uint32_t offRate,
			  uint32_t ftabChars,
			  bool useMm, // = false,
			  bool useShmem, // = false,
			  bool mmSweep, // = false,
			  bool loadNames, // = false,
			  bool loadSASamp, // = true,
			  bool loadFtab, // = true,
			  bool loadRstarts, // = true,
			  bool verbose, // = false,
			  bool startVerbose, // = false,
			  bool passMemExc, // = false,
			  bool sanityCheck) : // = false) :
	Ebwt<index_t>(in,
				  color,
				  needEntireReverse,
				  fw,
				  overrideOffRate,
				  offRatePlus,
				  useMm,
				  useShmem,
				  mmSweep,
				  loadNames,
				  loadSASamp,
				  loadFtab,
				  loadRstarts,
				  verbose,
				  startVerbose,
				  passMemExc,
				  sanityCheck,
				  true)
	{
		this->_in1Str = in + ".5." + gEbwt_ext;
		this->_in2Str = in + ".5." + gEbwt_ext;
		readIntoMemory(
					   in5,
					   in6,
					   mmFile5,
					   mmFile6,
					   tidx,
					   localOffset,
					   switchEndian,
					   bytesRead,
					   color,
					   needEntireReverse,
					   loadSASamp,
					   loadFtab,
					   loadRstarts,
					   false,              //justHeader
					   lineRate,
					   offRate,
					   ftabChars,
					   mmSweep,
					   loadNames,
					   startVerbose);
		
		_tidx = tidx;
		_localOffset = localOffset;
		
		// If the offRate has been overridden, reflect that in the
		// _eh._offRate field
		if(offRatePlus > 0 && this->_overrideOffRate == -1) {
			this->_overrideOffRate = this->_eh._offRate + offRatePlus;
		}
		if(this->_overrideOffRate > this->_eh._offRate) {
			this->_eh.setOffRate(this->_overrideOffRate);
			assert_eq(this->_overrideOffRate, this->_eh._offRate);
		}
		assert(this->repOk());
	}


	/// Construct an Ebwt from the given header parameters and string
	/// vector, optionally using a blockwise suffix sorter with the
	/// given 'bmax' and 'dcv' parameters.  The string vector is
	/// ultimately joined and the joined string is passed to buildToDisk().
	template<typename TStr>
	LocalEbwt(
			  TStr& s,
			  full_index_t tidx,
			  full_index_t local_offset,
			  index_t local_size,
			  bool packed,
			  int color,
			  int needEntireReverse,
			  int32_t lineRate,
			  int32_t offRate,
			  int32_t ftabChars,
			  const string& file,   // base filename for EBWT files
			  bool fw,
			  int dcv,
			  EList<RefRecord>& szs,
			  index_t sztot,
			  const RefReadInParams& refparams,
			  uint32_t seed,
			  ostream& out5,
			  ostream& out6,
			  int32_t overrideOffRate = -1,
			  bool verbose = false,
			  bool passMemExc = false,
			  bool sanityCheck = false) :
	Ebwt<index_t>(packed,
				  color,
				  needEntireReverse,
				  lineRate,
				  offRate,
				  ftabChars,
				  file,
				  fw,
				  dcv,
				  szs,
				  sztot,
				  refparams,
				  seed,
				  overrideOffRate,
				  verbose,
				  passMemExc,
				  sanityCheck)
	{
		const EbwtParams<index_t>& eh = this->_eh;
		assert(eh.repOk());
		uint32_t be = this->toBe();
		assert(out5.good());
		assert(out6.good());
		writeIndex<full_index_t>(out5, tidx, be);
		writeIndex<full_index_t>(out5, local_offset, be);
		writeU32(out5, eh._len,      be); // length of string (and bwt and suffix array)
		if(eh._len > 0) {
			assert_gt(szs.size(), 0);
			assert_gt(sztot, 0);
			// Not every fragment represents a distinct sequence - many
			// fragments may correspond to a single sequence.  Count the
			// number of sequences here by counting the number of "first"
			// fragments.
			this->_nPat = 0;
			this->_nFrag = 0;
			for(size_t i = 0; i < szs.size(); i++) {
				if(szs[i].len > 0) this->_nFrag++;
				if(szs[i].first && szs[i].len > 0) this->_nPat++;
			}
			assert_eq(this->_nPat, 1);
			assert_geq(this->_nFrag, this->_nPat);
			this->_rstarts.reset();
			writeIndex(out5, this->_nPat, be);
			assert_eq(this->_nPat, 1);
			this->_plen.init(new index_t[this->_nPat], this->_nPat);
			// For each pattern, set plen
			int npat = -1;
			for(size_t i = 0; i < szs.size(); i++) {
				if(szs[i].first && szs[i].len > 0) {
					if(npat >= 0) {
						writeIndex(out5, this->plen()[npat], be);
					}
					npat++;
					this->plen()[npat] = (szs[i].len + szs[i].off);
				} else {
					this->plen()[npat] += (szs[i].len + szs[i].off);
				}
			}
			assert_eq((index_t)npat, this->_nPat-1);
			writeIndex(out5, this->plen()[npat], be);
			// Write the number of fragments
			writeIndex(out5, this->_nFrag, be);
			
			if(refparams.reverse == REF_READ_REVERSE) {
				EList<RefRecord> tmp(EBWT_CAT);
                reverseRefRecords(szs, tmp, false, verbose);
				this->szsToDisk(tmp, out5, refparams.reverse);
			} else {
				this->szsToDisk(szs, out5, refparams.reverse);
			}
			
			VMSG_NL("Constructing suffix-array element generator");
			KarkkainenBlockwiseSA<TStr> bsa(s, s.length()+1, dcv, seed, this->_sanity, this->_passMemExc, this->_verbose);
			assert(bsa.suffixItrIsReset());
			assert_eq(bsa.size(), s.length()+1);
			VMSG_NL("Converting suffix-array elements to index image");
			buildToDisk(bsa, s, out5, out6);
		}
		
		out5.flush(); out6.flush();
		if(out5.fail() || out6.fail()) {
			cerr << "An error occurred writing the index to disk.  Please check if the disk is full." << endl;
			throw 1;
		}
	}
	
	template <typename TStr> void buildToDisk(
											  InorderBlockwiseSA<TStr>& sa,
											  const TStr& s,
											  ostream& out1, 
											  ostream& out2);
	
	// I/O
	void readIntoMemory(
						FILE *in5,
						FILE *in6,
						char *mmFile5,
						char *mmFile6,
						full_index_t& tidx,
						full_index_t& localOffset,
						bool switchEndian,
						size_t bytesRead,
						int color,
						int needEntireRev, 
						bool loadSASamp, 
						bool loadFtab,
						bool loadRstarts, 
						bool justHeader, 
						int32_t lineRate,
						int32_t offRate,
						int32_t ftabChars,
						bool mmSweep, 
						bool loadNames, 
						bool startVerbose);
	
	/**
	 * Sanity-check various pieces of the Ebwt
	 */
	void sanityCheckAll(int reverse) const {
		if(this->_eh._len > 0) {
			PARENT_CLASS::sanityCheckAll(reverse);
		}
	}
    
    bool empty() const { return this->_eh._len == 0; }
	
public:
	full_index_t _tidx;
	full_index_t _localOffset;
};

/**
 * Build an Ebwt from a string 's' and its suffix array 'sa' (which
 * might actually be a suffix array *builder* that builds blocks of the
 * array on demand).  The bulk of the Ebwt, i.e. the ebwt and offs
 * arrays, is written directly to disk.  This is by design: keeping
 * those arrays in memory needlessly increases the footprint of the
 * building process.  Instead, we prefer to build the Ebwt directly
 * "to disk" and then read it back into memory later as necessary.
 *
 * It is assumed that the header values and join-related values (nPat,
 * plen) have already been written to 'out1' before this function
 * is called.  When this function is finished, it will have
 * additionally written ebwt, zOff, fchr, ftab and eftab to the primary
 * file and offs to the secondary file.
 *
 * Assume DNA/RNA/any alphabet with 4 or fewer elements.
 * Assume occ array entries are 32 bits each.
 *
 * @param sa            the suffix array to convert to a Ebwt
 * @param s             the original string
 * @param out
 */
template <typename index_t, typename full_index_t>
template <typename TStr>
void LocalEbwt<index_t, full_index_t>::buildToDisk(
									 InorderBlockwiseSA<TStr>& sa,
									 const TStr& s,
									 ostream& out5,
									 ostream& out6)
{
	assert_leq(s.length(), std::numeric_limits<index_t>::max());
	const EbwtParams<index_t>& eh = this->_eh;
	
	assert(eh.repOk());
	assert_eq(s.length()+1, sa.size());
	assert_eq(s.length(), eh._len);
	assert_gt(eh._lineRate, 3);
	assert(sa.suffixItrIsReset());
	
	index_t len = eh._len;
	index_t ftabLen = eh._ftabLen;
	index_t sideSz = eh._sideSz;
	index_t ebwtTotSz = eh._ebwtTotSz;
	index_t fchr[] = {0, 0, 0, 0, 0};
	EList<index_t> ftab(EBWT_CAT);
	index_t zOff = (index_t)OFF_MASK;
	
	// Save # of occurrences of each character as we walk along the bwt
	index_t occ[4] = {0, 0, 0, 0};
	index_t occSave[4] = {0, 0, 0, 0};
	
	// Record rows that should "absorb" adjacent rows in the ftab.
	// The absorbed rows represent suffixes shorter than the ftabChars
	// cutoff.
	uint8_t absorbCnt = 0;
	EList<uint8_t> absorbFtab(EBWT_CAT);
	try {
		VMSG_NL("Allocating ftab, absorbFtab");
		ftab.resize(ftabLen);
		ftab.fillZero();
		absorbFtab.resize(ftabLen);
		absorbFtab.fillZero();
	} catch(bad_alloc &e) {
		cerr << "Out of memory allocating ftab[] or absorbFtab[] "
		<< "in Ebwt::buildToDisk() at " << __FILE__ << ":"
		<< __LINE__ << endl;
		throw e;
	}
	
	// Allocate the side buffer; holds a single side as its being
	// constructed and then written to disk.  Reused across all sides.
#ifdef SIXTY4_FORMAT
	EList<uint64_t> ebwtSide(EBWT_CAT);
#else
	EList<uint8_t> ebwtSide(EBWT_CAT);
#endif
	try {
#ifdef SIXTY4_FORMAT
		ebwtSide.resize(sideSz >> 3);
#else
		ebwtSide.resize(sideSz);
#endif
	} catch(bad_alloc &e) {
		cerr << "Out of memory allocating ebwtSide[] in "
		<< "Ebwt::buildToDisk() at " << __FILE__ << ":"
		<< __LINE__ << endl;
		throw e;
	}
	
	// Points to the base offset within ebwt for the side currently
	// being written
	index_t side = 0;
	
	// Whether we're assembling a forward or a reverse bucket
	bool fw;
	int sideCur = 0;
	fw = true;
	
	// Have we skipped the '$' in the last column yet?
	ASSERT_ONLY(bool dollarSkipped = false);

	index_t si = 0;   // string offset (chars)
	ASSERT_ONLY(uint32_t lastSufInt = 0);
	ASSERT_ONLY(bool inSA = true); // true iff saI still points inside suffix
	// array (as opposed to the padding at the
	// end)
	// Iterate over packed bwt bytes
	VMSG_NL("Entering Ebwt loop");
	ASSERT_ONLY(uint32_t beforeEbwtOff = (uint32_t)out5.tellp());
	while(side < ebwtTotSz) {
		// Sanity-check our cursor into the side buffer
		assert_geq(sideCur, 0);
		assert_lt(sideCur, (int)eh._sideBwtSz);
		assert_eq(0, side % sideSz); // 'side' must be on side boundary
		ebwtSide[sideCur] = 0; // clear
		assert_lt(side + sideCur, ebwtTotSz);
		// Iterate over bit-pairs in the si'th character of the BWT
#ifdef SIXTY4_FORMAT
		for(int bpi = 0; bpi < 32; bpi++, si++) {
#else
		for(int bpi = 0; bpi < 4; bpi++, si++) {
#endif
			int bwtChar;
			bool count = true;
			if(si <= len) {
				// Still in the SA; extract the bwtChar
				index_t saElt = (index_t)sa.nextSuffix();
				// (that might have triggered sa to calc next suf block)
				if(saElt == 0) {
					// Don't add the '$' in the last column to the BWT
					// transform; we can't encode a $ (only A C T or G)
					// and counting it as, say, an A, will mess up the
					// LR mapping
					bwtChar = 0; count = false;
					ASSERT_ONLY(dollarSkipped = true);
					zOff = si; // remember the SA row that
					// corresponds to the 0th suffix
				} else {
					bwtChar = (int)(s[saElt-1]);
					assert_lt(bwtChar, 4);
					// Update the fchr
					fchr[bwtChar]++;
				}
				// Update ftab
				if((len-saElt) >= (index_t)eh._ftabChars) {
					// Turn the first ftabChars characters of the
					// suffix into an integer index into ftab.  The
					// leftmost (lowest index) character of the suffix
					// goes in the most significant bit pair if the
					// integer.
					uint32_t sufInt = 0;
					for(int i = 0; i < eh._ftabChars; i++) {
						sufInt <<= 2;
						assert_lt((index_t)i, len-saElt);
						sufInt |= (unsigned char)(s[saElt+i]);
					}
					// Assert that this prefix-of-suffix is greater
					// than or equal to the last one (true b/c the
					// suffix array is sorted)
#ifndef NDEBUG
					if(lastSufInt > 0) assert_geq(sufInt, lastSufInt);
					lastSufInt = sufInt;
#endif
					// Update ftab
					assert_lt(sufInt+1, ftabLen);
					ftab[sufInt+1]++;
					if(absorbCnt > 0) {
						// Absorb all short suffixes since the last
						// transition into this transition
						absorbFtab[sufInt] = absorbCnt;
						absorbCnt = 0;
					}
				} else {
					// Otherwise if suffix is fewer than ftabChars
					// characters long, then add it to the 'absorbCnt';
					// it will be absorbed into the next transition
					assert_lt(absorbCnt, 255);
					absorbCnt++;
				}
				// Suffix array offset boundary? - update offset array
				if((si & eh._offMask) == si) {
					assert_lt((si >> eh._offRate), eh._offsLen);
					// Write offsets directly to the secondary output
					// stream, thereby avoiding keeping them in memory
					writeIndex(out6, saElt, this->toBe());
				}
			} else {
				// Strayed off the end of the SA, now we're just
				// padding out a bucket
#ifndef NDEBUG
				if(inSA) {
					// Assert that we wrote all the characters in the
					// string before now
					assert_eq(si, len+1);
					inSA = false;
				}
#endif
				// 'A' used for padding; important that padding be
				// counted in the occ[] array
				bwtChar = 0;
			}
			if(count) occ[bwtChar]++;
			// Append BWT char to bwt section of current side
			if(fw) {
				// Forward bucket: fill from least to most
#ifdef SIXTY4_FORMAT
				ebwtSide[sideCur] |= ((uint64_t)bwtChar << (bpi << 1));
				if(bwtChar > 0) assert_gt(ebwtSide[sideCur], 0);
#else
				pack_2b_in_8b(bwtChar, ebwtSide[sideCur], bpi);
				assert_eq((ebwtSide[sideCur] >> (bpi*2)) & 3, bwtChar);
#endif
			} else {
				// Backward bucket: fill from most to least
#ifdef SIXTY4_FORMAT
				ebwtSide[sideCur] |= ((uint64_t)bwtChar << ((31 - bpi) << 1));
				if(bwtChar > 0) assert_gt(ebwtSide[sideCur], 0);
#else
				pack_2b_in_8b(bwtChar, ebwtSide[sideCur], 3-bpi);
				assert_eq((ebwtSide[sideCur] >> ((3-bpi)*2)) & 3, bwtChar);
#endif
			}
		} // end loop over bit-pairs
		assert_eq(dollarSkipped ? 3 : 0, (occ[0] + occ[1] + occ[2] + occ[3]) & 3);
#ifdef SIXTY4_FORMAT
		assert_eq(0, si & 31);
#else
		assert_eq(0, si & 3);
#endif
		
		sideCur++;
		if(sideCur == (int)eh._sideBwtSz) {
			sideCur = 0;
			index_t *uside = reinterpret_cast<index_t*>(ebwtSide.ptr());
			// Write 'A', 'C', 'G' and 'T' tallies
			side += sideSz;
			assert_leq(side, eh._ebwtTotSz);
			uside[(sideSz / sizeof(index_t))-4] = endianizeIndex(occSave[0], this->toBe());
			uside[(sideSz / sizeof(index_t))-3] = endianizeIndex(occSave[1], this->toBe());
			uside[(sideSz / sizeof(index_t))-2] = endianizeIndex(occSave[2], this->toBe());
			uside[(sideSz / sizeof(index_t))-1] = endianizeIndex(occSave[3], this->toBe());
			occSave[0] = occ[0];
			occSave[1] = occ[1];
			occSave[2] = occ[2];
			occSave[3] = occ[3];
			// Write backward side to primary file
			out5.write((const char *)ebwtSide.ptr(), sideSz);
		}
	}
	VMSG_NL("Exited Ebwt loop");
	assert_neq(zOff, (index_t)OFF_MASK);
	if(absorbCnt > 0) {
		// Absorb any trailing, as-yet-unabsorbed short suffixes into
		// the last element of ftab
		absorbFtab[ftabLen-1] = absorbCnt;
	}
	// Assert that our loop counter got incremented right to the end
	assert_eq(side, eh._ebwtTotSz);
	// Assert that we wrote the expected amount to out1
	assert_eq(((uint32_t)out5.tellp() - beforeEbwtOff), eh._ebwtTotSz);
	// assert that the last thing we did was write a forward bucket
	
	//
	// Write zOff to primary stream
	//
	writeIndex(out5, zOff, this->toBe());
	
	//
	// Finish building fchr
	//
	// Exclusive prefix sum on fchr
	for(int i = 1; i < 4; i++) {
		fchr[i] += fchr[i-1];
	}
	assert_eq(fchr[3], len);
	// Shift everybody up by one
	for(int i = 4; i >= 1; i--) {
		fchr[i] = fchr[i-1];
	}
	fchr[0] = 0;
	if(this->_verbose) {
		for(int i = 0; i < 5; i++)
			cout << "fchr[" << "ACGT$"[i] << "]: " << fchr[i] << endl;
	}
	// Write fchr to primary file
	for(int i = 0; i < 5; i++) {
		writeIndex(out5, fchr[i], this->toBe());
	}
	
	//
	// Finish building ftab and build eftab
	//
	// Prefix sum on ftable
	index_t eftabLen = 0;
	assert_eq(0, absorbFtab[0]);
	for(index_t i = 1; i < ftabLen; i++) {
		if(absorbFtab[i] > 0) eftabLen += 2;
	}
	assert_leq(eftabLen, (index_t)eh._ftabChars*2);
	eftabLen = eh._ftabChars*2;
	EList<index_t> eftab(EBWT_CAT);
	try {
		eftab.resize(eftabLen);
		eftab.fillZero();
	} catch(bad_alloc &e) {
		cerr << "Out of memory allocating eftab[] "
		<< "in Ebwt::buildToDisk() at " << __FILE__ << ":"
		<< __LINE__ << endl;
		throw e;
	}
	index_t eftabCur = 0;
	for(index_t i = 1; i < ftabLen; i++) {
		index_t lo = ftab[i] + Ebwt<index_t>::ftabHi(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, i-1);
		if(absorbFtab[i] > 0) {
			// Skip a number of short pattern indicated by absorbFtab[i]
			index_t hi = lo + absorbFtab[i];
			assert_lt(eftabCur*2+1, eftabLen);
			eftab[eftabCur*2] = lo;
			eftab[eftabCur*2+1] = hi;
			ftab[i] = (eftabCur++) ^ (index_t)OFF_MASK; // insert pointer into eftab
			assert_eq(lo, Ebwt<index_t>::ftabLo(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, i));
			assert_eq(hi, Ebwt<index_t>::ftabHi(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, i));
		} else {
			ftab[i] = lo;
		}
	}
	assert_eq(Ebwt<index_t>::ftabHi(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, ftabLen-1), len+1);
	// Write ftab to primary file
	for(index_t i = 0; i < ftabLen; i++) {
		writeIndex(out5, ftab[i], this->toBe());
	}
	// Write eftab to primary file
	for(index_t i = 0; i < eftabLen; i++) {
		writeIndex(out5, eftab[i], this->toBe());
	}
	
	// Note: if you'd like to sanity-check the Ebwt, you'll have to
	// read it back into memory first!
	assert(!this->isInMemory());
	VMSG_NL("Exiting Ebwt::buildToDisk()");
}

/**
 * Read an Ebwt from file with given filename.
 */
template <typename index_t, typename full_index_t>
void LocalEbwt<index_t, full_index_t>::readIntoMemory(
										FILE *in5,
										FILE *in6,
										char *mmFile5,
										char *mmFile6,
										full_index_t& tidx,
										full_index_t& localOffset,
										bool switchEndian,
										size_t bytesRead,
										int color,
										int entireRev,
										bool loadSASamp,
										bool loadFtab,
										bool loadRstarts,
										bool justHeader,
										int32_t lineRate,
										int32_t offRate,
										int32_t ftabChars,
										bool mmSweep,
										bool loadNames,
										bool startVerbose)
{
#ifdef BOWTIE_MM
	char *mmFile[] = { mmFile5, mmFile6 };
#endif
	
	// Reads header entries one by one from primary stream
	tidx = readIndex<full_index_t>(in5, switchEndian); bytesRead += sizeof(full_index_t);
	localOffset = readIndex<full_index_t>(in5, switchEndian); bytesRead += sizeof(full_index_t);
	uint32_t len = readU32(in5, switchEndian); bytesRead += 4;
	
	// Create a new EbwtParams from the entries read from primary stream
	this->_eh.init(len, lineRate, offRate, ftabChars, color, entireRev);
	
	if(len <= 0) {
		return;
	}
	
	// Set up overridden suffix-array-sample parameters
	uint32_t offsLen = this->_eh._offsLen;
	uint32_t offRateDiff = 0;
	uint32_t offsLenSampled = offsLen;
	if(this->_overrideOffRate > offRate) {
		offRateDiff = this->_overrideOffRate - offRate;
	}
	if(offRateDiff > 0) {
		offsLenSampled >>= offRateDiff;
		if((offsLen & ~((index_t)OFF_MASK << offRateDiff)) != 0) {
			offsLenSampled++;
		}
	}
	
	// Can't override the offrate or isarate and use memory-mapped
	// files; ultimately, all processes need to copy the sparser sample
	// into their own memory spaces.
	if(this->_useMm && (offRateDiff)) {
		cerr << "Error: Can't use memory-mapped files when the offrate is overridden" << endl;
		throw 1;
	}
	
	// Read nPat from primary stream
	this->_nPat = readIndex<index_t>(in5, switchEndian);
	assert_eq(this->_nPat, 1);
	bytesRead += sizeof(index_t);
	this->_plen.reset();
	
	// Read plen from primary stream
	if(this->_useMm) {
#ifdef BOWTIE_MM
		this->_plen.init((index_t*)(mmFile[0] + bytesRead), this->_nPat, false);
		bytesRead += this->_nPat*sizeof(index_t);
		fseek(in5, this->_nPat*sizeof(index_t), SEEK_CUR);
#endif
	} else {
		try {
			if(this->_verbose || startVerbose) {
				cerr << "Reading plen (" << this->_nPat << "): ";
				logTime(cerr);
			}
			this->_plen.init(new index_t[this->_nPat], this->_nPat, true);
			if(switchEndian) {
				for(index_t i = 0; i < this->_nPat; i++) {
					this->plen()[i] = readIndex<index_t>(in5, switchEndian);
				}
			} else {
				size_t r = MM_READ(in5, (void*)(this->plen()), this->_nPat*sizeof(index_t));
				if(r != (size_t)(this->_nPat*sizeof(index_t))) {
					cerr << "Error reading _plen[] array: " << r << ", " << this->_nPat*sizeof(index_t) << endl;
					throw 1;
				}
			}
		} catch(bad_alloc& e) {
			cerr << "Out of memory allocating plen[] in Ebwt::read()"
			<< " at " << __FILE__ << ":" << __LINE__ << endl;
			throw e;
		}
	}

	bool shmemLeader;
	
	// TODO: I'm not consistent on what "header" means.  Here I'm using
	// "header" to mean everything that would exist in memory if we
	// started to build the Ebwt but stopped short of the build*() step
	// (i.e. everything up to and including join()).
	if(justHeader) return;
	
	this->_nFrag = readIndex<index_t>(in5, switchEndian);
	bytesRead += sizeof(index_t);
	if(this->_verbose || startVerbose) {
		cerr << "Reading rstarts (" << this->_nFrag*3 << "): ";
		logTime(cerr);
	}
	assert_geq(this->_nFrag, this->_nPat);
	this->_rstarts.reset();
	if(loadRstarts) {
		if(this->_useMm) {
#ifdef BOWTIE_MM
			this->_rstarts.init((index_t*)(mmFile[0] + bytesRead), this->_nFrag*3, false);
			bytesRead += this->_nFrag*sizeof(index_t)*3;
			fseek(in5, this->_nFrag*sizeof(index_t)*3, SEEK_CUR);
#endif
		} else {
			this->_rstarts.init(new index_t[this->_nFrag*3], this->_nFrag*3, true);
			if(switchEndian) {
				for(index_t i = 0; i < this->_nFrag*3; i += 3) {
					// fragment starting position in joined reference
					// string, text id, and fragment offset within text
					this->rstarts()[i]   = readIndex<index_t>(in5, switchEndian);
					this->rstarts()[i+1] = readIndex<index_t>(in5, switchEndian);
					this->rstarts()[i+2] = readIndex<index_t>(in5, switchEndian);
				}
			} else {
				size_t r = MM_READ(in5, (void *)this->rstarts(), this->_nFrag*sizeof(index_t)*3);
				if(r != (size_t)(this->_nFrag*sizeof(index_t)*3)) {
					cerr << "Error reading _rstarts[] array: " << r << ", " << (this->_nFrag*sizeof(index_t)*3) << endl;
					throw 1;
				}
			}
		}
	} else {
		// Skip em
		assert(this->rstarts() == NULL);
		bytesRead += this->_nFrag*sizeof(index_t)*3;
		fseek(in5, this->_nFrag*sizeof(index_t)*3, SEEK_CUR);
	}
	
	this->_ebwt.reset();
	if(this->_useMm) {
#ifdef BOWTIE_MM
		this->_ebwt.init((uint8_t*)(mmFile[0] + bytesRead), this->_eh._ebwtTotLen, false);
		bytesRead += this->_eh._ebwtTotLen;
		fseek(in5, this->_eh._ebwtTotLen, SEEK_CUR);
#endif
	} else {
		// Allocate ebwt (big allocation)
		if(this->_verbose || startVerbose) {
			cerr << "Reading ebwt (" << this->_eh._ebwtTotLen << "): ";
			logTime(cerr);
		}
		bool shmemLeader = true;
		if(this->useShmem_) {
			uint8_t *tmp = NULL;
			shmemLeader = ALLOC_SHARED_U8(
										  (this->_in1Str + "[ebwt]"), this->_eh._ebwtTotLen, &tmp,
										  "ebwt[]", (this->_verbose || startVerbose));
			assert(tmp != NULL);
			this->_ebwt.init(tmp, this->_eh._ebwtTotLen, false);
			if(this->_verbose || startVerbose) {
				cerr << "  shared-mem " << (shmemLeader ? "leader" : "follower") << endl;
			}
		} else {
			try {
				this->_ebwt.init(new uint8_t[this->_eh._ebwtTotLen], this->_eh._ebwtTotLen, true);
			} catch(bad_alloc& e) {
				cerr << "Out of memory allocating the ebwt[] array for the Bowtie index.  Please try" << endl
				<< "again on a computer with more memory." << endl;
				throw 1;
			}
		}
		if(shmemLeader) {
			// Read ebwt from primary stream
			uint64_t bytesLeft = this->_eh._ebwtTotLen;
			char *pebwt = (char*)this->ebwt();
            
			while (bytesLeft>0){
				size_t r = MM_READ(in5, (void *)pebwt, bytesLeft);
				if(MM_IS_IO_ERR(in5, r, bytesLeft)) {
					cerr << "Error reading _ebwt[] array: " << r << ", "
                    << bytesLeft << endl;
					throw 1;
				}
				pebwt += r;
				bytesLeft -= r;
			}
			if(switchEndian) {
				uint8_t *side = this->ebwt();
				for(size_t i = 0; i < this->_eh._numSides; i++) {
					index_t *cums = reinterpret_cast<index_t*>(side + this->_eh._sideSz - sizeof(index_t)*2);
					cums[0] = endianSwapIndex(cums[0]);
					cums[1] = endianSwapIndex(cums[1]);
					side += this->_eh._sideSz;
				}
			}
#ifdef BOWTIE_SHARED_MEM
			if(useShmem_) NOTIFY_SHARED(this->ebwt(), this->_eh._ebwtTotLen);
#endif
		} else {
			// Seek past the data and wait until master is finished
			fseek(in5, this->_eh._ebwtTotLen, SEEK_CUR);
#ifdef BOWTIE_SHARED_MEM
			if(useShmem_) WAIT_SHARED(this->ebwt(), this->_eh._ebwtTotLen);
#endif
		}
	}
	
	// Read zOff from primary stream
	this->_zOff = readIndex<index_t>(in5, switchEndian);
	bytesRead += sizeof(index_t);
	assert_lt(this->_zOff, len);
	
	try {
		// Read fchr from primary stream
		if(this->_verbose || startVerbose) cerr << "Reading fchr (5)" << endl;
		this->_fchr.reset();
		if(this->_useMm) {
#ifdef BOWTIE_MM
			this->_fchr.init((index_t*)(mmFile[0] + bytesRead), 5, false);
			bytesRead += 5*sizeof(index_t);
			fseek(in5, 5*sizeof(index_t), SEEK_CUR);
#endif
		} else {
			this->_fchr.init(new index_t[5], 5, true);
			for(index_t i = 0; i < 5; i++) {
				this->fchr()[i] = readIndex<index_t>(in5, switchEndian);
				assert_leq(this->fchr()[i], len);
				assert(i <= 0 || this->fchr()[i] >= this->fchr()[i-1]);
			}
		}
		assert_gt(this->fchr()[4], this->fchr()[0]);
		// Read ftab from primary stream
		if(this->_verbose || startVerbose) {
			if(loadFtab) {
				cerr << "Reading ftab (" << this->_eh._ftabLen << "): ";
				logTime(cerr);
			} else {
				cerr << "Skipping ftab (" << this->_eh._ftabLen << "): ";
			}
		}
		this->_ftab.reset();
		if(loadFtab) {
			if(this->_useMm) {
#ifdef BOWTIE_MM
				this->_ftab.init((index_t*)(mmFile[0] + bytesRead), this->_eh._ftabLen, false);
				bytesRead += this->_eh._ftabLen*sizeof(index_t);
				fseek(in5, this->_eh._ftabLen*sizeof(index_t), SEEK_CUR);
#endif
			} else {
				this->_ftab.init(new index_t[this->_eh._ftabLen], this->_eh._ftabLen, true);
				if(switchEndian) {
					for(uint32_t i = 0; i < this->_eh._ftabLen; i++)
						this->ftab()[i] = readIndex<index_t>(in5, switchEndian);
				} else {
					size_t r = MM_READ(in5, (void *)this->ftab(), this->_eh._ftabLen*sizeof(index_t));
					if(r != (size_t)(this->_eh._ftabLen*sizeof(index_t))) {
						cerr << "Error reading _ftab[] array: " << r << ", " << (this->_eh._ftabLen*sizeof(index_t)) << endl;
						throw 1;
					}
				}
			}
			// Read etab from primary stream
			if(this->_verbose || startVerbose) {
				if(loadFtab) {
					cerr << "Reading eftab (" << this->_eh._eftabLen << "): ";
					logTime(cerr);
				} else {
					cerr << "Skipping eftab (" << this->_eh._eftabLen << "): ";
				}
				
			}
			this->_eftab.reset();
			if(this->_useMm) {
#ifdef BOWTIE_MM
				this->_eftab.init((index_t*)(mmFile[0] + bytesRead), this->_eh._eftabLen, false);
				bytesRead += this->_eh._eftabLen*sizeof(index_t);
				fseek(in5, this->_eh._eftabLen*sizeof(index_t), SEEK_CUR);
#endif
			} else {
				this->_eftab.init(new index_t[this->_eh._eftabLen], this->_eh._eftabLen, true);
				if(switchEndian) {
					for(uint32_t i = 0; i < this->_eh._eftabLen; i++)
						this->eftab()[i] = readIndex<index_t>(in5, switchEndian);
				} else {
					size_t r = MM_READ(in5, (void *)this->eftab(), this->_eh._eftabLen*sizeof(index_t));
					if(r != (size_t)(this->_eh._eftabLen*sizeof(index_t))) {
						cerr << "Error reading _eftab[] array: " << r << ", " << (this->_eh._eftabLen*sizeof(index_t)) << endl;
						throw 1;
					}
				}
			}
			for(uint32_t i = 0; i < this->_eh._eftabLen; i++) {
				if(i > 0 && this->eftab()[i] > 0) {
					assert_geq(this->eftab()[i], this->eftab()[i-1]);
				} else if(i > 0 && this->eftab()[i-1] == 0) {
					assert_eq(0, this->eftab()[i]);
				}
			}
		} else {
			assert(this->ftab() == NULL);
			assert(this->eftab() == NULL);
			// Skip ftab
			bytesRead += this->_eh._ftabLen*sizeof(index_t);
			fseek(in5, this->_eh._ftabLen*sizeof(index_t), SEEK_CUR);
			// Skip eftab
			bytesRead += this->_eh._eftabLen*sizeof(index_t);
			fseek(in5, this->_eh._eftabLen*sizeof(index_t), SEEK_CUR);
		}
	} catch(bad_alloc& e) {
		cerr << "Out of memory allocating fchr[], ftab[] or eftab[] arrays for the Bowtie index." << endl
		<< "Please try again on a computer with more memory." << endl;
		throw 1;
	}
	
	this->_offs.reset();
	if(loadSASamp) {
		bytesRead = 4; // reset for secondary index file (already read 1-sentinel)		
		shmemLeader = true;
		if(this->_verbose || startVerbose) {
			cerr << "Reading offs (" << offsLenSampled << " " << std::setw(2) << sizeof(index_t)*8 << "-bit words): ";
			logTime(cerr);
		}
		
		if(!this->_useMm) {
			if(!this->useShmem_) {
				// Allocate offs_
				try {
					this->_offs.init(new index_t[offsLenSampled], offsLenSampled, true);
				} catch(bad_alloc& e) {
					cerr << "Out of memory allocating the offs[] array  for the Bowtie index." << endl
					<< "Please try again on a computer with more memory." << endl;
					throw 1;
				}
			} else {
				index_t *tmp = NULL;
				shmemLeader = ALLOC_SHARED_U32(
											   (this->_in2Str + "[offs]"), offsLenSampled*2, &tmp,
											   "offs", (this->_verbose || startVerbose));
				this->_offs.init((index_t*)tmp, offsLenSampled, false);
			}
		}
		
		if(this->_overrideOffRate < 32) {
			if(shmemLeader) {
				// Allocate offs (big allocation)
				if(switchEndian || offRateDiff > 0) {
					assert(!this->_useMm);
					const uint32_t blockMaxSz = (2 * 1024 * 1024); // 2 MB block size
					const uint32_t blockMaxSzUIndex = (blockMaxSz / sizeof(index_t)); // # UIndexs per block
					char *buf;
					try {
						buf = new char[blockMaxSz];
					} catch(std::bad_alloc& e) {
						cerr << "Error: Out of memory allocating part of _offs array: '" << e.what() << "'" << endl;
						throw e;
					}
					for(index_t i = 0; i < offsLen; i += blockMaxSzUIndex) {
					  index_t block = min<index_t>((index_t)blockMaxSzUIndex, (index_t)(offsLen - i));
						size_t r = MM_READ(in6, (void *)buf, block * sizeof(index_t));
						if(r != (size_t)(block * sizeof(index_t))) {
							cerr << "Error reading block of _offs[] array: " << r << ", " << (block * sizeof(index_t)) << endl;
							throw 1;
						}
						index_t idx = i >> offRateDiff;
						for(index_t j = 0; j < block; j += (1 << offRateDiff)) {
							assert_lt(idx, offsLenSampled);
							this->offs()[idx] = ((index_t*)buf)[j];
							if(switchEndian) {
								this->offs()[idx] = endianSwapIndex(this->offs()[idx]);
							}
							idx++;
						}
					}
					delete[] buf;
				} else {
					if(this->_useMm) {
#ifdef BOWTIE_MM
						this->_offs.init((index_t*)(mmFile[1] + bytesRead), offsLen, false);
						bytesRead += (offsLen * sizeof(index_t));
						fseek(in6, (offsLen * sizeof(index_t)), SEEK_CUR);
#endif
					} else {
						// If any of the high two bits are set
						if((offsLen & 0xc0000000) != 0) {
							if(sizeof(char *) <= 4) {
								cerr << "Sanity error: sizeof(char *) <= 4 but offsLen is " << hex << offsLen << endl;
								throw 1;
							}
							// offsLen << 2 overflows, so do it in four reads
							char *offs = (char *)this->offs();
							for(size_t i = 0; i < sizeof(index_t); i++) {
								size_t r = MM_READ(in6, (void*)offs, offsLen);
								if(r != (size_t)(offsLen)) {
									cerr << "Error reading block of _offs[] array: " << r << ", " << offsLen << endl;
									throw 1;
								}
								offs += offsLen;
							}
						} else {
							// Do it all in one read
							size_t r = MM_READ(in6, (void*)this->offs(), offsLen * sizeof(index_t));
							if(r != (size_t)(offsLen * sizeof(index_t))) {
								cerr << "Error reading _offs[] array: " << r << ", " << (offsLen * sizeof(index_t)) << endl;
								throw 1;
							}
						}
					}
				}
#ifdef BOWTIE_SHARED_MEM				
				if(this->useShmem_) NOTIFY_SHARED(this->offs(), offsLenSampled*sizeof(index_t));
#endif
			} else {
				// Not the shmem leader
				fseek(in6, offsLenSampled*sizeof(index_t), SEEK_CUR);
#ifdef BOWTIE_SHARED_MEM				
				if(this->useShmem_) WAIT_SHARED(this->offs(), offsLenSampled*sizeof(index_t));
#endif
			}
		}
	}
	
	this->postReadInit(this->_eh); // Initialize fields of Ebwt not read from file
	if(this->_verbose || startVerbose) this->print(cerr, this->_eh);
}

/**
 * Extended Burrows-Wheeler transform data.
 * HierEbwt is a specialized Ebwt index that represents one global index and a large set of local indexes.
 *
 */
template <typename index_t = uint32_t, typename local_index_t = uint16_t>
class HierEbwt : public Ebwt<index_t> {
	typedef Ebwt<index_t> PARENT_CLASS;
public:
	/// Construct an Ebwt from the given input file
	HierEbwt(const string& in,
			 int color,
			 int needEntireReverse,
			 bool fw,
			 int32_t overrideOffRate, // = -1,
			 int32_t offRatePlus, // = -1,
			 bool useMm, // = false,
			 bool useShmem, // = false,
			 bool mmSweep, // = false,
			 bool loadNames, // = false,
			 bool loadSASamp, // = true,
			 bool loadFtab, // = true,
			 bool loadRstarts, // = true,
			 bool verbose, // = false,
			 bool startVerbose, // = false,
			 bool passMemExc, // = false,
			 bool sanityCheck, // = false
             bool skipLoading = false) :
	         Ebwt<index_t>(in,
						   color,
						   needEntireReverse,
						   fw,
						   overrideOffRate,
						   offRatePlus,
						   useMm,
						   useShmem,
						   mmSweep,
						   loadNames,
						   loadSASamp,
						   loadFtab,
						   loadRstarts,
						   verbose,
						   startVerbose,
						   passMemExc,
						   sanityCheck,
						   skipLoading),
	         _in5(NULL),
	         _in6(NULL)
	{
		_in5Str = in + ".5." + gEbwt_ext;
		_in6Str = in + ".6." + gEbwt_ext;
        
        if(!skipLoading && false) {
            readIntoMemory(
                           color,       // expect index to be colorspace?
                           fw ? -1 : needEntireReverse, // need REF_READ_REVERSE
                           loadSASamp,  // load the SA sample portion?
                           loadFtab,    // load the ftab & eftab?
                           loadRstarts, // load the rstarts array?
                           true,        // stop after loading the header portion?
                           &(this->_eh),
                           mmSweep,     // mmSweep
                           loadNames,   // loadNames
                           startVerbose); // startVerbose
            // If the offRate has been overridden, reflect that in the
            // _eh._offRate field
            if(offRatePlus > 0 && this->_overrideOffRate == -1) {
                this->_overrideOffRate = this->_eh._offRate + offRatePlus;
            }
            if(this->_overrideOffRate > this->_eh._offRate) {
                this->_eh.setOffRate(this->_overrideOffRate);
                assert_eq(this->_overrideOffRate, this->_eh._offRate);
            }
            assert(this->repOk());
        }
	}
	
	/// Construct an Ebwt from the given header parameters and string
	/// vector, optionally using a blockwise suffix sorter with the
	/// given 'bmax' and 'dcv' parameters.  The string vector is
	/// ultimately joined and the joined string is passed to buildToDisk().
	template<typename TStr>
	HierEbwt(
			 TStr& s,
			 bool packed,
			 int color,
			 int needEntireReverse,
			 int32_t lineRate,
			 int32_t offRate,
			 int32_t ftabChars,
             int32_t localOffRate,
             int32_t localFtabChars,
			 const string& file,   // base filename for EBWT files
			 bool fw,
			 bool useBlockwise,
			 TIndexOffU bmax,
			 TIndexOffU bmaxSqrtMult,
			 TIndexOffU bmaxDivN,
			 int dcv,
			 EList<FileBuf*>& is,
			 EList<RefRecord>& szs,
			 index_t sztot,
			 const RefReadInParams& refparams,
			 uint32_t seed,
			 int32_t overrideOffRate = -1,
			 bool verbose = false,
			 bool passMemExc = false,
			 bool sanityCheck = false);
	        	
	~HierEbwt() {
		clearLocalEbwts();
	}
    
    /**
	 * Load this Ebwt into memory by reading it in from the _in1 and
	 * _in2 streams.
	 */
	void loadIntoMemory(
                        int color,
                        int needEntireReverse,
                        bool loadSASamp,
                        bool loadFtab,
                        bool loadRstarts,
                        bool loadNames,
                        bool verbose)
	{
		readIntoMemory(
                       color,       // expect index to be colorspace?
                       needEntireReverse, // require reverse index to be concatenated reference reversed
                       loadSASamp,  // load the SA sample portion?
                       loadFtab,    // load the ftab (_ftab[] and _eftab[])?
                       loadRstarts, // load the r-starts (_rstarts[])?
                       false,       // stop after loading the header portion?
                       NULL,        // params
                       false,       // mmSweep
                       loadNames,   // loadNames
                       verbose);    // startVerbose
	}
	
	// I/O
	void readIntoMemory(
                        int color,
                        int needEntireRev,
                        bool loadSASamp,
                        bool loadFtab,
                        bool loadRstarts,
                        bool justHeader,
                        EbwtParams<index_t> *params,
                        bool mmSweep,
                        bool loadNames,
                        bool startVerbose);
	
	/**
	 * Frees memory associated with the Ebwt.
	 */
	void evictFromMemory() {
		assert(PARENT_CLASS::isInMemory());
		clearLocalEbwts();
		PARENT_CLASS::evictFromMemory();		
	}
	
	/**
	 * Sanity-check various pieces of the Ebwt
	 */
	void sanityCheckAll(int reverse) const {
		PARENT_CLASS::sanityCheckAll(reverse);
		for(size_t tidx = 0; tidx < _localEbwts.size(); tidx++) {
			for(size_t local_idx = 0; local_idx < _localEbwts[tidx].size(); local_idx++) {
				assert(_localEbwts[tidx][local_idx] != NULL);
				_localEbwts[tidx][local_idx]->sanityCheckAll(reverse);
			}
		}
	}
    
    const LocalEbwt<local_index_t, index_t>* getLocalEbwt(index_t tidx, index_t offset) const {
        assert_lt(tidx, _localEbwts.size());
        const EList<LocalEbwt<local_index_t, index_t>*>& localEbwts = _localEbwts[tidx];
        index_t offsetidx = offset / local_index_interval;
        if(offsetidx >= localEbwts.size()) {
            return NULL;
        } else {
            return localEbwts[offsetidx];
        }
    }
    
    const LocalEbwt<local_index_t, index_t>* prevLocalEbwt(const LocalEbwt<local_index_t, index_t>* currLocalEbwt) const {
        assert(currLocalEbwt != NULL);
        index_t tidx = currLocalEbwt->_tidx;
        index_t offset = currLocalEbwt->_localOffset;
        if(offset < local_index_interval) {
            return NULL;
        } else {
            return getLocalEbwt(tidx, offset - local_index_interval);
        }
    }
    
    const LocalEbwt<local_index_t, index_t>* nextLocalEbwt(const LocalEbwt<local_index_t, index_t>* currLocalEbwt) const {
        assert(currLocalEbwt != NULL);
        index_t tidx = currLocalEbwt->_tidx;
        index_t offset = currLocalEbwt->_localOffset;
        return getLocalEbwt(tidx, offset + local_index_interval);
    }
	
	void clearLocalEbwts() {
		for(size_t tidx = 0; tidx < _localEbwts.size(); tidx++) {
			for(size_t local_idx = 0; local_idx < _localEbwts[tidx].size(); local_idx++) {
				assert(_localEbwts[tidx][local_idx] != NULL);
				delete _localEbwts[tidx][local_idx];
			}
			
			_localEbwts[tidx].clear();
		}
		
		_localEbwts.clear();
	}
	

public:
	index_t                                  _nrefs;      /// the number of reference sequences
	EList<index_t>                           _refLens;    /// approx lens of ref seqs (excludes trailing ambig chars)
	
	EList<EList<LocalEbwt<local_index_t, index_t>*> > _localEbwts;
	index_t                                  _nlocalEbwts;
	
	FILE                                     *_in5;    // input fd for primary index file
	FILE                                     *_in6;    // input fd for secondary index file
	string                                   _in5Str;
	string                                   _in6Str;
	
	char                                     *mmFile5_;
	char                                     *mmFile6_;
};
    
/// Construct an Ebwt from the given header parameters and string
/// vector, optionally using a blockwise suffix sorter with the
/// given 'bmax' and 'dcv' parameters.  The string vector is
/// ultimately joined and the joined string is passed to buildToDisk().
template <typename index_t, typename local_index_t>
template <typename TStr>
HierEbwt<index_t, local_index_t>::HierEbwt(
                                           TStr& s,
                                           bool packed,
                                           int color,
                                           int needEntireReverse,
                                           int32_t lineRate,
                                           int32_t offRate,
                                           int32_t ftabChars,
                                           int32_t localOffRate,
                                           int32_t localFtabChars,
                                           const string& file,   // base filename for EBWT files
                                           bool fw,
                                           bool useBlockwise,
                                           TIndexOffU bmax,
                                           TIndexOffU bmaxSqrtMult,
                                           TIndexOffU bmaxDivN,
                                           int dcv,
                                           EList<FileBuf*>& is,
                                           EList<RefRecord>& szs,
                                           index_t sztot,
                                           const RefReadInParams& refparams,
                                           uint32_t seed,
                                           int32_t overrideOffRate,
                                           bool verbose,
                                           bool passMemExc,
                                           bool sanityCheck) :
    Ebwt<index_t>(s,
                  packed,
                  color,
                  needEntireReverse,
                  lineRate,
                  offRate,
                  ftabChars,
                  file,
                  fw,
                  useBlockwise,
                  bmax,
                  bmaxSqrtMult,
                  bmaxDivN,
                  dcv,
                  is,
                  szs,
                  sztot,
                  refparams,
                  seed,
                  overrideOffRate,
                  verbose,
                  passMemExc,
                  sanityCheck),
    _in5(NULL),
    _in6(NULL)
{
    _in5Str = file + ".5." + gEbwt_ext;
    _in6Str = file + ".6." + gEbwt_ext;
    
    // Open output files
    ofstream fout5(_in5Str.c_str(), ios::binary);
    if(!fout5.good()) {
        cerr << "Could not open index file for writing: \"" << _in5Str.c_str() << "\"" << endl
        << "Please make sure the directory exists and that permissions allow writing by" << endl
        << "Bowtie." << endl;
        throw 1;
    }
    ofstream fout6(_in6Str.c_str(), ios::binary);
    if(!fout6.good()) {
        cerr << "Could not open index file for writing: \"" << _in6Str.c_str() << "\"" << endl
        << "Please make sure the directory exists and that permissions allow writing by" << endl
        << "Bowtie." << endl;
        throw 1;
    }
    
    // split the whole genome into a set of local indexes
    _nrefs = 0;
    _nlocalEbwts = 0;
    
    index_t cumlen = 0;
    typedef EList<RefRecord, 1> EList_RefRecord;
    EList<EList<EList_RefRecord> > all_local_recs;
    // For each unambiguous stretch...
    for(index_t i = 0; i < szs.size(); i++) {
        const RefRecord& rec = szs[i];
        if(rec.first) {
            if(_nrefs > 0) {
                // refLens_ links each reference sequence with the total number
                // of ambiguous and unambiguous characters in it.
                _refLens.push_back(cumlen);
            }
            cumlen = 0;
            _nrefs++;
            all_local_recs.expand();
            assert_eq(_nrefs, all_local_recs.size());
        } else if(i == 0) {
            cerr << "First record in reference index file was not marked as "
            << "'first'" << endl;
            throw 1;
        }
        
        assert_gt(_nrefs, 0);
        assert_eq(_nrefs, all_local_recs.size());
        EList<EList_RefRecord>& ref_local_recs = all_local_recs[_nrefs-1];
        index_t next_cumlen = cumlen + rec.off + rec.len;
        index_t local_off = (cumlen / local_index_interval) * local_index_interval;
        if(local_off >= local_index_interval) {
            local_off -= local_index_interval;
        }
        for(;local_off < next_cumlen; local_off += local_index_interval) {
            if(local_off + local_index_size < cumlen) {
                continue;
            }
            index_t local_idx = local_off / local_index_interval;
            
            if(local_idx >= ref_local_recs.size()) {
                assert_eq(local_idx, ref_local_recs.size());
                ref_local_recs.expand();
                _nlocalEbwts++;
            }
            assert_lt(local_idx, ref_local_recs.size());
            EList_RefRecord& local_recs = ref_local_recs[local_idx];
            assert_gt(local_off + local_index_size, cumlen);
            local_recs.expand();
            if(local_off + local_index_size <= cumlen + rec.off) {
                local_recs.back().off = local_off + local_index_size - std::max(local_off, cumlen);
                local_recs.back().len = 0;
            } else {
                if(local_off < cumlen + rec.off) {
                    local_recs.back().off = rec.off - (local_off > cumlen ? local_off - cumlen : 0);
                } else {
                    local_recs.back().off = 0;
                }
                local_recs.back().len = std::min(next_cumlen, local_off + local_index_size) - std::max(local_off, cumlen + rec.off);
            }
            local_recs.back().first = (local_recs.size() == 1);
        }
        cumlen = next_cumlen;
    }
    
    // Store a cap entry for the end of the last reference seq
    _refLens.push_back(cumlen);
    
#ifndef NDEBUG
    EList<RefRecord> temp_szs;
    index_t temp_sztot = 0;
    index_t temp_nlocalEbwts = 0;
    for(size_t tidx = 0; tidx < all_local_recs.size(); tidx++) {
        assert_lt(tidx, _refLens.size());
        EList<EList_RefRecord>& ref_local_recs = all_local_recs[tidx];
        assert_eq((_refLens[tidx] + local_index_interval - 1) / local_index_interval, ref_local_recs.size());
        temp_szs.expand();
        temp_szs.back().off = 0;
        temp_szs.back().len = 0;
        temp_szs.back().first = true;
        index_t temp_ref_len = 0;
        index_t temp_ref_sztot = 0;
        temp_nlocalEbwts += ref_local_recs.size();
        for(size_t i = 0; i < ref_local_recs.size(); i++) {
            EList_RefRecord& local_recs = ref_local_recs[i];
            index_t local_len = 0;
            for(size_t j = 0; j < local_recs.size(); j++) {
                assert(local_recs[j].off != 0 || local_recs[j].len != 0);
                assert(j != 0 || local_recs[j].first);
                RefRecord local_rec = local_recs[j];
                if(local_len < local_index_interval && local_recs[j].off > 0){
                    if(local_len + local_recs[j].off > local_index_interval) {
                        temp_ref_len += (local_index_interval - local_len);
                        local_rec.off = local_index_interval - local_len;
                    } else {
                        temp_ref_len += local_recs[j].off;
                    }
                } else {
                    local_rec.off = 0;
                }
                local_len += local_recs[j].off;
                if(local_len < local_index_interval && local_recs[j].len > 0) {
                    if(local_len + local_recs[j].len > local_index_interval) {
                        temp_ref_len += (local_index_interval - local_len);
                        temp_ref_sztot += (local_index_interval - local_len);
                        local_rec.len = local_index_interval - local_len;
                    } else {
                        temp_ref_len += local_recs[j].len;
                        temp_ref_sztot += local_recs[j].len;
                    }
                } else {
                    local_rec.len = 0;
                }
                local_len += local_recs[j].len;
                if(local_rec.off > 0) {
                    if(temp_szs.back().len > 0) {
                        temp_szs.expand();
                        temp_szs.back().off = local_rec.off;
                        temp_szs.back().len = local_rec.len;
                        temp_szs.back().first = false;
                    } else {
                        temp_szs.back().off += local_rec.off;
                        temp_szs.back().len = local_rec.len;
                    }
                } else if(local_rec.len > 0) {
                    temp_szs.back().len += local_rec.len;
                }
            }
            if(i + 1 < ref_local_recs.size()) {
                assert_eq(local_len, local_index_size);
                assert_eq(temp_ref_len % local_index_interval, 0);
            } else {
                assert_eq(local_len, _refLens[tidx] % local_index_interval);
            }
        }
        assert_eq(temp_ref_len, _refLens[tidx]);
        temp_sztot += temp_ref_sztot;
    }
    assert_eq(temp_sztot, sztot);
    for(size_t i = 0; i < temp_szs.size(); i++) {
        assert_lt(i, szs.size());
        assert_eq(temp_szs[i].off, szs[i].off);
        assert_eq(temp_szs[i].len, szs[i].len);
        assert_eq(temp_szs[i].first, szs[i].first);
    }
    assert_eq(temp_szs.size(), szs.size());
    assert_eq(_nlocalEbwts, temp_nlocalEbwts);
#endif
    
    uint32_t be = this->toBe();
    assert(fout5.good());
    assert(fout6.good());
    
    // When building an Ebwt, these header parameters are known
    // "up-front", i.e., they can be written to disk immediately,
    // before we join() or buildToDisk()
    writeI32(fout5, 1, be); // endian hint for priamry stream
    writeI32(fout6, 1, be); // endian hint for secondary stream
    writeIndex<index_t>(fout5, _nlocalEbwts, be); // number of local Ebwts
    writeI32(fout5, local_lineRate,  be); // 2^lineRate = size in bytes of 1 line
    writeI32(fout5, 2, be); // not used
    writeI32(fout5, (int32_t)localOffRate,   be); // every 2^offRate chars is "marked"
    writeI32(fout5, (int32_t)localFtabChars, be); // number of 2-bit chars used to address ftab
    int32_t flags = 1;
    if(this->_eh._color) flags |= EBWT_COLOR;
    if(this->_eh._entireReverse) flags |= EBWT_ENTIRE_REV;
    writeI32(fout5, -flags, be); // BTL: chunkRate is now deprecated
    
    // build local FM indexes
    index_t curr_sztot = 0;
    bool firstIndex = true;
    for(size_t tidx = 0; tidx < _refLens.size(); tidx++) {
        index_t refLen = _refLens[tidx];
        index_t local_offset = 0;
        _localEbwts.expand();
        assert_lt(tidx, _localEbwts.size());
        while(local_offset < refLen) {
            index_t index_size = std::min<index_t>(refLen - local_offset, local_index_size);
            assert_lt(tidx, all_local_recs.size());
            assert_lt(local_offset / local_index_interval, all_local_recs[tidx].size());
            EList_RefRecord& local_szs = all_local_recs[tidx][local_offset / local_index_interval];
            
            EList<RefRecord> conv_local_szs;
            index_t local_len = 0, local_sztot = 0, local_sztot_interval = 0;
            for(size_t i = 0; i < local_szs.size(); i++) {
                assert(local_szs[i].off != 0 || local_szs[i].len != 0);
                assert(i != 0 || local_szs[i].first);
                conv_local_szs.push_back(local_szs[i]);
                local_len += local_szs[i].off;
                if(local_len < local_index_interval && local_szs[i].len > 0) {
                    if(local_len + local_szs[i].len > local_index_interval) {
                        local_sztot_interval += (local_index_interval - local_len);
                    } else {
                        local_sztot_interval += local_szs[i].len;
                    }
                }
                local_sztot += local_szs[i].len;
                local_len += local_szs[i].len;
            }
            TStr local_s;
            local_s.resize(local_sztot);
            if(refparams.reverse == REF_READ_REVERSE) {
                local_s.install(s.buf() + s.length() - curr_sztot - local_sztot, local_sztot);
            } else {
                local_s.install(s.buf() + curr_sztot, local_sztot);
            }
            LocalEbwt<local_index_t, index_t>* localEbwt = new LocalEbwt<local_index_t, index_t>(
                                                                                                 local_s,
                                                                                                 tidx,
                                                                                                 local_offset,
                                                                                                 index_size,
                                                                                                 packed,
                                                                                                 color,
                                                                                                 needEntireReverse,
                                                                                                 local_lineRate,
                                                                                                 localOffRate,      // suffix-array sampling rate
                                                                                                 localFtabChars,    // number of chars in initial arrow-pair calc
                                                                                                 file,               // basename for .?.ebwt files
                                                                                                 fw,                 // fw
                                                                                                 dcv,                // difference-cover period
                                                                                                 conv_local_szs,     // list of reference sizes
                                                                                                 local_sztot,        // total size of all unambiguous ref chars
                                                                                                 refparams,          // reference read-in parameters
                                                                                                 seed,               // pseudo-random number generator seed
                                                                                                 fout5,
                                                                                                 fout6,
                                                                                                 -1,                 // override offRate
                                                                                                 false,              // be silent
                                                                                                 passMemExc,         // pass exceptions up to the toplevel so that we can adjust memory settings automatically
                                                                                                 sanityCheck);       // verify results and internal consistency
            firstIndex = false;
            _localEbwts[tidx].push_back(localEbwt);
            curr_sztot += local_sztot_interval;
            local_offset += local_index_interval;
        }
    }
    assert_eq(curr_sztot, sztot);
    
    
    fout5 << '\0';
    fout5.flush(); fout6.flush();
    if(fout5.fail() || fout6.fail()) {
        cerr << "An error occurred writing the index to disk.  Please check if the disk is full." << endl;
        throw 1;
    }
    VMSG_NL("Returning from initFromVector");
    
    // Close output files
    fout5.flush();
    int64_t tellpSz5 = (int64_t)fout5.tellp();
    VMSG_NL("Wrote " << fout5.tellp() << " bytes to primary EBWT file: " << _in5Str.c_str());
    fout5.close();
    bool err = false;
    if(tellpSz5 > fileSize(_in5Str.c_str())) {
        err = true;
        cerr << "Index is corrupt: File size for " << _in5Str.c_str() << " should have been " << tellpSz5
        << " but is actually " << fileSize(_in5Str.c_str()) << "." << endl;
    }
    fout6.flush();
    int64_t tellpSz6 = (int64_t)fout6.tellp();
    VMSG_NL("Wrote " << fout6.tellp() << " bytes to secondary EBWT file: " << _in6Str.c_str());
    fout6.close();
    if(tellpSz6 > fileSize(_in6Str.c_str())) {
        err = true;
        cerr << "Index is corrupt: File size for " << _in6Str.c_str() << " should have been " << tellpSz6
        << " but is actually " << fileSize(_in6Str.c_str()) << "." << endl;
    }
    if(err) {
        cerr << "Please check if there is a problem with the disk or if disk is full." << endl;
        throw 1;
    }
    // Reopen as input streams
    VMSG_NL("Re-opening _in5 and _in5 as input streams");
    if(this->_sanity) {
        VMSG_NL("Sanity-checking Bt2");
        assert(!this->isInMemory());
        readIntoMemory(
                       color,                       // colorspace?
                       fw ? -1 : needEntireReverse, // 1 -> need the reverse to be reverse-of-concat
                       true,                        // load SA sample (_offs[])?
                       true,                        // load ftab (_ftab[] & _eftab[])?
                       true,                        // load r-starts (_rstarts[])?
                       false,                       // just load header?
                       NULL,                        // Params object to fill
                       false,                       // mm sweep?
                       true,                        // load names?
                       false);                      // verbose startup?
        sanityCheckAll(refparams.reverse);
        evictFromMemory();
        assert(!this->isInMemory());
    }
    VMSG_NL("Returning from HierEbwt constructor");
}

    
/**
 * Read an Ebwt from file with given filename.
 */
template <typename index_t, typename local_index_t>
void HierEbwt<index_t, local_index_t>::readIntoMemory(
													  int color,
													  int needEntireRev,
													  bool loadSASamp,
													  bool loadFtab,
													  bool loadRstarts,
													  bool justHeader,
													  EbwtParams<index_t> *params,
													  bool mmSweep,
													  bool loadNames,
													  bool startVerbose)
{
    PARENT_CLASS::readIntoMemory(color,
                                 needEntireRev,
                                 loadSASamp,
                                 loadFtab,
                                 loadRstarts,
                                 justHeader || needEntireRev == 1,
                                 params,
                                 mmSweep,
                                 loadNames,
                                 startVerbose);

	bool switchEndian; // dummy; caller doesn't care
#ifdef BOWTIE_MM
	char *mmFile[] = { NULL, NULL };
#endif
	if(_in5Str.length() > 0) {
		if(this->_verbose || startVerbose) {
			cerr << "  About to open input files: ";
			logTime(cerr);
		}
        // Initialize our primary and secondary input-stream fields
		if(_in5 != NULL) fclose(_in5);
		if(this->_verbose || startVerbose) cerr << "Opening \"" << _in5Str.c_str() << "\"" << endl;
		if((_in5 = fopen(_in5Str.c_str(), "rb")) == NULL) {
			cerr << "Could not open index file " << _in5Str.c_str() << endl;
		}
		if(loadSASamp) {
			if(_in6 != NULL) fclose(_in6);
			if(this->_verbose || startVerbose) cerr << "Opening \"" << _in6Str.c_str() << "\"" << endl;
			if((_in6 = fopen(_in6Str.c_str(), "rb")) == NULL) {
				cerr << "Could not open index file " << _in6Str.c_str() << endl;
			}
		}
		if(this->_verbose || startVerbose) {
			cerr << "  Finished opening input files: ";
			logTime(cerr);
		}
		
#ifdef BOWTIE_MM
		if(this->_useMm /*&& !justHeader*/) {
			const char *names[] = {_in5Str.c_str(), _in6Str.c_str()};
            int fds[] = { fileno(_in5), fileno(_in6) };
			for(int i = 0; i < (loadSASamp ? 2 : 1); i++) {
				if(this->_verbose || startVerbose) {
					cerr << "  ¯ " << (i+1) << ": ";
					logTime(cerr);
				}
				struct stat sbuf;
				if (stat(names[i], &sbuf) == -1) {
					perror("stat");
					cerr << "Error: Could not stat index file " << names[i] << " prior to memory-mapping" << endl;
					throw 1;
				}
                mmFile[i] = (char*)mmap((void *)0, (size_t)sbuf.st_size,
										PROT_READ, MAP_SHARED, fds[(size_t)i], 0);
				if(mmFile[i] == (void *)(-1)) {
					perror("mmap");
					cerr << "Error: Could not memory-map the index file " << names[i] << endl;
					throw 1;
				}
				if(mmSweep) {
					int sum = 0;
					for(off_t j = 0; j < sbuf.st_size; j += 1024) {
						sum += (int) mmFile[i][j];
					}
					if(startVerbose) {
						cerr << "  Swept the memory-mapped ebwt index file 1; checksum: " << sum << ": ";
						logTime(cerr);
					}
				}
			}
			mmFile5_ = mmFile[0];
			mmFile6_ = loadSASamp ? mmFile[1] : NULL;
		}
#endif
	}
#ifdef BOWTIE_MM
	else if(this->_useMm && !justHeader) {
		mmFile[0] = mmFile5_;
		mmFile[1] = mmFile6_;
	}
	if(this->_useMm && !justHeader) {
		assert(mmFile[0] == mmFile5_);
		assert(mmFile[1] == mmFile6_);
	}
#endif
	
	if(this->_verbose || startVerbose) {
		cerr << "  Reading header: ";
		logTime(cerr);
	}
	
	// Read endianness hints from both streams
	size_t bytesRead = 0;
	switchEndian = false;
	uint32_t one = readU32(_in5, switchEndian); // 1st word of primary stream
	bytesRead += 4;
	if(loadSASamp) {
#ifndef NDEBUG
		assert_eq(one, readU32(_in6, switchEndian)); // should match!
#else
		readU32(_in6, switchEndian);
#endif
	}
	if(one != 1) {
		assert_eq((1u<<24), one);
		assert_eq(1, endianSwapU32(one));
		switchEndian = true;
	}
	
	// Can't switch endianness and use memory-mapped files; in order to
	// support this, someone has to modify the file to switch
	// endiannesses appropriately, and we can't do this inside Bowtie
	// or we might be setting up a race condition with other processes.
	if(switchEndian && this->_useMm) {
		cerr << "Error: Can't use memory-mapped files when the index is the opposite endianness" << endl;
		throw 1;
	}	
	
	_nlocalEbwts      = readIndex<index_t>(_in5, switchEndian); bytesRead += sizeof(index_t);
	int32_t lineRate  = readI32(_in5, switchEndian); bytesRead += 4;
	readI32(_in5, switchEndian); bytesRead += 4;
	int32_t offRate   = readI32(_in5, switchEndian); bytesRead += 4;
	// TODO: add isaRate to the actual file format (right now, the
	// user has to tell us whether there's an ISA sample and what the
	// sampling rate is.
	int32_t ftabChars = readI32(_in5, switchEndian); bytesRead += 4;
	/*int32_t flag  =*/ readI32(_in5, switchEndian); bytesRead += 4;
    
    if(this->_verbose || startVerbose) {
        cerr << "    number of local indexes: " << _nlocalEbwts << endl
             << "    local offRate: " << offRate << endl
             << "    local ftabLen: " << (1 << (2 * ftabChars)) << endl
             << "    local ftabSz: "  << (2 << (2 * ftabChars)) << endl
        ;
    }
	
	clearLocalEbwts();
	
	index_t tidx = 0, localOffset = 0;
	string base = "";
	for(size_t i = 0; i < _nlocalEbwts; i++) {
		LocalEbwt<local_index_t, index_t> *localEbwt = new LocalEbwt<local_index_t, index_t>(base,
                                                                                             _in5,
                                                                                             _in6,
                                                                                             mmFile5_,
                                                                                             mmFile6_,
                                                                                             tidx,
                                                                                             localOffset,
                                                                                             switchEndian,
                                                                                             bytesRead,
                                                                                             color,
                                                                                             needEntireRev,
                                                                                             this->fw_,
                                                                                             -1, // overrideOffRate
                                                                                             -1, // offRatePlus
                                                                                             (uint32_t)lineRate,
                                                                                             (uint32_t)offRate,
                                                                                             (uint32_t)ftabChars,
                                                                                             this->_useMm,
                                                                                             this->useShmem_,
                                                                                             mmSweep,
                                                                                             loadNames,
                                                                                             loadSASamp,
                                                                                             loadFtab,
                                                                                             loadRstarts,
                                                                                             false,  // _verbose
                                                                                             false,
                                                                                             this->_passMemExc,
                                                                                             this->_sanity);
		
		if(tidx >= _localEbwts.size()) {
			assert_eq(tidx, _localEbwts.size());
			_localEbwts.expand();
		}
		assert_eq(tidx + 1, _localEbwts.size());
		_localEbwts.back().push_back(localEbwt);
	}	
		
#ifdef BOWTIE_MM
    fseek(_in5, 0, SEEK_SET);
	fseek(_in6, 0, SEEK_SET);
#else
	rewind(_in5); rewind(_in6);
#endif
}

#endif /*HIEREBWT_H_*/
