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

#ifndef HGFM_H_
#define HGFM_H_

#include "hier_idx_common.h"
#include "gfm.h"

/**
 * Extended Burrows-Wheeler transform data.
 * LocalEbwt is a specialized Ebwt index that represents ~64K bps
 * and therefore uses two bytes as offsets within 64K bps.
 * This class has only two additional member variables to denote the genomic sequenuce it represents:
 * (1) the contig index and (2) the offset within the contig.
 *
 */
template <typename index_t = uint16_t, typename full_index_t = uint32_t>
class LocalGFM : public GFM<index_t> {
	typedef GFM<index_t> PARENT_CLASS;
public:
	/// Construct an Ebwt from the given input file
	LocalGFM(const string& in,
             ALTDB<index_t>* altdb,
             FILE *in5,
             FILE *in6,
             char *mmFile5,
             char *mmFile6,
             full_index_t& tidx,
             full_index_t& localOffset,
             full_index_t& joinedOffset,
             bool switchEndian,
             size_t& bytesRead,
             size_t& bytesRead2,
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
	GFM<index_t>(in,
                 altdb,
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
                 true, // load Splice Sites
                 verbose,
                 startVerbose,
                 passMemExc,
                 sanityCheck,
                 true)
	{
		this->_in1Str = in + ".5." + gfm_ext;
		this->_in2Str = in + ".5." + gfm_ext;
		readIntoMemory(
					   in5,
					   in6,
					   mmFile5,
					   mmFile6,
					   tidx,
					   localOffset,
                       joinedOffset,
					   switchEndian,
					   bytesRead,
                       bytesRead2,
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
        _joinedOffset = joinedOffset;
		
		// If the offRate has been overridden, reflect that in the
		// _eh._offRate field
		if(offRatePlus > 0 && this->_overrideOffRate == -1) {
			this->_overrideOffRate = this->_gh._offRate + offRatePlus;
		}
		if(this->_overrideOffRate > this->_gh._offRate) {
			this->_gh.setOffRate(this->_overrideOffRate);
			assert_eq(this->_overrideOffRate, this->_gh._offRate);
		}
		assert(this->repOk());
	}


	/// Construct an Ebwt from the given header parameters and string
	/// vector, optionally using a blockwise suffix sorter with the
	/// given 'bmax' and 'dcv' parameters.  The string vector is
	/// ultimately joined and the joined string is passed to buildToDisk().
	template<typename TStr>
	LocalGFM(
             TStr& s,
             const EList<full_index_t>& sa,
             PathGraph<full_index_t>* pg,
             full_index_t tidx,
             full_index_t localOffset,
             full_index_t joinedOffset,
             EList<ALT<full_index_t> >& alts,
             index_t local_size,
             bool packed,
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
	GFM<index_t>(packed,
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
		const GFMParams<index_t>& gh = this->_gh;
		assert(gh.repOk());
		uint32_t be = this->toBe();
		assert(out5.good());
		assert(out6.good());
        _tidx = tidx;
        _localOffset = localOffset;
        _joinedOffset = joinedOffset;
		writeIndex<full_index_t>(out5, tidx, be);
		writeIndex<full_index_t>(out5, localOffset, be);
        writeIndex<full_index_t>(out5, joinedOffset, be);
        writeIndex<index_t>(out5, gh._len, be); // length of string (and bwt and suffix array)
        streampos headerPos = out5.tellp();
        writeIndex<index_t>(out5, 0, be); // gbwtLen
        writeIndex<index_t>(out5, 0, be); // num of nodes
        writeIndex<index_t>(out5, 0, be); // eftabLen        
		if(gh._len > 0) {
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
            
            if(alts.empty()) {
                assert(pg == NULL);
                buildToDisk(sa, s, out5, out6, headerPos);
            } else {
                assert(pg != NULL);
                // Re-initialize GFM parameters to reflect real number of edges (gbwt string)
                this->_gh.init(
                               this->_gh.len(),
                               pg->getNumEdges(),
                               pg->getNumNodes(),
                               this->_gh.lineRate(),
                               this->_gh.offRate(),
                               this->_gh.ftabChars(),
                               0,
                               this->_gh.entireReverse());
                buildToDisk(*pg, s, out5, out6, headerPos);
            }
		}
		
		out5.flush(); out6.flush();
		if(out5.fail() || out6.fail()) {
			cerr << "An error occurred writing the index to disk.  Please check if the disk is full." << endl;
			throw 1;
		}
	}
	
	template <typename TStr> void buildToDisk(
											  PathGraph<full_index_t>& gbwt,
											  const TStr& s,
											  ostream& out1, 
											  ostream& out2,
                                              streampos& headerPos);
    
    template <typename TStr> void buildToDisk(
                                              const EList<full_index_t>& sa,
                                              const TStr& s,
                                              ostream& out1,
                                              ostream& out2,
                                              streampos& headerPos);
	
	// I/O
	void readIntoMemory(
						FILE *in5,
						FILE *in6,
						char *mmFile5,
						char *mmFile6,
						full_index_t& tidx,
						full_index_t& localOffset,
                        full_index_t& joinedOffset,
						bool switchEndian,
						size_t& bytesRead,
                        size_t& bytesRead2,
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
		if(this->_gh._len > 0) {
			PARENT_CLASS::sanityCheckAll(reverse);
		}
	}
    
    bool empty() const { return this->_gh._len == 0; }
	
public:
	full_index_t _tidx;
	full_index_t _localOffset;
    full_index_t _joinedOffset;
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
void LocalGFM<index_t, full_index_t>::buildToDisk(
                                                  PathGraph<full_index_t>& gbwt,
                                                  const TStr& s,
                                                  ostream& out5,
                                                  ostream& out6,
                                                  streampos& headerPos)
{
	assert_leq(s.length(), std::numeric_limits<index_t>::max());
	const GFMParams<index_t>& gh = this->_gh;
	
	assert(gh.repOk());
    assert_lt(s.length(), gh.gbwtLen());
    assert_eq(s.length(), gh._len);
	assert_gt(gh._lineRate, 3);
    
    index_t  gbwtLen = gh._gbwtLen;
    streampos out5pos = out5.tellp();
    out5.seekp(headerPos);
    writeIndex<index_t>(out5, gbwtLen, this->toBe());
    writeIndex<index_t>(out5, gh._numNodes, this->toBe());
    headerPos = out5.tellp();
    out5.seekp(out5pos);
	index_t ftabLen = gh._ftabLen;
	index_t sideSz = gh._sideSz;
	index_t gbwtTotSz = gh._gbwtTotSz;
	index_t fchr[] = {0, 0, 0, 0, 0};
	EList<index_t> ftab(EBWT_CAT);
	EList<index_t> zOffs;
	
	// Save # of occurrences of each character as we walk along the bwt
	index_t occ[4] = {0, 0, 0, 0};
	index_t occSave[4] = {0, 0, 0, 0};
    // # of occurrences of 1 in M arrays
    index_t M_occ = 0, M_occSave = 0;
    // Location in F that corresponds to 1 in M
    index_t F_loc = 0, F_locSave = 0;
	
	try {
		VMSG_NL("Allocating ftab, absorbFtab");
		ftab.resize(ftabLen);
		ftab.fillZero();
    } catch(bad_alloc &e) {
		cerr << "Out of memory allocating ftab[] or absorbFtab[] "
		<< "in LocalGFM::buildToDisk() at " << __FILE__ << ":"
		<< __LINE__ << endl;
		throw e;
	}
	
	// Allocate the side buffer; holds a single side as its being
	// constructed and then written to disk.  Reused across all sides.
#ifdef SIXTY4_FORMAT
	EList<uint64_t> gfmSide(EBWT_CAT);
#else
	EList<uint8_t> gfmSide(EBWT_CAT);
#endif
	try {
        // Used to calculate ftab and eftab, but having gfm costs a lot of memory
        this->_gfm.init(new uint8_t[gh._gbwtTotLen], gh._gbwtTotLen, true);
#ifdef SIXTY4_FORMAT
		gfmSide.resize(sideSz >> 3);
#else
		gfmSide.resize(sideSz);
#endif
	} catch(bad_alloc &e) {
		cerr << "Out of memory allocating ebwtSide[] in "
		<< "LocalGFM::buildToDisk() at " << __FILE__ << ":"
		<< __LINE__ << endl;
		throw e;
	}
	
	// Points to the base offset within ebwt for the side currently
	// being written
	index_t side = 0;
	
	// Whether we're assembling a forward or a reverse bucket
    bool fw = true;
	int sideCur = 0;
	
	index_t si = 0;   // string offset (chars)
	ASSERT_ONLY(bool inSA = true); // true iff saI still points inside suffix
	// array (as opposed to the padding at the
	// end)
	// Iterate over packed bwt bytes
	VMSG_NL("Entering LocalGFM loop");
	ASSERT_ONLY(uint32_t beforeGbwtOff = (uint32_t)out5.tellp());
	while(side < gbwtTotSz) {
        // Sanity-check our cursor into the side buffer
		assert_geq(sideCur, 0);
		assert_lt(sideCur, (int)gh._sideGbwtSz);
		assert_eq(0, side % sideSz); // 'side' must be on side boundary
		gfmSide[sideCur] = 0; // clear
        if(sideCur == 0) {
            memset(gfmSide.ptr(), 0, gh._sideGbwtSz);
            gfmSide[sideCur] = 0; // clear
        }
        assert_lt(side + sideCur, gbwtTotSz);
		// Iterate over bit-pairs in the si'th character of the BWT
#ifdef SIXTY4_FORMAT
		for(int bpi = 0; bpi < 32; bpi++, si++) {
#else
		for(int bpi = 0; bpi < 4; bpi++, si++) {
#endif
			int gbwtChar = 0;
            int F = 0, M = 0;
            full_index_t pos = 0;
            bool count = true;
			if(si < gbwtLen) {
                gbwt.nextRow(gbwtChar, F, M, pos);

				// (that might have triggered sa to calc next suf block)
				if(gbwtChar == 'Z') {
					// Don't add the '$' in the last column to the BWT
					// transform; we can't encode a $ (only A C T or G)
					// and counting it as, say, an A, will mess up the
					// LR mapping
					gbwtChar = 0; count = false;
#ifndef NDEBUG
                    if(zOffs.size() > 0) {
                        assert_gt(si, zOffs.back());
                    }
#endif
                    zOffs.push_back(si); // remember GBWT row that corresponds to the 0th suffix
				} else {
                    gbwtChar = asc2dna[gbwtChar];
                    assert_lt(gbwtChar, 4);
                    // Update the fchr
                    fchr[gbwtChar]++;
				}
                assert_lt(F, 2);
                assert_lt(M, 2);
                if(M == 1) {
                    assert_neq(F_loc, numeric_limits<index_t>::max());
                    F_loc = gbwt.nextFLocation();
#ifndef NDEBUG
                    if(F_loc > 0) {
                        assert_gt(F_loc, F_locSave);
                    }
#endif
                }
                // Suffix array offset boundary? - update offset array
                if(M == 1 && (M_occ & gh._offMask) == M_occ) {
                    assert_lt((M_occ >> gh._offRate), gh._offsLen);
                    // Write offsets directly to the secondary output
                    // stream, thereby avoiding keeping them in memory
                    writeIndex<index_t>(out6, pos, this->toBe());
                }
            } else {
				// Strayed off the end of the SA, now we're just
				// padding out a bucket
#ifndef NDEBUG
				if(inSA) {
					// Assert that we wrote all the characters in the
					// string before now
					assert_eq(si, gbwtLen);
					inSA = false;
				}
#endif
				// 'A' used for padding; important that padding be
				// counted in the occ[] array
				gbwtChar = 0;
                F = M = 0;
			}
            if(count) occ[gbwtChar]++;
            if(M) M_occ++;
            // Append BWT char to bwt section of current side
            if(fw) {
                // Forward bucket: fill from least to most
#ifdef SIXTY4_FORMAT
                gfmSide[sideCur] |= ((uint64_t)gbwtChar << (bpi << 1));
                if(gbwtChar > 0) assert_gt(gfmSide[sideCur], 0);
                assert(false);
                cerr << "Not implemented" << endl;
                exit(1);
#else
                pack_2b_in_8b(gbwtChar, gfmSide[sideCur], bpi);
                assert_eq((gfmSide[sideCur] >> (bpi*2)) & 3, gbwtChar);
                
                int F_sideCur = (gh._sideGbwtSz + sideCur) >> 1;
                int F_bpi = bpi + ((sideCur & 0x1) << 2); // Can be used as M_bpi as well
                pack_1b_in_8b(F, gfmSide[F_sideCur], F_bpi);
                assert_eq((gfmSide[F_sideCur] >> F_bpi) & 1, F);
                
                int M_sideCur = F_sideCur + (gh._sideGbwtSz >> 2);
                pack_1b_in_8b(M, gfmSide[M_sideCur], F_bpi);
                assert_eq((gfmSide[M_sideCur] >> F_bpi) & 1, M);
#endif
            } else {
				// Backward bucket: fill from most to least
#ifdef SIXTY4_FORMAT
                gfmSide[sideCur] |= ((uint64_t)gbwtChar << ((31 - bpi) << 1));
                if(gbwtChar > 0) assert_gt(gfmSide[sideCur], 0);
                // To be implemented ...
                assert(false);
                cerr << "Not implemented" << endl;
                exit(1);
#else
                pack_2b_in_8b(gbwtChar, gfmSide[sideCur], 3-bpi);
				assert_eq((gfmSide[sideCur] >> ((3-bpi)*2)) & 3, gbwtChar);
                // To be implemented ...
                assert(false);
                cerr << "Not implemented" << endl;
                exit(1);
#endif
            }
        } // end loop over bit-pairs
        assert_eq(0, (occ[0] + occ[1] + occ[2] + occ[3] + zOffs.size()) & 3);
#ifdef SIXTY4_FORMAT
		assert_eq(0, si & 31);
#else
		assert_eq(0, si & 3);
#endif
		
		sideCur++;
		if((sideCur << 1) == (int)gh._sideGbwtSz) {
			sideCur = 0;
			index_t *uside = reinterpret_cast<index_t*>(gfmSide.ptr());
			// Write 'A', 'C', 'G' and 'T' tallies
			side += sideSz;
			assert_leq(side, gh._gbwtTotSz);
            uside[(sideSz / sizeof(index_t))-6] = endianizeIndex(F_locSave, this->toBe());
            uside[(sideSz / sizeof(index_t))-5] = endianizeIndex(M_occSave, this->toBe());
			uside[(sideSz / sizeof(index_t))-4] = endianizeIndex(occSave[0], this->toBe());
			uside[(sideSz / sizeof(index_t))-3] = endianizeIndex(occSave[1], this->toBe());
			uside[(sideSz / sizeof(index_t))-2] = endianizeIndex(occSave[2], this->toBe());
			uside[(sideSz / sizeof(index_t))-1] = endianizeIndex(occSave[3], this->toBe());
            F_locSave = F_loc;
            M_occSave = M_occ;
			occSave[0] = occ[0];
			occSave[1] = occ[1];
			occSave[2] = occ[2];
			occSave[3] = occ[3];
			// Write backward side to primary file
			out5.write((const char *)gfmSide.ptr(), sideSz);
            
            //
            memcpy(((char*)this->_gfm.get()) + side - sideSz, (const char *)gfmSide.ptr(), sideSz);
		}
	}
	VMSG_NL("Exited LocalGFM loop");
	// Assert that our loop counter got incremented right to the end
	assert_eq(side, gh._gbwtTotSz);
	// Assert that we wrote the expected amount to out5
	assert_eq(((uint32_t)out5.tellp() - beforeGbwtOff), gh._gbwtTotSz);
	// assert that the last thing we did was write a forward bucket
	
    //
    // Write zOffs to primary stream
    //
    assert_gt(zOffs.size(), 0);
    writeIndex<index_t>(out5, zOffs.size(), this->toBe());
    for(size_t i = 0; i < zOffs.size(); i++) {
        writeIndex<index_t>(out5, zOffs[i], this->toBe());
    }
        
    //
    // Finish building fchr
    //
    // Exclusive prefix sum on fchr
    for(int i = 1; i < 4; i++) {
        fchr[i] += fchr[i-1];
    }
    assert_lt(fchr[3], gbwtLen);
    // Shift everybody up by one
    for(int i = 4; i >= 1; i--) {
        fchr[i] = fchr[i-1];
    }
    fchr[0] = 0;
    // Write fchr to primary file
    for(int i = 0; i < 5; i++) {
        writeIndex<index_t>(out5, fchr[i], this->toBe());
    }
    this->_fchr.init(new index_t[5], 5, true);
    memcpy(this->_fchr.get(), fchr, sizeof(index_t) * 5);
        
    // Initialize _zGbwtByteOffs and _zGbwtBpOffs
    this->_zOffs = zOffs;
    this->postReadInit(gh);
	
    // Build ftab and eftab
    EList<pair<index_t, index_t> > tFtab;
    tFtab.resizeExact(ftabLen - 1);
    for(index_t i = 0; i + 1 < ftabLen; i++) {
        index_t q = i;
        pair<index_t, index_t> range(0, gh._gbwtLen);
        SideLocus<index_t> tloc, bloc;
        SideLocus<index_t>::initFromTopBot(range.first, range.second, gh, this->gfm(), tloc, bloc);
        index_t j = 0;
        for(; j < gh._ftabChars; j++) {
            int nt = q & 0x3; q >>= 2;
            if(bloc.valid()) {
                range = this->mapGLF(tloc, bloc, nt);
            } else {
                range = this->mapGLF1(range.first, tloc, nt);
            }
            if(range.first == (index_t)INDEX_MAX || range.first >= range.second) {
                break;
            }
            if(range.first + 1 == range.second) {
                tloc.initFromRow(range.first, gh, this->gfm());
                bloc.invalidate();
            } else {
                SideLocus<index_t>::initFromTopBot(range.first, range.second, gh, this->gfm(), tloc, bloc);
            }
        }
            
        if(range.first >= range.second || j < gh._ftabChars) {
            if(i == 0) {
                tFtab[i].first = tFtab[i].second = 0;
            } else {
                tFtab[i].first = tFtab[i].second = tFtab[i-1].second;
            }
        } else {
            tFtab[i].first = range.first;
            tFtab[i].second = range.second;
        }
            
#ifndef NDEBUG
        if(gbwt.ftab.size() > i) {
            assert_eq(tFtab[i].first, gbwt.ftab[i].first);
            assert_eq(tFtab[i].second, gbwt.ftab[i].second);
        }
#endif
    }
        
    // Clear memory
    this->_gfm.reset();
    this->_fchr.reset();
    this->_zOffs.clear();
    this->_zGbwtByteOffs.clear();
    this->_zGbwtBpOffs.clear();
    
    //
    // Finish building ftab and build eftab
    //
    // Prefix sum on ftable
    index_t eftabLen = 0;
    for(index_t i = 1; i + 1 < ftabLen; i++) {
        if(tFtab[i-1].second != tFtab[i].first) {
            eftabLen += 2;
        }
    }
    if(gh._gbwtLen + (eftabLen >> 1) < gh._gbwtLen) {
        cerr << "Too many eftab entries: "
             << gh._gbwtLen << " + " << (eftabLen >> 1)
             << " > " << (index_t)INDEX_MAX << endl;
        throw 1;
    }
    EList<index_t> eftab(EBWT_CAT);
    try {
        eftab.resize(eftabLen);
        eftab.fillZero();
    } catch(bad_alloc &e) {
        cerr << "Out of memory allocating eftab[] "
        << "in LocalGFM::buildToDisk() at " << __FILE__ << ":"
        << __LINE__ << endl;
        throw e;
    }
    index_t eftabCur = 0;
    ftab[0] = tFtab[0].first;
    ftab[1] = tFtab[0].second;
    for(index_t i = 1; i + 1 < ftabLen; i++) {
        if(ftab[i] != tFtab[i].first) {
            index_t lo = ftab[i];
            index_t hi = tFtab[i].first;
            assert_lt(eftabCur*2+1, eftabLen);
            eftab[eftabCur*2] = lo;
            eftab[eftabCur*2+1] = hi;
            assert_leq(lo, hi + 4);
            ftab[i] = (eftabCur++) ^ (index_t)INDEX_MAX; // insert pointer into eftab
            assert_eq(lo, GFM<index_t>::ftabLo(ftab.ptr(), eftab.ptr(), gbwtLen, ftabLen, eftabLen, i));
            assert_eq(hi, GFM<index_t>::ftabHi(ftab.ptr(), eftab.ptr(), gbwtLen, ftabLen, eftabLen, i));
        }
        ftab[i+1] = tFtab[i].second;
    }
#ifndef NDEBUG
    for(index_t i = 0; i + 1 < ftabLen; i++ ){
        assert_eq(tFtab[i].first, GFM<index_t>::ftabHi(ftab.ptr(), eftab.ptr(), gbwtLen, ftabLen, eftabLen, i));
        assert_eq(tFtab[i].second, GFM<index_t>::ftabLo(ftab.ptr(), eftab.ptr(), gbwtLen, ftabLen, eftabLen, i+1));
    }
#endif
    // Write ftab to primary file
    for(index_t i = 0; i < ftabLen; i++) {
        writeIndex<index_t>(out5, ftab[i], this->toBe());
    }
    // Write eftab to primary file
    out5pos = out5.tellp();
    out5.seekp(headerPos);
    writeIndex<index_t>(out5, eftabLen, this->toBe());
    out5.seekp(out5pos);
    for(index_t i = 0; i < eftabLen; i++) {
        writeIndex<index_t>(out5, eftab[i], this->toBe());
    }
    // Note: if you'd like to sanity-check the Ebwt, you'll have to
    // read it back into memory first!
    assert(!this->isInMemory());
    VMSG_NL("Exiting LocalGFM::buildToDisk()");
}
    
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
void LocalGFM<index_t, full_index_t>::buildToDisk(
                                                  const EList<full_index_t>& sa,
                                                  const TStr& s,
                                                  ostream& out5,
                                                  ostream& out6,
                                                  streampos& headerPos)
{
    assert_leq(s.length(), std::numeric_limits<index_t>::max());
    const GFMParams<index_t>& gh = this->_gh;
    assert(gh.repOk());
    assert(gh.linearFM());
    assert_lt(s.length(), gh.gbwtLen());
    assert_eq(s.length(), gh._len);
    assert_gt(gh._lineRate, 3);
    
    index_t  len = gh._len;
    index_t  gbwtLen = gh._gbwtLen;
    assert_eq(len + 1, gbwtLen);
    streampos out5pos = out5.tellp();
    out5.seekp(headerPos);
    writeIndex<index_t>(out5, gbwtLen, this->toBe());
    writeIndex<index_t>(out5, gh._numNodes, this->toBe());
    headerPos = out5.tellp();
    out5.seekp(out5pos);
    
    index_t ftabLen = gh._ftabLen;
    index_t sideSz = gh._sideSz;
    index_t gbwtTotSz = gh._gbwtTotSz;
    index_t fchr[] = {0, 0, 0, 0, 0};
    EList<index_t> ftab(EBWT_CAT);
    EList<index_t> zOffs;
    
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
        << "in LocalGFM::buildToDisk() at " << __FILE__ << ":"
        << __LINE__ << endl;
        throw e;
    }
    
    // Allocate the side buffer; holds a single side as its being
    // constructed and then written to disk.  Reused across all sides.
#ifdef SIXTY4_FORMAT
    EList<uint64_t> gfmSide(EBWT_CAT);
#else
    EList<uint8_t> gfmSide(EBWT_CAT);
#endif
    try {
#ifdef SIXTY4_FORMAT
        gfmSide.resize(sideSz >> 3);
#else
        gfmSide.resize(sideSz);
#endif
    } catch(bad_alloc &e) {
        cerr << "Out of memory allocating gfmSide[] in "
        << "LocalGFM::buildToDisk() at " << __FILE__ << ":"
        << __LINE__ << endl;
        throw e;
    }
    
    // Points to the base offset within ebwt for the side currently
    // being written
    index_t side = 0;
    
    // Whether we're assembling a forward or a reverse bucket
    bool fw = true;
    int sideCur = 0;
    
    // Have we skipped the '$' in the last column yet?
    ASSERT_ONLY(bool dollarSkipped = false);
    
    index_t si = 0;   // string offset (chars)
    ASSERT_ONLY(uint32_t lastSufInt = 0);
    ASSERT_ONLY(bool inSA = true); // true iff saI still points inside suffix
    // array (as opposed to the padding at the
    // end)
    // Iterate over packed bwt bytes
    VMSG_NL("Entering LocalGFM loop");
    ASSERT_ONLY(uint32_t beforeGbwtOff = (uint32_t)out5.tellp());
    while(side < gbwtTotSz) {
        // Sanity-check our cursor into the side buffer
        assert_geq(sideCur, 0);
        assert_lt(sideCur, (int)gh._sideGbwtSz);
        assert_eq(0, side % sideSz); // 'side' must be on side boundary
        gfmSide[sideCur] = 0; // clear
        assert_lt(side + sideCur, gbwtTotSz);
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
                index_t saElt = (index_t)sa[si];
                // (that might have triggered sa to calc next suf block)
                if(saElt == 0) {
                    // Don't add the '$' in the last column to the BWT
                    // transform; we can't encode a $ (only A C T or G)
                    // and counting it as, say, an A, will mess up the
                    // LR mapping
                    bwtChar = 0; count = false;
                    ASSERT_ONLY(dollarSkipped = true);
                    zOffs.push_back(si); // remember the SA row that
                    // corresponds to the 0th suffix
                } else {
                    bwtChar = (int)(s[saElt-1]);
                    assert_lt(bwtChar, 4);
                    // Update the fchr
                    fchr[bwtChar]++;
                }
                // Update ftab
                if((len-saElt) >= (index_t)gh._ftabChars) {
                    // Turn the first ftabChars characters of the
                    // suffix into an integer index into ftab.  The
                    // leftmost (lowest index) character of the suffix
                    // goes in the most significant bit pair if the
                    // integer.
                    uint32_t sufInt = 0;
                    for(int i = 0; i < gh._ftabChars; i++) {
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
                if((si & gh._offMask) == si) {
                    assert_lt((si >> gh._offRate), gh._offsLen);
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
                gfmSide[sideCur] |= ((uint64_t)bwtChar << (bpi << 1));
                if(bwtChar > 0) assert_gt(gfmSide[sideCur], 0);
#else
                pack_2b_in_8b(bwtChar, gfmSide[sideCur], bpi);
                assert_eq((gfmSide[sideCur] >> (bpi*2)) & 3, bwtChar);
#endif
            } else {
                // Backward bucket: fill from most to least
#ifdef SIXTY4_FORMAT
                gfmSide[sideCur] |= ((uint64_t)bwtChar << ((31 - bpi) << 1));
                if(bwtChar > 0) assert_gt(gfmSide[sideCur], 0);
#else
                pack_2b_in_8b(bwtChar, gfmSide[sideCur], 3-bpi);
                assert_eq((gfmSide[sideCur] >> ((3-bpi)*2)) & 3, bwtChar);
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
        if(sideCur == (int)gh._sideGbwtSz) {
            sideCur = 0;
            index_t *uside = reinterpret_cast<index_t*>(gfmSide.ptr());
            // Write 'A', 'C', 'G' and 'T' tallies
            side += sideSz;
            assert_leq(side, gh._gbwtTotSz);
            uside[(sideSz / sizeof(index_t))-4] = endianizeIndex(occSave[0], this->toBe());
            uside[(sideSz / sizeof(index_t))-3] = endianizeIndex(occSave[1], this->toBe());
            uside[(sideSz / sizeof(index_t))-2] = endianizeIndex(occSave[2], this->toBe());
            uside[(sideSz / sizeof(index_t))-1] = endianizeIndex(occSave[3], this->toBe());
            occSave[0] = occ[0];
            occSave[1] = occ[1];
            occSave[2] = occ[2];
            occSave[3] = occ[3];
            // Write backward side to primary file
            out5.write((const char *)gfmSide.ptr(), sideSz);
        }
    }
    VMSG_NL("Exited LocalGFM loop");
    if(absorbCnt > 0) {
        // Absorb any trailing, as-yet-unabsorbed short suffixes into
        // the last element of ftab
        absorbFtab[ftabLen-1] = absorbCnt;
    }
    // Assert that our loop counter got incremented right to the end
    assert_eq(side, gh._gbwtTotSz);
    // Assert that we wrote the expected amount to out5
    assert_eq(((uint32_t)out5.tellp() - beforeGbwtOff), gh._gbwtTotSz);
    // assert that the last thing we did was write a forward bucket
        
    //
    // Write zOffs to primary stream
    //
    assert_eq(zOffs.size(), 1);
    writeIndex<index_t>(out5, zOffs.size(), this->toBe());
    for(size_t i = 0; i < zOffs.size(); i++) {
        assert_neq(zOffs[i], (index_t)OFF_MASK);
        writeIndex<index_t>(out5, zOffs[i], this->toBe());
    }
        
    //
    // Finish building fchr
    //
    // Exclusive prefix sum on fchr
    for(int i = 1; i < 4; i++) {
        fchr[i] += fchr[i-1];
    }
    assert_lt(fchr[3], gbwtLen);
    // Shift everybody up by one
    for(int i = 4; i >= 1; i--) {
        fchr[i] = fchr[i-1];
    }
    fchr[0] = 0;
    // Write fchr to primary file
    for(int i = 0; i < 5; i++) {
        writeIndex<index_t>(out5, fchr[i], this->toBe());
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
    assert_leq(eftabLen, (index_t)gh._ftabChars*2);
    eftabLen = gh._ftabChars*2;
    EList<index_t> eftab(EBWT_CAT);
    try {
        eftab.resize(eftabLen);
        eftab.fillZero();
    } catch(bad_alloc &e) {
        cerr << "Out of memory allocating eftab[] "
        << "in LocalGFM::buildToDisk() at " << __FILE__ << ":"
        << __LINE__ << endl;
        throw e;
    }
    index_t eftabCur = 0;
    for(index_t i = 1; i < ftabLen; i++) {
        index_t lo = ftab[i] + GFM<index_t>::ftabHi(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, i-1);
        if(absorbFtab[i] > 0) {
            // Skip a number of short pattern indicated by absorbFtab[i]
            index_t hi = lo + absorbFtab[i];
            assert_lt(eftabCur*2+1, eftabLen);
            eftab[eftabCur*2] = lo;
            eftab[eftabCur*2+1] = hi;
            ftab[i] = (eftabCur++) ^ (index_t)OFF_MASK; // insert pointer into eftab
            assert_eq(lo, GFM<index_t>::ftabLo(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, i));
            assert_eq(hi, GFM<index_t>::ftabHi(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, i));
        } else {
            ftab[i] = lo;
        }
    }
    assert_eq(GFM<index_t>::ftabHi(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, ftabLen-1), len+1);
    // Write ftab to primary file
    for(index_t i = 0; i < ftabLen; i++) {
        writeIndex(out5, ftab[i], this->toBe());
    }
    // Write eftab to primary file
    out5pos = out5.tellp();
    out5.seekp(headerPos);
    writeIndex<index_t>(out5, eftabLen, this->toBe());
    out5.seekp(out5pos);
    for(index_t i = 0; i < eftabLen; i++) {
        writeIndex(out5, eftab[i], this->toBe());
    }
        
    // Note: if you'd like to sanity-check the Ebwt, you'll have to
    // read it back into memory first!
    assert(!this->isInMemory());
    VMSG_NL("Exiting LocalGFM::buildToDisk()");
}
    
/**
 * Read an Ebwt from file with given filename.
 */
template <typename index_t, typename full_index_t>
void LocalGFM<index_t, full_index_t>::readIntoMemory(
                                                     FILE *in5,
                                                     FILE *in6,
                                                     char *mmFile5,
                                                     char *mmFile6,
                                                     full_index_t& tidx,
                                                     full_index_t& localOffset,
                                                     full_index_t& joinedOffset,
                                                     bool switchEndian,
                                                     size_t& bytesRead,
                                                     size_t& bytesRead2,
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
	tidx             = readIndex<full_index_t>(in5, switchEndian); bytesRead += sizeof(full_index_t);
	localOffset      = readIndex<full_index_t>(in5, switchEndian); bytesRead += sizeof(full_index_t);
    joinedOffset     = readIndex<full_index_t>(in5, switchEndian); bytesRead += sizeof(full_index_t);
    index_t len      = readIndex<index_t>(in5, switchEndian); bytesRead += sizeof(index_t);
    index_t gbwtLen  = readIndex<index_t>(in5, switchEndian); bytesRead += sizeof(index_t);
    index_t numNodes = readIndex<index_t>(in5, switchEndian); bytesRead += sizeof(index_t);
    index_t eftabLen = readIndex<index_t>(in5, switchEndian); bytesRead += sizeof(index_t);
	
	// Create a new EbwtParams from the entries read from primary stream
	this->_gh.init(len, gbwtLen, numNodes, lineRate, offRate, ftabChars, eftabLen, entireRev);
	
	if(len <= 0) {
		return;
	}
	
	// Set up overridden suffix-array-sample parameters
	uint32_t offsLen = this->_gh._offsLen;
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
	
	this->_gfm.reset();
	if(this->_useMm) {
#ifdef BOWTIE_MM
		this->_gfm.init((uint8_t*)(mmFile[0] + bytesRead), this->_gh._gbwtTotLen, false);
		bytesRead += this->_gh._gbwtTotLen;
		fseek(in5, this->_gh._gbwtTotLen, SEEK_CUR);
#endif
	} else {
		// Allocate ebwt (big allocation)
		if(this->_verbose || startVerbose) {
			cerr << "Reading ebwt (" << this->_gh._gbwtTotLen << "): ";
			logTime(cerr);
		}
		bool shmemLeader = true;
		if(this->useShmem_) {
			uint8_t *tmp = NULL;
			shmemLeader = ALLOC_SHARED_U8(
										  (this->_in1Str + "[gfm]"), this->_gh._gbwtTotLen, &tmp,
										  "gfm[]", (this->_verbose || startVerbose));
			assert(tmp != NULL);
			this->_gfm.init(tmp, this->_gh._gbwtTotLen, false);
			if(this->_verbose || startVerbose) {
				cerr << "  shared-mem " << (shmemLeader ? "leader" : "follower") << endl;
			}
		} else {
			try {
				this->_gfm.init(new uint8_t[this->_gh._gbwtTotLen], this->_gh._gbwtTotLen, true);
			} catch(bad_alloc& e) {
				cerr << "Out of memory allocating the ebwt[] array for the Bowtie index.  Please try" << endl
				<< "again on a computer with more memory." << endl;
				throw 1;
			}
		}
		if(shmemLeader) {
			// Read ebwt from primary stream
			uint64_t bytesLeft = this->_gh._gbwtTotLen;
			char *pgbwt = (char*)this->gfm();
            
			while (bytesLeft>0){
				size_t r = MM_READ(in5, (void *)pgbwt, bytesLeft);
				if(MM_IS_IO_ERR(in5, r, bytesLeft)) {
					cerr << "Error reading _ebwt[] array: " << r << ", "
                    << bytesLeft << endl;
					throw 1;
				}
				pgbwt += r;
				bytesLeft -= r;
			}
			if(switchEndian) {
				uint8_t *side = this->gfm();
				for(size_t i = 0; i < this->_gh._numSides; i++) {
					index_t *cums = reinterpret_cast<index_t*>(side + this->_gh._sideSz - sizeof(index_t)*2);
					cums[0] = endianSwapIndex(cums[0]);
					cums[1] = endianSwapIndex(cums[1]);
					side += this->_gh._sideSz;
				}
			}
#ifdef BOWTIE_SHARED_MEM
			if(useShmem_) NOTIFY_SHARED(this->gfm(), this->_gh._gbwtTotLen);
#endif
		} else {
			// Seek past the data and wait until master is finished
			fseek(in5, this->_gh._gbwtTotLen, SEEK_CUR);
#ifdef BOWTIE_SHARED_MEM
			if(useShmem_) WAIT_SHARED(this->gfm(), this->_gh._gbwtTotLen);
#endif
		}
	}
	
	// Read zOff from primary stream
    this->_zOffs.clear();
    index_t num_zOffs = readIndex<index_t>(in5, switchEndian);
    bytesRead += sizeof(index_t);
    for(index_t i = 0; i < num_zOffs; i++) {
        index_t zOff = readIndex<index_t>(in5, switchEndian);
        bytesRead += sizeof(index_t);
        assert_lt(zOff, gbwtLen);
        this->_zOffs.push_back(zOff);
    }	
	
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
				assert_leq(this->fchr()[i], gbwtLen);
				assert(i <= 0 || this->fchr()[i] >= this->fchr()[i-1]);
			}
		}
		assert_gt(this->fchr()[4], this->fchr()[0]);
		// Read ftab from primary stream
		if(this->_verbose || startVerbose) {
			if(loadFtab) {
				cerr << "Reading ftab (" << this->_gh._ftabLen << "): ";
				logTime(cerr);
			} else {
				cerr << "Skipping ftab (" << this->_gh._ftabLen << "): ";
			}
		}
		this->_ftab.reset();
		if(loadFtab) {
			if(this->_useMm) {
#ifdef BOWTIE_MM
				this->_ftab.init((index_t*)(mmFile[0] + bytesRead), this->_gh._ftabLen, false);
				bytesRead += this->_gh._ftabLen*sizeof(index_t);
				fseek(in5, this->_gh._ftabLen*sizeof(index_t), SEEK_CUR);
#endif
			} else {
				this->_ftab.init(new index_t[this->_gh._ftabLen], this->_gh._ftabLen, true);
				if(switchEndian) {
					for(uint32_t i = 0; i < this->_gh._ftabLen; i++)
						this->ftab()[i] = readIndex<index_t>(in5, switchEndian);
				} else {
					size_t r = MM_READ(in5, (void *)this->ftab(), this->_gh._ftabLen*sizeof(index_t));
					if(r != (size_t)(this->_gh._ftabLen*sizeof(index_t))) {
						cerr << "Error reading _ftab[] array: " << r << ", " << (this->_gh._ftabLen*sizeof(index_t)) << endl;
						throw 1;
					}
				}
			}
			// Read etab from primary stream
			if(this->_verbose || startVerbose) {
				if(loadFtab) {
					cerr << "Reading eftab (" << this->_gh._eftabLen << "): ";
					logTime(cerr);
				} else {
					cerr << "Skipping eftab (" << this->_gh._eftabLen << "): ";
				}
				
			}
			this->_eftab.reset();
			if(this->_useMm) {
#ifdef BOWTIE_MM
				this->_eftab.init((index_t*)(mmFile[0] + bytesRead), this->_gh._eftabLen, false);
				bytesRead += this->_gh._eftabLen*sizeof(index_t);
				fseek(in5, this->_gh._eftabLen*sizeof(index_t), SEEK_CUR);
#endif
			} else {
				this->_eftab.init(new index_t[this->_gh._eftabLen], this->_gh._eftabLen, true);
				if(switchEndian) {
					for(uint32_t i = 0; i < this->_gh._eftabLen; i++)
						this->eftab()[i] = readIndex<index_t>(in5, switchEndian);
				} else {
					size_t r = MM_READ(in5, (void *)this->eftab(), this->_gh._eftabLen*sizeof(index_t));
					if(r != (size_t)(this->_gh._eftabLen*sizeof(index_t))) {
						cerr << "Error reading _eftab[] array: " << r << ", " << (this->_gh._eftabLen*sizeof(index_t)) << endl;
						throw 1;
					}
				}
			}
			for(uint32_t i = 0; i < this->_gh._eftabLen; i++) {
				if(i > 0 && this->eftab()[i] > 0) {
					assert_geq(this->eftab()[i] + 4, this->eftab()[i-1]);
				} else if(i > 0 && this->eftab()[i-1] == 0) {
					assert_eq(0, this->eftab()[i]);
				}
			}
		} else {
			assert(this->ftab() == NULL);
			assert(this->eftab() == NULL);
			// Skip ftab
			bytesRead += this->_gh._ftabLen*sizeof(index_t);
			fseek(in5, this->_gh._ftabLen*sizeof(index_t), SEEK_CUR);
			// Skip eftab
			bytesRead += this->_gh._eftabLen*sizeof(index_t);
			fseek(in5, this->_gh._eftabLen*sizeof(index_t), SEEK_CUR);
		}
	} catch(bad_alloc& e) {
		cerr << "Out of memory allocating fchr[], ftab[] or eftab[] arrays for the Bowtie index." << endl
		<< "Please try again on a computer with more memory." << endl;
		throw 1;
	}
	
	this->_offs.reset();
	if(loadSASamp) {
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
						this->_offs.init((index_t*)(mmFile[1] + bytesRead2), offsLen, false);
						bytesRead2 += (offsLen * sizeof(index_t));
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
	
	this->postReadInit(this->_gh); // Initialize fields of Ebwt not read from file
	if(this->_verbose || startVerbose) this->print(cerr, this->_gh);
}

/**
 * Extended Burrows-Wheeler transform data.
 * HierEbwt is a specialized Ebwt index that represents one global index and a large set of local indexes.
 *
 */
template <typename index_t = uint32_t, typename local_index_t = uint16_t>
class HGFM : public GFM<index_t> {
	typedef GFM<index_t> PARENT_CLASS;
public:
	/// Construct an Ebwt from the given input file
	HGFM(const string& in,
         ALTDB<index_t>* altdb,
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
         bool loadSpliceSites, // = true,
         bool verbose, // = false,
         bool startVerbose, // = false,
         bool passMemExc, // = false,
         bool sanityCheck, // = false
         bool skipLoading = false) :
    GFM<index_t>(in,
                 altdb,
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
                 loadSpliceSites,
                 verbose,
                 startVerbose,
                 passMemExc,
                 sanityCheck,
                 skipLoading),
    _in5(NULL),
    _in6(NULL)
    {
        _in5Str = in + ".5." + gfm_ext;
        _in6Str = in + ".6." + gfm_ext;
    }
	
	/// Construct an Ebwt from the given header parameters and string
	/// vector, optionally using a blockwise suffix sorter with the
	/// given 'bmax' and 'dcv' parameters.  The string vector is
	/// ultimately joined and the joined string is passed to buildToDisk().
	template<typename TStr>
	HGFM(
         TStr& s,
         bool packed,
         int needEntireReverse,
         int32_t lineRate,
         int32_t offRate,
         int32_t ftabChars,
         int32_t localOffRate,
         int32_t localFtabChars,
         int nthreads,
         const string& snpfile,
         const string& htfile,
         const string& ssfile,
         const string& exonfile,
         const string& svfile,
         const string& outfile,   // base filename for GFM files
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
    
	~HGFM() {
        clearLocalGFMs();
	}
    
    /**
	 * Load this Ebwt into memory by reading it in from the _in1 and
	 * _in2 streams.
	 */
	void loadIntoMemory(
                        int needEntireReverse,
                        bool loadSASamp,
                        bool loadFtab,
                        bool loadRstarts,
                        bool loadNames,
                        bool verbose)
	{
		readIntoMemory(
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
                        int needEntireRev,
                        bool loadSASamp,
                        bool loadFtab,
                        bool loadRstarts,
                        bool justHeader,
                        GFMParams<index_t> *params,
                        bool mmSweep,
                        bool loadNames,
                        bool startVerbose);
	
	/**
	 * Frees memory associated with the Ebwt.
	 */
	void evictFromMemory() {
		assert(PARENT_CLASS::isInMemory());
		clearLocalGFMs();
		PARENT_CLASS::evictFromMemory();		
	}
	
	/**
	 * Sanity-check various pieces of the Ebwt
	 */
	void sanityCheckAll(int reverse) const {
		PARENT_CLASS::sanityCheckAll(reverse);
		for(size_t tidx = 0; tidx < _localGFMs.size(); tidx++) {
			for(size_t local_idx = 0; local_idx < _localGFMs[tidx].size(); local_idx++) {
				assert(_localGFMs[tidx][local_idx] != NULL);
				_localGFMs[tidx][local_idx]->sanityCheckAll(reverse);
			}
		}
	}
    
    const LocalGFM<local_index_t, index_t>* getLocalGFM(index_t tidx, index_t offset) const {
        assert_lt(tidx, _localGFMs.size());
        const EList<LocalGFM<local_index_t, index_t>*>& localGFMs = _localGFMs[tidx];
        index_t offsetidx = offset / local_index_interval;
        if(offsetidx >= localGFMs.size()) {
            return NULL;
        } else {
            return localGFMs[offsetidx];
        }
    }
    
    const LocalGFM<local_index_t, index_t>* prevLocalGFM(const LocalGFM<local_index_t, index_t>* currLocalGFM) const {
        assert(currLocalGFM != NULL);
        index_t tidx = currLocalGFM->_tidx;
        index_t offset = currLocalGFM->_localOffset;
        if(offset < local_index_interval) {
            return NULL;
        } else {
            return getLocalGFM(tidx, offset - local_index_interval);
        }
    }
    
    const LocalGFM<local_index_t, index_t>* nextLocalGFM(const LocalGFM<local_index_t, index_t>* currLocalGFM) const {
        assert(currLocalGFM != NULL);
        index_t tidx = currLocalGFM->_tidx;
        index_t offset = currLocalGFM->_localOffset;
        return getLocalGFM(tidx, offset + local_index_interval);
    }
	
	void clearLocalGFMs() {
		for(size_t tidx = 0; tidx < _localGFMs.size(); tidx++) {
			for(size_t local_idx = 0; local_idx < _localGFMs[tidx].size(); local_idx++) {
				assert(_localGFMs[tidx][local_idx] != NULL);
				delete _localGFMs[tidx][local_idx];
			}
			
			_localGFMs[tidx].clear();
		}
		
		_localGFMs.clear();
	}
	

public:
	index_t                                  _nrefs;      /// the number of reference sequences
	EList<index_t>                           _refLens;    /// approx lens of ref seqs (excludes trailing ambig chars)
	
	EList<EList<LocalGFM<local_index_t, index_t>*> > _localGFMs;
	index_t                                  _nlocalGFMs;
    //index_t                                  _local_index_interval;
	
	FILE                                     *_in5;    // input fd for primary index file
	FILE                                     *_in6;    // input fd for secondary index file
	string                                   _in5Str;
	string                                   _in6Str;
	
	char                                     *mmFile5_;
	char                                     *mmFile6_;
    
private:
    struct ThreadParam {
        // input
        SString<char>                s;
        EList<ALT<index_t> >         alts;
        EList<Haplotype<index_t> >   haplotypes;
        bool                         bigEndian;
        index_t                      local_offset;
        index_t                      curr_sztot;
        EList<RefRecord>             conv_local_szs;
        index_t                      local_sztot;
        index_t                      index_size;
        string                       file;
        EList<index_t>               sa;
        index_t                      dcv;
        index_t                      seed;
        
        // output
        RefGraph<index_t>*           rg;
        PathGraph<index_t>*          pg;
        
        // communication
        bool                         done;
        bool                         last;
        bool                         mainThread;
    };
    static void gbwt_worker(void* vp);
};

    
template <typename index_t, typename local_index_t>
void HGFM<index_t, local_index_t>::gbwt_worker(void* vp)
{
    ThreadParam& tParam = *(ThreadParam*)vp;
    while(!tParam.last) {
        if(tParam.mainThread) {
            assert(!tParam.done);
            if(tParam.s.length() <= 0) {
                tParam.done = true;
                return;
            }
        } else {
            while(tParam.done) {
                if(tParam.last) return;
#if defined(_TTHREAD_WIN32_)
                Sleep(1);
#elif defined(_TTHREAD_POSIX_)
                const static timespec ts = {0, 1000000};  // 1 millisecond
                nanosleep(&ts, NULL);
#endif
            }
            if(tParam.s.length() <= 0) {
                tParam.done = true;
                continue;
            }
        }
        if(tParam.alts.empty()) {
            KarkkainenBlockwiseSA<SString<char> > bsa(
                                                      tParam.s,
                                                      (index_t)(tParam.s.length()+1),
                                                      1,
                                                      tParam.dcv,
                                                      tParam.seed,
                                                      false,  /* this->_sanity */
                                                      false,  /* this->_passMemExc */
                                                      false); /* this->_verbose */
            assert(bsa.suffixItrIsReset());
            assert_eq(bsa.size(), tParam.s.length()+1);
            tParam.sa.clear();
            for(index_t i = 0; i < bsa.size(); i++) {
                tParam.sa.push_back(bsa.nextSuffix());
            }
        } else {
            while(true) {
                tParam.rg = new RefGraph<index_t>(
                                                  tParam.s,
                                                  tParam.conv_local_szs,
                                                  tParam.alts,
                                                  tParam.haplotypes,
                                                  tParam.file,
                                                  1,        /* num threads */
                                                  false);   /* verbose? */
                tParam.pg = new PathGraph<index_t>(
                                                   *tParam.rg,
                                                   tParam.file,
                                                   1,         /* num threads */
                                                   false);    /* verbose? */
                
                if(!tParam.pg->generateEdges(*tParam.rg)) {
                    cerr << "An error occurred - generateEdges" << endl;
                    throw 1;
                }
                if(tParam.pg->getNumEdges() > local_max_gbwt) {
                    delete tParam.pg; tParam.pg = NULL;
                    delete tParam.rg; tParam.rg = NULL;
                    if(tParam.alts.size() <= 1) {
                        tParam.alts.clear();
                    } else {
                        for(index_t s = 2; s < tParam.alts.size(); s += 2) {
                            tParam.alts[s >> 1] = tParam.alts[s];
                        }
                        tParam.alts.resize(tParam.alts.size() >> 1);
                        tParam.haplotypes.clear();
                        for(index_t a = 0; a < tParam.alts.size(); a++) {
                            const ALT<index_t>& alt = tParam.alts[a];
                            if(!alt.snp()) continue;
                            tParam.haplotypes.expand();
                            tParam.haplotypes.back().left = alt.pos;
                            if(alt.deletion()) {
                                tParam.haplotypes.back().right = alt.pos + alt.len - 1;
                            } else {
                                tParam.haplotypes.back().right = alt.pos;
                            }
                            tParam.haplotypes.back().alts.clear();
                            tParam.haplotypes.back().alts.push_back(a);
                        }
                    }
                    continue;
                }
                break;
            }
        }
        tParam.done = true;
        if(tParam.mainThread) break;
    }
}
    
/// Construct an Ebwt from the given header parameters and string
/// vector, optionally using a blockwise suffix sorter with the
/// given 'bmax' and 'dcv' parameters.  The string vector is
/// ultimately joined and the joined string is passed to buildToDisk().
template <typename index_t, typename local_index_t>
template <typename TStr>
HGFM<index_t, local_index_t>::HGFM(
                                   TStr& s,
                                   bool packed,
                                   int needEntireReverse,
                                   int32_t lineRate,
                                   int32_t offRate,
                                   int32_t ftabChars,
                                   int32_t localOffRate,
                                   int32_t localFtabChars,
                                   int nthreads,
                                   const string& snpfile,
                                   const string& htfile,
                                   const string& ssfile,
                                   const string& exonfile,
                                   const string& svfile,
                                   const string& outfile,   // base filename for EBWT files
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
    GFM<index_t>(s,
                 packed,
                 needEntireReverse,
                 lineRate,
                 offRate,
                 ftabChars,
                 nthreads,
                 snpfile,
                 htfile,
                 ssfile,
                 exonfile,
                 svfile,
                 outfile,
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
    _in5Str = outfile + ".5." + gfm_ext;
    _in6Str = outfile + ".6." + gfm_ext;
    
    int32_t local_lineRate;
    if(snpfile == "" && ssfile == "" && exonfile == "") {
        local_lineRate = local_lineRate_fm;
    } else {
        local_lineRate = local_lineRate_gfm;
    }
    
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
    
    // Split the whole genome into a set of local indexes
    _nrefs = 0;
    _nlocalGFMs = 0;
    
    index_t cumlen = 0;
    typedef EList<RefRecord, 1> EList_RefRecord;
    ELList<EList_RefRecord> all_local_recs;
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
                _nlocalGFMs++;
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
    index_t temp_nlocalGFMs = 0;
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
        temp_nlocalGFMs += ref_local_recs.size();
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
    assert_eq(_nlocalGFMs, temp_nlocalGFMs);
#endif
    
    uint32_t be = this->toBe();
    assert(fout5.good());
    assert(fout6.good());
    
    // When building an Ebwt, these header parameters are known
    // "up-front", i.e., they can be written to disk immediately,
    // before we join() or buildToDisk()
    writeI32(fout5, 1, be); // endian hint for priamry stream
    writeI32(fout6, 1, be); // endian hint for secondary stream
    writeIndex<index_t>(fout5, _nlocalGFMs, be); // number of local Ebwts
    writeI32(fout5, local_lineRate,  be); // 2^lineRate = size in bytes of 1 line
    writeI32(fout5, 2, be); // not used
    writeI32(fout5, (int32_t)localOffRate,   be); // every 2^offRate chars is "marked"
    writeI32(fout5, (int32_t)localFtabChars, be); // number of 2-bit chars used to address ftab
    int32_t flags = 1;
    if(this->_gh._entireReverse) flags |= GFM_ENTIRE_REV;
    writeI32(fout5, -flags, be); // BTL: chunkRate is now deprecated
    
    assert_gt(this->_nthreads, 0);
    AutoArray<tthread::thread*> threads(this->_nthreads - 1);    
    EList<ThreadParam> tParams;
    for(index_t t = 0; t < (index_t)this->_nthreads; t++) {
        tParams.expand();
        tParams.back().s.clear();
        tParams.back().rg = NULL;
        tParams.back().pg = NULL;
        tParams.back().file = outfile;
        tParams.back().done = true;
        tParams.back().last = false;
        tParams.back().dcv = 1024;
        tParams.back().seed = seed;
        if(t + 1 < (index_t)this->_nthreads) {
            tParams.back().mainThread = false;
            threads[t] = new tthread::thread(gbwt_worker, (void*)&tParams.back());
        } else {
            tParams.back().mainThread = true;
        }
    }
    
    // build local FM indexes
    index_t curr_sztot = 0;
    EList<ALT<index_t> > alts;
    for(size_t tidx = 0; tidx < _refLens.size(); tidx++) {
        index_t refLen = _refLens[tidx];
        index_t local_offset = 0;
        _localGFMs.expand();
        assert_lt(tidx, _localGFMs.size());
        while(local_offset < refLen) {
            index_t t = 0;
            while(local_offset < refLen && t < (index_t)this->_nthreads) {
                assert_lt(t, tParams.size());
                ThreadParam& tParam = tParams[t];
                
                tParam.index_size = std::min<index_t>(refLen - local_offset, local_index_size);
                assert_lt(tidx, all_local_recs.size());
                assert_lt(local_offset / local_index_interval, all_local_recs[tidx].size());
                EList_RefRecord& local_szs = all_local_recs[tidx][local_offset / local_index_interval];
                
                tParam.conv_local_szs.clear();
                index_t local_len = 0, local_sztot = 0, local_sztot_interval = 0;
                for(size_t i = 0; i < local_szs.size(); i++) {
                    assert(local_szs[i].off != 0 || local_szs[i].len != 0);
                    assert(i != 0 || local_szs[i].first);
                    tParam.conv_local_szs.push_back(local_szs[i]);
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
                
                // Extract sequence corresponding to this local index
                tParam.s.resize(local_sztot);
                if(refparams.reverse == REF_READ_REVERSE) {
                    tParam.s.install(s.buf() + s.length() - curr_sztot - local_sztot, local_sztot);
                } else {
                    tParam.s.install(s.buf() + curr_sztot, local_sztot);
                }
                
                // Extract ALTs corresponding to this local index
                index_t firstSNP = (index_t)INDEX_MAX;
                tParam.alts.clear();
                ALT<index_t> alt;
                alt.pos = curr_sztot;
                index_t alt_i = (index_t)this->_alts.bsearchLoBound(alt);
                for(; alt_i < this->_alts.size(); alt_i++) {
                    const ALT<index_t>& alt = this->_alts[alt_i];
                    if(alt.snp()) {
                        if(alt.mismatch()) {
                            if(curr_sztot + local_sztot <= alt.pos) break;
                        } else if(alt.insertion()) {
                            if(curr_sztot + local_sztot < alt.pos) break;
                        } else {
                            assert(alt.deletion());
                            if(curr_sztot + local_sztot < alt.pos + alt.len) break;
                        }
                        if(curr_sztot <= alt.pos) {
                            tParam.alts.push_back(alt);
                            tParam.alts.back().pos -= curr_sztot;
                            if(firstSNP == (index_t)OFF_MASK) {
                                firstSNP = alt_i;
                            }
                        }
                    } else if(alt.splicesite()) {
                        if(alt.excluded) continue;
                        if(curr_sztot + local_sztot <= alt.right + 1) continue;
                        if(curr_sztot <= alt.left) {
                            tParam.alts.push_back(alt);
                            tParam.alts.back().left -= curr_sztot;
                            tParam.alts.back().right -= curr_sztot;
                        }
                    } else {
                        assert(false);
                    }
                }
                
                // Extract haplotypes
                tParam.haplotypes.clear();
                Haplotype<index_t> haplotype;
                haplotype.left = curr_sztot;
                index_t haplotpye_i = (index_t)this->_haplotypes.bsearchLoBound(haplotype);
                for(; haplotpye_i < this->_haplotypes.size(); haplotpye_i++) {
                    const Haplotype<index_t>& haplotype = this->_haplotypes[haplotpye_i];
                    if(curr_sztot + local_sztot <= haplotype.right) continue;
                    if(curr_sztot <= haplotype.left) {
                        tParam.haplotypes.push_back(haplotype);
                        tParam.haplotypes.back().left -= curr_sztot;
                        tParam.haplotypes.back().right -= curr_sztot;
                        assert_neq(firstSNP, (index_t)INDEX_MAX);
                        for(index_t a = 0; a < tParam.haplotypes.back().alts.size(); a++) {
                            tParam.haplotypes.back().alts[a] -= firstSNP;
                        }
                    }
                }
                
                tParam.local_offset = local_offset;
                tParam.curr_sztot = curr_sztot;
                tParam.local_sztot = local_sztot;
                
                assert(tParam.rg == NULL);
                assert(tParam.pg == NULL);
                tParam.done = false;
                
                curr_sztot += local_sztot_interval;
                local_offset += local_index_interval;
                
                t++;
            }
            
            if(!tParams.back().done) {
                gbwt_worker((void*)&tParams.back());
            }
            
            for(index_t t2 = 0; t2 < t; t2++) {
                ThreadParam& tParam = tParams[t2];
                while(!tParam.done) {
#if defined(_TTHREAD_WIN32_)
                    Sleep(1);
#elif defined(_TTHREAD_POSIX_)
                    const static timespec ts = {0, 1000000};  // 1 millisecond
                    nanosleep(&ts, NULL);
#endif
                }
                
                LocalGFM<local_index_t, index_t>(
                                                 tParam.s,
                                                 tParam.sa,
                                                 tParam.pg,
                                                 (index_t)tidx,
                                                 tParam.local_offset,
                                                 tParam.curr_sztot,
                                                 tParam.alts,
                                                 tParam.index_size,
                                                 packed,
                                                 needEntireReverse,
                                                 local_lineRate,
                                                 localOffRate,      // suffix-array sampling rate
                                                 localFtabChars,    // number of chars in initial arrow-pair calc
                                                 outfile,           // basename for .?.ebwt files
                                                 fw,                 // fw
                                                 dcv,                // difference-cover period
                                                 tParam.conv_local_szs,     // list of reference sizes
                                                 tParam.local_sztot,        // total size of all unambiguous ref chars
                                                 refparams,          // reference read-in parameters
                                                 seed,               // pseudo-random number generator seed
                                                 fout5,
                                                 fout6,
                                                 -1,                 // override offRate
                                                 false,              // be silent
                                                 passMemExc,         // pass exceptions up to the toplevel so that we can adjust memory settings automatically
                                                 sanityCheck);       // verify results and internal consistency
                tParam.s.clear();
                if(tParam.rg != NULL) {
                    assert(tParam.pg != NULL);
                    delete tParam.rg; tParam.rg = NULL;
                    delete tParam.pg; tParam.pg = NULL;
                }
            }
        }
    }
    assert_eq(curr_sztot, sztot);
    if(this->_nthreads > 1) {
        for(index_t i = 0; i + 1 < (index_t)this->_nthreads; i++) {
            tParams[i].last = true;
            threads[i]->join();
        }
    }
    
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
    VMSG_NL("Wrote " << fout5.tellp() << " bytes to primary GFM file: " << _in5Str.c_str());
    fout5.close();
    bool err = false;
    if(tellpSz5 > fileSize(_in5Str.c_str())) {
        err = true;
        cerr << "Index is corrupt: File size for " << _in5Str.c_str() << " should have been " << tellpSz5
        << " but is actually " << fileSize(_in5Str.c_str()) << "." << endl;
    }
    fout6.flush();
    int64_t tellpSz6 = (int64_t)fout6.tellp();
    VMSG_NL("Wrote " << fout6.tellp() << " bytes to secondary GFM file: " << _in6Str.c_str());
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
void HGFM<index_t, local_index_t>::readIntoMemory(
                                                  int needEntireRev,
                                                  bool loadSASamp,
                                                  bool loadFtab,
                                                  bool loadRstarts,
                                                  bool justHeader,
                                                  GFMParams<index_t> *params,
                                                  bool mmSweep,
                                                  bool loadNames,
                                                  bool startVerbose)
{
    PARENT_CLASS::readIntoMemory(needEntireRev,
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
    size_t bytesRead = 0, bytesRead2 = 4;
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
	
	_nlocalGFMs        = readIndex<index_t>(_in5, switchEndian); bytesRead += sizeof(index_t);
	int32_t lineRate  = readI32(_in5, switchEndian); bytesRead += 4;
	readI32(_in5, switchEndian); bytesRead += 4;
	int32_t offRate   = readI32(_in5, switchEndian); bytesRead += 4;
	// TODO: add isaRate to the actual file format (right now, the
	// user has to tell us whether there's an ISA sample and what the
	// sampling rate is.
	int32_t ftabChars = readI32(_in5, switchEndian); bytesRead += 4;
	/*int32_t flag  =*/ readI32(_in5, switchEndian); bytesRead += 4;
    
    if(this->_verbose || startVerbose) {
        cerr << "    number of local indexes: " << _nlocalGFMs << endl
             << "    local offRate: " << offRate << endl
             << "    local ftabLen: " << (1 << (2 * ftabChars)) << endl
             << "    local ftabSz: "  << (2 << (2 * ftabChars)) << endl
        ;
    }
	
	clearLocalGFMs();
	
    index_t tidx = 0, localOffset = 0, joinedOffset = 0;
    string base = "";
	for(size_t i = 0; i < _nlocalGFMs; i++) {
		LocalGFM<local_index_t, index_t> *localGFM = new LocalGFM<local_index_t, index_t>(base,
                                                                                          NULL,
                                                                                          _in5,
                                                                                          _in6,
                                                                                          mmFile5_,
                                                                                          mmFile6_,
                                                                                          tidx,
                                                                                          localOffset,
                                                                                          joinedOffset,
                                                                                          switchEndian,
                                                                                          bytesRead,
                                                                                          bytesRead2,
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
        
		if(tidx >= _localGFMs.size()) {
			assert_eq(tidx, _localGFMs.size());
			_localGFMs.expand();
		}
		assert_eq(tidx + 1, _localGFMs.size());
		_localGFMs.back().push_back(localGFM);
	}	
		
#ifdef BOWTIE_MM
    fseek(_in5, 0, SEEK_SET);
	fseek(_in6, 0, SEEK_SET);
#else
	rewind(_in5); rewind(_in6);
#endif
}

#endif /*HGFM_H_*/
