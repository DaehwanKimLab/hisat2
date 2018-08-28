/*
 * Copyright 2018, Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>
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

#ifndef RFM_H_
#define RFM_H_

#include "hgfm.h"

/**
 * Extended Burrows-Wheeler transform data.
 * LocalEbwt is a specialized Ebwt index that represents ~64K bps
 * and therefore uses two bytes as offsets within 64K bps.
 * This class has only two additional member variables to denote the genomic sequenuce it represents:
 * (1) the contig index and (2) the offset within the contig.
 *
 */
template <typename index_t>
class LocalRFM : public GFM<index_t> {
	typedef GFM<index_t> PARENT_CLASS;
public:
	/// Construct an Ebwt from the given input file
	LocalRFM(const string& in,
             ALTDB<index_t>* altdb,
             FILE *in1,
             FILE *in2,
             char *mmFile1,
             char *mmFile2,
             bool switchEndian,
             size_t& bytesRead,
             size_t& bytesRead2,
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
             bool sanityCheck, // = false)
             bool useHaplotype) : // = false
	GFM<index_t>(in,
                 altdb,
                 NULL,
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
                 useHaplotype,
                 true)
	{
        this->_in1 = in1;
        this->_in2 = in2;
        
        this->_repeat = true;
        bool subIndex = true;
        PARENT_CLASS::readIntoMemory(
                                     needEntireReverse,
                                     loadSASamp,
                                     loadFtab,
                                     loadRstarts,
                                     false,              //justHeader
                                     &(this->_gh),
                                     mmSweep,
                                     loadNames,
                                     startVerbose,
                                     subIndex);
		
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
	LocalRFM(
             TStr& s,
             InorderBlockwiseSA<TStr>* sa,
             PathGraph<index_t>* pg,
             EList<ALT<index_t> >& alts,
             index_t local_size,
             const EList<string>& refnames,
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
             ostream& out1,
             ostream& out2,
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
		assert(out1.good());
		assert(out2.good());
        assert_eq(refnames.size(), 1);
        this->_refnames.push_back(refnames[0]);

		writeIndex<index_t>(out1, gh._len, be); // length of string (and bwt and suffix array)
        streampos headerPos = out1.tellp();
        writeIndex<index_t>(out1, 0, be); // gbwtLen
        writeIndex<index_t>(out1, 0, be); // num of nodes
        writeI32(out1, lineRate, be);
        writeI32(out1, 0, be);
        writeI32(out1, offRate, be);
        writeI32(out1, ftabChars, be);
        writeIndex<index_t>(out1, 0, be); // eftabLen
        writeI32(out1, 0, be); // flag
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
			writeIndex(out1, this->_nPat, be);
			assert_eq(this->_nPat, 1);
			this->_plen.init(new index_t[this->_nPat], this->_nPat);
			// For each pattern, set plen
			int npat = -1;
			for(size_t i = 0; i < szs.size(); i++) {
				if(szs[i].first && szs[i].len > 0) {
					if(npat >= 0) {
						writeIndex(out1, this->plen()[npat], be);
					}
					npat++;
					this->plen()[npat] = (szs[i].len + szs[i].off);
				} else {
					this->plen()[npat] += (szs[i].len + szs[i].off);
				}
			}
			assert_eq((index_t)npat, this->_nPat-1);
			writeIndex(out1, this->plen()[npat], be);
			// Write the number of fragments
			writeIndex(out1, this->_nFrag, be);
			
			if(refparams.reverse == REF_READ_REVERSE) {
				EList<RefRecord> tmp(EBWT_CAT);
                reverseRefRecords(szs, tmp, false, verbose);
				this->szsToDisk(tmp, out1, refparams.reverse);
			} else {
				this->szsToDisk(szs, out1, refparams.reverse);
			}
            
            if(alts.empty()) {
                assert(pg == NULL);
                PARENT_CLASS::buildToDisk(*sa, s, out1, out2, headerPos);
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
                PARENT_CLASS::buildToDisk(*pg, s, out1, out2, headerPos);
            }
		}
        
        // Now write reference sequence names on the end
        assert_eq(this->_refnames.size(), this->_nPat);
        for(index_t i = 0; i < this->_refnames.size(); i++) {
            out1 << this->_refnames[i].c_str();
            if(i + 1 < this->_refnames.size()) cout << endl;
            else out1 << '\0';
        }		
		out1.flush(); out2.flush();
		if(out1.fail() || out2.fail()) {
			cerr << "An error occurred writing the index to disk.  Please check if the disk is full." << endl;
			throw 1;
		}
	}
	
	/**
	 * Sanity-check various pieces of the Ebwt
	 */
	void sanityCheckAll(int reverse) const {
		if(this->_gh._len > 0) {
			PARENT_CLASS::sanityCheckAll(reverse);
		}
	}
    
    bool empty() const { return this->_gh._len == 0; }
};

/**
 * Extended Burrows-Wheeler transform data.
 * RFM is a specialized Ebwt index that represents one global index and a large set of local indexes.
 *
 */
template <typename index_t = uint32_t>
class RFM : public GFM<index_t> {
	typedef GFM<index_t> PARENT_CLASS;
public:
	/// Construct a GFM from the given input file
	RFM(const string& in,
         ALTDB<index_t>* altdb,
         RepeatDB<index_t>* repeatdb,
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
         bool useHaplotype, // = false
         bool skipLoading = false) :
    GFM<index_t>(in,
                 altdb,
                 repeatdb,
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
                 useHaplotype,
                 skipLoading),
    _in1(NULL),
    _in2(NULL)
    {
        _in1Str = in + ".1." + gfm_ext;
        _in2Str = in + ".2." + gfm_ext;
    }
	
	/// Construct a HGFM from the given header parameters and string
	/// vector, optionally using a blockwise suffix sorter with the
	/// given 'bmax' and 'dcv' parameters.  The string vector is
	/// ultimately joined and the joined string is passed to buildToDisk().
	template<typename TStr>
	RFM(
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
         const string& repeatfile,
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
         bool localIndex, // create local indexes?
         EList<RefRecord>* parent_szs,
         EList<string>* parent_refnames,
         uint32_t seed,
         int32_t overrideOffRate = -1,
         bool verbose = false,
         bool passMemExc = false,
         bool sanityCheck = false);
    
	RFM() {
        clearLocalRFMs();
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
		clearLocalRFMs();
		PARENT_CLASS::evictFromMemory();		
	}
	
	/**
	 * Sanity-check various pieces of the Ebwt
	 */
	void sanityCheckAll(int reverse) const {
		PARENT_CLASS::sanityCheckAll(reverse);
		for(size_t i = 0; i < _localRFMs.size(); i++) {
            assert(_localRFMs[i] != NULL);
            _localRFMs[i]->sanityCheckAll(reverse);
		}
	}
    
    void getReferenceNames(EList<string>& refnames) {
        for(size_t i = 0; i < _localRFMs.size(); i++) {
            assert_eq(_localRFMs[i]->refnames().size(), 1);
            refnames.push_back(_localRFMs[i]->refnames()[0]);
        }
    }
    
    void getReferenceLens(EList<size_t>& reflens) {
        for(size_t i = 0; i < _localRFMs.size(); i++) {
            reflens.push_back(_localRFMs[i]->plen()[0]);
        }
    }
    
    index_t getLocalRFM_idx(index_t readLen) {
        for(size_t i = 0; i < _localRFMLens.size(); i++) {
            if(_localRFMLens[i].first >= readLen) {
                return i;
            }
        }
        return _localRFMLens.size() - 1;
    }
    
    LocalRFM<index_t>& getLocalRFM(index_t idx) {
        assert_lt(idx, _localRFMLens.size());
        return *_localRFMs[idx];
    }
    
    RB_KmerTable& getKmertable(index_t idx) {
        assert_lt(idx, this->_repeat_kmertables.size());
        return this->_repeat_kmertables[idx];
    }

    void clearLocalRFMs() {
		for(size_t i = 0; i < _localRFMs.size(); i++) {
            assert(_localRFMs[i] != NULL);
            delete _localRFMs[i];
        }
		_localRFMs.clear();
	}
	

public:
	EList<index_t>              _refLens;    /// approx lens of ref seqs (excludes trailing ambig chars)
	EList<LocalRFM<index_t>*>   _localRFMs;
    
    EList<pair<index_t, index_t> >      _localRFMLens;
    EList<pair<streampos, streampos> >  _localRFMFilePos;
	
	FILE                        *_in1;    // input fd for primary index file
	FILE                        *_in2;    // input fd for secondary index file
	string                      _in1Str;
	string                      _in2Str;
	
	char                        *mmFile1_;
	char                        *mmFile2_;
    
private:
    struct WorkerParam {
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
        InorderBlockwiseSA<SString<char> >* sa;
        index_t                      dcv;
        index_t                      seed;
        EList<string>                refnames;

        // output
        RefGraph<index_t>*           rg;
        PathGraph<index_t>*          pg;
        
        int                          threads;
    };
    static void build_worker(void* vp);
};

    
template <typename index_t>
void RFM<index_t>::build_worker(void* vp)
{
    WorkerParam& tParam = *(WorkerParam*)vp;
    if(tParam.alts.empty()) {
        tParam.sa = NULL;
        tParam.sa = new KarkkainenBlockwiseSA<SString<char> >(
                                                              tParam.s,
                                                              (index_t)(tParam.s.length()+1),
                                                              tParam.threads,
                                                              tParam.dcv,
                                                              tParam.seed,
                                                              false,  /* this->_sanity */
                                                              false,  /* this->_passMemExc */
                                                              false); /* this->_verbose */
        assert(tParam.sa->suffixItrIsReset());
        assert_eq(tParam.sa->size(), tParam.s.length()+1);
    } else {
        while(true) {
            tParam.rg = NULL, tParam.pg = NULL;
            bool exploded = false;
            try {
                tParam.rg = new RefGraph<index_t>(
                                                  tParam.s,
                                                  tParam.conv_local_szs,
                                                  tParam.alts,
                                                  tParam.haplotypes,
                                                  tParam.file,
                                                  tParam.threads,        /* num threads */
                                                  false);   /* verbose? */
                tParam.pg = new PathGraph<index_t>(
                                                   *tParam.rg,
                                                   tParam.file,
                                                   local_max_gbwt,
                                                   1,         /* num threads */
                                                   false);    /* verbose? */
            } catch (const NongraphException& err) {
                cerr << "Warning: no variants or splice sites in this graph (" << tParam.curr_sztot << ")" << endl;
                delete tParam.rg;
                delete tParam.pg;
                tParam.alts.clear();
                continue;
            } catch (const ExplosionException& err) {
                exploded = true;
            }
            
            if(!exploded) {
                if(!tParam.pg->generateEdges(*tParam.rg)) {
                    cerr << "An error occurred - generateEdges" << endl;
                    throw 1;
                }
                exploded = tParam.pg->getNumEdges() > local_max_gbwt;
            }
            if(exploded) {
                cerr << "Warning: a local graph exploded (offset: " << tParam.curr_sztot << ", length: " << tParam.local_sztot << ")" << endl;
                
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
}
    
/// Construct a GFM from the given header parameters and string
/// vector, optionally using a blockwise suffix sorter with the
/// given 'bmax' and 'dcv' parameters.  The string vector is
/// ultimately joined and the joined string is passed to buildToDisk().
template <typename index_t>
template <typename TStr>
RFM<index_t>::RFM(
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
                  const string& repeatfile,
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
                  bool localIndex,
                  EList<RefRecord>* parent_szs,
                  EList<string>* parent_refnames,
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
                 repeatfile,
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
                 parent_szs,
                 parent_refnames,
                 seed,
                 overrideOffRate,
                 verbose,
                 passMemExc,
                 sanityCheck),
    _in1(NULL),
    _in2(NULL)
{
    _in1Str = outfile + ".1." + gfm_ext;
    _in2Str = outfile + ".2." + gfm_ext;
    
    // Open output files
    ofstream fout1(_in1Str.c_str(), ios::binary);
    if(!fout1.good()) {
        cerr << "Could not open index file for writing: \"" << _in1Str.c_str() << "\"" << endl
        << "Please make sure the directory exists and that permissions allow writing by" << endl
        << "HISAT2." << endl;
        throw 1;
    }
    ofstream fout2(_in2Str.c_str(), ios::binary);
    if(!fout2.good()) {
        cerr << "Could not open index file for writing: \"" << _in2Str.c_str() << "\"" << endl
        << "Please make sure the directory exists and that permissions allow writing by" << endl
        << "HISAT2." << endl;
        throw 1;
    }
    
    _localRFMLens.resizeExact(szs.size());
    for(size_t i = 0; i < _localRFMLens.size(); i++) {
        _localRFMLens[i].first = numeric_limits<index_t>::max();
        _localRFMLens[i].second = 0;
    }
    const EList<Repeat<index_t> >& repeats = this->_repeatdb.repeats();
    for(size_t i = 0; i < repeats.size(); i++) {
        index_t id = repeats[i].repID;
        index_t len = repeats[i].repLen;
        assert_lt(id, _localRFMLens.size());
        if(_localRFMLens[id].first > len) {
            _localRFMLens[id].first = len;
        }
        if(_localRFMLens[id].second < len) {
            _localRFMLens[id].second = len;
        }
    }
    
    // Split the whole genome into a set of local indexes
    uint32_t be = this->toBe();
    assert(fout1.good());
    assert(fout2.good());
    
    // When building an Ebwt, these header parameters are known
    // "up-front", i.e., they can be written to disk immediately,
    // before we join() or buildToDisk()
    writeI32(fout1, 1, be); // endian hint for priamry stream
    writeI32(fout2, 1, be); // endian hint for secondary stream
    int version = GFM<index_t>::getIndexVersion();
    writeI32(fout1, version, be); // version
    index_t nLocalRFMs = szs.size();
    writeIndex<index_t>(fout1, nLocalRFMs, be); // number of local Ebwts
    for(size_t i = 0; i < _localRFMLens.size(); i++) {
        writeIndex<index_t>(fout1, _localRFMLens[i].first, be);
        writeIndex<index_t>(fout1, _localRFMLens[i].second, be);
    }
    streampos filepos = fout1.tellp();
    _localRFMFilePos.resizeExact(szs.size());
    for(size_t i = 0; i < _localRFMFilePos.size(); i++) {
        writeIndex<uint64_t>(fout1, 0, be);
        writeIndex<uint64_t>(fout1, 0, be);
    }
    
    assert_gt(this->_nthreads, 0);
    WorkerParam tParam;
    tParam.rg = NULL;
    tParam.pg = NULL;
    tParam.sa = NULL;
    tParam.file = outfile;
    tParam.dcv = 1024;
    tParam.seed = seed;
    tParam.threads = this->_nthreads;
    
    // build local FM indexes
    assert_eq(szs.size(), this->_refnames.size());
    index_t curr_sztot = 0;
    EList<ALT<index_t> > alts;
    for(size_t tidx = 0; tidx < szs.size(); tidx++) {
        assert(szs[tidx].first);
        index_t refLen = szs[tidx].len;
        _localRFMs.expand();
        assert_lt(tidx, _localRFMs.size());
        
        tParam.index_size = refLen;
        tParam.conv_local_szs.clear();
        tParam.conv_local_szs.push_back(szs[tidx]);
        tParam.refnames.clear();
        tParam.refnames.push_back(this->_refnames[tidx]);
        
        // Extract sequence corresponding to this local index
        tParam.s.resize(refLen);
        if(refparams.reverse == REF_READ_REVERSE) {
            tParam.s.install(s.buf() + s.length() - curr_sztot - refLen, refLen);
        } else {
            tParam.s.install(s.buf() + curr_sztot, refLen);
        }
        
#if 0
        // Extract ALTs corresponding to this local index
        map<index_t, index_t> alt_map;
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
                    alt_map[alt_i] = (index_t)tParam.alts.size();
                    tParam.alts.push_back(alt);
                    tParam.alts.back().pos -= curr_sztot;
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
                assert(alt.exon());
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
                for(index_t a = 0; a < tParam.haplotypes.back().alts.size(); a++) {
                    index_t alt_i = tParam.haplotypes.back().alts[a];
                    if(alt_map.find(alt_i) == alt_map.end()) {
                        assert(false);
                        tParam.haplotypes.pop_back();
                        break;
                    }
                    tParam.haplotypes.back().alts[a] = alt_map[alt_i];
                }
            }
        }
#endif
        
        tParam.local_offset = 0;
        tParam.curr_sztot = curr_sztot;
        tParam.local_sztot = refLen;
        
        assert(tParam.rg == NULL);
        assert(tParam.pg == NULL);
        assert(tParam.sa == NULL);
        curr_sztot += refLen;
        
        build_worker(&tParam);
        
        _localRFMFilePos[tidx].first = fout1.tellp();
        _localRFMFilePos[tidx].second = fout2.tellp();
        
        LocalRFM<index_t>(
                          tParam.s,
                          tParam.sa,
                          tParam.pg,
                          tParam.alts,
                          tParam.index_size,
                          tParam.refnames,
                          packed,
                          needEntireReverse,
                          lineRate,
                          offRate,               // suffix-array sampling rate
                          ftabChars,             // number of chars in initial arrow-pair calc
                          outfile,               // basename for .?.ebwt files
                          fw,                    // fw
                          dcv,                   // difference-cover period
                          tParam.conv_local_szs, // list of reference sizes
                          tParam.local_sztot,    // total size of all unambiguous ref chars
                          refparams,             // reference read-in parameters
                          seed,                  // pseudo-random number generator seed
                          fout1,
                          fout2,
                          -1,                    // override offRate
                          false,                 // be silent
                          passMemExc,            // pass exceptions up to the toplevel so that we can adjust memory settings automatically
                          sanityCheck);          // verify results and internal consistency
        if(tParam.rg != NULL) {
            assert(tParam.pg != NULL);
            delete tParam.rg; tParam.rg = NULL;
            delete tParam.pg; tParam.pg = NULL;
        }
        if(tParam.sa != NULL) {
            delete tParam.sa; tParam.sa = NULL;
        }
    }
    assert_eq(curr_sztot, sztot);
    
    streampos origpos = fout1.tellp();
    fout1.seekp(filepos);
    for(size_t i = 0; i < _localRFMFilePos.size(); i++) {
        writeIndex<uint64_t>(fout1, _localRFMFilePos[i].first, be);
        writeIndex<uint64_t>(fout1, _localRFMFilePos[i].second, be);
    }
    fout1.seekp(origpos);
    
    fout1 << '\0';
    fout1.flush(); fout2.flush();
    if(fout1.fail() || fout2.fail()) {
        cerr << "An error occurred writing the index to disk.  Please check if the disk is full." << endl;
        throw 1;
    }
    VMSG_NL("Returning from initFromVector");
    
    // Close output files
    fout1.flush();
    int64_t tellpSz1 = (int64_t)fout1.tellp();
    VMSG_NL("Wrote " << fout1.tellp() << " bytes to primary GFM file: " << _in1Str.c_str());
    fout1.close();
    bool err = false;
    if(tellpSz1 > fileSize(_in1Str.c_str())) {
        err = true;
        cerr << "Index is corrupt: File size for " << _in1Str.c_str() << " should have been " << tellpSz1
        << " but is actually " << fileSize(_in1Str.c_str()) << "." << endl;
    }
    fout2.flush();
    int64_t tellpSz2 = (int64_t)fout2.tellp();
    VMSG_NL("Wrote " << fout2.tellp() << " bytes to secondary GFM file: " << _in2Str.c_str());
    fout2.close();
    if(tellpSz2 > fileSize(_in2Str.c_str())) {
        err = true;
        cerr << "Index is corrupt: File size for " << _in2Str.c_str() << " should have been " << tellpSz2
        << " but is actually " << fileSize(_in2Str.c_str()) << "." << endl;
    }
    if(err) {
        cerr << "Please check if there is a problem with the disk or if disk is full." << endl;
        throw 1;
    }
    // Reopen as input streams
    VMSG_NL("Re-opening _in1 and _in2 as input streams");
    if(this->_sanity) {
        VMSG_NL("Sanity-checking ht2");
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
    VMSG_NL("Returning from HGFM constructor");
}

    
/**
 * Read an Ebwt from file with given filename.
 */
template <typename index_t>
void RFM<index_t>::readIntoMemory(
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
    bool switchEndian; // dummy; caller doesn't care
#ifdef BOWTIE_MM
	char *mmFile[] = { NULL, NULL };
#endif
	if(_in1Str.length() > 0) {
		if(this->_verbose || startVerbose) {
			cerr << "  About to open input files: ";
			logTime(cerr);
		}
        // Initialize our primary and secondary input-stream fields
		if(_in1 != NULL) fclose(_in1);
		if(this->_verbose || startVerbose) cerr << "Opening \"" << _in1Str.c_str() << "\"" << endl;
		if((_in1 = fopen(_in1Str.c_str(), "rb")) == NULL) {
			cerr << "Could not open index file " << _in1Str.c_str() << endl;
		}
		if(loadSASamp) {
			if(_in2 != NULL) fclose(_in2);
			if(this->_verbose || startVerbose) cerr << "Opening \"" << _in2Str.c_str() << "\"" << endl;
			if((_in2 = fopen(_in2Str.c_str(), "rb")) == NULL) {
				cerr << "Could not open index file " << _in2Str.c_str() << endl;
			}
		}
		if(this->_verbose || startVerbose) {
			cerr << "  Finished opening input files: ";
			logTime(cerr);
		}
		
#ifdef BOWTIE_MM
		if(this->_useMm /*&& !justHeader*/) {
			const char *names[] = {_in1Str.c_str(), _in2Str.c_str()};
            int fds[] = { fileno(_in1), fileno(_in2) };
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
			mmFile1_ = mmFile[0];
			mmFile2_ = loadSASamp ? mmFile[1] : NULL;
		}
#endif
	}
#ifdef BOWTIE_MM
	else if(this->_useMm && !justHeader) {
		mmFile[0] = mmFile1_;
		mmFile[1] = mmFile2_;
	}
	if(this->_useMm && !justHeader) {
		assert(mmFile[0] == mmFile1_);
		assert(mmFile[1] == mmFile2_);
	}
#endif
	
	if(this->_verbose || startVerbose) {
		cerr << "  Reading header: ";
		logTime(cerr);
	}
	
	// Read endianness hints from both streams
    size_t bytesRead = 0, bytesRead2 = 4;
	switchEndian = false;
	uint32_t one = readU32(_in1, switchEndian); // 1st word of primary stream
	bytesRead += 4;
	if(loadSASamp) {
#ifndef NDEBUG
		assert_eq(one, readU32(_in2, switchEndian)); // should match!
#else
		readU32(_in2, switchEndian);
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
	
    readI32(_in1, switchEndian); bytesRead += 4;
	index_t nlocalRFMs = readIndex<index_t>(_in1, switchEndian); bytesRead += sizeof(index_t);
    
    _localRFMLens.resizeExact(nlocalRFMs);
    for(size_t i = 0; i < _localRFMLens.size(); i++) {
        _localRFMLens[i].first = readIndex<index_t>(_in1, switchEndian);
        _localRFMLens[i].second = readIndex<index_t>(_in1, switchEndian);
    }
    
    _localRFMFilePos.resizeExact(nlocalRFMs);
    for(size_t i = 0; i < _localRFMFilePos.size(); i++) {
        _localRFMFilePos[i].first = readIndex<uint64_t>(_in1, switchEndian);
        _localRFMFilePos[i].second = readIndex<uint64_t>(_in1, switchEndian);
    }
	
	clearLocalRFMs();
	
    string base = "";
	for(size_t i = 0; i < nlocalRFMs; i++) {
		LocalRFM<index_t> *localRFM = new LocalRFM<index_t>(base,
                                                            NULL,
                                                            _in1,
                                                            _in2,
                                                            mmFile1_,
                                                            mmFile2_,
                                                            switchEndian,
                                                            bytesRead,
                                                            bytesRead2,
                                                            needEntireRev,
                                                            this->fw_,
                                                            -1, // overrideOffRate
                                                            -1, // offRatePlus
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
                                                            this->_sanity,
                                                            false); // use haplotypes?
        
		_localRFMs.push_back(localRFM);
	}	
		
#ifdef BOWTIE_MM
    fseek(_in1, 0, SEEK_SET);
	fseek(_in2, 0, SEEK_SET);
#else
	rewind(_in1);
    rewind(_in2);
#endif
}

#endif /*RFM_H_*/
