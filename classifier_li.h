/*
 * Copyright 2014, Daehwan Kim <infphilo@gmail.com>
 *
 * This file is part of HISAT.
 *
 * HISAT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CLASSIFIER_H_
#define CLASSIFIER_H_

#include "hi_aligner.h"

struct IDCount {
	uint32_t id;
	uint32_t count;
	uint32_t weightedCount;
	uint32_t timeStamp ;
};


/**
 * With a hierarchical indexing, SplicedAligner provides several alignment strategies
 * , which enable effective alignment of RNA-seq reads
 */
template <typename index_t, typename local_index_t>
class Classifier : public HI_Aligner<index_t, local_index_t> {
    
public:
	
	/**
	 * Initialize with index.
	 */
	Classifier(const Ebwt<index_t>& ebwt,
               const EList<string>& refnames) :
    HI_Aligner<index_t, local_index_t>(
                                       ebwt,
                                       0,    // don't make use of splice sites found by earlier reads
                                       true), // no spliced alignment
    _refnames(refnames)
    {
    }
    
    ~Classifier() {
    }
    
    /**
     * Aligns a read or a pair
     * This funcion is called per read or pair
     */
    virtual
    int go(
           const Scoring&           sc,
           const Ebwt<index_t>&     ebwtFw,
           const Ebwt<index_t>&     ebwtBw,
           const BitPairReference&  ref,
           SwAligner&               swa,
           SpliceSiteDB&            ssdb,
           WalkMetrics&             wlm,
           PerReadMetrics&          prm,
           SwMetrics&               swm,
           HIMetrics&               him,
           RandomSource&            rnd,
           AlnSinkWrap<index_t>&    sink)
    {
        _speciesMap.clear();
        _genusMap.clear();
        
        const index_t increment = 10;
        const index_t minPartialLen = 22;
	size_t bestScore = 0, secondBestScore = 0 ;
	for(index_t rdi = 0; rdi < (this->_paired ? 2 : 1); rdi++) {
		assert(this->_rds[rdi] != NULL);
		const Read& rd = *(this->_rds[rdi]);
		index_t rdlen = rd.length();
		index_t fwi;

		bool done[2] = {false, false};
		size_t cur[2] = {0, 0 } ;
		const size_t maxDiff = ( rdlen / 2 > 2 * minPartialLen ) ? rdlen / 2 : ( 2 * minPartialLen ) ;
		while(!done[0] || !done[1]) {
			for(fwi = 0; fwi < 2; fwi++) {
				if(done[fwi]) continue;
				size_t mineFw = 0, mineRc = 0;
				bool fw = (fwi == 0);
				ReadBWTHit<index_t>& hit = this->_hits[rdi][fwi];
				bool pseudogeneStop = false, anchorStop = false;
				this->partialSearch(
						ebwtFw,
						rd,
						sc,
						fw,
						0,
						mineFw,
						mineRc,
						hit,
						rnd,
						pseudogeneStop,
						anchorStop);
				if(hit.done()) {
					done[fwi] = true;
					cur[fwi] = rdlen ;
					continue;
				}
				cur[fwi] = hit.cur() ;

				BWTHit<index_t>& lastHit = hit.getPartialHit(hit.offsetSize() - 1);
				if(lastHit.len() > increment) {
					if ( lastHit.len() < minPartialLen )
						hit.setOffset(hit.cur() - increment);
					else
						hit.setOffset( hit.cur() + 1 ) ;
				}

				if(hit.cur() + minPartialLen >= rdlen) {
					hit.done(true);
					done[fwi] = true;
					continue;
				}
			}

			//ReadBWTHit<index_t>& lastHitForward = (this->_hits[rdi][0] ).getPartialHit(hit.offsetSize() - 1).cur() ;
			//ReadBWTHit<index_t>& lastHitBackward = (this->_hits[rdi][1] ).getPartialHit(hit.offsetSize() - 1).cur() ;
			//size_t curForward = this->_hits[rdi][0].cur()  ;
			//size_t curBackward = this->_hits[rdi][1].cur() ;
			//cout<< cur[0] << " " << cur[1] << " " << done[0] << " " << done[1] << endl ;

			if ( cur[0] > cur[1] + maxDiff )
			{
				this->_hits[rdi][1].done( true ) ;	
				done[1] = true ;
			}
			else if ( cur[1] > cur[0] + maxDiff )
			{
				this->_hits[rdi][0].done( true ) ;	
				done[0] = true ;
			}
		}

		// choose fw or rc of a read
		index_t avgHitLength[2] = {0, 0};
		index_t totalHitLength[2] = {0, 0};
		for(fwi = 0; fwi < 2; fwi++) {
			ReadBWTHit<index_t>& hit = this->_hits[rdi][fwi];
			index_t numHits = 0;
			for(size_t i = 0; i < hit.offsetSize(); i++) {
				if(hit.getPartialHit(i).len() < minPartialLen) continue;
				totalHitLength[fwi] += hit.getPartialHit(i).len();
				numHits++;
			}

			if(numHits > 0) {
				avgHitLength[fwi] = totalHitLength[fwi] / numHits;
			}
		}

		if(avgHitLength[0] > avgHitLength[1]) {
			fwi = 0;
		} else {
			fwi = 1;
		}

		bool fw = (fwi == 0);
		const ReportingParams& rp = sink.reportingParams();
		ReadBWTHit<index_t>& hit = this->_hits[rdi][fwi];
		assert(hit.done());
		// choose candidate partial alignments for further alignment
		const index_t maxGenomeHitSize = rp.khits;
		index_t offsetSize = hit.offsetSize();
		this->_genomeHits.clear();

		int hitLen[100] ;
		int hitSize[100] ;
		size_t hiMap[100] ;
		for ( size_t hi = 0 ; hi < offsetSize ; ++hi )
		{
			hitLen[hi] = hit.getPartialHit(hi).len() ;
			hitSize[hi] = hit.getPartialHit(hi).size() ;
			hiMap[hi] = hi ;
		}

		// Change to quicksort in future
		for ( size_t hi = 0 ; hi < offsetSize ; ++hi )
		{
			for ( size_t hj = hi + 1 ; hj < offsetSize ; ++hj )
			{
				//if ( hitLen[ hiMap[hi] ] < hitLen[ hiMap[hj] ] )
				if ( hitSize[ hiMap[hi] ] > hitSize[ hiMap[hj] ] ) // When use size()
				{
					size_t tmp = hiMap[hi] ;
					hiMap[hi] = hiMap[hj] ;
					hiMap[hj] = tmp ;
				}
				else if ( hitSize[ hiMap[hi] ] == hitSize[ hiMap[hj] ] && hitLen[ hiMap[hi] ] < hitLen[ hiMap[hj] ] )
				{
					size_t tmp = hiMap[hi] ;
					hiMap[hi] = hiMap[hj] ;
					hiMap[hj] = tmp ;
				}
			}
		}

		size_t usedPortion = 0 ;
		size_t genomeHitCnt = 0 ;
		for(size_t hi = 0; hi < offsetSize; hi++) {
			/*index_t hj = 0;
			  for(; hj < offsetSize; hj++) {
			  BWTHit<index_t>& partialHit_j = hit.getPartialHit(hj);
			  if(partialHit_j.empty() ||
			  partialHit_j.hasGenomeCoords() ||
			  partialHit_j.len() < minPartialLen) continue;
			  else break;
			  }
			  if(hj >= offsetSize) break;
			  for(index_t hk = hj + 1; hk < offsetSize; hk++) {
			  BWTHit<index_t>& partialHit_j = hit.getPartialHit(hj);
			  BWTHit<index_t>& partialHit_k = hit.getPartialHit(hk);
			  if(partialHit_k.empty() ||
			  partialHit_k.hasGenomeCoords() ||
			  partialHit_k.len() < minPartialLen) continue;

			  if(partialHit_j.size() > partialHit_k.size() ||
			  (partialHit_j.size() == partialHit_k.size() && partialHit_j.len() < partialHit_k.len())) {
			  hj = hk;
			  }
			  }*/
			BWTHit<index_t>& partialHit = hit.getPartialHit( hiMap[ hi ] );

			if ( partialHit.len() < minPartialLen )
				continue ;
				//break ;

			assert(!partialHit.hasGenomeCoords());
			usedPortion += partialHit.len() ;
			bool straddled = false;
			this->getGenomeIdx(
					ebwtFw,
					ref,
					rnd,
					partialHit._top,
					partialHit._bot,
					fw,
					maxGenomeHitSize - genomeHitCnt, // this->_genomeHits.size(),
					hit._len - partialHit._bwoff - partialHit._len,
					partialHit._len,
					partialHit._coords,
					wlm,
					prm,
					him,
					false, // reject straddled
					straddled);
			if(!partialHit.hasGenomeCoords()) continue;
			EList<Coord>& coords = partialHit._coords;
			assert_gt(coords.size(), 0);
			//const index_t genomeHit_size = this->_genomeHits.size();
			if(genomeHitCnt + coords.size() >= maxGenomeHitSize) {
				coords.shufflePortion(0, coords.size(), rnd);
			}
			for(index_t k = 0; k < coords.size(); k++, ++genomeHitCnt) {
				//cout << genomeHitCnt << " " << maxGenomeHitSize << endl ;
				if ( genomeHitCnt >= maxGenomeHitSize )
					break ;
				const Coord& coord = coords[k] ;

				assert_lt(coord.ref(), _refnames.size());
				const string& refName = _refnames[coord.ref()];
				uint64_t id = 0;
				for(size_t ni = 0; ni < refName.length(); ni++) {
					if(refName[ni] < '0' || refName[ni] > '9') break;
					id *= 10;
					id += (refName[ni] - '0');
				}

				uint32_t speciesID = (uint32_t)(id >> 32);
				uint32_t genusID = (uint32_t)(id & 0xffffffff);

				assert_gt(partialHit.len(), 15);
				uint32_t addWeight = (uint32_t)((partialHit.len() - 15) * (partialHit.len() - 15));

				bool found = false;
				uint32_t newScore = 0 ;
				for(size_t mi = 0; mi < _speciesMap.size(); mi++) {
					if(_speciesMap[mi].id == speciesID) {
						found = true;
						if ( _speciesMap[mi].timeStamp != hi )
						{
							_speciesMap[mi].count += 1;
							_speciesMap[mi].weightedCount += addWeight;
							_speciesMap[mi].timeStamp = hi;
							//newScore = _speciesMap[mi].weightedCount ;
						}
						break;
					}
				}

				if(!found) {
					_speciesMap.expand();
					_speciesMap.back().id = speciesID;
					_speciesMap.back().count = 1;
					_speciesMap.back().timeStamp = hi ;
					_speciesMap.back().weightedCount = addWeight;
					//newScore = addWeight ;
				}

				found = false;
				for(size_t mi = 0; mi < _genusMap.size(); mi++) {
					if(_genusMap[mi].id == genusID) {
						found = true;
						if(_genusMap[mi].timeStamp != hi ) {
							_genusMap[mi].count += 1;
							_genusMap[mi].weightedCount += addWeight;
							_genusMap[mi].timeStamp = hi ;
							newScore = _genusMap[mi].weightedCount ;
						}
						break;
					}
				}

				if(!found) {
					_genusMap.expand();
					_genusMap.back().id = genusID;
					_genusMap.back().count = 1;
					_genusMap.back().weightedCount = addWeight;
					_genusMap.back().timeStamp = hi ;
					newScore = addWeight ;
				}

				// classification of bacterial sequences
#ifndef NDEBUG
				cout << this->_rds[rdi]->name << "\t"
					// << refName << "\t"
					<< speciesID << "\t"
					<< genusID << "\t"
					<< genomeHit.refoff() << "\t"
					<< genomeHit.len() << "M" << endl;
#endif
				if ( newScore > bestScore )
				{
					secondBestScore = bestScore ; 
					bestScore = newScore ;
				}
				else if ( newScore > secondBestScore )
				{
					secondBestScore = newScore ;
				}

			} // for k
			if ( genomeHitCnt >= maxGenomeHitSize )
				break ;

			//cout<< bestScore << " " << secondBestScore << " " << totalHitLength[fwi] << " " << usedPortion << endl ;
			if ( rdi == (this->_paired ? 2: 1) - 1 && bestScore > secondBestScore + 
					( totalHitLength[fwi] - usedPortion - 15 ) * ( totalHitLength[fwi] - usedPortion - 15 ) )
			{
				//cout << "saved\n" ;
				break ;
			}
		} // for hi
	} // for rdi 
        
        uint32_t speciesID = 0xffffffff, genusID = 0xffffffff;
        uint32_t speciesWeightedCount = 0, genusWeightedCount = 0;
        
        for(size_t mi = 0; mi < _speciesMap.size(); mi++) {
            if(_speciesMap[mi].weightedCount > speciesWeightedCount) {
                speciesID = _speciesMap[mi].id;
                speciesWeightedCount = _speciesMap[mi].weightedCount;
            }
        }
        
        for(size_t mi = 0; mi < _genusMap.size(); mi++) {
            if(_genusMap[mi].weightedCount > genusWeightedCount) {
                genusID = _genusMap[mi].id;
                genusWeightedCount = _genusMap[mi].weightedCount;
            }
        }
        
        if(genusID != 0xffffffff) {
            cout << this->_rds[0]->name << "\t"
            << speciesID << "\t"
            << genusID << endl;
        }
        
        return EXTEND_POLICY_FULFILLED;
    }
    
    bool getGenomeIdx(
                      const Ebwt<index_t>&       ebwt,
                      const BitPairReference&    ref,
                      RandomSource&              rnd,
                      index_t                    top,
                      index_t                    bot,
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
        index_t nelt = bot - top;
        nelt = min<index_t>(nelt, maxelt);
        coords.clear();
        him.globalgenomecoords += (bot - top);
        this->_offs.resize(nelt);
        this->_offs.fill(std::numeric_limits<index_t>::max());
        this->_sas.init(top, rdlen, EListSlice<index_t, 16>(this->_offs, 0, nelt));
        this->_gws.init(ebwt, ref, this->_sas, rnd, met);
        for(index_t off = 0; off < nelt; off++) {
            WalkResult<index_t> wr;
            this->_gws.advanceElement(
                                off,
                                ebwt,         // forward Bowtie index for walking left
                                ref,          // bitpair-encoded reference
                                this->_sas,   // SA range with offsets
                                this->_gwstate,     // GroupWalk state; scratch space
                                wr,           // put the result here
                                met,          // metrics
                                prm);         // per-read metrics
            // Coordinate of the seed hit w/r/t the pasted reference string
            coords.expand();
            coords.back().init(wr.toff, 0, fw);
        }
        
        return true;
    }
 
private:
        
    EList<string>   _refnames;
    EList<IDCount>  _speciesMap;
    EList<IDCount>  _genusMap;
};



#endif /*CLASSIFIER_H_*/
