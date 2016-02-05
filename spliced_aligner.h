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

#ifndef SPLICED_ALIGNER_H_
#define SPLICED_ALIGNER_H_

#include "hi_aligner.h"

/**
 * With a hierarchical indexing, SplicedAligner provides several alignment strategies
 * , which enable effective alignment of RNA-seq reads
 */
template <typename index_t, typename local_index_t>
class SplicedAligner : public HI_Aligner<index_t, local_index_t> {

public:
	/**
	 * Initialize with index.
	 */
	SplicedAligner(
                   const GFM<index_t>& gfm,
                   const TranscriptomePolicy& tpol,
                   bool anchorStop,
                   size_t minIntronLen,
                   size_t maxIntronLen,
                   bool secondary = false,
                   bool local = false,
                   uint64_t threads_rids_mindist = 0) :
    HI_Aligner<index_t, local_index_t>(gfm,
                                       tpol,
                                       anchorStop,
                                       minIntronLen,
                                       maxIntronLen,
                                       secondary,
                                       local,
                                       threads_rids_mindist)
    {
    }
    
    ~SplicedAligner() {
    }
    
    /**
     * Given a partial alignment of a read, try to further extend
     * the alignment bidirectionally using a combination of
     * local search, extension, and global search
     */
    virtual
    void hybridSearch(
                      const Scoring&                     sc,
                      const GFM<index_t>&                gfm,
                      const ALTDB<index_t>&              altdb,
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
                      AlnSinkWrap<index_t>&              sink);
    
    /**
     * Given a partial alignment of a read, try to further extend
     * the alignment bidirectionally using a combination of
     * local search, extension, and global search
     */
    virtual
    int64_t hybridSearch_recur(
                               const Scoring&                   sc,
                               const GFM<index_t>&              gfm,
                               const ALTDB<index_t>&            altdb,
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
                               index_t                          dep = 0);
};

/**
 * Given a partial alignment of a read, try to further extend
 * the alignment bidirectionally using a combination of
 * local search, extension, and global search
 */
template <typename index_t, typename local_index_t>
void SplicedAligner<index_t, local_index_t>::hybridSearch(
                                                          const Scoring&                 sc,
                                                          const GFM<index_t>&            gfm,
                                                          const ALTDB<index_t>&          altdb,
                                                          const BitPairReference&        ref,
                                                          SwAligner&                     swa,
                                                          SpliceSiteDB&                  ssdb,
                                                          index_t                        rdi,
                                                          bool                           fw,
                                                          WalkMetrics&                   wlm,
                                                          PerReadMetrics&                prm,
                                                          SwMetrics&                     swm,
                                                          HIMetrics&                     him,
                                                          RandomSource&                  rnd,
                                                          AlnSinkWrap<index_t>&          sink)
{
    assert_lt(rdi, 2);
    assert(this->_rds[rdi] != NULL);
    him.localatts++;
    
    // before further alignment using local search, extend the partial alignments directly
    // by comparing with the corresponding genomic sequences
    // this extension is performed without any mismatches allowed
    for(index_t hi = 0; hi < this->_genomeHits.size(); hi++) {
        GenomeHit<index_t>& genomeHit = this->_genomeHits[hi];
        index_t leftext = (index_t)INDEX_MAX, rightext = (index_t)INDEX_MAX;
        genomeHit.extend(
                         *(this->_rds[rdi]),
                         gfm,
                         ref,
                         altdb,
                         ssdb,
                         swa,
                         swm,
                         prm,
                         sc,
                         this->_minsc[rdi],
                         rnd,
                         INDEX_MAX,
                         (index_t)this->_minIntronLen,
                         (index_t)this->_maxIntronLen,
                         this->_minAnchorLen,
                         this->_minAnchorLen_noncan,
                         leftext,
                         rightext);
    }
    
    // for the candidate alignments, examine the longest (best) one first
    this->_genomeHits_done.resize(this->_genomeHits.size());
    this->_genomeHits_done.fill(false);
    for(size_t hi = 0; hi < this->_genomeHits.size(); hi++) {
        index_t hj = 0;
        for(; hj < this->_genomeHits.size(); hj++) {
            if(!this->_genomeHits_done[hj]) break;
        }
        if(hj >= this->_genomeHits.size()) break;
        for(index_t hk = hj + 1; hk < this->_genomeHits.size(); hk++) {
            if(this->_genomeHits_done[hk]) continue;
            GenomeHit<index_t>& genomeHit_j = this->_genomeHits[hj];
            GenomeHit<index_t>& genomeHit_k = this->_genomeHits[hk];
            if(genomeHit_k.hitcount() > genomeHit_j.hitcount() ||
               (genomeHit_k.hitcount() == genomeHit_j.hitcount() && genomeHit_k.len() > genomeHit_j.len())) {
                hj = hk;
            }
        }
        
        // given a candidate partial alignment, extend it bidirectionally
        him.anchoratts++;
        GenomeHit<index_t>& genomeHit = this->_genomeHits[hj];
        hybridSearch_recur(
                           sc,
                           gfm,
                           altdb,
                           ref,
                           swa,
                           ssdb,
                           rdi,
                           genomeHit,
                           genomeHit.rdoff(),
                           genomeHit.len(),
                           wlm,
                           prm,
                           swm,
                           him,
                           rnd,
                           sink);
        this->_genomeHits_done[hj] = true;
    }
}


/**
 * Given a partial alignment of a read, try to further extend
 * the alignment bidirectionally using a combination of
 * local search, extension, and global search
 */
template <typename index_t, typename local_index_t>
int64_t SplicedAligner<index_t, local_index_t>::hybridSearch_recur(
                                                                   const Scoring&                   sc,
                                                                   const GFM<index_t>&              gfm,
                                                                   const ALTDB<index_t>&            altdb,
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
                                                                   index_t                          dep)
{
    int64_t maxsc = numeric_limits<int64_t>::min();
    him.localsearchrecur++;
    assert_lt(rdi, 2);
    assert(this->_rds[rdi] != NULL);
    const Read& rd = *(this->_rds[rdi]);
    index_t rdlen = (index_t)rd.length();
    if(hit.score() < this->_minsc[rdi]) return maxsc;
    
    // if it's already examined, just return
    if(hitoff == hit.rdoff() - hit.trim5() && hitlen == hit.len() + hit.trim5() + hit.trim3()) {
        if(this->isSearched(hit, rdi)) return maxsc;
        this->addSearched(hit, rdi);
    }
    
    // for effective use of memory allocation and deallocation
    if(this->_coords.size() <= dep) {
        this->_coords.expand();
        assert_leq(this->_local_genomeHits.size(), dep);
        this->_local_genomeHits.expand();
        assert_leq(this->_spliceSites.size(), dep);
        this->_spliceSites.expand();
    }
    EList<Coord>& coords = this->_coords[dep];
    EList<GenomeHit<index_t> >& local_genomeHits = this->_local_genomeHits[dep];
    EList<SpliceSite>& spliceSites = this->_spliceSites[dep];
    
    // daehwan - for debugging purposes
#if 0
    cout << rd.name << "\t"
    << (hit.fw() ? "+" : "-") << "\t"
    << hitoff << "\t"
    << hitoff + hitlen << "\t"
    << "( " << hit.rdoff() << "\t"
    << hit.rdoff() + hit.len() << " )" << "\t"
    << hit.refoff() << "\t"
    << hit.getRightOff() << "\t"
    << hit.score() << "\t"
    << "dep: " << dep << "\t";
    Edit::print(cout, hit.edits());
    cout << endl;
#endif
    
    assert_leq(hitoff + hitlen, rdlen);
    // if this is a full alignment, report it
    if(hitoff == 0 && hitlen == rdlen) {
        if(!this->redundant(sink, rdi, hit)) {
            bool another_spliced = false;
            if(!ssdb.empty()) {
                int64_t best_score = hit.score();
                local_genomeHits.clear();
                this->_anchors_added.clear();
                
                local_genomeHits.expand();
                local_genomeHits.back() = hit;
                this->_anchors_added.push_back(0);
                
                index_t fragoff = 0, fraglen = 0, left = 0, right = 0;
                hit.getLeft(fragoff, fraglen, left);
                const index_t minMatchLen = (index_t)this->_minK;
                index_t min_left_anchor = rdlen, min_right_anchor = rdlen;
                // make use of a list of known or novel splice sites to further align the read
                if(fraglen >= minMatchLen &&
                   left >= minMatchLen &&
                   hit.trim5() == 0 &&
                   !this->_tpol.no_spliced_alignment()) {
                    spliceSites.clear();
                    ssdb.getLeftSpliceSites(hit.ref(), left + minMatchLen, minMatchLen, spliceSites);
                    for(size_t si = 0; si < spliceSites.size(); si++) {
                        const SpliceSite& ss = spliceSites[si];
                        if(!ss._fromfile && ss._readid + this->_thread_rids_mindist > rd.rdid) continue;
                        if(left + fraglen - 1 < ss.right()) continue;
                        index_t frag2off = ss.left() -  (ss.right() - left);
                        if(frag2off + 1 < hitoff) continue;
                        GenomeHit<index_t> tempHit;
                        tempHit.init(hit.fw(),
                                     0,
                                     hitoff,
                                     0, // trim5
                                     0, // trim3
                                     hit.ref(),
                                     frag2off + 1,
                                     (index_t)INDEX_MAX,
                                     this->_sharedVars);
                        if(!tempHit.compatibleWith(
                                                   hit,
                                                   (index_t)this->_minIntronLen,
                                                   (index_t)this->_maxIntronLen,
                                                   this->_tpol.no_spliced_alignment()))
                            continue;
                        int64_t minsc = max<int64_t>(this->_minsc[rdi], best_score);
                        bool combined = tempHit.combineWith(
                                                            hit,
                                                            rd,
                                                            gfm,
                                                            ref,
                                                            altdb,
                                                            ssdb,
                                                            swa,
                                                            swm,
                                                            sc,
                                                            minsc,
                                                            rnd,
                                                            (index_t)this->_minK_local,
                                                            (index_t)this->_minIntronLen,
                                                            (index_t)this->_maxIntronLen,
                                                            1,
                                                            1,
                                                            &ss);
                        if(rdi == 0) minsc = max(minsc, sink.bestUnp1());
                        else         minsc = max(minsc, sink.bestUnp2());
                        index_t leftAnchorLen = 0, nedits = 0;
                        tempHit.getLeftAnchor(leftAnchorLen, nedits);
                        if(combined &&
                           tempHit.score() >= minsc &&
                           nedits <= leftAnchorLen / 4) { // prevent (short) anchors from having many mismatches
                            if(!this->redundant(sink, rdi, tempHit)) {
                                another_spliced = true;
                                if(tempHit.score() > best_score)
                                    best_score = tempHit.score();
                                local_genomeHits.expand();
                                local_genomeHits.back() = tempHit;
                                this->_anchors_added.push_back(1);
                                index_t temp_fragoff = 0, temp_fraglen = 0, temp_left = 0;
                                tempHit.getLeft(temp_fragoff, temp_fraglen, temp_left);
                                if(temp_fraglen < min_left_anchor)
                                    min_left_anchor = temp_fraglen;
                            }
                        }
                    }
                }
                
                size_t num_local_genomeHits = local_genomeHits.size();
                for(size_t i = 0; i < num_local_genomeHits; i++) {
                    local_genomeHits[i].getRight(fragoff, fraglen, right);
                    if(local_genomeHits[i].score() < best_score) continue;
                    // make use of a list of known or novel splice sites to further align the read
                    if(fraglen >= minMatchLen &&
                       local_genomeHits[i].trim3() == 0 &&
                       !this->_tpol.no_spliced_alignment()) {
                        spliceSites.clear();
                        assert_gt(fraglen, 0);
                        ssdb.getRightSpliceSites(local_genomeHits[i].ref(), right + fraglen - minMatchLen, minMatchLen, spliceSites);
                        for(size_t si = 0; si < spliceSites.size(); si++) {
                            const GenomeHit<index_t>& canHit = local_genomeHits[i];
                            const SpliceSite& ss = spliceSites[si];
                            if(!ss._fromfile && ss._readid + this->_thread_rids_mindist > rd.rdid) continue;
                            if(right > ss.left()) continue;
                            index_t frag2off = ss.right() - ss.left() + right + fraglen - 1;
                            GenomeHit<index_t> tempHit;
                            tempHit.init(canHit.fw(),
                                         fragoff + fraglen,
                                         rdlen - fragoff - fraglen,
                                         0, // trim5
                                         0, // trim3
                                         canHit.ref(),
                                         frag2off,
                                         (index_t)INDEX_MAX,
                                         this->_sharedVars);
                            if(!canHit.compatibleWith(tempHit, (index_t)this->_minIntronLen, (index_t)this->_maxIntronLen, this->_tpol.no_spliced_alignment())) continue;
                            GenomeHit<index_t> combinedHit = canHit;
                            int64_t minsc = max<int64_t>(this->_minsc[rdi], best_score);
                            bool combined = combinedHit.combineWith(
                                                                    tempHit,
                                                                    rd,
                                                                    gfm,
                                                                    ref,
                                                                    altdb,
                                                                    ssdb,
                                                                    swa,
                                                                    swm,
                                                                    sc,
                                                                    minsc,
                                                                    rnd,
                                                                    (index_t)this->_minK_local,
                                                                    (index_t)this->_minIntronLen,
                                                                    (index_t)this->_maxIntronLen,
                                                                    1,
                                                                    1,
                                                                    &ss);
                            if(rdi == 0) minsc = max(minsc, sink.bestUnp1());
                            else         minsc = max(minsc, sink.bestUnp2());
                            index_t rightAnchorLen = 0, nedits = 0;
                            combinedHit.getRightAnchor(rightAnchorLen, nedits);
                            if(combined &&
                               combinedHit.score() >= minsc &&
                               nedits <= rightAnchorLen / 4) { // prevent (short) anchors from having many mismatches
                                if(!this->redundant(sink, rdi, combinedHit)) {
                                    another_spliced = true;
                                    if(combinedHit.score() > best_score)
                                        best_score = tempHit.score();
                                    local_genomeHits.expand();
                                    local_genomeHits.back() = combinedHit;
                                    this->_anchors_added.push_back(this->_anchors_added[i] + 1);
                                    
                                    index_t temp_fragoff = 0, temp_fraglen = 0, temp_right = 0;
                                    combinedHit.getLeft(temp_fragoff, temp_fraglen, temp_right);
                                    if(temp_fraglen < min_right_anchor)
                                        min_right_anchor = temp_fraglen;
                                }
                            }
                        }
                    }
                }
                
                assert_eq(local_genomeHits.size(), this->_anchors_added.size());
                for(size_t i = 0; i < local_genomeHits.size(); i++) {
                    const GenomeHit<index_t>& canHit = local_genomeHits[i];
                    if(!this->_secondary && canHit.score() < best_score) continue;
                    // if(min(min_left_anchor, min_right_anchor) <= this->_minK_local) {
                    
                    // daehwan - for debugging purposes
                    // if(this->_anchors_added[i] < this->_anchors_added.back()) continue;
                    
                    //}
                    if(!this->redundant(sink, rdi, canHit)) {
                        this->reportHit(sc, gfm, altdb, ref, ssdb, sink, rdi, canHit);
                        maxsc = max<int64_t>(maxsc, canHit.score());
                    }
                }
            }
            else {
                this->reportHit(sc, gfm, altdb, ref, ssdb, sink, rdi, hit);
                maxsc = max<int64_t>(maxsc, hit.score());
            }
            return maxsc;
        }
    } else if(hitoff > 0 && (hitoff + hitlen == rdlen || hitoff + hitoff < rdlen - hitlen)) {
        // Decide which side to extend first (left or right)
        if(!ssdb.empty()) {
            // extend the partial alignment in the left direction
            index_t fragoff = 0, fraglen = 0, left = 0;
            hit.getLeft(fragoff, fraglen, left);
            const index_t minMatchLen = (index_t)this->_minK_local;
            // make use of a list of known or novel splice sites to further align the read
            if(fraglen >= minMatchLen && left >= minMatchLen && !this->_tpol.no_spliced_alignment()) {
                spliceSites.clear();
                ssdb.getLeftSpliceSites(hit.ref(), left + minMatchLen, minMatchLen + min<index_t>(minMatchLen, fragoff), spliceSites);
                for(size_t si = 0; si < spliceSites.size(); si++) {
                    const SpliceSite& ss = spliceSites[si];
                    if(!ss._fromfile && ss._readid + this->_thread_rids_mindist > rd.rdid) continue;
                    if(left + fraglen - 1 < ss.right()) continue;
                    index_t frag2off = ss.left() -  (ss.right() - left);
                    if(frag2off + 1 < hit.rdoff()) continue;
                    index_t toff = frag2off + 1 - fragoff;
                    index_t joinedOff = 0;
                    bool success = gfm.textOffToJoined(hit.ref(), toff, joinedOff);
                    if(!success) {
                        continue;
                    }
#ifndef NDEBUG
                    index_t debug_tid = 0, debug_toff = 0, debug_tlen = 0;
                    bool debug_straddled = false;
                    gfm.joinedToTextOff(1, // qlen
                                        joinedOff,
                                        debug_tid,
                                        debug_toff,
                                        debug_tlen,
                                        false,
                                        debug_straddled);
                    assert_eq(hit.ref(), debug_tid);
                    assert_eq(toff, debug_toff);
#endif
                    GenomeHit<index_t> tempHit;
                    tempHit.init(hit.fw(),
                                 0, // rdoff
                                 fragoff, // len
                                 0, // trim5
                                 0, // trim3
                                 hit.ref(),
                                 toff,
                                 joinedOff,
                                 this->_sharedVars);
                    if(!tempHit.compatibleWith(hit, (index_t)this->_minIntronLen, (index_t)this->_maxIntronLen, this->_tpol.no_spliced_alignment())) continue;
                    int64_t minsc = this->_minsc[rdi];
                    bool combined = tempHit.combineWith(
                                                        hit,
                                                        rd,
                                                        gfm,
                                                        ref,
                                                        altdb,
                                                        ssdb,
                                                        swa,
                                                        swm,
                                                        sc,
                                                        minsc,
                                                        rnd,
                                                        (index_t)this->_minK_local,
                                                        (index_t)this->_minIntronLen,
                                                        (index_t)this->_maxIntronLen,
                                                        1,
                                                        1,
                                                        &ss);
                    if(!this->_secondary) {
                        if(rdi == 0) minsc = max(minsc, sink.bestUnp1());
                        else         minsc = max(minsc, sink.bestUnp2());
                    }
                    if(combined &&
                       tempHit.score() >= minsc &&
                       // soft-clipping might be better
                       tempHit.score() + sc.sc(0) * hit.rdoff() >= hit.score()) {
                        assert_eq(tempHit.trim5(), 0);
                        assert_leq(tempHit.rdoff() + tempHit.len() + tempHit.trim3(), rdlen);
                        int64_t tmp_maxsc = hybridSearch_recur(
                                                               sc,
                                                               gfm,
                                                               altdb,
                                                               ref,
                                                               swa,
                                                               ssdb,
                                                               rdi,
                                                               tempHit,
                                                               tempHit.rdoff(),
                                                               tempHit.len() + tempHit.trim3(),
                                                               wlm,
                                                               prm,
                                                               swm,
                                                               him,
                                                               rnd,
                                                               sink,
                                                               dep + 1);
                        maxsc = max<int64_t>(maxsc, tmp_maxsc);
                    }
                }
            }
        }
        
        bool use_localindex = true;
        if(hitoff == hit.rdoff() && hitoff <= this->_minK) {
            index_t leftext = (index_t)INDEX_MAX, rightext = (index_t)0;
            GenomeHit<index_t> tempHit = hit;
            tempHit.extend(
                           rd,
                           gfm,
                           ref,
                           altdb,
                           ssdb,
                           swa,
                           swm,
                           prm,
                           sc,
                           this->_minsc[rdi],
                           rnd,
                           (index_t)this->_minK_local,
                           (index_t)this->_minIntronLen,
                           (index_t)this->_maxIntronLen,
                           this->_minAnchorLen,
                           this->_minAnchorLen_noncan,
                           leftext,
                           rightext,
                           1);
            if(tempHit.rdoff() == 0) {
                use_localindex = false;
            }
        }
        
        // Choose a local index based on the genomic location of the partial alignment
        const HGFM<index_t, local_index_t>* hGFM = (const HGFM<index_t, local_index_t>*)(&gfm);
        const LocalGFM<local_index_t, index_t>* lGFM = hGFM->getLocalGFM(hit.ref(), hit.refoff());
        assert_leq(lGFM->_localOffset, hit.refoff());
        bool success = false, first = true;
        index_t count = 0;
        // Use at most two local indexes
        const index_t max_count = 2;
        int64_t prev_score = hit.score();
        local_genomeHits.clear();
        while(!success && count++ < max_count && use_localindex) {
            if(him.localindexatts >= this->max_localindexatts) break;
            if(first) {
                first = false;
            } else {
                lGFM = hGFM->prevLocalGFM(lGFM);
                if(lGFM == NULL || lGFM->empty()) break;
            }
            // local index search
            index_t extlen = 0;
            local_index_t top = (local_index_t)INDEX_MAX, bot = (local_index_t)INDEX_MAX;
            local_index_t node_top = (local_index_t)INDEX_MAX, node_bot = (local_index_t)INDEX_MAX;
            index_t extoff = hitoff - 1;
            if(extoff > 0) extoff -= 1;
            if(extoff < this->_minAnchorLen) {
                extoff = this->_minAnchorLen;
            }
            index_t nelt = (index_t)INDEX_MAX;
            index_t max_nelt = std::max<index_t>(5, extlen);
            bool no_extension = false;
            bool uniqueStop= false;
            index_t minUniqueLen = (index_t)this->_minK_local;
            for(; extoff < rdlen; extoff++) {
                extlen = 0;
                uniqueStop = true;
                him.localindexatts++;
                this->_local_node_iedge_count.clear();
                nelt = this->localGFMSearch(
                                            *lGFM,    // BWT index
                                            rd,       // read to align
                                            sc,       // scoring scheme
                                            sink.reportingParams(),
                                            hit.fw(),
                                            extoff,
                                            extlen,
                                            top,
                                            bot,
                                            node_top,
                                            node_bot,
                                            this->_local_node_iedge_count,
                                            rnd,
                                            uniqueStop,
                                            minUniqueLen);
                if(extoff + 1 - extlen >= hitoff) {
                    no_extension = true;
                    break;
                }
                if(nelt <= max_nelt) break;
            }
            assert_leq(node_top, node_bot);
            assert_eq(nelt, (index_t)(node_bot - node_top));
            assert_leq(extlen, extoff + 1);
            if(nelt > 0 &&
               nelt <= max_nelt &&
               extlen >= this->_minAnchorLen &&
               !no_extension) {
                assert_leq(nelt, max_nelt);
                coords.clear();
                bool straddled = false;
                // get genomic locations for this local search
                this->getGenomeCoords_local(
                                            *lGFM,
                                            altdb,
                                            ref,
                                            rnd,
                                            top,
                                            bot,
                                            node_top,
                                            node_bot,
                                            this->_local_node_iedge_count,
                                            hit.fw(),
                                            extoff + 1 - extlen,
                                            extlen,
                                            coords,
                                            wlm,
                                            prm,
                                            him,
                                            true, // reject straddled?
                                            straddled);
                assert_leq(coords.size(), nelt);
                coords.sort();
                for(int ri = (int)coords.size() - 1; ri >= 0; ri--) {
                    const Coord& coord = coords[ri];
                    GenomeHit<index_t> tempHit;
                    tempHit.init(coord.orient(),
                                 extoff + 1 - extlen,
                                 extlen,
                                 0, // trim5
                                 0, // trim3
                                 (index_t)coord.ref(),
                                 (index_t)coord.off(),
                                 (index_t)coord.joinedOff(),
                                 this->_sharedVars);
                    if(!tempHit.adjustWithALT(*this->_rds[rdi], gfm, altdb, ref)) continue;
                    // check if the partial alignment is compatible with the new alignment using the local index
                    if(!tempHit.compatibleWith(hit, (index_t)this->_minIntronLen, (index_t)this->_maxIntronLen, this->_tpol.no_spliced_alignment())) {
                        if(count == 1) continue;
                        else break;
                    }
                    if(uniqueStop) {
                        assert_eq(coords.size(), 1);
                        index_t leftext = (index_t)INDEX_MAX, rightext = (index_t)0;
                        tempHit.extend(
                                       rd,
                                       gfm,
                                       ref,
                                       altdb,
                                       ssdb,
                                       swa,
                                       swm,
                                       prm,
                                       sc,
                                       this->_minsc[rdi],
                                       rnd,
                                       (index_t)this->_minK_local,
                                       (index_t)this->_minIntronLen,
                                       (index_t)this->_maxIntronLen,
                                       this->_minAnchorLen,
                                       this->_minAnchorLen_noncan,
                                       leftext,
                                       rightext);
                    }
                    // combine the partial alignment and the new alignment
                    int64_t minsc = this->_minsc[rdi];
                    bool combined = tempHit.combineWith(
                                                        hit,
                                                        rd,
                                                        gfm,
                                                        ref,
                                                        altdb,
                                                        ssdb,
                                                        swa,
                                                        swm,
                                                        sc,
                                                        minsc,
                                                        rnd,
                                                        (index_t)this->_minK_local,
                                                        (index_t)this->_minIntronLen,
                                                        (index_t)this->_maxIntronLen,
                                                        this->_minAnchorLen,
                                                        this->_minAnchorLen_noncan);
                    if(!this->_secondary) {
                        if(rdi == 0) minsc = max(minsc, sink.bestUnp1());
                        else         minsc = max(minsc, sink.bestUnp2());
                    }
                    if(combined && tempHit.score() >= minsc) {
                        assert_eq(tempHit.trim5(), 0);
                        assert_leq(tempHit.rdoff() + tempHit.len() + tempHit.trim3(), rdlen);
                        if(tempHit.score() >= prev_score - sc.mmpMax) {
                            // extend the new partial alignment recursively
                            int64_t tmp_maxsc = hybridSearch_recur(
                                                                   sc,
                                                                   gfm,
                                                                   altdb,
                                                                   ref,
                                                                   swa,
                                                                   ssdb,
                                                                   rdi,
                                                                   tempHit,
                                                                   tempHit.rdoff(),
                                                                   tempHit.len() + tempHit.trim3(),
                                                                   wlm,
                                                                   prm,
                                                                   swm,
                                                                   him,
                                                                   rnd,
                                                                   sink,
                                                                   dep + 1);
                            maxsc = max<int64_t>(maxsc, tmp_maxsc);
                        } else {
                            local_genomeHits.push_back(tempHit);
                        }
                    }
                }
            }
            if(maxsc >= prev_score - sc.mmpMax) success = true;
            if(!success &&
               (him.localindexatts >= this->max_localindexatts || count == max_count || hGFM->prevLocalGFM(lGFM) == NULL)) {
                for(index_t ti = 0; ti < local_genomeHits.size(); ti++) {
                    GenomeHit<index_t>& tempHit = local_genomeHits[ti];
                    int64_t minsc = this->_minsc[rdi];
                    if(!this->_secondary) {
                        if(rdi == 0) minsc = max(minsc, sink.bestUnp1());
                        else         minsc = max(minsc, sink.bestUnp2());
                    }
                    if(tempHit.score() >= minsc) {
                        int64_t tmp_maxsc = hybridSearch_recur(
                                                               sc,
                                                               gfm,
                                                               altdb,
                                                               ref,
                                                               swa,
                                                               ssdb,
                                                               rdi,
                                                               tempHit,
                                                               tempHit.rdoff(),
                                                               tempHit.len() + tempHit.trim3(),
                                                               wlm,
                                                               prm,
                                                               swm,
                                                               him,
                                                               rnd,
                                                               sink,
                                                               dep + 1);
                        maxsc = max<int64_t>(maxsc, tmp_maxsc);
                    }
                }
            }
        } // while(!success && count++ < 2)
        
        if(!success) {
            if(hitoff > this->_minK &&
               him.localindexatts < this->max_localindexatts) {
                index_t extlen = 0;
                index_t top = (index_t)INDEX_MAX, bot = (index_t)INDEX_MAX;
                index_t node_top = (index_t)INDEX_MAX, node_bot = (index_t)INDEX_MAX;
                this->_node_iedge_count.clear();
                index_t extoff = hitoff - 1;
                bool uniqueStop = true;
                // perform global search for long introns
                index_t nelt = this->globalGFMSearch(
                                                     gfm,    // GFM index
                                                     rd,     // read to align
                                                     sc,     // scoring scheme
                                                     sink.reportingParams(),
                                                     hit.fw(),
                                                     extoff,
                                                     extlen,
                                                     top,
                                                     bot,
                                                     node_top,
                                                     node_bot,
                                                     this->_node_iedge_count,
                                                     rnd,
                                                     uniqueStop);
                if(nelt > 0 && nelt <= 5 && extlen >= this->_minK) {
                    coords.clear();
                    bool straddled = false;
                    this->getGenomeCoords(
                                          gfm,
                                          altdb,
                                          ref,
                                          rnd,
                                          top,
                                          bot,
                                          node_top,
                                          node_bot,
                                          this->_node_iedge_count,
                                          hit.fw(),
                                          bot - top,
                                          extoff + 1 - extlen,
                                          extlen,
                                          coords,
                                          wlm,
                                          prm,
                                          him,
                                          true, // reject straddled?
                                          straddled);
                    assert_leq(coords.size(), nelt);
                    if(coords.size() > 1) coords.sort();
                    for(int ri = (int)coords.size() - 1; ri >= 0; ri--) {
                        const Coord& coord = coords[ri];
                        GenomeHit<index_t> tempHit;
                        tempHit.init(coord.orient(),
                                     extoff + 1 - extlen,
                                     extlen,
                                     0, // trim5
                                     0, // trim3
                                     (index_t)coord.ref(),
                                     (index_t)coord.off(),
                                     (index_t)coord.joinedOff(),
                                     this->_sharedVars);
                        if(!tempHit.adjustWithALT(*this->_rds[rdi], gfm, altdb, ref)) continue;
                        if(!tempHit.compatibleWith(hit, (index_t)this->_minIntronLen, (index_t)this->_maxIntronLen, this->_tpol.no_spliced_alignment())) continue;
                        if(uniqueStop) {
                            assert_eq(coords.size(), 1);
                            index_t leftext = (index_t)INDEX_MAX, rightext = (index_t)0;
                            tempHit.extend(
                                           rd,
                                           gfm,
                                           ref,
                                           altdb,
                                           ssdb,
                                           swa,
                                           swm,
                                           prm,
                                           sc,
                                           this->_minsc[rdi],
                                           rnd,
                                           (index_t)this->_minK_local,
                                           (index_t)this->_minIntronLen,
                                           (index_t)this->_maxIntronLen,
                                           this->_minAnchorLen,
                                           this->_minAnchorLen_noncan,
                                           leftext,
                                           rightext);
                        }
                        int64_t minsc = this->_minsc[rdi];
                        bool combined = tempHit.combineWith(
                                                            hit,
                                                            rd,
                                                            gfm,
                                                            ref,
                                                            altdb,
                                                            ssdb,
                                                            swa,
                                                            swm,
                                                            sc,
                                                            minsc,
                                                            rnd,
                                                            (index_t)this->_minK_local,
                                                            (index_t)this->_minIntronLen,
                                                            (index_t)this->_maxIntronLen,
                                                            this->_minAnchorLen,
                                                            this->_minAnchorLen_noncan);
                        if(!this->_secondary) {
                            if(rdi == 0) minsc = max(minsc, sink.bestUnp1());
                            else         minsc = max(minsc, sink.bestUnp2());
                        }
                        if(combined && tempHit.score() >= minsc) {
                            assert_eq(tempHit.trim5(), 0);
                            assert_leq(tempHit.rdoff() + tempHit.len() + tempHit.trim3(), rdlen);
                            int64_t tmp_maxsc = hybridSearch_recur(
                                                                   sc,
                                                                   gfm,
                                                                   altdb,
                                                                   ref,
                                                                   swa,
                                                                   ssdb,
                                                                   rdi,
                                                                   tempHit,
                                                                   tempHit.rdoff(),
                                                                   tempHit.len() + tempHit.trim3(),
                                                                   wlm,
                                                                   prm,
                                                                   swm,
                                                                   him,
                                                                   rnd,
                                                                   sink,
                                                                   dep + 1);
                            maxsc = max<int64_t>(maxsc, tmp_maxsc);
                        }
                    }
                }
            }
            GenomeHit<index_t> tempHit = hit;
            index_t trimMax = (index_t)((tempHit.score() - max<int64_t>(maxsc, this->_minsc[rdi])) / sc.sc(0));
            if(tempHit.rdoff() < trimMax) {
                index_t trim5 = tempHit.rdoff();
                GenomeHit<index_t> trimedHit = tempHit;
                trimedHit.trim5(trim5,
                                rd,
                                ssdb,
                                sc,
                                (index_t)this->_minK_local,
                                (index_t)this->_minIntronLen,
                                (index_t)this->_maxIntronLen,
                                this->_minAnchorLen,
                                this->_minAnchorLen_noncan,
                                ref);
                assert_leq(trimedHit.len() + trimedHit.trim5() + trimedHit.trim3(), rdlen);
                int64_t tmp_score = trimedHit.score();
                if(tmp_score > max<int64_t>(maxsc, this->_minsc[rdi])) {
                    int64_t tmp_maxsc = hybridSearch_recur(
                                                           sc,
                                                           gfm,
                                                           altdb,
                                                           ref,
                                                           swa,
                                                           ssdb,
                                                           rdi,
                                                           trimedHit,
                                                           0,
                                                           trimedHit.len() + trimedHit.trim5() + trimedHit.trim3(),
                                                           wlm,
                                                           prm,
                                                           swm,
                                                           him,
                                                           rnd,
                                                           sink,
                                                           dep + 1);
                    maxsc = max<int64_t>(maxsc, tmp_maxsc);
                    // return maxsc;
                }
            }
            // extend the partial alignment directly comparing with the corresponding genomic sequence
            // with mismatches or a gap allowed
            int64_t minsc = this->_minsc[rdi];
            assert_geq(tempHit.score(), minsc);
            index_t mm = (index_t)((tempHit.score() - minsc) / sc.mmpMax);
            index_t leftext = (index_t)INDEX_MAX, rightext = (index_t)0;
            index_t num_mismatch_allowed = 1;
            if(hitoff <= this->_minK_local) {
                num_mismatch_allowed = min<index_t>(tempHit.rdoff(), mm);
            }
            him.localextatts++;
            tempHit.extend(
                           rd,
                           gfm,
                           ref,
                           altdb,
                           ssdb,
                           swa,
                           swm,
                           prm,
                           sc,
                           this->_minsc[rdi],
                           rnd,
                           (index_t)this->_minK_local,
                           (index_t)this->_minIntronLen,
                           (index_t)this->_maxIntronLen,
                           this->_minAnchorLen,
                           this->_minAnchorLen_noncan,
                           leftext,
                           rightext,
                           num_mismatch_allowed);
            if(!this->_secondary) {
                if(rdi == 0) minsc = max(minsc, sink.bestUnp1());
                else         minsc = max(minsc, sink.bestUnp2());
            }
            if(tempHit.score() >= minsc && leftext >= min<index_t>((index_t)this->_minK_local, hit.rdoff())) {
                assert_eq(tempHit.trim5(), 0);
                assert_leq(tempHit.rdoff() + tempHit.len() + tempHit.trim3(), rdlen);
                int64_t tmp_maxsc = hybridSearch_recur(
                                                       sc,
                                                       gfm,
                                                       altdb,
                                                       ref,
                                                       swa,
                                                       ssdb,
                                                       rdi,
                                                       tempHit,
                                                       tempHit.rdoff(),
                                                       tempHit.len() + tempHit.trim3(),
                                                       wlm,
                                                       prm,
                                                       swm,
                                                       him,
                                                       rnd,
                                                       sink,
                                                       dep + 1);
                maxsc = max<int64_t>(maxsc, tmp_maxsc);
            } else if(hitoff > this->_minK_local) {
                // skip some bases of a read
                index_t jumplen = hitoff > this->_minK ? (index_t)this->_minK : (index_t)this->_minK_local;
                assert_leq(hitoff, hit.rdoff());
                int64_t expected_score = hit.score() - (hit.rdoff() - hitoff) / jumplen * sc.mmpMax - sc.mmpMax;
                if(expected_score >= minsc) {
                    assert_lt(hitlen + jumplen, rdlen);
                    assert_eq(hit.trim5(), 0);
                    assert_leq(hitoff + hitlen, rdlen);
                    int64_t tmp_maxsc = hybridSearch_recur(
                                                           sc,
                                                           gfm,
                                                           altdb,
                                                           ref,
                                                           swa,
                                                           ssdb,
                                                           rdi,
                                                           hit,
                                                           hitoff - jumplen,
                                                           hitlen + jumplen,
                                                           wlm,
                                                           prm,
                                                           swm,
                                                           him,
                                                           rnd,
                                                           sink,
                                                           dep + 1);
                    maxsc = max<int64_t>(maxsc, tmp_maxsc);
                }
            }
        }
    } else {
        // extend the partial alignment in the right direction
        assert_lt(hitoff + hitlen, rdlen);
        if(!ssdb.empty()) {
            index_t fragoff = 0, fraglen = 0, right = 0;
            hit.getRight(fragoff, fraglen, right);
            const index_t minMatchLen = (index_t)this->_minK_local;
            // make use of a list of known or novel splice sites to further align the read
            if(fraglen >= minMatchLen && !this->_tpol.no_spliced_alignment()) {
                spliceSites.clear();
                assert_gt(fraglen, 0);
                assert_leq(fragoff + fraglen, rdlen);
                index_t right_unmapped_len = rdlen - fragoff - fraglen;
                ssdb.getRightSpliceSites(hit.ref(), right + fraglen - minMatchLen, minMatchLen + min<index_t>(minMatchLen, right_unmapped_len), spliceSites);
                for(size_t si = 0; si < spliceSites.size(); si++) {
                    const SpliceSite& ss = spliceSites[si];
                    if(!ss._fromfile && ss._readid + this->_thread_rids_mindist > rd.rdid) continue;
                    if(right > ss.left()) continue;
                    index_t frag2off = ss.right() - ss.left() + right + fraglen - 1;
                    GenomeHit<index_t> tempHit;
                    tempHit.init(hit.fw(),
                                 fragoff + fraglen,
                                 rdlen - fragoff - fraglen,
                                 0, // trim5
                                 0, // trim3
                                 hit.ref(),
                                 frag2off,
                                 (index_t)INDEX_MAX,
                                 this->_sharedVars);
                    if(!hit.compatibleWith(tempHit, (index_t)this->_minIntronLen, (index_t)this->_maxIntronLen, this->_tpol.no_spliced_alignment())) continue;
                    GenomeHit<index_t> combinedHit = hit;
                    int64_t minsc = this->_minsc[rdi];
                    bool combined = combinedHit.combineWith(
                                                            tempHit,
                                                            rd,
                                                            gfm,
                                                            ref,
                                                            altdb,
                                                            ssdb,
                                                            swa,
                                                            swm,
                                                            sc,
                                                            minsc,
                                                            rnd,
                                                            (index_t)this->_minK_local,
                                                            (index_t)this->_minIntronLen,
                                                            (index_t)this->_maxIntronLen,
                                                            1,
                                                            1,
                                                            &ss);
                    if(!this->_secondary) {
                        if(rdi == 0) minsc = max(minsc, sink.bestUnp1());
                        else         minsc = max(minsc, sink.bestUnp2());
                    }
                    if(combined && combinedHit.score() >= minsc &&
                       // soft-clipping might be better
                       combinedHit.score() + sc.sc(0) * (rdlen - hit.rdoff() - hit.len() - hit.trim5()) >= hit.score()) {
                        assert_leq(combinedHit.trim5(), combinedHit.rdoff());
                        assert_eq(combinedHit.rdoff() + combinedHit.len(), rdlen);
                        int64_t tmp_maxsc = hybridSearch_recur(
                                                               sc,
                                                               gfm,
                                                               altdb,
                                                               ref,
                                                               swa,
                                                               ssdb,
                                                               rdi,
                                                               combinedHit,
                                                               combinedHit.rdoff() - combinedHit.trim5(),
                                                               combinedHit.len() + combinedHit.trim5(),
                                                               wlm,
                                                               prm,
                                                               swm,
                                                               him,
                                                               rnd,
                                                               sink,
                                                               dep + 1);
                        maxsc = max<int64_t>(maxsc, tmp_maxsc);
                    }
                }
            }
        }
        
        bool use_localindex = true;
        if(hit.len() == hitlen && hitoff + hitlen + this->_minK > rdlen) {
            index_t leftext = (index_t)0, rightext = (index_t)INDEX_MAX;
            GenomeHit<index_t> tempHit = hit;
            tempHit.extend(
                           rd,
                           gfm,
                           ref,
                           altdb,
                           ssdb,
                           swa,
                           swm,
                           prm,
                           sc,
                           this->_minsc[rdi],
                           rnd,
                           (index_t)this->_minK_local,
                           (index_t)this->_minIntronLen,
                           (index_t)this->_maxIntronLen,
                           this->_minAnchorLen,
                           this->_minAnchorLen_noncan,
                           leftext,
                           rightext,
                           1);
            if(tempHit.rdoff() + tempHit.len()== rdlen) {
                use_localindex = false;
            }
        }
        
        // Choose a local index based on the genomic location of the partial alignment
        const HGFM<index_t, local_index_t>* hGFM = (const HGFM<index_t, local_index_t>*)(&gfm);
        const LocalGFM<local_index_t, index_t>* lGFM = hGFM->getLocalGFM(hit.ref(), hit.refoff());
        bool success = false, first = true;
        index_t count = 0;
        // Use at most two local indexes
        const index_t max_count = 2;
        int64_t prev_score = hit.score();
        local_genomeHits.clear();
        while(!success && count++ < max_count && use_localindex) {
            if(him.localindexatts >= this->max_localindexatts) break;
            if(first) {
                first = false;
            } else {
                lGFM = hGFM->nextLocalGFM(lGFM);
                if(lGFM == NULL || lGFM->empty()) break;
            }
            // local index search
            index_t extlen = 0;
            local_index_t top = (local_index_t)INDEX_MAX, bot = (local_index_t)INDEX_MAX;
            local_index_t node_top = (local_index_t)INDEX_MAX, node_bot = (local_index_t)INDEX_MAX;
            index_t extoff = hitoff + hitlen + (index_t)this->_minK_local;
            if(extoff + 1 < rdlen) extoff += 1;
            if(extoff >= rdlen) {
                extoff = rdlen - 1;
            }
            index_t nelt = (index_t)INDEX_MAX;
            index_t max_nelt = std::max<index_t>(5, extlen);
            bool no_extension = false;
            bool uniqueStop;
            index_t minUniqueLen = (index_t)this->_minK_local;
            index_t maxHitLen = max<index_t>(extoff - hitoff - hitlen, (index_t)this->_minK_local);
            for(; maxHitLen < extoff + 1 && extoff < rdlen;) {
                extlen = 0;
                uniqueStop = false;
                him.localindexatts++;
                this->_local_node_iedge_count.clear();
                nelt = this->localGFMSearch(
                                            *lGFM,    // GFM index
                                            rd,       // read to align
                                            sc,       // scoring scheme
                                            sink.reportingParams(),
                                            hit.fw(),
                                            extoff,
                                            extlen,
                                            top,
                                            bot,
                                            node_top,
                                            node_bot,
                                            this->_local_node_iedge_count,
                                            rnd,
                                            uniqueStop,
                                            minUniqueLen,
                                            maxHitLen);
                if(extoff < hitoff + hitlen) {
                    no_extension = true;
                    break;
                }
                if(nelt <= max_nelt) break;
                if(extoff + 1 < rdlen) extoff++;
                else {
                    if(extlen < maxHitLen) break;
                    else maxHitLen++;
                }
            }
            assert_leq(node_top, node_bot);
            assert_eq(nelt, (index_t)(node_bot - node_top));
            assert_leq(extlen, extoff + 1);
            assert_leq(extoff, rdlen);
            if(nelt > 0 &&
               nelt <= max_nelt &&
               extlen >= this->_minAnchorLen &&
               !no_extension) {
                assert_leq(nelt, max_nelt);
                coords.clear();
                bool straddled = false;
                // get genomic locations for this local search
                this->getGenomeCoords_local(
                                            *lGFM,
                                            altdb,
                                            ref,
                                            rnd,
                                            top,
                                            bot,
                                            node_top,
                                            node_bot,
                                            this->_local_node_iedge_count,
                                            hit.fw(),
                                            extoff + 1 - extlen,
                                            extlen,
                                            coords,
                                            wlm,
                                            prm,
                                            him,
                                            true, // reject straddled?
                                            straddled);
                assert_leq(coords.size(), nelt);
                if(coords.size() > 1) coords.sort();
                for(index_t ri = 0; ri < coords.size(); ri++) {
                    const Coord& coord = coords[ri];
                    GenomeHit<index_t> tempHit;
                    tempHit.init(coord.orient(),
                                 extoff + 1 - extlen,
                                 extlen,
                                 0, // trim5
                                 0, // trim3
                                 (index_t)coord.ref(),
                                 (index_t)coord.off(),
                                 (index_t)coord.joinedOff(),
                                 this->_sharedVars);
                    if(!tempHit.adjustWithALT(*this->_rds[rdi], gfm, altdb, ref)) continue;
                    // check if the partial alignment is compatible with the new alignment using the local index
                    if(!hit.compatibleWith(tempHit, (index_t)this->_minIntronLen, (index_t)this->_maxIntronLen, this->_tpol.no_spliced_alignment())) {
                        if(count == 1) continue;
                        else break;
                    }
                    index_t leftext = (index_t)0, rightext = (index_t)INDEX_MAX;
                    tempHit.extend(
                                   rd,
                                   gfm,
                                   ref,
                                   altdb,
                                   ssdb,
                                   swa,
                                   swm,
                                   prm,
                                   sc,
                                   this->_minsc[rdi],
                                   rnd,
                                   (index_t)this->_minK_local,
                                   (index_t)this->_minIntronLen,
                                   (index_t)this->_maxIntronLen,
                                   this->_minAnchorLen,
                                   this->_minAnchorLen_noncan,
                                   leftext,
                                   rightext);
                    GenomeHit<index_t> combinedHit = hit;
                    int64_t minsc = this->_minsc[rdi];
                    // combine the partial alignment and the new alignment
                    bool combined = combinedHit.combineWith(
                                                            tempHit,
                                                            rd,
                                                            gfm,
                                                            ref,
                                                            altdb,
                                                            ssdb,
                                                            swa,
                                                            swm,
                                                            sc,
                                                            minsc,
                                                            rnd,
                                                            (index_t)this->_minK_local,
                                                            (index_t)this->_minIntronLen,
                                                            (index_t)this->_maxIntronLen,
                                                            this->_minAnchorLen,
                                                            this->_minAnchorLen_noncan);
                    if(!this->_secondary) {
                        if(rdi == 0) minsc = max(minsc, sink.bestUnp1());
                        else         minsc = max(minsc, sink.bestUnp2());
                    }
                    if(combined && combinedHit.score() >= minsc) {
                        assert_leq(combinedHit.trim5(), combinedHit.rdoff());
                        if(combinedHit.score() >= prev_score - sc.mmpMax) {
                            // extend the new partial alignment recursively
                            int64_t tmp_maxsc = hybridSearch_recur(
                                                                   sc,
                                                                   gfm,
                                                                   altdb,
                                                                   ref,
                                                                   swa,
                                                                   ssdb,
                                                                   rdi,
                                                                   combinedHit,
                                                                   combinedHit.rdoff() - combinedHit.trim5(),
                                                                   combinedHit.len() + combinedHit.trim5(),
                                                                   wlm,
                                                                   prm,
                                                                   swm,
                                                                   him,
                                                                   rnd,
                                                                   sink,
                                                                   dep + 1);
                            maxsc = max<int64_t>(maxsc, tmp_maxsc);
                        } else {
                            local_genomeHits.push_back(combinedHit);
                        }
                    }
                }
            }
            // int64_t minsc = (rdi == 0 ? sink.bestUnp1() : sink.bestUnp2());
            if(maxsc >= prev_score - sc.mmpMax) success = true;
            if(!success &&
               (him.localindexatts >= this->max_localindexatts || count == max_count || hGFM->nextLocalGFM(lGFM) == NULL) ) {
                for(index_t ti = 0; ti < local_genomeHits.size(); ti++) {
                    GenomeHit<index_t>& tempHit = local_genomeHits[ti];
                    int64_t minsc = this->_minsc[rdi];
                    if(!this->_secondary) {
                        if(rdi == 0) minsc = max(minsc, sink.bestUnp1());
                        else         minsc = max(minsc, sink.bestUnp2());
                    }
                    if(tempHit.score() >= minsc) {
                        int64_t tmp_maxsc = hybridSearch_recur(
                                                               sc,
                                                               gfm,
                                                               altdb,
                                                               ref,
                                                               swa,
                                                               ssdb,
                                                               rdi,
                                                               tempHit,
                                                               tempHit.rdoff() - tempHit.trim5(),
                                                               tempHit.len() + tempHit.trim5(),
                                                               wlm,
                                                               prm,
                                                               swm,
                                                               him,
                                                               rnd,
                                                               sink,
                                                               dep + 1);
                        maxsc = max<int64_t>(maxsc, tmp_maxsc);
                    }
                }
            }
        } // while(!success && count++ < 2)
        
        if(!success) {
            // perform global search for long introns
            if(hitoff + hitlen + this->_minK + 1 < rdlen &&
               him.localindexatts < this->max_localindexatts) {
                index_t extlen = 0;
                index_t top = (index_t)INDEX_MAX, bot = (index_t)INDEX_MAX;
                index_t node_top = (index_t)INDEX_MAX, node_bot = (index_t)INDEX_MAX;
                this->_node_iedge_count.clear();
                index_t extoff = hitoff + hitlen + (index_t)this->_minK + 1;
                bool uniqueStop = true;
                index_t nelt = this->globalGFMSearch(
                                                     gfm,    // GFM index
                                                     rd,     // read to align
                                                     sc,     // scoring scheme
                                                     sink.reportingParams(),
                                                     hit.fw(),
                                                     extoff,
                                                     extlen,
                                                     top,
                                                     bot,
                                                     node_top,
                                                     node_bot,
                                                     this->_node_iedge_count,
                                                     rnd,
                                                     uniqueStop);
                if(nelt > 0 && nelt <= 5 && extlen >= this->_minK) {
                    coords.clear();
                    bool straddled = false;
                    this->getGenomeCoords(
                                          gfm,
                                          altdb,
                                          ref,
                                          rnd,
                                          top,
                                          bot,
                                          node_top,
                                          node_bot,
                                          this->_node_iedge_count,
                                          hit.fw(),
                                          bot - top,
                                          extoff + 1 - extlen,
                                          extlen,
                                          coords,
                                          wlm,
                                          prm,
                                          him,
                                          true, // reject straddled
                                          straddled);
                    assert_leq(coords.size(), nelt);
                    coords.sort();
                    for(index_t ri = 0; ri < coords.size(); ri++) {
                        const Coord& coord = coords[ri];
                        GenomeHit<index_t> tempHit;
                        tempHit.init(coord.orient(),
                                     extoff + 1 - extlen,
                                     extlen,
                                     0, // trim5
                                     0, // trim3
                                     (index_t)coord.ref(),
                                     (index_t)coord.off(),
                                     (index_t)coord.joinedOff(),
                                     this->_sharedVars);
                        if(!tempHit.adjustWithALT(*this->_rds[rdi], gfm, altdb, ref)) continue;
                        if(!hit.compatibleWith(tempHit, (index_t)this->_minIntronLen, (index_t)this->_maxIntronLen, this->_tpol.no_spliced_alignment())) continue;
                        index_t leftext = (index_t)0, rightext = (index_t)INDEX_MAX;
                        tempHit.extend(
                                       rd,
                                       gfm,
                                       ref,
                                       altdb,
                                       ssdb,
                                       swa,
                                       swm,
                                       prm,
                                       sc,
                                       this->_minsc[rdi],
                                       rnd,
                                       (index_t)this->_minK_local,
                                       (index_t)this->_minIntronLen,
                                       (index_t)this->_maxIntronLen,
                                       this->_minAnchorLen,
                                       this->_minAnchorLen_noncan,
                                       leftext,
                                       rightext);
                        GenomeHit<index_t> combinedHit = hit;
                        int64_t minsc = this->_minsc[rdi];
                        bool combined = combinedHit.combineWith(
                                                                tempHit,
                                                                rd,
                                                                gfm,
                                                                ref,
                                                                altdb,
                                                                ssdb,
                                                                swa,
                                                                swm,
                                                                sc,
                                                                minsc,
                                                                rnd,
                                                                (index_t)this->_minK_local,
                                                                (index_t)this->_minIntronLen,
                                                                (index_t)this->_maxIntronLen,
                                                                this->_minAnchorLen,
                                                                this->_minAnchorLen_noncan);
                        if(!this->_secondary) {
                            if(rdi == 0) minsc = max(minsc, sink.bestUnp1());
                            else         minsc = max(minsc, sink.bestUnp2());
                        }
                        if(combined && combinedHit.score() >= minsc) {
                            assert_leq(combinedHit.trim5(), combinedHit.rdoff());
                            int64_t tmp_maxsc = hybridSearch_recur(
                                                                   sc,
                                                                   gfm,
                                                                   altdb,
                                                                   ref,
                                                                   swa,
                                                                   ssdb,
                                                                   rdi,
                                                                   combinedHit,
                                                                   combinedHit.rdoff() - combinedHit.trim5(),
                                                                   combinedHit.len() + combinedHit.trim5(),
                                                                   wlm,
                                                                   prm,
                                                                   swm,
                                                                   him,
                                                                   rnd,
                                                                   sink,
                                                                   dep + 1);
                            maxsc = max<int64_t>(maxsc, tmp_maxsc);
                        }
                    }
                }
            }
            GenomeHit<index_t> tempHit = hit;
            assert(tempHit.trim5() == 0 || hitoff == 0);
            index_t trimLen = rdlen - hitoff - tempHit.len() - tempHit.trim5();
            index_t trimMax = (index_t)((tempHit.score() - max<int64_t>(maxsc, this->_minsc[rdi])) / sc.sc(0));
            if(trimLen < trimMax) {
                index_t trim3 = rdlen - hitoff - tempHit.len() - tempHit.trim5();
                GenomeHit<index_t> trimedHit = tempHit;
                trimedHit.trim3(trim3,
                                rd,
                                ssdb,
                                sc,
                                (index_t)this->_minK_local,
                                (index_t)this->_minIntronLen,
                                (index_t)this->_maxIntronLen,
                                this->_minAnchorLen,
                                this->_minAnchorLen_noncan,
                                ref);
                assert_leq(trimedHit.trim5(), trimedHit.rdoff());
                assert_leq(trimedHit.len() + trimedHit.trim5() + trimedHit.trim3(), rdlen);
                int64_t tmp_score = trimedHit.score();
                if(tmp_score > max<int64_t>(maxsc, this->_minsc[rdi])) {
                    int64_t tmp_maxsc = hybridSearch_recur(
                                                           sc,
                                                           gfm,
                                                           altdb,
                                                           ref,
                                                           swa,
                                                           ssdb,
                                                           rdi,
                                                           trimedHit,
                                                           trimedHit.rdoff() - trimedHit.trim5(),
                                                           trimedHit.len() + trimedHit.trim5() + trimedHit.trim3(),
                                                           wlm,
                                                           prm,
                                                           swm,
                                                           him,
                                                           rnd,
                                                           sink,
                                                           dep + 1);
                    maxsc = max<int64_t>(maxsc, tmp_maxsc);
                    // return maxsc;
                }
            }
            // extend the partial alignment directly comparing with the corresponding genomic sequence
            // with mismatches or a gap allowed
            int64_t minsc = this->_minsc[rdi];
            assert_geq(tempHit.score(), minsc);
            index_t leftext = (index_t)0, rightext = (index_t)INDEX_MAX;
            index_t mm = (index_t)((tempHit.score() - minsc) / sc.mmpMax);
            index_t num_mismatch_allowed = 1;
            if(rdlen - hitoff - hitlen <= this->_minK_local) {
                num_mismatch_allowed = min<index_t>(rdlen - tempHit.rdoff() - tempHit.len(), mm);
            }
            him.localextatts++;
            tempHit.extend(
                           rd,
                           gfm,
                           ref,
                           altdb,
                           ssdb,
                           swa,
                           swm,
                           prm,
                           sc,
                           this->_minsc[rdi],
                           rnd,
                           (index_t)this->_minK_local,
                           (index_t)this->_minIntronLen,
                           (index_t)this->_maxIntronLen,
                           this->_minAnchorLen,
                           this->_minAnchorLen_noncan,
                           leftext,
                           rightext,
                           num_mismatch_allowed);
            if(!this->_secondary) {
                if(rdi == 0) minsc = max(minsc, sink.bestUnp1());
                else         minsc = max(minsc, sink.bestUnp2());
            }
            if(tempHit.score() >= minsc && rightext >= min<index_t>((index_t)this->_minK_local, rdlen - hit.len() - hit.rdoff())) {
                assert_eq(tempHit.trim3(), 0);
                assert_leq(tempHit.trim5(), tempHit.rdoff());
                int64_t tmp_maxsc = hybridSearch_recur(
                                                       sc,
                                                       gfm,
                                                       altdb,
                                                       ref,
                                                       swa,
                                                       ssdb,
                                                       rdi,
                                                       tempHit,
                                                       tempHit.rdoff() - tempHit.trim5(),
                                                       tempHit.len() + tempHit.trim5(),
                                                       wlm,
                                                       prm,
                                                       swm,
                                                       him,
                                                       rnd,
                                                       sink,
                                                       dep + 1);
                maxsc = max<int64_t>(maxsc, tmp_maxsc);
            } else if(hitoff + hitlen + this->_minK_local < rdlen) {
                // skip some bases of a read
                index_t jumplen = hitoff + hitlen + this->_minK < rdlen ? (index_t)this->_minK : (index_t)this->_minK_local;
                assert_lt(hitoff + hitlen + jumplen, rdlen);
                assert_leq(hit.len(), hitlen);
                int64_t expected_score = hit.score() - (hitlen - hit.len()) / jumplen * sc.mmpMax - sc.mmpMax;
                if(expected_score >= minsc) {
                    assert_eq(hit.trim3(), 0);
                    int64_t tmp_maxsc = hybridSearch_recur(
                                                           sc,
                                                           gfm,
                                                           altdb,
                                                           ref,
                                                           swa,
                                                           ssdb,
                                                           rdi,
                                                           hit,
                                                           hitoff,
                                                           hitlen + jumplen,
                                                           wlm,
                                                           prm,
                                                           swm,
                                                           him,
                                                           rnd,
                                                           sink,
                                                           dep + 1);
                    maxsc = max<int64_t>(maxsc, tmp_maxsc);
                }
            }
        }
    }
    
    return maxsc;
}

#endif /*SPLICED_ALIGNER_H_*/
