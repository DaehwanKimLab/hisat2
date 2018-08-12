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

#ifndef REPEAT_H_
#define REPEAT_H_

#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include "assert_helpers.h"
#include "word_io.h"
#include "mem_ids.h"
#include "ref_coord.h"
#include "alt.h"

using namespace std;

template <typename index_t>
class RepeatCoord {
public:
    bool operator< (const RepeatCoord<index_t>& o) const {
        if(joinedOff != o.joinedOff)
            return joinedOff < o.joinedOff;
        if(fw != o.fw)
            return fw;
        if(alleleID != o.alleleID)
            return alleleID < o.alleleID;
        return false;
    }
        
public:
    RepeatCoord() {};
    RepeatCoord(index_t l_tid, index_t l_toff, index_t l_joinedOff, bool l_fw, index_t l_alleleID) :
        tid(l_tid), toff(l_toff), joinedOff(l_joinedOff), fw(l_fw) {};
    index_t tid;
    index_t toff;
    index_t joinedOff;
    bool    fw;
    index_t alleleID;
};

template <typename index_t>
class RepeatAllele {
public:
    RepeatAllele() {
        reset();
    }
    
    void init(index_t                             alleleID_,
              index_t                             allelePos_,
              index_t                             alleleLen_,
              const EList<index_t>&               snpIDs_) {
        alleleID = alleleID_;
        allelePos = allelePos_;
        alleleLen = alleleLen_;
        snpIDs = snpIDs_;
    }
    
    void reset() {
        alleleID = 0;
        allelePos = 0;
        alleleLen = 0;
        snpIDs.clear();
    }
    
    bool operator< (const RepeatAllele& o) const {
        if(alleleID != o.alleleID)
            return alleleID < o.alleleID;
        if(snpIDs.size() != o.snpIDs.size())
            return snpIDs.size() < o.snpIDs.size();
        for(size_t i = 0; i < snpIDs.size(); i++) {
            if(snpIDs[i] != o.snpIDs[i])
                return snpIDs[i] < o.snpIDs[i];
        }
        assert(false);
        return false;
    }
    
#ifndef NDEBUG
    bool repOk() const {
        return true;
    }
#endif
    
    bool write(ofstream& f_out, bool bigEndian) const {
        writeIndex<index_t>(f_out, alleleID, bigEndian);
        writeIndex<index_t>(f_out, allelePos, bigEndian);
        writeIndex<index_t>(f_out, alleleLen, bigEndian);
        writeIndex<index_t>(f_out, snpIDs.size(), bigEndian);
        for(index_t i = 0; i < snpIDs.size(); i++) {
            writeIndex<index_t>(f_out, snpIDs[i], bigEndian);
        }
        return true;
    }
    
    bool read(ifstream& f_in, bool bigEndian) {
        alleleID = readIndex<index_t>(f_in, bigEndian);
        allelePos = readIndex<index_t>(f_in, bigEndian);
        alleleLen = readIndex<index_t>(f_in, bigEndian);
        index_t numSNPs = readIndex<index_t>(f_in, bigEndian);
        snpIDs.resizeExact(numSNPs);
        for(index_t i = 0; i < numSNPs; i++) {
            snpIDs[i] = readIndex<index_t>(f_in, bigEndian);
        }
        return true;
    }
    
    bool compatible(index_t left,
                    index_t right,
                    const EList<index_t>& cmp_snpIDs,
                    pair<index_t, index_t> alt_range) const {
        if(left < allelePos || right > allelePos + alleleLen)
            return false;
        
        if(snpIDs.size() < cmp_snpIDs.size())
            return false;
        
        index_t i = 0;
        for(; i < snpIDs.size(); i++) {
            index_t snpID = snpIDs[i];
            if(snpID >= alt_range.first)
                break;
        }
        if(snpIDs.size() - i < cmp_snpIDs.size())
            return false;
        
        for(index_t j = 0; j < cmp_snpIDs.size(); j++) {
            assert_lt(i + j, snpIDs.size());
            if(snpIDs[i+j] != cmp_snpIDs[j])
                return false;
        }
        
        i += cmp_snpIDs.size();
        if(i < snpIDs.size() && snpIDs[i] < alt_range.second)
            return false;
        
        return true;
    }
    
public:
    index_t                      alleleID;
    index_t                      allelePos;
    index_t                      alleleLen;
    EList<index_t>               snpIDs;
};

// sorting functions
template <typename index_t>
struct sort_pair_loci {
    bool operator()(const pair<RepeatCoord<index_t>, index_t>& a, const pair<RepeatCoord<index_t>, index_t>& b) {
        return a.first.joinedOff < b.first.joinedOff;
    }
};

template <typename index_t>
struct sort_pair_loci_by_index {
    bool operator()(const pair<RepeatCoord<index_t>, index_t>& a, const pair<RepeatCoord<index_t>, index_t>& b) {
        return a.second < b.second;
    }
};

template <typename index_t>
class Repeat {
public:
    void init(const string&                       repName_,
              index_t                             repID_,
              index_t                             repPos_,
              index_t                             repLen_) {
        repName = repName_;
        repID = repID_;
        repPos = repPos_;
        repLen = repLen_;
    }
    
    bool write(ofstream& f_out, bool bigEndian) const {
        writeIndex<index_t>(f_out, repID, bigEndian);
        writeIndex<index_t>(f_out, repPos, bigEndian);
        writeIndex<index_t>(f_out, repLen, bigEndian);
        writeIndex<index_t>(f_out, alleles.size(), bigEndian);
        for(index_t i = 0; i < alleles.size(); i++) {
            alleles[i].write(f_out, bigEndian);
        }
        writeIndex<index_t>(f_out, positions.size(), bigEndian);
        for(index_t i = 0; i < positions.size(); i++) {
            writeIndex<index_t>(f_out, positions[i].joinedOff, bigEndian);
            writeU8(f_out, positions[i].fw);
            writeIndex<index_t>(f_out, positions[i].alleleID, bigEndian);
        }
        return true;
    }
    
    bool read(ifstream& f_in, bool bigEndian) {
        repID = readIndex<index_t>(f_in, bigEndian);
        repPos = readIndex<index_t>(f_in, bigEndian);
        repLen = readIndex<index_t>(f_in, bigEndian);
        index_t numAlleles = readIndex<index_t>(f_in, bigEndian);
        alleles.resizeExact(numAlleles);
        for(index_t i = 0; i < numAlleles; i++) {
            alleles[i].read(f_in, bigEndian);
        }
        index_t numPositions = readIndex<index_t>(f_in, bigEndian);
        positions.resizeExact(numPositions);
        for(index_t i = 0; i < numPositions; i++) {
            positions[i].tid = 0;
            positions[i].toff = 0;
            positions[i].joinedOff = readIndex<index_t>(f_in, bigEndian);
            positions[i].fw = readU8(f_in);
            positions[i].alleleID = readIndex<index_t>(f_in, bigEndian);
            assert_lt(positions[i].alleleID, alleles.size());
        }
        return true;
    }
    
public:
    string                         repName;
    index_t                        repID;
    index_t                        repPos;
    index_t                        repLen;
    EList<RepeatAllele<index_t> >  alleles;
    EList<RepeatCoord<index_t> >   positions;
};

template <typename index_t>
class RepeatDB {
public:
    RepeatDB() {}
    
    virtual ~RepeatDB() {}
    
    bool empty() const { return _repeats.size() == 0; }
    
    EList<Repeat<index_t> >&       repeats()       { return _repeats; }
    const EList<Repeat<index_t> >& repeats() const { return _repeats; }
    
    void write(ofstream& f_out, bool bigEndian) const {
        writeIndex<index_t>(f_out, (index_t)_repeats.size(), bigEndian);
        if(_repeats.size() > 0) {
            for(index_t i = 0; i < _repeats.size(); i++) {
                _repeats[i].write(f_out, bigEndian);
            }
        }
    }
    
    void read(ifstream& f_in, bool bigEndian) {
        index_t numRepeats = readIndex<index_t>(f_in, bigEndian);
        if(numRepeats > 0) {
            _repeats.resizeExact(numRepeats);
            for(index_t i = 0; i < _repeats.size(); i++) {
                _repeats[i].read(f_in, bigEndian);
            }
        }
    }
    
    // Build an internal table to enable rapid search of repeats
    //  and converts joined offsets to chromosome IDs (tid) and loci (toff)
    void construct(const index_t* rstarts, index_t rlen) {
        _repeatMap.clear();
        for(index_t r = 0; r < _repeats.size(); r++) {
            _repeatMap.expand();
            if(_repeatMap.size() == 1) {
                _repeatMap.back().first = _repeats[r].repLen;
            } else {
                _repeatMap.back().first = _repeatMap[_repeatMap.size() - 2].first + _repeats[r].repLen;
            }
            _repeatMap.back().second = r;
        }
        
        EList<pair<RepeatCoord<index_t>, index_t> > joinedOffList;
        for(index_t r = 0; r < _repeats.size(); r++) {
            Repeat<index_t>& repeat = _repeats[r];
            EList<RepeatCoord<index_t> >& positions = repeat.positions;
            for(index_t p = 0; p < positions.size(); p++) {
                RepeatAllele<index_t>& allele = repeat.alleles[positions[p].alleleID];
                if(positions[p].fw) {
                    positions[p].joinedOff -= allele.allelePos;
                } else {
                    assert_leq(allele.allelePos + allele.alleleLen, repeat.repLen);
                    positions[p].joinedOff -= (repeat.repLen - allele.allelePos - allele.alleleLen);
                }
                joinedOffList.expand();
                joinedOffList.back().first.joinedOff = positions[p].joinedOff;
                joinedOffList.back().first.tid = 0;
                joinedOffList.back().first.toff = 0;
                joinedOffList.back().first.fw = positions[p].fw;
                joinedOffList.back().first.alleleID = positions[p].alleleID;
                joinedOffList.back().second = joinedOffList.size() - 1;
            }
        }
        
        sort(joinedOffList.begin(), joinedOffList.end(), sort_pair_loci<index_t>());
        
        index_t j = 0, r = 0;
        while(j < joinedOffList.size() && r < rlen) {
            index_t off = joinedOffList[j].first.joinedOff;
            index_t lower = rstarts[r*3];
            index_t upper;
            if(r == rlen - 1) {
                upper = numeric_limits<index_t>::max();
            } else {
                upper = rstarts[(r+1)*3];
            }
            assert_gt(upper, lower);
            if(off > upper) {
                r++;
                continue;
            }
            assert_geq(off, lower);
            joinedOffList[j].first.tid = rstarts[(r*3)+1];
            joinedOffList[j].first.toff = off - lower + rstarts[(r*3)+2];
            j++;
        }
        
        sort(joinedOffList.begin(), joinedOffList.end(), sort_pair_loci_by_index<index_t>());
        
        index_t count = 0;
        for(index_t r = 0; r < _repeats.size(); r++) {
            EList<RepeatCoord<index_t> >& positions = _repeats[r].positions;
            for(index_t p = 0; p < positions.size(); p++) {
                assert_lt(count, joinedOffList.size());
                assert_eq(positions[p].joinedOff, joinedOffList[count].first.joinedOff);
                positions[p] = joinedOffList[count].first;
                count++;
            }
        }
    }
    
    bool findCoords(index_t               left,  // left offset in the repeat sequence
                    index_t               right, // right offset
                    const EList<index_t>& snpIDs, // SNP IDs
                    const ALTDB<index_t>& altdb,
                    EList<pair<RepeatCoord<index_t>, RepeatCoord<index_t> > >& near_positions,
                    index_t dist = 1000) const {
        near_positions.clear();
        
        // Find a repeat corresponding to a given location (left, right)
        pair<index_t, index_t> repeat(left, numeric_limits<index_t>::max());
        index_t repeatIdx = _repeatMap.bsearchLoBound(repeat);
        assert_lt(repeatIdx, _repeats.size());
        if(right > _repeatMap[repeatIdx].first)
            return false;
        
        const EList<RepeatAllele<index_t> >& alleles = _repeats[repeatIdx].alleles;
        index_t adjLeft = left, adjRight = right;
        if(repeatIdx > 0) {
            adjLeft -= _repeatMap[repeatIdx-1].first;
            adjRight -= _repeatMap[repeatIdx-1].first;
        }
        pair<index_t, index_t> alt_range = get_alt_range(altdb, left, right);
        
        const EList<RepeatCoord<index_t> >& positions = _repeats[repeatIdx].positions;
        for(index_t p = 0; p < positions.size(); p++) {
            const RepeatCoord<index_t>& position = positions[p];
            assert_lt(position.alleleID, alleles.size());
            const RepeatAllele<index_t>& allele = alleles[position.alleleID];
            if(!allele.compatible(adjLeft, adjRight, snpIDs, alt_range))
                continue;
            
            near_positions.expand();
            near_positions.back().first = position;
            if(positions[p].fw) {
                near_positions.back().first.joinedOff += adjLeft;
                near_positions.back().first.toff += adjLeft;
            } else {
                const index_t len = right - left;
                assert_leq(adjLeft + len, _repeats[repeatIdx].repLen);
                index_t rc_adjLeft = _repeats[repeatIdx].repLen - adjLeft - len;
                near_positions.back().first.joinedOff += rc_adjLeft;
                near_positions.back().first.toff += rc_adjLeft;
            }
        }
        
        return near_positions.size() > 0;
    }
    
    bool findCommonCoords(index_t               left,    // left offset in the repeat sequence
                          index_t               right,   // right offset
                          const EList<index_t>& snpIDs,  // SNP IDs
                          index_t               left2,   // left offset 2 in the repeat sequence
                          index_t               right2,  // right offset 2
                          const EList<index_t>& snpIDs2, // SNP IDs
                          const ALTDB<index_t>& altdb,
                          EList<pair<RepeatCoord<index_t>, RepeatCoord<index_t> > >& common_positions,
                          index_t dist = 1000) const {
        common_positions.clear();
        
        // Find a repeat corresponding to a given location (left, right)
        assert_lt(left, right);
        pair<index_t, index_t> repeat(left, numeric_limits<index_t>::max());
        index_t repeatIdx = _repeatMap.bsearchLoBound(repeat);
        assert_lt(repeatIdx, _repeats.size());
        if(right > _repeatMap[repeatIdx].first)
            return false;
        const EList<RepeatAllele<index_t> >& alleles = _repeats[repeatIdx].alleles;
        index_t adjLeft = left, adjRight = right;
        if(repeatIdx > 0) {
            adjLeft -= _repeatMap[repeatIdx-1].first;
            adjRight -= _repeatMap[repeatIdx-1].first;
        }
        pair<index_t, index_t> alt_range = get_alt_range(altdb, left, right);
        
        // Find a repeat cooresponding to a given location (left2, right2)
        assert_lt(left2, right2);
        pair<index_t, index_t> repeat2(left2, numeric_limits<index_t>::max());
        index_t repeatIdx2 = _repeatMap.bsearchLoBound(repeat2);
        assert_lt(repeatIdx2, _repeats.size());
        if(right2 > _repeatMap[repeatIdx2].first)
            return false;
        const EList<RepeatAllele<index_t> >& alleles2 = _repeats[repeatIdx2].alleles;
        index_t adjLeft2 = left2, adjRight2 = right2;
        if(repeatIdx2 > 0) {
            adjLeft2 -= _repeatMap[repeatIdx2-1].first;
            adjRight2 -= _repeatMap[repeatIdx2-1].first;
        }
        pair<index_t, index_t> alt_range2 = get_alt_range(altdb, left2, right2);
        
        const EList<RepeatCoord<index_t> >& positions = _repeats[repeatIdx].positions;
        const EList<RepeatCoord<index_t> >& positions2 = _repeats[repeatIdx2].positions;
        index_t jsave = 0;
        for(index_t i = 0; i < positions.size(); i++) {
            const RepeatAllele<index_t>& allele = alleles[positions[i].alleleID];
            if(!allele.compatible(adjLeft, adjRight, snpIDs, alt_range))
                continue;
            index_t i_pos = positions[i].joinedOff;
            for(index_t j = jsave; j < positions2.size(); j++) {
                index_t j_pos = positions2[j].joinedOff;
                if(j_pos + dist < i_pos) {
                    jsave = j + 1;
                    continue;
                }
                if(i_pos + dist < j_pos)
                    break;
                
                const RepeatAllele<index_t>& allele2 = alleles2[positions2[j].alleleID];
                if(!allele2.compatible(adjLeft2, adjRight2, snpIDs2, alt_range2))
                    continue;
            
                // DK - debugging purposes
                if(positions[i].toff >= 14074044 - 500 && positions[i].toff <= 14074044 + 500 &&
                   positions2[j].toff >= 14073894 - 500 && positions2[j].toff <= 14073894 + 500) {
                    int dk = 0;
                    dk += 1;
                }
            
                common_positions.expand();
                common_positions.back().first = positions[i];
                if(positions[i].fw) {
                    common_positions.back().first.joinedOff += adjLeft;
                    common_positions.back().first.toff += adjLeft;
                } else {
                    const index_t len = right - left;
                    assert_leq(adjLeft + len, _repeats[repeatIdx].repLen);
                    index_t rc_adjLeft = _repeats[repeatIdx].repLen - adjLeft - len;
                    common_positions.back().first.joinedOff += rc_adjLeft;
                    common_positions.back().first.toff += rc_adjLeft;
                }
                common_positions.back().second = positions2[j];
                if(positions2[j].fw) {
                    common_positions.back().second.toff += adjLeft2;
                    common_positions.back().second.joinedOff += adjLeft2;
                } else {
                    const index_t len = right2 - left2;
                    assert_leq(adjLeft2 + len, _repeats[repeatIdx2].repLen);
                    index_t rc_adjLeft2 = _repeats[repeatIdx2].repLen - adjLeft2 - len;
                    common_positions.back().second.joinedOff += rc_adjLeft2;
                    common_positions.back().second.toff += rc_adjLeft2;
                }
            }
        }
        
        return common_positions.size() > 0;
    }
    
private:
    pair<index_t, index_t> get_alt_range(const ALTDB<index_t>& altdb,
                                         index_t left,
                                         index_t right) const {
        pair<index_t, index_t> alt_range;
        ALT<index_t> cmp_alt;
        cmp_alt.pos = left;
        alt_range.first = alt_range.second = (index_t)altdb.alts().bsearchLoBound(cmp_alt);
        for(; alt_range.second < altdb.alts().size(); alt_range.second++) {
            const ALT<index_t>& alt = altdb.alts()[alt_range.second];
            if(alt.left > right) break;
        }
        return alt_range;
    }
    
private:
    EList<Repeat<index_t> >        _repeats;
    EList<pair<index_t, index_t> > _repeatMap; // pos to repeat id
};

#endif /*ifndef REPEAT_H_*/
