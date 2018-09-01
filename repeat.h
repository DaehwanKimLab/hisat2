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
    
    void init(index_t allelePos_,
              index_t alleleLen_) {
        allelePos = allelePos_;
        alleleLen = alleleLen_;
    }
    
    void reset() {
        allelePos = 0;
        alleleLen = 0;
    }
    
    bool operator< (const RepeatAllele& o) const {
        if(allelePos != o.allelePos)
            return allelePos < o.allelePos;
        return alleleLen < o.alleleLen;
    }
    
#ifndef NDEBUG
    bool repOk() const {
        return true;
    }
#endif
    
    bool write(ofstream& f_out, bool bigEndian) const {
        writeU16(f_out, allelePos, bigEndian);
        writeU16(f_out, alleleLen, bigEndian);
        return true;
    }
    
    bool read(ifstream& f_in, bool bigEndian) {
        allelePos = readU16(f_in, bigEndian);
        alleleLen = readU16(f_in, bigEndian);
        return true;
    }
    
    bool compatible(index_t left, index_t right) const {
        if(left < allelePos || right > allelePos + alleleLen)
            return false;

        return true;
    }
    
public:
    uint16_t allelePos;
    uint16_t alleleLen;
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
    void init(const string& repName_,
              index_t       repID_,
              index_t       repPos_,
              index_t       repLen_) {
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
        if(_repeats.size() <= 0) {
            writeIndex<index_t>(f_out, 0, bigEndian);
            return;
        }
        EList<index_t> repeatGroup;
        for(index_t i = 0; i < _repeats.size(); i++) {
#ifndef NDEBUG
            if(i + 1 < _repeats.size()) {
                assert_leq(_repeats[i].repID, _repeats[i+1].repID);
            }
#endif
            if(_repeats[i].repID > repeatGroup.size()) {
                repeatGroup.push_back(i);
                assert_eq(_repeats[i].repID, repeatGroup.size());
            }
        }
        repeatGroup.push_back(_repeats.size());
        assert_eq(_repeats.back().repID + 1, repeatGroup.size());
        writeIndex<index_t>(f_out, repeatGroup.size(), bigEndian);
        streampos filepos = f_out.tellp();
        EList<streampos> repeatFilePos;
        for(index_t i = 0; i < repeatGroup.size(); i++) {
            writeIndex<uint64_t>(f_out, 0, bigEndian);
        }
        
        for(index_t i = 0; i < repeatGroup.size(); i++) {
            index_t begin = (i == 0 ? 0 : repeatGroup[i-1]);
            index_t end = repeatGroup[i];
            writeIndex<index_t>(f_out, end - begin, bigEndian);
            for(index_t j = begin; j < end; j++) {
                _repeats[j].write(f_out, bigEndian);
            }
            repeatFilePos.push_back(f_out.tellp());
        }
        assert_eq(repeatFilePos.size(), repeatGroup.size());
        
        streampos origpos = f_out.tellp();
        f_out.seekp(filepos);
        for(index_t i = 0; i < repeatFilePos.size(); i++) {
            writeIndex<uint64_t>(f_out, repeatFilePos[i], bigEndian);
        }
        f_out.seekp(origpos);
    }
    
    void read(ifstream& f_in, bool bigEndian, const EList<uint8_t>& includeRepeat) {
        index_t numRepeatGroup = readIndex<index_t>(f_in, bigEndian);
        EList<streampos> filePos; filePos.resizeExact(numRepeatGroup);
        for(index_t i = 0; i < numRepeatGroup; i++) {
            filePos[i] = readIndex<uint64_t>(f_in, bigEndian);
        }
        for(index_t i = 0, repID = 0; i < numRepeatGroup; i++) {
            if(!includeRepeat[i])
                continue;
            if(i > 0) {
                f_in.seekg(filePos[i-1]);
            }
            index_t numRepeats = readIndex<index_t>(f_in, bigEndian);
            index_t repeat_size = _repeats.size();
            _repeats.resizeExact(repeat_size + numRepeats);
            for(index_t j = 0; j < numRepeats; j++) {
                _repeats[repeat_size+j].read(f_in, bigEndian);
                _repeats[repeat_size+j].repID = repID;
            }
            repID++;
        }
        f_in.seekg(filePos.back());
    }
    
    // Build an internal table to enable rapid search of repeats
    //  and converts joined offsets to chromosome IDs (tid) and loci (toff)
    void construct(const index_t* rstarts, index_t rlen) {
        _repeatMap.clear();
        if(_repeats.empty())
            return;
        for(index_t r = 0; r < _repeats.size(); r++) {
            if(_repeats[r].repID >= _repeatMap.size()) {
                _repeatMap.expand();
                _repeatMap.back().clear();
            }
            EList<pair<index_t, index_t> >& repeatMap = _repeatMap.back();
            repeatMap.expand();
            if(repeatMap.size() == 1) {
                repeatMap.back().first = _repeats[r].repLen;
            } else {
                repeatMap.back().first = repeatMap[repeatMap.size() - 2].first + _repeats[r].repLen;
            }
            repeatMap.back().second = r;
        }
        
        EList<pair<RepeatCoord<index_t>, index_t> > joinedOffList;
        for(index_t r = 0; r < _repeats.size(); r++) {
            Repeat<index_t>& repeat = _repeats[r];
            EList<RepeatCoord<index_t> >& positions = repeat.positions;
            for(index_t p = 0; p < positions.size(); p++) {
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
            if(off >= upper) {
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
            Repeat<index_t>& repeat = _repeats[r];
            EList<RepeatCoord<index_t> >& positions = _repeats[r].positions;
            for(index_t p = 0; p < positions.size(); p++) {
                assert_lt(count, joinedOffList.size());
                assert_eq(positions[p].joinedOff, joinedOffList[count].first.joinedOff);
                positions[p] = joinedOffList[count].first;
                
                RepeatAllele<index_t>& allele = repeat.alleles[positions[p].alleleID];
                if(positions[p].fw) {
                    positions[p].joinedOff -= allele.allelePos;
                    positions[p].toff -= allele.allelePos;
                } else {
                    assert_leq(allele.allelePos + allele.alleleLen, repeat.repLen);
                    index_t subLen = repeat.repLen - allele.allelePos - allele.alleleLen;
                    positions[p].joinedOff -= subLen;
                    positions[p].toff -= subLen;
                }
                
                count++;
            }
        }
    }
    
    bool repeatExist(index_t repID, index_t left, index_t right) const {
        if(repID >= _repeatMap.size())
            return false;
        
        // Find a repeat corresponding to a given location (left, right)
        const EList<pair<index_t, index_t> >& repeatMap = _repeatMap[repID];
        pair<index_t, index_t> repeat(left, numeric_limits<index_t>::max());
        index_t repeatIdx = repeatMap.bsearchLoBound(repeat);
        assert_lt(repeatIdx, repeatMap.size());
        if(right > repeatMap[repeatIdx].first)
            return false;
        return true;
    }
    
    bool getCoords(index_t               repID,
                   index_t               left,  // left offset in the repeat sequence
                   index_t               right, // right offset
                   const EList<index_t>& snpIDs, // SNP IDs
                   const ALTDB<index_t>& altdb,
                   EList<pair<RepeatCoord<index_t>, RepeatCoord<index_t> > >& near_positions,
                   index_t max_positions = numeric_limits<index_t>::max()) const {
        near_positions.clear();
        
        if(repID >= _repeatMap.size())
            return false;
        
        // Find a repeat corresponding to a given location (left, right)
        const EList<pair<index_t, index_t> >& repeatMap = _repeatMap[repID];
        pair<index_t, index_t> repeat(left, numeric_limits<index_t>::max());
        index_t repeatIdx = repeatMap.bsearchLoBound(repeat);
        assert_lt(repeatIdx, repeatMap.size());
        if(right > repeatMap[repeatIdx].first)
            return false;
        
        index_t repeatIdx_ = repeatMap[repeatIdx].second;
        assert_lt(repeatIdx_, _repeats.size());
        
        const EList<RepeatAllele<index_t> >& alleles = _repeats[repeatIdx_].alleles;
        index_t adjLeft = left, adjRight = right;
        if(repeatIdx > 0) {
            adjLeft -= repeatMap[repeatIdx-1].first;
            adjRight -= repeatMap[repeatIdx-1].first;
        }
        const EList<RepeatCoord<index_t> >& positions = _repeats[repeatIdx_].positions;
        for(index_t p = 0; p < positions.size(); p++) {
            const RepeatCoord<index_t>& position = positions[p];
            assert_lt(position.alleleID, alleles.size());
            const RepeatAllele<index_t>& allele = alleles[position.alleleID];
            if(!allele.compatible(adjLeft, adjRight))
                continue;
            
            near_positions.expand();
            near_positions.back().first = position;
            if(positions[p].fw) {
                near_positions.back().first.joinedOff += adjLeft;
                near_positions.back().first.toff += adjLeft;
            } else {
                const index_t len = right - left;
                assert_leq(adjLeft + len, _repeats[repeatIdx_].repLen);
                index_t rc_adjLeft = _repeats[repeatIdx_].repLen - adjLeft - len;
                near_positions.back().first.joinedOff += rc_adjLeft;
                near_positions.back().first.toff += rc_adjLeft;
            }
            
            if(near_positions.size() >= max_positions)
                break;
        }
        
        return near_positions.size() > 0;
    }
    
    bool findCoords(index_t               anchor_left,
                    index_t               anchor_right,
                    index_t               repID,
                    index_t               left,  // left offset in the repeat sequence
                    index_t               right, // right offset
                    const EList<index_t>& snpIDs, // SNP IDs
                    const ALTDB<index_t>& altdb,
                    EList<pair<RepeatCoord<index_t>, RepeatCoord<index_t> > >& near_positions,
                    index_t max_positions = numeric_limits<index_t>::max(),
                    index_t dist = 1000) const {
        near_positions.clear();
        
        if(repID >= _repeatMap.size())
            return false;

        // Find a repeat corresponding to a given location (left, right)
        const EList<pair<index_t, index_t> >& repeatMap = _repeatMap[repID];
        pair<index_t, index_t> repeat(left, numeric_limits<index_t>::max());
        index_t repeatIdx = repeatMap.bsearchLoBound(repeat);
        assert_lt(repeatIdx, repeatMap.size());
        if(right > repeatMap[repeatIdx].first)
            return false;
        
        index_t repeatIdx_ = repeatMap[repeatIdx].second;
        assert_lt(repeatIdx_, _repeats.size());
        
        const EList<RepeatAllele<index_t> >& alleles = _repeats[repeatIdx_].alleles;
        index_t adjLeft = left, adjRight = right;
        if(repeatIdx > 0) {
            adjLeft -= repeatMap[repeatIdx-1].first;
            adjRight -= repeatMap[repeatIdx-1].first;
        }
        const EList<RepeatCoord<index_t> >& positions = _repeats[repeatIdx_].positions;
        
        RepeatCoord<index_t> cmp;
        cmp.joinedOff = (anchor_left >= dist ? anchor_left - dist : 0);
        index_t p = positions.bsearchLoBound(cmp);
        for(; p < positions.size(); p++) {
            const RepeatCoord<index_t>& position = positions[p];
            index_t pos = positions[p].joinedOff + adjLeft;
            if(pos + dist < anchor_left)
                continue;
            if(anchor_right + dist < pos)
                break;
            
            assert_lt(position.alleleID, alleles.size());
            const RepeatAllele<index_t>& allele = alleles[position.alleleID];
            if(!allele.compatible(adjLeft, adjRight))
                continue;
            
            near_positions.expand();
            near_positions.back().first = position;
            if(positions[p].fw) {
                near_positions.back().first.joinedOff += adjLeft;
                near_positions.back().first.toff += adjLeft;
            } else {
                const index_t len = right - left;
                assert_leq(adjLeft + len, _repeats[repeatIdx_].repLen);
                index_t rc_adjLeft = _repeats[repeatIdx_].repLen - adjLeft - len;
                near_positions.back().first.joinedOff += rc_adjLeft;
                near_positions.back().first.toff += rc_adjLeft;
            }
            
            if(near_positions.size() >= max_positions)
                break;
        }
        
        return near_positions.size() > 0;
    }
    
    bool findCommonCoords(index_t               repID,
                          index_t               left,    // left offset in the repeat sequence
                          index_t               right,   // right offset
                          const EList<index_t>& snpIDs,  // SNP IDs
                          index_t               repID2,
                          index_t               left2,   // left offset 2 in the repeat sequence
                          index_t               right2,  // right offset 2
                          const EList<index_t>& snpIDs2, // SNP IDs
                          const ALTDB<index_t>& altdb,
                          EList<pair<RepeatCoord<index_t>, RepeatCoord<index_t> > >& common_positions,
                          index_t max_positions = numeric_limits<index_t>::max(),
                          index_t dist = 1000) const {
        common_positions.clear();
        
        if(repID >= _repeatMap.size() || repID2 >= _repeatMap.size())
            return false;
        
        // Find a repeat corresponding to a given location (left, right)
        const EList<pair<index_t, index_t> >& repeatMap = _repeatMap[repID];
        assert_lt(left, right);
        pair<index_t, index_t> repeat(left, numeric_limits<index_t>::max());
        index_t repeatIdx = repeatMap.bsearchLoBound(repeat);
        assert_lt(repeatIdx, repeatMap.size());
        if(right > repeatMap[repeatIdx].first)
            return false;
        index_t repeatIdx_ = repeatMap[repeatIdx].second;
        assert_lt(repeatIdx_, _repeats.size());
        const EList<RepeatAllele<index_t> >& alleles = _repeats[repeatIdx_].alleles;
        index_t adjLeft = left, adjRight = right;
        if(repeatIdx > 0) {
            adjLeft -= repeatMap[repeatIdx-1].first;
            adjRight -= repeatMap[repeatIdx-1].first;
        }
        
        // Find a repeat cooresponding to a given location (left2, right2)
        const EList<pair<index_t, index_t> >& repeatMap2 = _repeatMap[repID2];
        assert_lt(left2, right2);
        pair<index_t, index_t> repeat2(left2, numeric_limits<index_t>::max());
        index_t repeatIdx2 = repeatMap2.bsearchLoBound(repeat2);
        assert_lt(repeatIdx2, repeatMap2.size());
        if(right2 > repeatMap2[repeatIdx2].first)
            return false;
        index_t repeatIdx2_ = repeatMap2[repeatIdx2].second;
        assert_lt(repeatIdx2_, _repeats.size());
        const EList<RepeatAllele<index_t> >& alleles2 = _repeats[repeatIdx2_].alleles;
        index_t adjLeft2 = left2, adjRight2 = right2;
        if(repeatIdx2 > 0) {
            adjLeft2 -= repeatMap2[repeatIdx2-1].first;
            adjRight2 -= repeatMap2[repeatIdx2-1].first;
        }
        
        const EList<RepeatCoord<index_t> >& positions = _repeats[repeatIdx_].positions;
        const EList<RepeatCoord<index_t> >& positions2 = _repeats[repeatIdx2_].positions;
        index_t jsave = 0;
        for(index_t i = 0; i < positions.size(); i++) {
            const RepeatAllele<index_t>& allele = alleles[positions[i].alleleID];
            if(!allele.compatible(adjLeft, adjRight))
                continue;
            index_t i_pos = positions[i].joinedOff + adjLeft;
            for(index_t j = jsave; j < positions2.size(); j++) {
                index_t j_pos = positions2[j].joinedOff + adjLeft2;
                if(j_pos + dist < i_pos) {
                    jsave = j + 1;
                    continue;
                }
                if(i_pos + dist < j_pos)
                    break;
                
                const RepeatAllele<index_t>& allele2 = alleles2[positions2[j].alleleID];
                if(!allele2.compatible(adjLeft2, adjRight2))
                    continue;
            
                common_positions.expand();
                common_positions.back().first = positions[i];
                if(positions[i].fw) {
                    common_positions.back().first.joinedOff += adjLeft;
                    common_positions.back().first.toff += adjLeft;
                } else {
                    const index_t len = right - left;
                    assert_leq(adjLeft + len, _repeats[repeatIdx_].repLen);
                    index_t rc_adjLeft = _repeats[repeatIdx_].repLen - adjLeft - len;
                    common_positions.back().first.joinedOff += rc_adjLeft;
                    common_positions.back().first.toff += rc_adjLeft;
                }
                common_positions.back().second = positions2[j];
                if(positions2[j].fw) {
                    common_positions.back().second.toff += adjLeft2;
                    common_positions.back().second.joinedOff += adjLeft2;
                } else {
                    const index_t len = right2 - left2;
                    assert_leq(adjLeft2 + len, _repeats[repeatIdx2_].repLen);
                    index_t rc_adjLeft2 = _repeats[repeatIdx2_].repLen - adjLeft2 - len;
                    common_positions.back().second.joinedOff += rc_adjLeft2;
                    common_positions.back().second.toff += rc_adjLeft2;
                }
                
                if(common_positions.size() >= max_positions)
                    break;
            }
            
            if(common_positions.size() >= max_positions)
                break;
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
    EList<Repeat<index_t> >         _repeats;
    ELList<pair<index_t, index_t> > _repeatMap; // pos to repeat id
};

#endif /*ifndef REPEAT_H_*/
