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

using namespace std;

template <typename index_t>
class RepeatCoord {
public:
    index_t tid;
    index_t toff;
    index_t joinedOff;
    bool    fw;
};

template <typename index_t>
class RepeatAllele {
public:
    RepeatAllele() {
        reset();
    }
    
    void init(const string&                       repName_,
              index_t                             alleleID_,
              index_t                             repID_,
              index_t                             repPos_,
              index_t                             repLen_,
              const EList<index_t>&               snpIDs_,
              const EList<RepeatCoord<index_t> >& positions_) {
        repName = repName_;
        alleleID = alleleID_;
        repID = repID_;
        repPos = repPos_;
        repLen = repLen_;
        snpIDs = snpIDs_;
        positions = positions_;
    }
    
    void reset() {
        repName = "";
        repID = std::numeric_limits<index_t>::max();
        alleleID = 0;
        repPos = 0;
        repLen = 0;
        snpIDs.clear();
        positions.clear();
    }
    
    bool operator< (const RepeatAllele& o) const {
        if(repID != o.repID)
            return repID < o.repID;
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
        writeIndex<index_t>(f_out, repID, bigEndian);
        writeIndex<index_t>(f_out, alleleID, bigEndian);
        writeIndex<index_t>(f_out, repPos, bigEndian);
        writeIndex<index_t>(f_out, repLen, bigEndian);
        writeIndex<index_t>(f_out, snpIDs.size(), bigEndian);
        for(index_t i = 0; i < snpIDs.size(); i++) {
            writeIndex<index_t>(f_out, snpIDs[i], bigEndian);
        }
        writeIndex<index_t>(f_out, positions.size(), bigEndian);
        for(index_t i = 0; i < positions.size(); i++) {
            writeIndex<index_t>(f_out, positions[i].joinedOff, bigEndian);
            writeU8(f_out, positions[i].fw);
        }
        return true;
    }
    
    bool read(ifstream& f_in, bool bigEndian) {
        repID = readIndex<index_t>(f_in, bigEndian);
        alleleID = readIndex<index_t>(f_in, bigEndian);
        repPos = readIndex<index_t>(f_in, bigEndian);
        repLen = readIndex<index_t>(f_in, bigEndian);
        index_t numSNPs = readIndex<index_t>(f_in, bigEndian);
        snpIDs.resizeExact(numSNPs);
        for(index_t i = 0; i < numSNPs; i++) {
            snpIDs[i] = readIndex<index_t>(f_in, bigEndian);
        }
        index_t numPositions = readIndex<index_t>(f_in, bigEndian);
        positions.resizeExact(numPositions);
        for(index_t i = 0; i < numPositions; i++) {
            positions[i].tid = 0;
            positions[i].toff = 0;
            positions[i].joinedOff = readIndex<index_t>(f_in, bigEndian);
            positions[i].fw = readU8(f_in);
        }
        return true;
    }
    
public:
    string                       repName;
    index_t                      alleleID;
    index_t                      repID;  //
    index_t                      repPos;
    index_t                      repLen;
    EList<index_t>               snpIDs;
    EList<RepeatCoord<index_t> > positions;
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
class RepeatDB {
public:
    RepeatDB() {}
    
    virtual ~RepeatDB() {}
    
    bool empty() const { return _repeatAlleles.size() == 0; }
    
    EList<RepeatAllele<index_t> >&       repeatAlleles()       { return _repeatAlleles; }
    const EList<RepeatAllele<index_t> >& repeatAlleles() const { return _repeatAlleles; }
    
    // Build an internal table to allows rapid search of repeats
    //  and converts joined offsets to chromosome IDs (tid) and loci (toff)
    void construct(const index_t* rstarts, index_t rlen) {
        _repeatMap.clear();
        for(index_t i = 0; i < _repeatAlleles.size(); i++) {
            _repeatMap.expand();
            _repeatMap.back().first = _repeatAlleles[i].repPos + _repeatAlleles[i].repLen;
            _repeatMap.back().second = i;
        }
        
        EList<pair<RepeatCoord<index_t>, index_t> > joinedOffList;
        for(index_t i = 0; i < _repeatAlleles.size(); i++) {
            const EList<RepeatCoord<index_t> >& positions = _repeatAlleles[i].positions;
            for(index_t j = 0; j < positions.size(); j++) {
                joinedOffList.expand();
                joinedOffList.back().first.joinedOff = positions[j].joinedOff;
                joinedOffList.back().first.tid = 0;
                joinedOffList.back().first.toff = 0;
                joinedOffList.back().first.fw = positions[j].fw;
                joinedOffList.back().second = joinedOffList.size() - 1;
            }
        }
        
        sort(&joinedOffList[0], &joinedOffList[0] + joinedOffList.size(), sort_pair_loci<index_t>());
        
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
        
        sort(&joinedOffList[0], &joinedOffList[0] + joinedOffList.size(), sort_pair_loci_by_index<index_t>());
        
        index_t count = 0;
        for(index_t i = 0; i < _repeatAlleles.size(); i++) {
            EList<RepeatCoord<index_t> >& positions = _repeatAlleles[i].positions;
            for(index_t j = 0; j < positions.size(); j++) {
                assert_lt(count, joinedOffList.size());
                assert_eq(positions[j].joinedOff, joinedOffList[count].first.joinedOff);
                positions[j] = joinedOffList[count].first;
                count++;
            }
        }
    }
    
    bool findCommonCoords(index_t pos,  // offset in the repeat sequence
                          index_t pos2, // offset in the repeat sequence
                          EList<pair<RepeatCoord<index_t>, RepeatCoord<index_t> > >& common_positions,
                          index_t dist = 1000) const {
        pair<index_t, index_t> repeat1(pos, 0);
        index_t repeatIdx = _repeatMap.bsearchLoBound(repeat1);
        assert_lt(repeatIdx, _repeatAlleles.size());
        const EList<RepeatCoord<index_t> >& positions = _repeatAlleles[repeatIdx].positions;
        index_t adjustedPos = pos;
        if(repeatIdx > 0) {
            adjustedPos -= _repeatMap[repeatIdx-1].first;
        }
        
        pair<index_t, index_t> repeat2(pos2, 0);
        index_t repeatIdx2 = _repeatMap.bsearchLoBound(repeat2);
        assert_lt(repeatIdx2, _repeatAlleles.size());
        const EList<RepeatCoord<index_t> >& positions2 = _repeatAlleles[repeatIdx2].positions;
        index_t adjustedPos2 = pos2;
        if(repeatIdx2 > 0) {
            adjustedPos2 -= _repeatMap[repeatIdx2-1].first;
        }
        
        index_t i = 0, j = 0;
        while(i < positions.size() && j < positions2.size()) {
            index_t i_pos = positions[i].joinedOff, j_pos = positions2[j].joinedOff;
            if(i_pos <= j_pos) {
                if(i_pos + dist >= j_pos) {
                    common_positions.expand();
                    common_positions.back().first = positions[i];
                    common_positions.back().first.toff += adjustedPos;
                    common_positions.back().first.joinedOff += adjustedPos;
                    common_positions.back().second = positions2[j];
                    common_positions.back().second.toff += adjustedPos2;
                    common_positions.back().second.joinedOff += adjustedPos2;
                }
                i += 1;
            } else {
                if(i_pos <= j_pos + dist) {
                    common_positions.expand();
                    common_positions.back().first = positions[i];
                    common_positions.back().first.toff += adjustedPos;
                    common_positions.back().first.joinedOff += adjustedPos;
                    common_positions.back().second = positions2[j];
                    common_positions.back().second.toff += adjustedPos2;
                    common_positions.back().second.joinedOff += adjustedPos2;
                }
                j += 1;
            }
        }

        return common_positions.size() > 0;
    }

private:
    void getRepeatLoci(EList<pair<index_t, index_t> >& pos_list) {
        for(size_t i = 0; i < _repeatAlleles.size(); i++) {
            const RepeatAllele<index_t>& repeatAllele = _repeatAlleles[i];
            const EList<RepeatCoord<index_t> >& repeatCoords = repeatAllele.positions;
            for(size_t j = 0; j < repeatCoords.size(); j++) {
                pos_list.expand();
                pos_list.back().first = repeatCoords[j].joinedOff;
                pos_list.back().second = repeatAllele.repLen;
            }
        }
        pos_list.sort();
        
#ifndef NDEBUG
        for(size_t i = 0; i + 1 < pos_list.size(); i++) {
            if(pos_list[i].first + pos_list[i].second > pos_list[i+1].first) {
                assert(false);
            }
        }
#endif
    }
    
private:
    EList<RepeatAllele<index_t> >  _repeatAlleles;
    EList<pair<index_t, index_t> > _repeatMap; // pos to repeat id
};

#endif /*ifndef REPEAT_H_*/
