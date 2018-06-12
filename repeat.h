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
        index_t numSNPs = readIndex<index_t>(f_in, bigEndian);
        snpIDs.resizeExact(numSNPs);
        for(index_t i = 0; i < numSNPs; i++) {
            snpIDs[i] = readIndex<index_t>(f_in, bigEndian);
        }
        index_t numPositions = readIndex<index_t>(f_in, bigEndian);
        positions.resizeExact(numPositions);
        for(index_t i = 0; i < numPositions; i++) {
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

template <typename index_t>
class RepeatDB {
public:
    RepeatDB() {}
    
    virtual ~RepeatDB() {}
    
    bool empty() const { return _repeatAlleles.size() == 0; }
    
    EList<RepeatAllele<index_t> >&       repeatAlleles()       { return _repeatAlleles; }
    const EList<RepeatAllele<index_t> >& repeatAlleles() const { return _repeatAlleles; }
    
    void constructMap() {
        _repeatMap.clear();
        for(size_t i = 0; i < _repeatAlleles.size(); i++) {
            _repeatMap[_repeatAlleles[i].repID] = i;
        }
    }
    
    bool isRepeat(index_t tid) const {
        return _repeatMap.find(tid) != _repeatMap.end();
    }
    
    bool findCommonCoords(index_t tid, index_t tid2, EList<index_t>& common_positions, index_t dist = 1000) const {
        //assert_neq(_repeatMap.find(tid), _repeatMap.end());
        //assert_neq(_repeatMap.find(tid2), _repeatMap.end());
        
        index_t repeatIdx = _repeatMap.find(tid)->second;
        index_t repeatIdx2 = _repeatMap.find(tid2)->second;
        
        assert_lt(repeatIdx, _repeatAlleles.size());
        const EList<RepeatCoord<index_t> >& positions = _repeatAlleles[repeatIdx].positions;
        
        assert_lt(repeatIdx2, _repeatAlleles.size());
        const EList<RepeatCoord<index_t> >& positions2 = _repeatAlleles[repeatIdx2].positions;
        
        index_t i = 0, j = 0;
        while(i < positions.size() && j < positions2.size()) {
            index_t i_pos = positions[i].joinedOff, j_pos = positions2[j].joinedOff;
            if(i_pos <= j_pos) {
                if(i_pos + dist >= j_pos) {
                    common_positions.push_back(i_pos);
                }
                i += 1;
            } else {
                if(i_pos <= j_pos + dist) {
                    common_positions.push_back(j_pos);
                }
                j += 1;
            }
        }
        
        return common_positions.size() > 0;
    }
    
    
#if 0
    template<typename TStr>
    index_t mask(TStr* s, EList<RefRecord>& szs) {
        index_t szs2_total = 0;
        EList<RefRecord> szs2;
        
        EList<pair<index_t, index_t> > pos_list;
        getRepeatLoci(pos_list);

        ASSERT_ONLY(index_t szs_total = 0);
        index_t pos_i = 0;
        index_t joinedOff = 0;
        for(size_t szs_i = 0; szs_i < szs.size(); szs_i++) {
            RefRecord sz = szs[szs_i];
            index_t left = joinedOff;
            index_t right = left + sz.len;
            ASSERT_ONLY(szs_total += sz.off);
            ASSERT_ONLY(szs_total += sz.len);
            while(pos_i < pos_list.size()) {
                index_t repeatPos = pos_list[pos_i].first;
                index_t repeatLen = pos_list[pos_i].second;
                assert_geq(repeatPos, left);
                if(repeatPos >= right) break;
                assert_leq(repeatPos + repeatLen, right);
                
                // bool first = (repeatPos == left);
                if(szs2.size() > 0 && !sz.first) {
                    if(repeatPos == left)
                    szs2.back().off += sz.off;
                    szs2_total += sz.off;
                } else {
                    szs2.expand();
                    szs2.back().off = sz.off;
                    szs2.back().first = sz.first;
                    szs2.back().len = repeatPos - left;
                    szs2_total += szs2.back().off;
                    szs2_total += szs2.back().len;
                }
                
                sz.off = repeatLen;
                sz.len = sz.len - repeatLen - szs2.back().len;
                sz.first = false;
                left = repeatPos + repeatLen;
                
                pos_i++;
            }
            
            if(sz.off > 0 || sz.len > 0) {
                if(szs2.size() > 0 && szs2.back().len == 0 && !sz.first) {
                    szs2.back().off += sz.off;
                    szs2.back().len = sz.len;
                } else {
                    szs2.push_back(sz);
                    szs2_total += szs2.back().off;
                    szs2_total += szs2.back().len;
                }
            }
            assert_eq(szs_total, szs2_total);
            
            joinedOff = right;
        }
        
#ifndef NDEBUG
        EList<RefRecord> szs3 = szs2;
        unmask(s, szs3);
        assert_eq(szs.size(), szs3.size());
        for(size_t szs_i = 0; szs_i < szs.size(); szs_i++) {
            assert_eq(szs[szs_i].off, szs3[szs_i].off);
            assert_eq(szs[szs_i].len, szs3[szs_i].len);
            assert_eq(szs[szs_i].first, szs3[szs_i].first);
        }
#endif
        
        szs = szs2;
        return szs2_total;
    }
    
    template<typename TStr>
    index_t unmask(TStr* s, EList<RefRecord>& szs) {
        index_t szs2_total = 0;
        EList<RefRecord> szs2;
        
        EList<pair<index_t, index_t> > pos_list;
        getRepeatLoci(pos_list);
        
        ASSERT_ONLY(index_t szs_total = 0);
        index_t pos_i = 0;
        index_t joinedOff = 0;
        for(size_t szs_i = 0; szs_i < szs.size(); szs_i++) {
            RefRecord sz = szs[szs_i];
            index_t left = joinedOff;
            index_t right = left + sz.len;
            ASSERT_ONLY(szs_total += sz.off);
            ASSERT_ONLY(szs_total += sz.len);
            while(pos_i < pos_list.size()) {
                index_t repeatPos = pos_list[pos_i].first;
                index_t repeatLen = pos_list[pos_i].second;
                assert_geq(repeatPos, left);
                if(repeatPos >= left + sz.off) {
                    assert_geq(repeatPos, right);
                    break;
                }
                assert_leq(repeatPos + repeatLen, left + sz.off);
                
                if(szs2.size() > 0 && repeatPos == left && !sz.first) {
                    szs2.back().len += repeatLen;
                } else {
                    szs2.expand();
                    szs2.back().off = repeatPos - left;
                    szs2.back().len = repeatLen;
                    szs2.back().first = sz.first;
                    szs2_total += szs2.back().off;
                    szs2_total += szs2.back().len;
                }
                
                sz.off = sz.off - repeatLen - szs2.back().off;
                sz.first = false;
                left = repeatPos + repeatLen;
                
                pos_i++;
            }
            
            if(sz.off > 0 || sz.len > 0) {
                if(szs2.size() > 0 && szs2.back().len == 0 && !sz.first) {
                    szs2.back().off += sz.off;
                    szs2.back().len = sz.len;
                } else {
                    szs2.push_back(sz);
                }
                szs2_total += szs2.back().off;
                szs2_total += szs2.back().len;
            }
            assert_eq(szs_total, szs2_total);
            
            joinedOff = right;
        }
        
        szs = szs2;
        return szs2_total;
    }
    
#endif
    
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
    map<index_t, index_t>          _repeatMap; // tid to repeat id
};

#endif /*ifndef REPEAT_H_*/
