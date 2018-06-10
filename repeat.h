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
#include "assert_helpers.h"
#include "word_io.h"
#include "mem_ids.h"

using namespace std;

template <typename index_t>
class RepeatAllele {
public:
    RepeatAllele() {
        reset();
    }
    
    void init(index_t               repID_,
              index_t               alleleID_,
              const EList<index_t>& snpIDs_,
              const EList<index_t>& positions_) {
        repID = repID_;
        alleleID = alleleID_;
        snpIDs = snpIDs_;
        positions = positions_;
    }
    
    void reset() {
        repID = std::numeric_limits<index_t>::max();
        alleleID = 0;
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
            writeIndex<index_t>(f_out, positions[i], bigEndian);
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
            positions[i] = readIndex<index_t>(f_in, bigEndian);
        }
        return true;
    }
    
private:
    index_t        repID;
    index_t        alleleID;
    EList<index_t> snpIDs;
    EList<index_t> positions;
};

template <typename index_t>
class RepeatDB {
public:
    RepeatDB() {}
    
    virtual ~RepeatDB() {}
    
    bool empty() const { return _repeats.size() == 0; }
    
    EList<RepeatAllele<index_t> >&       repeats()       { return _repeats; }
    const EList<RepeatAllele<index_t> >& repeats() const { return _repeats; }
    
private:
    EList<RepeatAllele<index_t> >       _repeats;
};

#endif /*ifndef REPEAT_H_*/
