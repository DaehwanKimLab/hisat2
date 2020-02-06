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

#ifndef __HISAT2_BIT_PACKED_ARRAY_H
#define __HISAT2_BIT_PACKED_ARRAY_H

#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include "assert_helpers.h"
#include "word_io.h"
#include "mem_ids.h"
#include "ds.h"

using namespace std;

class BitPackedArray {
public:
    BitPackedArray () {}
    ~BitPackedArray();

    /**
     * Return true iff there are no items
     * @return
     */
    inline bool empty() const { return cur_ == 0; }
    inline size_t size() const { return cur_; }

    TIndexOffU get(size_t idx) const;

    inline TIndexOffU operator[](size_t i) const { return get(i); }
    void pushBack(TIndexOffU val);

    void init(size_t max_value);
    void reset();

    void writeFile(const char *filename);
    void writeFile(const string& filename);
    void writeFile(ofstream &fp);

    void readFile(const char *filename);
    void readFile(const string& filename);
    void readFile(ifstream &fp);

    void dump() const;

    size_t getMemUsage() const;

private:
    void init_by_log2(size_t ceil_log2);

    void put(size_t index, TIndexOffU val);
    inline uint64_t bitToMask(size_t bit) const
    {
        return (uint64_t) ((1ULL << bit) - 1);
    }

    TIndexOffU getItem(uint64_t *block, size_t idx, size_t offset) const;
    void setItem(uint64_t *block, size_t idx, size_t offset, TIndexOffU val);

    pair<size_t, size_t> indexToAddress(size_t index) const;
    pair<size_t, size_t> columnToPosition(size_t col) const;


    void expand(size_t count = 1);
    void allocSize(size_t sz);
    void allocItems(size_t count);


private:
    size_t item_bit_size_;      // item bit size(e.g. 33bit)

    size_t elm_bit_size_;       // 64bit
    size_t items_per_block_bit_;
    size_t items_per_block_bit_mask_;
    size_t items_per_block_;    // number of items in block

    size_t cur_;                // current item count
    size_t sz_;                 // maximum item count

    size_t block_size_;         // block size in byte

    // List of packed array
    EList<uint64_t *> blocks_;
};


#endif //__HISAT2_BIT_PACKED_ARRAY_H
