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

#include <iostream>
#include <vector>
#include <algorithm>
#include "timer.h"
#include "aligner_sw.h"
#include "aligner_result.h"
#include "scoring.h"
#include "sstring.h"

#include "bit_packed_array.h"

TIndexOffU BitPackedArray::get(size_t index) const
{
    assert_lt(index, cur_);

    pair<size_t, size_t> addr = indexToAddress(index);
    uint64_t *block = blocks_[addr.first];
    pair<size_t, size_t> pos = columnToPosition(addr.second);
    TIndexOffU val = getItem(block, pos.first, pos.second);

    return val;
}


#define write_fp(x) fp.write((const char *)&(x), sizeof((x)))

void BitPackedArray::writeFile(ofstream &fp)
{
    size_t sz = 0;

    write_fp(item_bit_size_);
    write_fp(elm_bit_size_);
    write_fp(items_per_block_bit_);
    write_fp(items_per_block_bit_mask_);
    write_fp(items_per_block_);

    write_fp(cur_);
    write_fp(sz_);

    write_fp(block_size_);

    // number of blocks
    sz = blocks_.size();
    write_fp(sz);
    for(size_t i = 0; i < sz; i++) {
        fp.write((const char *)blocks_[i], block_size_);
    }
}

void BitPackedArray::writeFile(const char *filename)
{
    ofstream fp(filename, std::ofstream::binary);
    writeFile(fp);
    fp.close();
}

void BitPackedArray::writeFile(const string &filename)
{
    writeFile(filename.c_str());
}


#define read_fp(x) fp.read((char *)&(x), sizeof((x)))

void BitPackedArray::readFile(ifstream &fp)
{
    size_t val_sz = 0;

    read_fp(val_sz);
    rt_assert_eq(val_sz, item_bit_size_);

    read_fp(val_sz);
    rt_assert_eq(val_sz, elm_bit_size_);

    read_fp(val_sz);
    rt_assert_eq(val_sz, items_per_block_bit_);

    read_fp(val_sz);
    rt_assert_eq(val_sz, items_per_block_bit_mask_);

    read_fp(val_sz);
    rt_assert_eq(val_sz, items_per_block_);

    // skip cur_
    size_t prev_cnt = 0;
    read_fp(prev_cnt);
    cur_ = 0;

    // skip sz_
    size_t prev_sz = 0;
    read_fp(prev_sz);
    sz_ = 0;

    // block_size_
    read_fp(val_sz);
    rt_assert_eq(val_sz, block_size_);

    // alloc blocks
    allocItems(prev_cnt);
    rt_assert_eq(prev_sz, sz_);

    // number of blocks
    read_fp(val_sz);
    rt_assert_eq(val_sz, blocks_.size());
    for(size_t i = 0; i < blocks_.size(); i++) {
        fp.read((char *)blocks_[i], block_size_);
    }
    cur_ = prev_cnt;
}

void BitPackedArray::readFile(const char *filename)
{
    ifstream fp(filename, std::ifstream::binary);
    readFile(fp);
    fp.close();
}

void BitPackedArray::readFile(const string &filename)
{
    readFile(filename.c_str());
}

void BitPackedArray::put(size_t index, TIndexOffU val)
{
    assert_lt(index, cur_);

    pair<size_t, size_t> addr = indexToAddress(index);
    uint64_t *block = blocks_[addr.first];
    pair<size_t, size_t> pos = columnToPosition(addr.second);

    setItem(block, pos.first, pos.second, val);
}

void BitPackedArray::pushBack(TIndexOffU val)
{
    if(cur_ == sz_) {
        allocItems(items_per_block_);
    }

    put(cur_++, val);

    assert_leq(cur_, sz_);
}

TIndexOffU BitPackedArray::getItem(uint64_t *block, size_t idx, size_t offset) const
{
    size_t remains = item_bit_size_;

    TIndexOffU val = 0;

    while(remains > 0) {
        size_t bits = min(elm_bit_size_ - offset, remains);
        uint64_t mask = bitToMask(bits);

        // get value from block
        TIndexOffU t = (block[idx] >> offset) & mask;
        val = val | (t << (item_bit_size_ - remains));

        remains -= bits;
        offset = 0;
        idx++;
    }

    return val;
}

void BitPackedArray::setItem(uint64_t *block, size_t idx, size_t offset, TIndexOffU val)
{
    size_t remains = item_bit_size_;

    while(remains > 0) {
        size_t bits = min(elm_bit_size_ - offset, remains);
        uint64_t mask = bitToMask(bits);
        uint64_t dest_mask = mask << offset;

        // get 'bits' lsb from val
        uint64_t t = val & mask;
        val >>= bits;

        // save 't' to block[idx]
        t <<= offset;
        block[idx] &= ~(dest_mask); // clear
        block[idx] |= t;

        idx++;
        remains -= bits;
        offset = 0;
    }
}

pair<size_t, size_t> BitPackedArray::indexToAddress(size_t index) const
{
    pair<size_t, size_t> addr;

    addr.first = index >> items_per_block_bit_;
    addr.second = index & items_per_block_bit_mask_;

    return addr;
}

pair<size_t, size_t> BitPackedArray::columnToPosition(size_t col) const {
    pair<size_t, size_t> pos;

    pos.first = (col * item_bit_size_) / elm_bit_size_;
    pos.second = (col * item_bit_size_) % elm_bit_size_;
    return pos;
}

void BitPackedArray::expand(size_t count)
{
    if((cur_ + count) > sz_) {
        allocItems(count);
    }

    cur_ += count;

    assert_leq(cur_, sz_);
}

void BitPackedArray::allocSize(size_t sz)
{
    size_t num_block = (sz * sizeof(uint64_t) + block_size_ - 1) / block_size_;

    for(size_t i = 0; i < num_block; i++) {
        uint64_t *ptr = new uint64_t[block_size_];
        blocks_.push_back(ptr);
        sz_ += items_per_block_;
    }
}

void BitPackedArray::allocItems(size_t count)
{
    size_t sz = (count * item_bit_size_ + elm_bit_size_ - 1) / elm_bit_size_;
    allocSize(sz);
}

void BitPackedArray::init(size_t max_value)
{
    item_bit_size_ = ceil(log2(max_value));
    elm_bit_size_ = sizeof(uint64_t) * 8;

    items_per_block_bit_ = 20;  // 1M
    items_per_block_ = 1ULL << (items_per_block_bit_);
    items_per_block_bit_mask_ = items_per_block_ - 1;

    block_size_ = (items_per_block_ * item_bit_size_ + elm_bit_size_ - 1) / elm_bit_size_ * sizeof(uint64_t);

    cur_ = 0;
    sz_ = 0;
}

void BitPackedArray::dump() const
{
    cerr << "item_bit_size_: " << item_bit_size_ << endl;
    cerr << "block_size_: " << block_size_ << endl;
    cerr << "items_per_block_: " << items_per_block_ << endl;
    cerr << "cur_: " << cur_ << endl;
    cerr << "sz_: " << sz_ << endl;
    cerr << "number of blocks: " << blocks_.size() << endl;
}

size_t BitPackedArray::getMemUsage() const
{
    size_t tot = blocks_.size() * block_size_;
    tot += blocks_.totalCapacityBytes();
    return tot;
}

BitPackedArray::~BitPackedArray()
{
    for(size_t i = 0; i < blocks_.size(); i++) {
        uint64_t *ptr = blocks_[i];
        delete [] ptr;
    }
}

