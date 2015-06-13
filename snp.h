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

#ifndef SNP_H_
#define SNP_H_

#include <iostream>
#include <fstream>
#include <limits>
#include "assert_helpers.h"
#include "word_io.h"
#include "mem_ids.h"

using namespace std;

enum SNP_TYPE {
    SNP_SGL = 0, // single nucleotide substitution
    SNP_INS,     // small insertion wrt reference genome
    SNP_DEL,     // small deletion wrt reference genome
    SNP_ALT,     // alternative sequence (to be implemented ...)
};

template <typename index_t>
struct SNP {
    index_t   pos;
    SNP_TYPE  type;
    uint32_t  len;
    uint64_t  seq;  // used to store 32 bp, but it can be used to store a pointer to EList<uint64_t>
    // in order to support a sequence longer than 32 bp
    
    bool operator< (const SNP& o) const {
        if(pos != o.pos) return pos < o.pos;
        if(type != o.type) return type < o.type;
        if(len != o.len) return len < o.len;
        if(seq != o.seq) return seq < o.seq;
        return false;
    }
    
    bool compatibleWith(const SNP& o) const {
        if(pos == o.pos) return false;
        
        // sort the two SNPs
        const SNP& a = (pos < o.pos ? *this : o);
        const SNP& b = (pos < o.pos ? o : *this);
        
        if(a.type == SNP_INS) {
            if(b.pos <= a.pos + a.len) {
                return false;
            }
        }

        return true;
    }
    
#ifndef NDEBUG
    bool repOk() const {
        if(type == SNP_SGL) {
            if(len != 1) {
                assert(false);
                return false;
            }
            
            if(seq > 3) {
                assert(false);
                return false;
            }
        } else if(type == SNP_DEL) {
            if(len <= 0) {
                assert(false);
                return false;
            }
            if(seq != 0) {
                assert(false);
                return false;
            }
        } else if(type == SNP_INS) {
            if(len <= 0) {
                assert(false);
                return false;
            }
        } else {
            assert(false);
            return false;
        }
        return true;
    }
#endif
    
    bool write(ofstream& f_out, bool bigEndian) const {
        writeIndex<index_t>(f_out, pos, bigEndian);
        writeU32(f_out, type, bigEndian);
        writeU32(f_out, len, bigEndian);
        writeIndex<uint64_t>(f_out, seq, bigEndian);
        return true;
    }
    
    bool read(ifstream& f_in, bool bigEndian) {
        pos = readIndex<index_t>(f_in, bigEndian);
        type = (SNP_TYPE)readU32(f_in, bigEndian);
        len = readU32(f_in, bigEndian);
        seq = readIndex<uint64_t>(f_in, bigEndian);
        return true;
    }
};


template <typename index_t>
class SNPDB {
public:
    SNPDB() {
        
    }
    
    virtual ~SNPDB() {
        
    }
    
    EList<SNP<index_t> >& snps()     { return _snps; }
    EList<string>&        snpnames() { return _snpnames; }
    
    const EList<SNP<index_t> >& snps() const     { return _snps; }
    const EList<string>&        snpnames() const { return _snpnames; }

private:
    EList<SNP<index_t> > _snps;
    EList<string>        _snpnames;
};

#endif /*ifndef SNP_H_*/
