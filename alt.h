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

#ifndef ALT_H_
#define ALT_H_

#include <iostream>
#include <fstream>
#include <limits>
#include "assert_helpers.h"
#include "word_io.h"
#include "mem_ids.h"

using namespace std;

enum ALT_TYPE {
    ALT_NONE = 0,
    ALT_SNP_SGL,     // single nucleotide substitution
    ALT_SNP_INS,     // small insertion wrt reference genome
    ALT_SNP_DEL,     // small deletion wrt reference genome
    ALT_SNP_ALT,     // alternative sequence (to be implemented ...)
    ALT_SPLICESITE,
    ALT_EXON
};

template <typename index_t>
struct ALT {
    ALT() {
        reset();
    }
    
    void reset() {
        type = ALT_NONE;
        pos = len = 0;
        seq = 0;
    }
    
    ALT_TYPE type;
    
    union {
        index_t pos;
        index_t left;
    };
    
    union {
        index_t len;
        index_t right;
    };
    
    union {
        uint64_t seq;  // used to store 32 bp, but it can be used to store a pointer to EList<uint64_t>
        struct {
            union {
                bool fw;
                bool reversed;
            };
            bool excluded;
        };
    };
        
public:
    // in order to support a sequence longer than 32 bp
    
    bool snp() const { return type == ALT_SNP_SGL || type == ALT_SNP_DEL || type == ALT_SNP_INS; }
    bool splicesite() const { return type == ALT_SPLICESITE; }
    bool mismatch() const { return type == ALT_SNP_SGL; }
    bool gap() const { return type == ALT_SNP_DEL || type == ALT_SNP_INS || type == ALT_SPLICESITE; }
    bool deletion() const { return type == ALT_SNP_DEL; }
    bool insertion() const { return type == ALT_SNP_INS; }
    bool exon() const { return type == ALT_EXON; }
    
    bool operator< (const ALT& o) const {
        if(pos != o.pos) return pos < o.pos;
        if(type != o.type) {
            if(type == ALT_NONE || o.type == ALT_NONE) {
                return type == ALT_NONE;
            }
            if(type == ALT_SNP_INS) return true;
            else if(o.type == ALT_SNP_INS) return false;
            return type < o.type;
        }
        if(len != o.len) return len < o.len;
        if(seq != o.seq) return seq < o.seq;
        return false;
    }
    
    bool compatibleWith(const ALT& o) const {
        if(pos == o.pos) return false;
        
        // sort the two SNPs
        const ALT& a = (pos < o.pos ? *this : o);
        const ALT& b = (pos < o.pos ? o : *this);
        
        if(a.snp()) {
            if(a.type == ALT_SNP_DEL || a.type == ALT_SNP_INS) {
                if(b.pos <= a.pos + a.len) {
                    return false;
                }
            }
        } else if(a.splicesite()) {
            if(b.pos <= a.right + 2) {
                return false;
            }
        } else {
            assert(false);
        }

        return true;
    }
    
#ifndef NDEBUG
    bool repOk() const {
        if(type == ALT_SNP_SGL) {
            if(len != 1) {
                assert(false);
                return false;
            }
            
            if(seq > 3) {
                assert(false);
                return false;
            }
        } else if(type == ALT_SNP_DEL) {
            if(len <= 0) {
                assert(false);
                return false;
            }
            if(seq != 0) {
                assert(false);
                return false;
            }
        } else if(type == ALT_SNP_INS) {
            if(len <= 0) {
                assert(false);
                return false;
            }
        } else if(type == ALT_SPLICESITE) {
            assert_lt(left, right);
            assert_leq(fw, 1);
        }else {
            assert(false);
            return false;
        }
        return true;
    }
#endif
    
    bool write(ofstream& f_out, bool bigEndian) const {
        writeIndex<index_t>(f_out, pos, bigEndian);
        writeU32(f_out, type, bigEndian);
        writeIndex<index_t>(f_out, len, bigEndian);
        writeIndex<uint64_t>(f_out, seq, bigEndian);
        return true;
    }
    
    bool read(ifstream& f_in, bool bigEndian) {
        pos = readIndex<index_t>(f_in, bigEndian);
        type = (ALT_TYPE)readU32(f_in, bigEndian);
        assert_neq(type, ALT_SNP_ALT);
        len = readIndex<index_t>(f_in, bigEndian);
        seq = readIndex<uint64_t>(f_in, bigEndian);
        return true;
    }
};


template <typename index_t>
class ALTDB {
public:
    ALTDB() :
    _snp(false),
    _ss(false),
    _exon(false)
    {}
    
    virtual ~ALTDB() {}
    
    bool hasSNPs() const { return _snp; }
    bool hasSpliceSites() const { return _ss; }
    bool hasExons() const { return _exon; }
    
    void setSNPs(bool snp) { _snp = snp; }
    void setSpliceSites(bool ss) { _ss = ss; }
    void setExons(bool exon) { _exon = exon; }
    
    EList<ALT<index_t> >& alts()     { return _alts; }
    EList<string>&        altnames() { return _altnames; }
    
    const EList<ALT<index_t> >& alts() const     { return _alts; }
    const EList<string>&        altnames() const { return _altnames; }

private:
    bool _snp;
    bool _ss;
    bool _exon;
    
    EList<ALT<index_t> > _alts;
    EList<string>        _altnames;
};

template <typename index_t>
struct Haplotype {
    Haplotype() {
        reset();
    }
    
    void reset() {
        left = right = 0;
        alts.clear();
    }
    
    bool operator< (const Haplotype& o) const {
        if(left != o.left) return left < o.left;
        if(right != o.right) return right < o.right;
        return false;
    }
    
    index_t left;
    index_t right;
    EList<index_t, 4> alts;
};

#endif /*ifndef ALT_H_*/
