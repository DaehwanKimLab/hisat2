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

#ifndef __REPEAT_KMER_H__
#define __REPEAT_KMER_H__

#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include "assert_helpers.h"
#include "word_io.h"
#include "mem_ids.h"
#include "ref_coord.h"
#include "ref_read.h"
#include "edit.h"
#include "ds.h"
#include "repeat.h"
#include "blockwise_sa.h"


class RB_Minimizer {
public:
    static pair<uint64_t, size_t>
    get_minimizer(const string& seq,
                  size_t off,
                  size_t window,
                  size_t k)
    {
        assert_leq(k, 32);
        assert_leq(off + window + k - 1, seq.length());
        pair<uint64_t, size_t> minimizer(get_kmer(seq, off, k), off);
        uint64_t kmer = minimizer.first;
        for(size_t i = off + 1; i < off + window; i++) {
            uint64_t next_kmer = get_next_kmer(kmer, seq[i+k-1], k);
            if(minimizer_leq(next_kmer, minimizer.first)) {
                minimizer.first = next_kmer;
                minimizer.second = i;
            }
            kmer = next_kmer;
        }
        return minimizer;
    }
    
    static void
    get_minimizer(const string& seq,
                  size_t window,
                  size_t k,
                  EList<pair<uint64_t, size_t> >& minimizers)
    {
        assert_leq(k, 32);
        assert_leq(window + k - 1, seq.length());
        
        minimizers.clear();
        pair<uint64_t, size_t> minimizer = get_minimizer(seq, 0, window, k);
        minimizers.push_back(minimizer);
        uint64_t kmer = get_kmer(seq, window - 1, k);
        for(size_t i = 1; i + window + k - 1 <= seq.length(); i++) {
            uint64_t next_kmer = get_next_kmer(kmer, seq[i+window+k-2], k);
            if(minimizer.second < i) {
               minimizer = get_minimizer(seq, i, window, k);
            } else if(minimizer_leq(next_kmer, minimizer.first)) {
                minimizer.first = next_kmer;
                minimizer.second = i + window - 1;
            }
            minimizers.push_back(minimizer);
            kmer = next_kmer;
        }
        
#ifndef NDEBUG
        assert_eq(minimizers.size() + window + k - 2, seq.length());
        for(size_t i = 0; i + window + k - 1 <= seq.length(); i++) {
            pair<uint64_t, size_t> minimizer = get_minimizer(seq, i, window, k);
            assert(minimizer == minimizers[i]);
        }
#endif
    }
    
protected:
    static bool
    minimizer_leq(uint64_t kmer, uint64_t kmer2)
    {
        return kmer <= kmer2;
    }
    
    static uint64_t
    get_kmer(const string& seq,
             size_t offset,
             size_t k)
    {
        assert_leq(offset + k, seq.length());
        uint64_t kmer = 0;
        for(size_t i = 0; i < k; i++) {
            kmer = (kmer << 2) | asc2dna[seq[offset + i]];
        }
        return kmer;
    }
    
    static uint64_t
    get_next_kmer(uint64_t kmer,
                  char base,
                  size_t k)
    {
        assert(base == 'A' || base == 'C' || base == 'G' || base == 'T');
        kmer &= (((uint64_t)1 << ((k-1)*2))) - 1;
        kmer = (kmer << 2) | asc2dna[base];
        return kmer;
    }
    
    static string get_string(uint64_t kmer, size_t k)
    {
        string seq = "";
        for(size_t i = 0; i < k; i++) {
            size_t nt = kmer & 0x3;
            seq.push_back("ACGT"[nt]);
            kmer >>= 2;
        }
        reverse(seq.begin(), seq.end());
        return seq;
    }
};

class RB_KmerTable {
public:
    RB_KmerTable()
    {
        
    }
    
    ~RB_KmerTable()
    {
        
    }


public:
    template<typename TStr>
    void build_kmer_table(const string& consensus, EList<pair<size_t, size_t> >& kmer_table, size_t k)
    {
#if 0
        kmer_table.clear();
        if(consensus.length() < k)
            return;
        size_t kmer = 0;
        for(size_t i = 0; i + k <= consensus.length(); i++) {
            if(i == 0) {
                kmer = extract_kmer(consensus, i, k);
            } else {
                kmer = next_kmer(kmer, consensus[i+k-1], k);
            }
            kmer_table.expand();
            kmer_table.back().first = kmer;
            kmer_table.back().second = i;
        }
        kmer_table.sort();
#endif
    }

private:
};


#endif /* __REPEAT_KMER_H__ */
