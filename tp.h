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

/*
 *  tp.h
 *
 */

#ifndef TP_H_
#define TP_H_

#include <iostream>
#include <stdint.h>

/**
 * Encapsulates alignment policy for transcriptome
 */
class TranscriptomePolicy {
    
public:
    
    TranscriptomePolicy() { reset(); }
    
    TranscriptomePolicy(
                        size_t minIntronLen,
                        size_t maxIntronLen,
                        uint32_t minAnchorLen = 7,
                        uint32_t minAnchorLen_noncan = 14,
                        bool no_spliced_alignment = false,
                        bool transcriptome_mapping_only = false,
                        bool transcriptome_assembly = false,
                        bool xs_only = false,
                        bool avoid_pseudogene = false)
    {
        init(minIntronLen,
             maxIntronLen,
             minAnchorLen,
             minAnchorLen_noncan,
             no_spliced_alignment,
             transcriptome_mapping_only,
             transcriptome_assembly,
             xs_only,
             avoid_pseudogene);
    }
    
    /**
     */
    void reset() {
        init(false, false, false);
    }
    
    /**
     */
    void init(
              size_t minIntronLen,
              size_t maxIntronLen,
              uint32_t minAnchorLen = 7,
              uint32_t minAnchorLen_noncan = 14,
              bool no_spliced_alignment = false,
              bool transcriptome_mapping_only = false,
              bool transcriptome_assembly = false,
              bool xs_only = false,
              bool avoid_pseudogene = false)
    {
        minIntronLen_ = minIntronLen;
        maxIntronLen_ = maxIntronLen;
        minAnchorLen_ = minAnchorLen;
        minAnchorLen_noncan_ = minAnchorLen_noncan;
        no_spliced_alignment_ = no_spliced_alignment;
        transcriptome_mapping_only_ = transcriptome_mapping_only;
        transcriptome_assembly_ = transcriptome_assembly;
        xs_only_ = xs_only;
        avoid_pseudogene_ = avoid_pseudogene;
    }
    
    size_t minIntronLen() const { return minIntronLen_; }
    size_t maxIntronLen() const { return maxIntronLen_; }
    uint32_t minAnchorLen() const { return minAnchorLen_; }
    uint32_t minAnchorLen_noncan() const { return minAnchorLen_noncan_; }
    bool no_spliced_alignment() const { return no_spliced_alignment_; }
    bool transcriptome_mapping_only() const { return transcriptome_mapping_only_; }
    bool transcriptome_assembly() const { return transcriptome_assembly_; }
    bool xs_only() const { return xs_only_; }
    bool avoid_pseudogene() const { return avoid_pseudogene_; }
    
private:
    size_t   minIntronLen_;
    size_t   maxIntronLen_;
    
    // Minimum anchor length required for canonical splice sites
    uint32_t minAnchorLen_;
    // Minimum anchor length required for non-canonical splice sites
    uint32_t minAnchorLen_noncan_;
    
    bool no_spliced_alignment_;
    bool transcriptome_mapping_only_;
    bool transcriptome_assembly_;
    bool xs_only_;
    bool avoid_pseudogene_;
};

#endif /*ndef TP_H_*/
