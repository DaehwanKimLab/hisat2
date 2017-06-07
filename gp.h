/*
 * Copyright 2016, Daehwan Kim <infphilo@gmail.com>
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
 *  gp.h
 *
 */

#ifndef GP_H_
#define GP_H_

#include <iostream>
#include <stdint.h>

/**
 * Encapsulates alignment policy for graph
 */
class GraphPolicy {
    
public:
    
    GraphPolicy() { reset(); }
    
    GraphPolicy(size_t maxAltsTried,
                bool useHaplotype,
                bool haplotypeOnly,
                bool enableCODIS)
    {
        init(maxAltsTried,
             useHaplotype,
             haplotypeOnly,
             enableCODIS);
    }
    
    /**
     */
    void reset() {
        init(0);
    }
    
    /**
     */
    void init(size_t maxAltsTried,
              bool useHaplotype = false,
              bool haplotypeOnly = false,
              bool enableCODIS = false)
    {
        maxAltsTried_ = maxAltsTried;
        useHaplotype_ = useHaplotype;
        haplotypeOnly_ = haplotypeOnly;
        enableCODIS_ = enableCODIS;
    }
    
    size_t maxAltsTried() const { return maxAltsTried_; }
    bool   useHaplotype() const { return useHaplotype_; }
    bool   haplotypeOnly() const { return haplotypeOnly_; }
    bool   enableCODIS() const { return enableCODIS_; }
    
    
private:
    size_t maxAltsTried_;
    bool   useHaplotype_;
    bool   haplotypeOnly_;
    bool   enableCODIS_;
};

#endif /*ndef GP_H_*/
