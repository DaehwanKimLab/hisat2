 /*
 * Copyright 2021, Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>
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

#ifndef __QUANT_H__
#define __QUANT_H__

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include "assert_helpers.h"
using namespace std;

class Transcript;
class Gene;

/**
 *
 */
class Quant
{
public:
    void init(vector<string>& infiles);
    
    void calculate();
    
    void showInfo() const;
    
private:
    size_t addGeneIfNotExist(const string& geneName, uint64_t geneLen);
    size_t addTranscriptIfNotExist(const string& geneName, const string& transcriptName, uint64_t transcriptLen);
    
private:
    map<string, uint64_t> seqLens;
    map<string, size_t> ID2Transcript;
    map<string, size_t> ID2Gene;
    
private:
    vector<Gene> genes;
    vector<Transcript> transcripts;
};

/**
 *
 */
class SeqElement
{
public:
    SeqElement() :
    name("unknown")
    , len(0)
    , count(0)
    {}
    
public:
    string name;
    uint64_t count;
    uint64_t len;
};


/**
 *
 */
class Transcript : public SeqElement
{
public:
};


class QuantCalc;
/**
 *
 */
class Gene : public SeqElement
{
public:
    vector<uint64_t> transcripts;
    
public:
    // compatibility matrix
    map<set<size_t>, size_t> compMat;
    vector<double> a;
    vector<double> n;
};


/**
 *
 */
class QuantCalc
{
public:
    static void calculate();
};


#endif /* __QUANT_H__ */

