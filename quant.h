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

typedef double   qfloat;
typedef uint64_t quint;
typedef int64_t  qint;

class Transcript;
class Gene;

/**
 *
 */
class Quant
{
public:
    void init(vector<string>& infiles, const bool sim = false);
    
    void calculate();
    
    void showInfo() const;
    
private:
    size_t addGeneIfNotExist(const string& geneName, quint geneLen);
    size_t addTranscriptIfNotExist(const string& geneName, const string& transcriptName, quint transcriptLen);
    
private:
    map<string, quint> seqLens;
    map<string, quint> ID2Transcript;
    map<string, quint> ID2Gene;
    
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
    quint  count;
    quint  len;
};


/**
 *
 */
class Transcript : public SeqElement
{
public:
};


class SNP;
/**
 *
 */
class Gene : public SeqElement
{
public:
    vector<quint> transcripts;
    vector<quint> transcriptLens;
    map<quint, quint> toLocalID;
    
public:
    vector<SNP> snps;
    
public:
    bool hasAlts(quint left, quint right) const;
    
public:
    // compatibility matrix
    map<set<quint>, quint> compMat;
    vector<qfloat>         a;
    vector<qfloat>         n;
};


/**
 *
 */
class SNP
{
public:
    enum TYPE
    {
        SINGLE = 0,
        DELETION,
        INSERTION
    };
    
public:
    string name;
    quint  pos;
    TYPE   type;
};


/**
 *
 */
class QuantCalc
{
public:
    static void calculate(const vector<quint>& tranLens,
                          const map<set<quint>, quint>& compMat,
                          vector<qfloat>& a,
                          vector<qfloat>& n);
};


#endif /* __QUANT_H__ */

