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

#include "quant.h"
#include "tokenize.h"
#include <iomanip>
#include <cmath>
#include <algorithm>


void Quant::init(vector<string>& infiles)
{
    // Adapt sequence files to ifstreams
    for(size_t i = 0; i < infiles.size(); i++) {
        ifstream fp(infiles[i], std::ifstream::in);
        if (!fp.is_open()) {
            cerr << "Error: could not open "<< infiles[i].c_str() << endl;
            throw 1;
        }
        
        string line;
        vector<string> fields;
        vector<string> transcriptFields;
        set<size_t> comp;
        while (getline(fp, line)) {
            if (line.empty()) {
                continue;
            }
            
            if (line[0] == '@') {
                
                if (line.length() < 3 || line[1] != 'S' || line[2] != 'Q') {
                    continue;
                }
                
                // @SQ     SN:ENSG00000100095      LN:7149
                fields.clear();
                tokenize(line, "\t", fields);
                
                assert (fields.size() == 3);
                string seqName = fields[1].substr(3);
                string lenStr = fields[2].substr(3);
                uint64_t len = atol(lenStr.c_str());
                seqLens[seqName] = len;
                continue;
            }
            
            fields.clear();
            tokenize(line, "\t", fields);
            
            auto geneName = fields[2];
            auto geneLen = seqLens[geneName];
            size_t gid = addGeneIfNotExist(geneName, geneLen);
            auto& gene = genes[gid];
            
            int32_t TI_i = -1, TO_i = -1, Zs_i = -1;
            for (size_t i = 12; i < fields.size(); ++i) {
                const string& field = fields[i];
                if (field[0] == 'T' && field[1] == 'I') {
                    TI_i = i;
                }
                else if (field[0] == 'T' && field[1] == 'O') {
                    TO_i = i;
                }
                else if (field[0] == 'Z' && field[1] == 's') {
                    Zs_i = i;
                }
            }
            assert (TI_i >= 0);
            transcriptFields.clear();
            tokenize(fields[TI_i].substr(5), "-", transcriptFields);
            assert (transcriptFields.size() == 2);
            string transcriptName = transcriptFields[0];
            
            if (Zs_i >= 0) {
                // DK - debugging purposes
                // transcriptName += "_alt";
            }
            
            auto transcriptLen = seqLens[transcriptName];
            size_t tid = addTranscriptIfNotExist(geneName, transcriptName, transcriptLen);
            transcripts[tid].count++;
            
            comp.clear();
            comp.insert(tid);
            
            if (TO_i >= 0) {
                vector<string> otherTranscripts;
                tokenize(fields[TO_i].substr(5), "|", otherTranscripts);
                for (const auto& transcriptStr : otherTranscripts) {
                    transcriptFields.clear();
                    tokenize(transcriptStr, "-", transcriptFields);
                    assert (transcriptFields.size() == 2);
                    
                    string otranscriptName = transcriptFields[0];
                    auto otranscriptLen = seqLens[otranscriptName];
                    if (Zs_i >= 0) {
                        // DK - debugging purposes
                        // otranscriptName += "_alt";
                    }
                    size_t oid = addTranscriptIfNotExist(geneName, otranscriptName, otranscriptLen);
                    assert (comp.find(oid) == comp.end());
                    comp.insert(oid);
                }
            }
            
            auto itr = gene.compMat.find(comp);
            if (itr == gene.compMat.end()) {
                gene.compMat[comp] = 1;
            } else {
                itr->second++;
            }
        }
    }
    
    calculate();
    
    showInfo();
}

void Quant::calculate()
{
    for (auto& gene : genes) {
        const size_t nt = gene.transcripts.size();
        assert (nt > 0);
        
        auto& a = gene.a;
        auto& n = gene.n;
        auto& compMat = gene.compMat;
        
        a = vector<double>(nt, 1.0 / nt);
        n = vector<double>(nt, 0);
        
        vector<double> an(nt, 0.0);
        while (true) {
            fill(n.begin(), n.end(), 0);
            for (auto compItem = compMat.begin(); compItem != compMat.end(); compItem++) {
                double total = 0;
                for (auto t : compItem->first) {
                    total += a[t];
                }
                
                for (auto t : compItem->first) {
                    auto nr = compItem->second; // number of reads
                    n[t] += (nr * a[t] / total);
                }
            }
            
            double total = 0.0;
            for (size_t t = 0; t < nt; ++t) {
                total += (n[t] / transcripts[t].len);
            }
            
            for (size_t t = 0; t < nt; ++t) {
                an[t] = n[t] / transcripts[t].len / total;
            }
            
            double diff = 0.0;
            for (size_t i = 0; i < a.size(); ++i) {
                diff += abs(an[i] - a[i]);
            }
            
            if (diff <= 0.0001) {
                break;
            }
            
            a = an;
        }
    }
}

void Quant::showInfo() const
{
    for (const auto& gene : genes) {
        const auto& a = gene.a;
        const auto& n = gene.n;
        const auto& compMat = gene.compMat;
        
        cout << endl << endl;
        cout << "Gene: " << gene.name << endl << endl;
        cout << "\t" << left << setw(20) << "Transcript Name" << right << setw(10) << "Length" << setw(10) << "Count" << endl;
        for (size_t t = 0; t < gene.transcripts.size(); ++t) {
            const auto tid = gene.transcripts[t];
            const auto& transcript = transcripts[tid];
            cout << "\t" << left << setw(20) << transcript.name << right << setw(10) << transcript.len << setw(10) << transcript.count << endl;
        }
        
        cout << endl << endl;
        cout << "\tCompatibility matrix:" << endl;
        for (auto compItem = compMat.begin(); compItem != compMat.end(); compItem++) {
            cout << "\t" << setw(10) << right << compItem->second << "\t" << left;
            for (auto id : compItem->first) {
                cout << "\t" << setw(20) << transcripts[id].name;
            }
            cout << endl;
        }
        
        cout << endl << endl;
        cout << "\tEstimated:" << endl;
        cout << "\t\t" << left << setw(20) << "Transcript Name" << right << setw(10) << "Count" << setw(14) << "Proportion" << endl;
        for (size_t t = 0; t < gene.transcripts.size(); ++t) {
            const auto tid = gene.transcripts[t];
            const auto& transcript = transcripts[tid];
            cout << "\t\t" << left << setw(20) << transcript.name << right << setw(10) << n[t] << setw(14) << a[t] << endl;
        }
    }
}

size_t Quant::addGeneIfNotExist(const string& geneName, uint64_t geneLen)
{
    auto itr = ID2Gene.find(geneName);
    if (itr != ID2Gene.end()) {
        return itr->second;
    }
    
    Gene gene;
    gene.name = geneName;
    gene.len = geneLen;
    gene.count = 0;
    
    genes.push_back(gene);
    size_t id = genes.size() - 1;
    ID2Gene[geneName] = id;
    return id;
}

size_t Quant::addTranscriptIfNotExist(const string& geneName, const string& transcriptName, uint64_t transcriptLen)
{
    auto itr = ID2Transcript.find(transcriptName);
    if (itr != ID2Transcript.end()) {
        return itr->second;
    }
    
    assert (ID2Gene.find(geneName) != ID2Gene.end());
    auto gene_id = ID2Gene[geneName];
    auto& gene = genes[gene_id];
    
    Transcript transcript;
    transcript.name = transcriptName;
    transcript.len = transcriptLen;
    transcript.count = 0;
    
    transcripts.push_back(transcript);
    auto id = transcripts.size() - 1;
    ID2Transcript[transcriptName] = id;
    
    gene.transcripts.push_back(id);
    
    return id;
}
