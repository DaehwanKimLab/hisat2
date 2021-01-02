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


void Quant::init(vector<string>& infiles, const bool sim)
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
        set<quint> comp;
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
                assert (fields.size() >= 3);
                string seqName = fields[1].substr(3);
                string lenStr = fields[2].substr(3);
                quint len = atol(lenStr.c_str());
                seqLens[seqName] = len;
                
                // DK - temporary parsing routine
                if (fields.size() > 3) {
                    if (fields[3][0] == 'A' && fields[3][1] == 'L') {
                        vector<string> snpstrs;
                        tokenize(fields[3].substr(3), "|", snpstrs);
                        
                        size_t gid = addGeneIfNotExist(seqName, len);
                        auto& gene = genes[gid];
                        
                        vector<string> snpfields;
                        for (const auto& snpstr : snpstrs) {
                            snpfields.clear();
                            tokenize(snpstr, "-", snpfields);
                            assert (snpfields.size() == 4);
                            
                            SNP snp;
                            snp.name = snpfields[0];
                            switch (snpfields[3][0]) {
                                case 'S': snp.type = SNP::SINGLE; break;
                                case 'D': snp.type = SNP::DELETION; break;
                                case 'I': snp.type = SNP::INSERTION; break;
                                default: break;
                            }
                            snp.pos = atol(snpfields[2].c_str()) - 1;
                            gene.snps.push_back(snp);
                        }
                    } else if (fields[3][0] == 'G' && fields[3][1] == 'E') {
                        auto geneName = fields[3].substr(3);
                        addTranscriptIfNotExist(geneName, seqName, len);
                        addTranscriptIfNotExist(geneName, seqName + "_alt", len);
                    }
                }
                
                continue;
            }
            
            fields.clear();
            tokenize(line, "\t", fields);
            
            auto geneName = fields[2];
            auto geneLen = seqLens[geneName];
            size_t gid = addGeneIfNotExist(geneName, geneLen);
            auto& gene = genes[gid];
            
            const auto& left_str = fields[3];
            quint left = atol(left_str.c_str()) - 1;
            const auto& cigar_str = fields[5];
            
            auto right = left;
            size_t cigar_i = 0;
            while (cigar_i < cigar_str.length()) {
                quint cigar_len = 0;
                while (cigar_i < cigar_str.length()) {
                    auto ch = cigar_str[cigar_i];
                    if (ch >= '0' && ch <= '9') {
                        cigar_len = cigar_len * 10 + ch - '0';
                        ++cigar_i;
                    } else {
                        break;
                    }
                }
                auto cigar_op = cigar_str[cigar_i];
                if (cigar_op == 'M' ||
                    cigar_op == 'D' ||
                    cigar_op == 'N') {
                    right += cigar_len;
                } else {
                    assert (cigar_op == 'I');
                }
                
                ++cigar_i;
            }
            
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
            
            bool ref = true, alt = true;
            if (gene.hasAlts(left, right)) {
                if (Zs_i >= 0) {
                    ref = false;
                    alt = true;
                } else {
                    ref = true;
                    alt = false;
                }
            } else {
                assert (Zs_i < 0);
            }
            
            const auto tid = ID2Transcript[transcriptName];
            transcripts[tid].count++;
            
            comp.clear();
            const auto local_tid = gene.toLocalID[tid];
            if (ref) {
                comp.insert(local_tid);
            }
            if (alt) {
                comp.insert(local_tid + 1);
            }
            
            if (TO_i >= 0) {
                vector<string> otherTranscripts;
                tokenize(fields[TO_i].substr(5), "|", otherTranscripts);
                for (const auto& transcriptStr : otherTranscripts) {
                    transcriptFields.clear();
                    tokenize(transcriptStr, "-", transcriptFields);
                    assert (transcriptFields.size() == 2);
                    
                    string otranscriptName = transcriptFields[0];
                    const auto otid = ID2Transcript[otranscriptName];

                    const auto local_otid = gene.toLocalID[otid];
                    assert (comp.find(otid) == comp.end());
                    if (ref) {
                        comp.insert(local_otid);
                    }
                    if (alt) {
                        comp.insert(local_otid + 1);
                    }
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
        QuantCalc::calculate(gene.transcriptLens,
                             gene.compMat,
                             gene.a,
                             gene.n);
    }
}

void Quant::showInfo() const
{
    for (const auto& gene : genes) {
        const auto& a = gene.a;
        const auto& n = gene.n;
        const auto& compMat = gene.compMat;
        
        cout << endl << endl;
        cout << "Gene: " << gene.name << endl;
        if (gene.snps.size() > 0) {
            for (const auto& snp : gene.snps) {
                cout << "\t" << setw(14) << snp.name << ": " << setw(12) << snp.pos << endl;
            }
        }
        cout << endl << endl;
        
        cout << "\t" << left << setw(20) << "Transcript Name" << right << setw(10) << "Length" << setw(10) << "Count" << endl;
        for (size_t t = 0; t < gene.transcripts.size(); ++t) {
            const auto tid = gene.transcripts[t];
            const auto& transcript = transcripts[tid];
            if (transcript.count <= 0) {
                continue;
            }
            
            cout << "\t" << left << setw(20) << transcript.name << right << setw(10) << transcript.len << setw(10) << transcript.count << endl;
        }
        
        cout << endl << endl;
        cout << "\tCompatibility matrix:" << endl;
        for (auto compItem = compMat.begin(); compItem != compMat.end(); compItem++) {
            cout << "\t" << setw(10) << right << compItem->second << "\t" << left;
            for (const auto local_tid : compItem->first) {
                const auto tid = gene.transcripts[local_tid];
                cout << "\t" << setw(20) << transcripts[tid].name;
            }
            cout << endl;
        }
        
        if (n.size() == 0) {
            continue;
        }
        
        cout << endl << endl;
        cout << "\tEstimated:" << endl;
        cout << "\t\t" << left << setw(20) << "Transcript Name" << right << setw(10) << "Count" << setw(14) << "Proportion" << endl;
        for (size_t t = 0; t < gene.transcripts.size(); ++t) {
            const auto tid = gene.transcripts[t];
            const auto& transcript = transcripts[tid];
            if (transcript.count <= 0 && n[t] < 1.0) {
                continue;
            }
            
            cout << "\t\t" << left << setw(20) << transcript.name
                 << right << setw(10) << (quint)n[t]
                 << setw(14) << quint(a[t] * 1000.0) / 1000.0;
            if (true) {
                cout << setw(10) << transcript.count;
            }
            cout << endl;
        }
    }
}

size_t Quant::addGeneIfNotExist(const string& geneName, quint geneLen)
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

size_t Quant::addTranscriptIfNotExist(const string& geneName, const string& transcriptName, quint transcriptLen)
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
    gene.transcriptLens.push_back(transcriptLen);
    gene.toLocalID[id] = gene.transcripts.size() - 1;
    
    return id;
}

bool Gene::hasAlts(quint left, quint right) const
{
    for (const auto& snp : snps) {
        if (snp.pos >= right) {
            break;
        } else if (snp.pos >= left){
            return true;
        }
    }
    return false;
}

void QuantCalc::calculate(const vector<quint>& tranLens,
                          const map<set<quint>, quint>& compMat,
                          vector<qfloat>& a,
                          vector<qfloat>& n)
{
    const quint nt = tranLens.size();
    assert (nt > 0);
    
    a = vector<qfloat>(nt, 1.0 / (qfloat)nt);
    n = vector<qfloat>(nt, 0.0);
    
    uint32_t iter = 0;
    vector<qfloat> an(nt, 0.0);
    while (true) {
        fill(n.begin(), n.end(), 0);
        for (auto compItem = compMat.begin(); compItem != compMat.end(); compItem++) {
            qfloat total = 0;
            for (auto t : compItem->first) {
                assert (t < a.size());
                total += a[t];
            }
            
            for (auto t : compItem->first) {
                auto nr = compItem->second; // number of reads
                n[t] += (nr * a[t] / total);
            }
        }
        
        qfloat total = 0.0;
        for (size_t t = 0; t < nt; ++t) {
            total += (n[t] / tranLens[t]);
        }
        
        for (size_t t = 0; t < nt; ++t) {
            an[t] = n[t] / tranLens[t] / total;
        }
        
        qfloat diff = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            diff += abs(an[i] - a[i]);
        }
        
        if (diff <= 0.00001) {
            break;
        }
        
#if 0
        cout << "iteration: " << iter << endl;
        for (size_t t = 0; t < nt; ++t) {
            cout
            << "\t" << setw(5) << t
            << "\t" << setw(8) << n[t]
            << "\t" << setw(10) << a[t] << " => " << setw(10) << an[t]
            << endl;
        }
        
        if (iter >= 5) {
            break;
        }
#endif
        
        a = an;
        ++iter;
    }
}
