/*
 * Copyright 2020, Yun (Leo) Zhang <imzhangyun@gmail.com>
 *
 * This file is part of HISAT-3N.
 *
 * HISAT-3N is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT-3N is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT-3N.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ALIGNMENT_3N_TABLE_H
#define ALIGNMENT_3N_TABLE_H

#include <string>
#include <regex>
#include "utility_3n_table.h"

extern bool uniqueOnly;
extern bool multipleOnly;
extern char convertFrom;
extern char convertTo;
extern char convertFromComplement;
extern char convertToComplement;

using namespace std;

regex MDExp("[0-9]+|[A-Z]|\\^[A-Z]+");
regex cigarExp("\\d+\\w");

/**
 * the class to store information from one SAM line
 */
class Alignment {
public:
    string chromosome;
    long long int location;
    int flag;
    bool mapped;
    string MD;
    string cigarString;
    char strand;
    string sequence;
    string quality;
    bool unique;
    string mapQ;
    int NH;
    vector<PosQuality> bases;

    void initialize() {
        chromosome.clear();
        location = -1;
        flag = -1;
        mapped = false;
        MD.clear();
        cigarString.clear();
        sequence.clear();
        quality.clear();
        unique = false;
        mapQ.clear();
        NH = -1;
        bases.clear();
    }

    /**
     * for start position in input Line, check if it contain the target information.
     */
    bool startWith(string* inputLine, int startPosition, string tag){
        for (int i = 0; i < tag.size(); i++){
            if (inputLine->at(startPosition+i) != tag[i]){
                return false;
            }
        }
        return true;
    }

    /**
     * extract the information from SAM line to Alignment.
     */
     void parseInfo(string* line) {
        int startPosition = 0;
        int endPosition = 0;
        int count = 0;

        while ((endPosition = line->find("\t", startPosition)) != string::npos) {
            if (count == 1) {
                flag = stoi(line->substr(startPosition, endPosition - startPosition));
                mapped = (flag & 4) == 0;
            } else if (count == 2) {
                chromosome = line->substr(startPosition, endPosition - startPosition);
            } else if (count == 3) {
                location = stoll(line->substr(startPosition, endPosition - startPosition));
            } else if (count == 4) {
                mapQ = line->substr(startPosition, endPosition - startPosition);
                if (mapQ == "1") {
                    unique = false;
                } else {
                    unique = true;
                }
            } else if (count == 5) {
                cigarString = line->substr(startPosition, endPosition - startPosition);
            } else if (count == 9) {
                sequence = line->substr(startPosition, endPosition - startPosition);
            } else if (count == 10) {
                quality = line->substr(startPosition, endPosition - startPosition);
            } else if (count > 10) {
                if (startWith(line, startPosition, "MD")) {
                    MD = line->substr(startPosition + 5, endPosition - startPosition - 5);
                } else if (startWith(line, startPosition, "NM")) {
                    NH = stoi(line->substr(startPosition + 5, endPosition - startPosition - 5));
                } else if (startWith(line, startPosition, "YZ")) {
                    strand = line->at(endPosition-1);
                }
            }
            startPosition = endPosition + 1;
            count++;
        }
     }

    /**
     * parse the sam line to alignment information
     */
    void parse(string* line) {
        initialize();
        parseInfo(line);
        if ((uniqueOnly && !unique) || (multipleOnly && unique)) {
            return;
        }
        appendBase();
        adjustPos();
    }

    /**
     * find all my target base, append to bases.
     */
    void appendBase() {
        smatch match;
        auto searchStart(MD.cbegin());
        int pos = 0;
        bases.reserve(sequence.size()/4);
        while (regex_search(searchStart, MD.cend(), match, MDExp)) {
            string test = match.str();
            if (isdigit(match.str().front())) { // the first char of match is digit this is match
                int len = stoi(match.str());
                for (int i = 0; i < len; i++) {
                    if ((strand == '+' && sequence[pos] == convertFrom) ||
                        (strand == '-' && sequence[pos] == convertFromComplement)) {
                        bases.emplace_back(pos, quality[pos], false);
                    }
                    pos ++;
                }
            } else if (isalpha(match.str().front())) { // this is mismatch or conversion
                char refBase = match.str().front();
                // for + strand, it should have C->T change
                // for - strand, it should have G->A change
                if ((strand == '+' && refBase == convertFrom && sequence[pos] == convertTo) ||
                    (strand == '-' && refBase == convertFromComplement && sequence[pos] == convertToComplement)){
                    bases.emplace_back(pos, quality[pos], true);
                }
                pos ++;
            } else { // deletion
                pos += match.str().size()-1;
            }
            searchStart = match.suffix().first;
        }
    }

    /**
     * adjust the read position in methylatedBases and unmethylatedBases to reference position
     */
    void adjustPos() {
        smatch match;
        auto searchStart(cigarString.cbegin());
        int readPos = 0;

        while (regex_search(searchStart, cigarString.cend(), match, cigarExp)) {
            string cigar = match.str();
            char cigarSymbol = cigar.back();
            int cigarLen = stoi(cigar.substr(0, cigar.size()-1));

            if (cigarSymbol == 'S') {
                if (readPos == 0) { // soft clip is at the begin of the read
                    for (int i = 0; i < bases.size(); i++) {
                        if (bases[i].readPos <= cigarLen) {
                            bases.erase(bases.begin()+i);
                        }
                    }
                } else { // soft clip is at the end of the read
                    for (int i = 0; i < bases.size(); i++) {
                        if (bases[i].readPos >= sequence.size()-cigarLen) {
                            bases.erase(bases.begin()+i);
                        }
                    }
                }
                readPos += cigarLen;
            } else if (cigarSymbol == 'N') {
                for (int i = 0; i < bases.size(); i++) {
                    if (bases[i].readPos >= readPos) {
                        bases[i].readPos += cigarLen;
                    }
                }
            } else if (cigarSymbol == 'M') {
                readPos += cigarLen;
            } else if (cigarSymbol == 'I') {
                readPos += cigarLen;
            } else if (cigarSymbol == 'D') {
                for (int i = 0; i < bases.size(); i++) {
                    if (bases[i].readPos >= readPos) {
                        bases[i].readPos += cigarLen;
                    }
                }
            }
            searchStart = match.suffix().first;
        }
    }

};

#endif //ALIGNMENT_3N_TABLE_H
