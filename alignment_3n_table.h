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
#include "utility_3n_table.h"

extern bool uniqueOnly;
extern bool multipleOnly;
extern char convertFrom;
extern char convertTo;
extern char convertFromComplement;
extern char convertToComplement;

using namespace std;

/**
 * the class to store information from one SAM line
 */
class Alignment {
public:
    string chromosome;
    long long int location;
    int flag;
    bool mapped;
    char strand;
    string sequence;
    string quality;
    bool unique;
    string mapQ;
    int NH;
    vector<PosQuality> bases;
    CIGAR cigarString;
    MD_tag MD;
    unsigned long long readNameID;

    void initialize() {
        chromosome.clear();
        location = -1;
        flag = -1;
        mapped = false;
        MD.initialize();
        cigarString.initialize();
        sequence.clear();
        quality.clear();
        unique = false;
        mapQ.clear();
        NH = -1;
        bases.clear();
        readNameID = 0;
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
     * generate a hash value for readName
     */
     void getNameHash(string& readName) {
         readNameID = 0;
         for (int i = 0; i < readName.size(); i++) {
             readNameID = (readNameID << 6) | int(readName[i]);
         }
     }

    /**
     * extract the information from SAM line to Alignment.
     */
     void parseInfo(string* line) {
        int startPosition = 0;
        int endPosition = 0;
        int count = 0;

        while ((endPosition = line->find("\t", startPosition)) != string::npos) {
            if (count == 0) {
                string readName = line->substr(startPosition, endPosition - startPosition);
                getNameHash(readName);
            } else if (count == 1) {
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
                cigarString.loadString(line->substr(startPosition, endPosition - startPosition));
            } else if (count == 9) {
                sequence = line->substr(startPosition, endPosition - startPosition);
            } else if (count == 10) {
                quality = line->substr(startPosition, endPosition - startPosition);
            } else if (count > 10) {
                if (startWith(line, startPosition, "MD")) {
                    MD.loadString(line->substr(startPosition + 5, endPosition - startPosition - 5));
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
    }

    /**
     *  scan all base in read sequence label them if they are qualified.
     */
    void appendBase() {
        if (!mapped) {
            return;
        }

        bases.reserve(sequence.size());
        for (int i = 0; i < sequence.size(); i++) {
            bases.emplace_back(i);
        }
        int pos = adjustPos();

        string match;
        while (MD.getNextSegment(match)) {
            if (isdigit(match.front())) { // the first char of match is digit this is match
                int len = stoi(match);
                for (int i = 0; i < len; i++) {
                    while (bases[pos].remove) {
                        pos++;
                    }
                    if ((strand == '+' && sequence[pos] == convertFrom) ||
                        (strand == '-' && sequence[pos] == convertFromComplement)) {
                        bases[pos].setQual(quality[pos], false);
                    } else {
                        bases[pos].remove = true;
                    }
                    pos ++;
                }
            } else if (isalpha(match.front())) { // this is mismatch or conversion
                char refBase = match.front();
                // for + strand, it should have C->T change
                // for - strand, it should have G->A change
                while (bases[pos].remove) {
                    pos++;
                }

                if ((strand == '+' && refBase == convertFrom && sequence[pos] == convertTo) ||
                    (strand == '-' && refBase == convertFromComplement && sequence[pos] == convertToComplement)){
                    bases[pos].setQual(quality[pos], true);
                } else {
                    bases[pos].remove = true;
                }
                pos ++;
            } else { // deletion. do nothing.

            }
        }
    }

    /**
     * adjust the reference position in bases
     */
    int  adjustPos() {

        int readPos = 0;
        int returnPos = 0;
        int seqLength = sequence.size();

        char cigarSymbol;
        int cigarLen;

        while (cigarString.getNextSegment(cigarLen, cigarSymbol)) {
            if (cigarSymbol == 'S') {
                if (readPos == 0) { // soft clip is at the begin of the read
                    returnPos = cigarLen;
                    for (int i = cigarLen; i < seqLength; i++) {
                        bases[i].refPos -= cigarLen;
                    }
                } else { // soft clip is at the end of the read
                    // do nothing
                }
                readPos += cigarLen;
            } else if (cigarSymbol == 'N') {
                for (int i = readPos; i < seqLength; i++) {
                    bases[i].refPos += cigarLen;
                }
            } else if (cigarSymbol == 'M') {
                for (int i = readPos; i < readPos+cigarLen; i++) {
                    bases[i].remove = false;
                }
                readPos += cigarLen;
            } else if (cigarSymbol == 'I') {
                for (int i = readPos + cigarLen; i < seqLength; i++) {
                    bases[i].refPos -= cigarLen;
                }
                readPos += cigarLen;
            } else if (cigarSymbol == 'D') {
                for (int i = readPos; i < seqLength; i++) {
                    bases[i].refPos += cigarLen;
                }
            }
        }
        return returnPos;
    }

};

#endif //ALIGNMENT_3N_TABLE_H
