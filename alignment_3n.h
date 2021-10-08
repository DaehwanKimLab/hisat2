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

#ifndef HISAT2_ALIGNMENT_3N_H
#define HISAT2_ALIGNMENT_3N_H

#include <string>
#include <vector>
#include <fstream>
#include <limits>
#include <algorithm>
#include "sstring.h"
#include "util.h"
#include "hisat2lib/ht2.h"
#include "read.h"
#include "outq.h"
#include "reference.h"
#include <unistd.h>
#include <queue>
#include "position_3n.h"
#include "utility_3n.h"


extern char usrInput_convertedFrom;
extern char usrInput_convertedTo;
extern char usrInput_convertedFromComplement;
extern char usrInput_convertedToComplement;

extern char hs3N_convertedFrom;
extern char hs3N_convertedTo;
extern char hs3N_convertedFromComplement;
extern char hs3N_convertedToComplement;

extern vector<ht2_handle_t> repeatHandles;
extern struct ht2_index_getrefnames_result *refNameMap;
extern int repeatLimit;
extern bool uniqueOutputOnly;

using namespace std;

struct ReportingMetrics;

/**
 * the data structure to store all information of one alignment result.
 */
class Alignment {
public:
    // basic information
    BTString readName;
    int flag;
    BTString chromosomeName;
    int chromosomeIndex; // the chromosome index to use getStretch() function
    long long int location;
    BTString MAPQ;
    BTString cigarString;
    vector<Cigar> cigarSegments;
    int cigarLength; // this is the length that read cover genome.
    BTString pairToChromosome;
    long long int pairToLocation;
    long long int pairingDistance;
    BTString readSequence;
    BTString readQuality;
    //tags
    int AS; // alignment score
    int NH; // number of alignment
    int XM; // number of mismatch
    int NM; // edit distance
    int YS; // mate's AS
    BTString MD;
    BTString YT; //"UU" for single-end. "CP" for concordant alignment, "DP" for disconcordant alignment, "UP" for else.
    // special tags in HISAT-3N
    int Yf; // number of conversion.
    int Zf; // number of unconverted base.
    char YZ;  // this tag shows alignment strand:
              // + for REF strand (conversionCount[0] is equal or smaller than conversionCount[1]),
              // - for REF-RC strand (conversionCount[1] is smaller)
    // unChanged tags
    BTString unChangedTags;
    BTString passThroughLine; // this is controlled by print_xr_ in SamConfig

    // for pairScore calculation
    static const int maxPairDistance = 500000;
    static const int penaltyFreeDistance_DNA = 1000;
    static const int penaltyFreeDistance_RNA = 100000;
    static const int distancePenaltyFraction_DNA = 100;
    static const int distancePenaltyFraction_RNA = 1000;
    static const int ASPenalty = 100;
    static const int concordantScoreBounce = 500000;

    // intermediate variable
    bool outputted = false; // whether the alignment is outputted.
    bool DNA = false;
    int cycle_3N; // indicate which cycle_3N make this alignment result. 0 or 3 for repeatHandles[0], else repeatHandles[1]
    bool paired;
    bool forward;
    bool mapped;
    bool concordant;
    int pairSegment; // 0 for first segment, 1 for second segment.
    struct ht2_repeat_expand_result *repeatResult = nullptr;
    int pairScore; // to identify the better pair
    bool mateMapped; // to adjust the YT tag
    bool repeat;
    bool pairToRepeat;
    RepeatMappingPositions repeatPositions; // to store the expanded repeat information
    int conversionCount[2] = {0}; // there are two type of conversion could happen, save the number of conversion separately.
    int unConversionCount[2] = {0}; // save the unconverted base count.
    string intToBase = "ACGTN";

    void initialize() {
        readName.clear();
        flag = -1;
        chromosomeName.clear();
        chromosomeIndex = -1;
        location = 0;
        MAPQ.clear();
        cigarString.clear();
        cigarSegments.clear();
        cigarLength = 0;
        pairToChromosome.clear();
        pairToLocation = 0;
        pairingDistance = 0;
        readSequence.clear();
        readQuality.clear();

        AS = numeric_limits<int>::min();
        NH = 0;
        XM = 0;
        NM = 0;
        YS = 0;
        MD.clear();
        YT.clear();
        Yf = 0;
        unChangedTags.clear();

        outputted = false;
        DNA = false;
        cycle_3N = -1;
        paired = false;
        forward = false;
        mapped = false;
        concordant = false;
        pairSegment = 0;
        if (repeatResult != nullptr) {
            free(repeatResult);
            repeatResult = nullptr;
        }
        pairScore = numeric_limits<int>::min();
        mateMapped = false;

        repeat = false;
        pairToRepeat = false;
        repeatPositions.initialize();
        conversionCount[0] = 0;
        conversionCount[1] = 0;
        unConversionCount[0] = 0;
        unConversionCount[1] = 0;
        passThroughLine.clear();
    }

    Alignment() {
        initialize();
    }

    ~Alignment() {
        if (repeatResult != nullptr) free(repeatResult);
    }

    /**
     * change YS tag for output.
     */
    void setYS (Alignment* input) {
        YS = input->AS;
    }

    void setYS(RepeatMappingPosition* input) {
        YS = input->AS;
    }

    /**
     * change concordant status and flag
     */
    void setConcordant(bool concordant_) {
        concordant = concordant_;
        if (concordant) {
            flag |= 2;
        } else {
            flag &= ~((int)2);
        }
    }

    /**
     * change mateMapped status and flag
     */
    void setMateMappingFlag(long long int *mateLocation) {
        if (mateLocation == NULL) { return; }
        mateMapped = *mateLocation != 0;
        if ((flag&8) && mateMapped) flag -= 8;
        else if (!(flag&8) && !mateMapped) flag += 8;
    }

    /**
     * change YT tag base on mateMapped and concordant information.
     */
    void setYT() {
        if (paired) {
            if (!mateMapped) {
                YT = "UP";
                return;
            }
            if (concordant) { YT = "CP"; }
            else { YT = "DP"; }
        } else {
            YT = "UU";
        }
    }

    /**
     * update NH and MAPQ by number of alignment.
     */
    void updateNH(int nAlignment) {
        if (!mapped) return;
        NH = nAlignment;
        if (nAlignment == 0) return;
        else if (nAlignment == 1) MAPQ = "60";
        else MAPQ = "1";
    }

    /**
     * extract information from flag, change flag to secondary alignment.
     */
    void extractFlagInfo() {
        paired = (flag & 1) != 0;
        forward = (flag & 16) == 0;
        if ((flag & 256) == 0) { // change all read to secondary alignment
            flag += 256;
        }
        mapped = (flag & 4) == 0;
        if (flag & 128) {
            pairSegment = 1;
        } else {
            pairSegment = 0; // it could be the first pair segment or it is single read.
        }
        concordant = (flag & 2) != 0;
        if (!mapped) {
            repeat = false;
        }
    }

    /**
     * calculate the pairScore for a pair of alignment result. Output pair Score and number of pair.
     * Do not update their pairScore.
     */
    int calculatePairScore(Alignment *inputAlignment, int &nPair);

    /**
     * make YZ tag.
     * if the conversion type 0 is less, the read is mapped to REF (+).
     * if the conversion type 1 is less, the read is mapped to REF-RC (-).
     */
    void makeYZ(char &YZ_string) {
        if (conversionCount[0] >= conversionCount[1]) {
            YZ_string = '+';
        } else {
            YZ_string = '-';
        }
    }

    /**
     * expand the repeat mapping location and construct MD for each location.
     * return ture if there is any mapping location pass the filter, else, false.
     */
    bool constructRepeatMD(BitPairReference* bitReference, MappingPositions &alignmentPositions) {

        if (!mapped) {
            return true;
        }

        // expand the repeat locations
        ht2_error_t err = ht2_repeat_expand((cycle_3N == 0 || cycle_3N == 3) ? repeatHandles[0] : repeatHandles[1],
                                            chromosomeName.toZBuf(),
                                            location - 1,
                                            readSequence.length(),
                                            &repeatResult);

        BTString chromosomeRepeat;
        long long int locationRepeat;
        for (int i = 0; i < repeatResult->count; i++) {
            struct ht2_position *pos = &repeatResult->positions[i];
            chromosomeRepeat = refNameMap->names[pos->chr_id];
            for (int j = 0; j < chromosomeRepeat.length(); j++) {
                if (chromosomeRepeat[j] == ' ') {
                    chromosomeRepeat.trimEnd(chromosomeRepeat.length() - j);
                    break;
                }
            }
            locationRepeat = (pos->pos) + 1;
            bool genomeForward = pos->direction == 0;
            if (!genomeForward) { continue; } // if the repeat mapping direction is different to the designed direction, ignore it.

            // if the mapping location is already exist, continue.
            if (alignmentPositions.positionExist(chromosomeRepeat, locationRepeat, pairSegment)){
                continue;
            }

            // get reference sequence
            ASSERT_ONLY(SStringExpandable<uint32_t> destU32);
            SStringExpandable<char> raw_refbuf;
            raw_refbuf.resize(cigarLength + 16);
            raw_refbuf.clear();
            int off = bitReference->getStretch(
                    reinterpret_cast<uint32_t*>(raw_refbuf.wbuf()),
                    (size_t)pos->chr_id,
                    (size_t)max<int>(locationRepeat-1, 0),
                    (size_t)cigarLength ASSERT_ONLY(, destU32));
            char* refSeq = raw_refbuf.wbuf() + off;
            BTString refSequence;
            refSequence.resize(cigarLength);
            for (int j = 0; j < cigarLength; j++) {
                refSequence.set(intToBase[*(refSeq + j)], j);
            }

            // check whether the refSequence is exist. if do, directly append the repeat.
            int repeatPositionsIndex;
            if (repeatPositions.sequenceExist(refSequence, repeatPositionsIndex)) {
                repeatPositions.append(chromosomeRepeat, locationRepeat, repeatPositionsIndex);
                continue;
            }

            BTString newMD;
            int newMismatch = 0;
            char repeatYZ;
            if (!constructRepeatMD(refSequence, newMD, newMismatch, repeatYZ)) {
                continue;
            }

            int newXM = XM + newMismatch;
            int newNM = NM + newMismatch;
            int newAS = AS - 6*newMismatch;
            repeatPositions.append(locationRepeat, chromosomeRepeat, refSequence,newAS, newMD, newXM, newNM, Yf, Zf, repeatYZ);

            // if there are too many mappingPosition exist return.
            if (repeatPositions.size() >= repeatLimit || alignmentPositions.size() > repeatLimit) {
                return true;
            }
        }
        if (repeatPositions.size() == 0) {
            return false;
        } else {
            return true;
        }
    }

    /**
     * for each repeat mapping position, construct its MD
     * return true if the mapping result does not have a lot of mismatch, else return false.
     */
    bool constructRepeatMD(BTString &refSeq, BTString &newMD_String, int &newMismatch, char &repeatYZ) {
        char buf[1024];

        conversionCount[0] = 0;
        conversionCount[1] = 0;
        unConversionCount[0] = 0;
        unConversionCount[1] = 0;

        int readPos = 0;
        long long int refPos = 0;
        int count = 0;
        int newXM = 0;

        char cigarSymbol;
        int cigarLen;
        for (int i = 0; i < cigarSegments.size(); i++) {
            cigarSymbol = cigarSegments[i].getLabel();
            cigarLen = cigarSegments[i].getLen();

            if (cigarSymbol == 'S') {
                readPos += cigarLen;
            } else if (cigarSymbol == 'N') {
                refPos += cigarLen;
            } else if (cigarSymbol == 'M') {
                for (int j = 0; j < cigarLen; j++) {
                    char readChar = readSequence[readPos];
                    char refChar = refSeq[refPos];
                    if (readChar == refChar) {
                        if (refChar == usrInput_convertedFrom)
                        {
                            unConversionCount[0]++;
                        }
                        else if (refChar == usrInput_convertedFromComplement)
                        {
                            unConversionCount[1]++;
                        }
                        count++;
                    } else {// mismatch
                        // output matched count
                        if (count != 0) {
                            itoa10<int>(count, buf);
                            newMD_String.append(buf);
                            count = 0;
                        }
                        // output mismatch
                        if (!newMD_String.empty() && isalpha(newMD_String[newMD_String.length()-1])) {
                            newMD_String.append('0');
                        }
                        if ((readChar == usrInput_convertedTo) && (refChar == usrInput_convertedFrom)) {
                            conversionCount[0]++;
                        } else if ((readChar == usrInput_convertedToComplement) && (refChar == usrInput_convertedFromComplement)) {
                            conversionCount[1]++;
                        } else {
                            // real mismatch
                            newXM++;
                        }
                        newMD_String.append(refChar);
                    }
                    if ((conversionCount[0] >= conversionCount[1] ? conversionCount[1]:conversionCount[0])> readSequence.length()/25) {
                        return false;
                    }
                    readPos++;
                    refPos++;
                }
            } else if (cigarSymbol == 'I') {
                readPos += cigarLen;
            } else if (cigarSymbol == 'D') {
                newMD_String.append('^');
                for (int j = 0; j < cigarLen; j++) {
                    newMD_String.append(refSeq[refPos]);
                    refPos++;
                }
            }
        }

        if (count != 0) {
            itoa10<int>(count, buf);
            newMD_String.append(buf);
        }
        if (isalpha(newMD_String[0])) { newMD_String.insert('0', 0); }
        if (isalpha(newMD_String[newMD_String.length()-1])) { newMD_String.append('0'); }

        int badConversion = 0;
        // choose the smaller conversionCount as bad count (will give penalty). choose the bigger conversionCount as Yf.
        if (conversionCount[0] >= conversionCount[1]) {
            badConversion = conversionCount[1];
            Yf = conversionCount[0];
            Zf = unConversionCount[0];
        } else {
            badConversion = conversionCount[0];
            Yf = conversionCount[1];
            Zf = unConversionCount[1];
        }

        newXM += badConversion;
        newMismatch = newXM - XM;

        if (newMismatch < 0){
            newMismatch = 0;
        }

        makeYZ(repeatYZ);
        return true;
    }

    /**
     * for each non-repeat mapping position, construct its MD
     * return true if the mapping result does not have a lot of mismatch, else return false.
     */
    bool constructMD(BitPairReference* bitReference) {
        if (!mapped) {
            return true;
        }
        char buf[1024];
        MD.clear();

        ASSERT_ONLY(SStringExpandable<uint32_t> destU32);
        SStringExpandable<char> raw_refbuf;
        raw_refbuf.resize(cigarLength + 16);
        raw_refbuf.clear();
        int off = bitReference->getStretch(
                reinterpret_cast<uint32_t*>(raw_refbuf.wbuf()),
                (size_t)chromosomeIndex,
                (size_t)max<int>(location-1, 0),
                (size_t)cigarLength ASSERT_ONLY(, destU32));
        char* refSeq = raw_refbuf.wbuf() + off;

        int readPos = 0;
        long long int refPos = 0;
        int count = 0;
        int newXM = 0;

        char cigarSymbol;
        int cigarLen;
        for (int i = 0; i < cigarSegments.size(); i++) {
            cigarSymbol = cigarSegments[i].getLabel();
            cigarLen = cigarSegments[i].getLen();
            if (cigarSymbol == 'S') {
                readPos += cigarLen;
            } else if (cigarSymbol == 'N') {
                refPos += cigarLen;
            } else if (cigarSymbol == 'M') {
                for (int j = 0; j < cigarLen; j++) {
                    char readChar = readSequence[readPos];
                    char refChar = intToBase[*(refSeq + refPos)];
                    if (readChar == refChar) {
                        if (refChar == usrInput_convertedFrom)
                        {
                            unConversionCount[0]++;
                        }
                        else if (refChar == usrInput_convertedFromComplement)
                        {
                            unConversionCount[1]++;
                        }
                        count++;
                    } else {// mismatch
                        // output matched count
                        if (count != 0) {
                            itoa10<int>(count, buf);
                            MD.append(buf);
                            count = 0;
                        }
                        // output mismatch
                        if (!MD.empty() && isalpha(MD[MD.length()-1])) {
                            MD.append('0');
                        }

                        if ((readChar == usrInput_convertedTo) && (refChar == usrInput_convertedFrom)) {
                            conversionCount[0]++;
                        } else if ((readChar == usrInput_convertedToComplement) && (refChar == usrInput_convertedFromComplement)) {
                            conversionCount[1]++;
                        } else {
                            // real mismatch
                            newXM++;
                        }
                        MD.append(refChar);
                    }
                    if ((conversionCount[0] >= conversionCount[1] ? conversionCount[1]:conversionCount[0])> readSequence.length()/25) {
                        return false;
                    }
                    readPos++;
                    refPos++;
                }
            } else if (cigarSymbol == 'I') {
                readPos += cigarLen;
            } else if (cigarSymbol == 'D') {
                if (count != 0) {
                    itoa10<int>(count, buf);
                    MD.append(buf);
                    count = 0;
                }
                MD.append('^');
                for (int j = 0; j < cigarLen; j++) {
                    MD.append(intToBase[*(refSeq + refPos)]);
                    refPos++;
                }
            }
        }

        if (count != 0) {
            itoa10<int>(count, buf);
            MD.append(buf);
        }
        if (isalpha(MD[0])) { MD.insert('0', 0); }
        if (isalpha(MD[MD.length()-1])) { MD.append('0'); }

        int badConversion = 0;
        // choose the smaller conversionCount as bad count (will give penalty). choose the bigger conversionCount as Yf.
        if (conversionCount[0] >= conversionCount[1]) {
            badConversion = conversionCount[1];
            Yf = conversionCount[0];
            Zf = unConversionCount[0];
        } else {
            badConversion = conversionCount[0];
            Yf = conversionCount[1];
            Zf = unConversionCount[1];
        }

        newXM += badConversion;
        newXM -= XM;

        if (newXM < 0){
            newXM = 0;
        }

        makeYZ(YZ);
        NM += newXM;
        XM += newXM;
        AS = AS - 6*newXM;
        BTString tmp;
        if (pairToRepeat) {
            repeatPositions.append(location, chromosomeName, tmp, AS, MD, XM, NM, Yf, Zf, YZ);
        }
        return true;
    }

    /**
     * output the tags for non-repeat alignment.
     */
    void outputTags(BTString& o) {
        char buf[1024];
        if (mapped) {
            o.append('\t');
            // AS
            assert(AS <= 0);
            o.append("AS:i:");
            itoa10<int>(AS, buf);
            o.append(buf);
            o.append('\t');
            // NH
            assert(NH > 0);
            o.append("NH:i:");
            itoa10<int>(NH, buf);
            o.append(buf);
            o.append('\t');
            // XM
            assert(XM >= 0);
            o.append("XM:i:");
            itoa10<int>(XM, buf);
            o.append(buf);
            o.append('\t');
            // NM
            assert(NM >= 0);
            o.append("NM:i:");
            itoa10<int>(NM, buf);
            o.append(buf);
            o.append('\t');
            // MD
            assert(!MD.empty());
            o.append("MD:Z:");
            o.append(MD.toZBuf());
            o.append('\t');
            // YS
            if (paired && mateMapped) {
                o.append("YS:i:");
                itoa10<int>(YS, buf);
                o.append(buf);
                o.append('\t');
            }
            // YZ
            o.append("YZ:A:");
            o.append(YZ);
            o.append('\t');
            // Yf
            o.append("Yf:i:");
            itoa10<int>(Yf, buf);
            o.append(buf);
            o.append('\t');
            //Zf
            o.append("Zf:i:");
            itoa10<int>(Zf, buf);
            o.append(buf);
        }
        // unchanged Tags
        if (!unChangedTags.empty()) {
            o.append('\t');
            o.append(unChangedTags.toZBuf());
        }
        o.append(passThroughLine.toZBuf());
    }

    /**
     * output the tags for repeat alignment.
     */
    void outputTags(BTString& o, RepeatMappingPosition* repeatInfo){
        // this function is for repeat alignment output.
        char buf[1024];
        o.append('\t');
        // AS
        assert(AS <= 0);
        o.append("AS:i:");
        itoa10<int>(repeatInfo->AS, buf);
        o.append(buf);
        o.append('\t');
        // NH
        assert(NH > 0);
        o.append("NH:i:");
        itoa10<int>(NH, buf);
        o.append(buf);
        o.append('\t');
        // XM
        assert(XM >= 0);
        o.append("XM:i:");
        itoa10<int>(repeatInfo->XM, buf);
        o.append(buf);
        o.append('\t');
        // NM
        assert(NM >= 0);
        o.append("NM:i:");
        itoa10<int>(repeatInfo->NM, buf);
        o.append(buf);
        o.append('\t');
        // MD
        assert(!MD.empty());
        o.append("MD:Z:");
        o.append(repeatInfo->MD.toZBuf());
        o.append('\t');
        // YS
        if (paired) {
            o.append("YS:i:");
            itoa10<int>(YS, buf);
            o.append(buf);
            o.append('\t');
        }
        //YT
        o.append("YT:Z:");
        o.append(YT.toZBuf());
        o.append('\t');
        // YS
        if (paired && mateMapped) {
            o.append("YS:i:");
            itoa10<int>(YS, buf);
            o.append(buf);
            o.append('\t');
        }
        // YZ
        o.append("YZ:A:");
        o.append(repeatInfo->YZ);
        o.append('\t');
        // Yf
        o.append("Yf:i:");
        itoa10<int>(repeatInfo->Yf, buf);
        o.append(buf);
        o.append('\t');
        // Zf
        o.append("Zf:i:");
        itoa10<int>(repeatInfo->Zf, buf);
        o.append(buf);
        o.append('\t');
        // unchanged Tags
        if (!unChangedTags.empty()) {
            o.append('\t');
            o.append(unChangedTags.toZBuf());
        }
        o.append(passThroughLine.toZBuf());
    }

    /**
     * output alignment. this function is for both repeat and non-repeat alignment.
     */
    void outputAlignment (BTString& o, RepeatMappingPosition* repeatInfo, long long int* oppoLocation, bool& primaryAlignment) {
        BTString* outputChromosome;
        long long int* outputLocation;

        if (repeatInfo == NULL) {
            if (outputted) { return; }
            outputted = true;
            outputChromosome = &chromosomeName;
            outputLocation = &location;
        } else {
            if (repeatInfo->outputted) { return; }
            repeatInfo->outputted = true;
            outputChromosome = &repeatInfo->repeatChromosome;
            outputLocation = &repeatInfo->repeatLocation;
        }

        //setMateMappingFlag(oppoLocation);
        setYT();

        char buf[1024];
        // readName
        o.append(readName.toZBuf());
        o.append('\t');
        // flag, if it is primary alignment, -256
        assert(flag >=0);
        itoa10<int>(flag-primaryAlignment*256, buf);
        o.append(buf);
        o.append('\t');
        // chromosome
        assert(!outputChromosome->empty());
        o.append(outputChromosome->toZBuf());
        o.append('\t');
        // location
        assert(*outputLocation >= 0);
        itoa10<int>(*outputLocation, buf);
        o.append(buf);
        o.append('\t');

        //MAPQ
        o.append(MAPQ.toZBuf());
        o.append('\t');
        // cigar
        o.append(cigarString.toZBuf());
        o.append('\t');
        // pair to chromosome
        if (paired && *oppoLocation!=0) {
            o.append("=");
            o.append('\t');
        } else {
            o.append("*");
            o.append('\t');
        }
        // pair to location
        if (paired) {
            itoa10<int>(*oppoLocation, buf);
            o.append(buf);
            o.append('\t');
        } else {
            o.append('0');
            o.append('\t');
        }
        // pairing distance
        if (paired) {
            itoa10<int>(*oppoLocation - *outputLocation, buf);
            o.append(buf);
            o.append('\t');
        } else {
            o.append('0');
            o.append('\t');
        }
        // read sequence
        o.append(readSequence.toZBuf());
        o.append('\t');
        // read quality
        o.append(readQuality.toZBuf());

        // make sure there is no '\t' at the beginning of unChangedTags
        while (!unChangedTags.empty() && unChangedTags[0] == '\t') {
            unChangedTags.remove(0);
        }

        // tags
        if (repeatInfo == NULL) {
            outputTags(o);
            o.append('\n');
        } else {
            outputTags(o, repeatInfo->flagInfoIndex == -1?repeatInfo:&repeatPositions.positions[repeatInfo->flagInfoIndex]);
            o.append('\n');
        }
    }

    /**
     * return true if two location is concordant.
     * return false, if there are not concordant or too far (>maxPairDistance).
     */
    static bool isConcordant(long long int &location1, bool &forward1, long long int &location2, bool &forward2);

    /**
     * this is the basic function to calculate DNA pair score.
     * if the distance between 2 alignments is more than penaltyFreeDistance_DNA, we reduce the score by the distance/100.
     * if two alignment is concordant we add concordantScoreBounce to make sure to select the concordant pair as best pair.
     */
    static int calculatePairScore_DNA (long long int &location0, int& AS0, bool& forward0, long long int &location1, int &AS1, bool &forward1, bool& concordant);

    /**
     * this is the basic function to calculate RNA pair score.
     * if the distance between 2 alignments is more than penaltyFreeDistance_RNA, we reduce the score by the distance/1000.
     * if two alignment is concordant we add concordantScoreBounce to make sure to select the concordant pair as best pair.
     */
    static int calculatePairScore_RNA (long long int &location0, int& XM0, bool& forward0, long long int &location1, int &XM1, bool &forward1, bool& concordant);
};

/**
 * the data structure to store, process, and output all Alignment
 */
class Alignments {
public:
    vector<Alignment*> alignments; // pool to store current alignment result.
    vector<Alignment*> freeAlignments; // free pointer pool for new alignment result. after output a alignment, return the pointer back to this pool.

    TReadId previousReadID;
    MappingPositions alignmentPositions; // the pool to save all alignment position

    BTString readName[2]; // the read name could be different for segment 1 and segment 2.
    BTDnaString readSequence[2]; // save the read sequence for output.
    BTString qualityScore[2]; // save the quality score for output.

    bool paired;
    const int repeatPoolLimit = 20; // this is the maximum number of repeat alignment we allowed.
    bool multipleAligned; // check whether we have multiple alignment, it is work unique mode.

    const int maxPairScore = 500000; // maximum pair score, if pairScore == maxPairScore, both math are perfect match and the pairDistance is small.

    BitPairReference* bitReference; // bit pair reference sequence
    bool DNA;
    int nRepeatAlignment; // count number of repeat alignment we received, for short sequence we could receive a lot of repeat alignment result.

    BTString passThroughLines[2];

    void initialize() {
        alignmentPositions.initialize();
        paired = false;
        multipleAligned = false;
        nRepeatAlignment = 0;

        for (int i = 0; i < 2; i++) {
            readName[i].clear();
            readSequence[i].clear();
            qualityScore[i].clear();
            passThroughLines[i].clear();
        }
        for (int i = 0; i < alignments.size(); i++) {
            alignments[i]->initialize();
            freeAlignments.push_back(alignments[i]);
        }
        alignments.clear();
    }

    Alignments(BitPairReference* ref, bool inputDNA): bitReference(ref), DNA(inputDNA) {
        initialize();
    }

    ~Alignments() {
        while (!freeAlignments.empty()) {
            delete freeAlignments.back();
            freeAlignments.pop_back();
        }
        for (int i = 0; i < alignments.size(); i++) {
            delete alignments[i];
        }
    }

    /**
     * get sequence for rd. if it already exist, ignore it.
     */
    void getSequence(const Read& rd) {
        int pairSegment = rd.mate == 0? rd.mate : rd.mate-1;
        if (readName[pairSegment].empty()) { readName[pairSegment] = rd.name; }
        if (readSequence[pairSegment].empty()) { readSequence[pairSegment] = rd.originalFw; }
        if (qualityScore[pairSegment].empty()) { qualityScore[pairSegment] = rd.qual; }
    }

    /**
     * return true if we want to receive more new alignment.
     */
    bool acceptNewAlignment() {
        if (uniqueOutputOnly && multipleAligned ||
            alignmentPositions.nBestSingle >= repeatLimit ||
            nRepeatAlignment > repeatPoolLimit ||
            alignmentPositions.nBestPair >= repeatLimit) {
            return false;
        }
        return true;
    }

    /**
     * return the alignment back to freeAlignment pool.
     */
    void returnToFreeAlignments (Alignment*& currentAlignment) {
        currentAlignment->initialize();
        freeAlignments.push_back(currentAlignment);
    }

    /**
     * get a Alignment pointer from freeAlignments, if freeAlignment is empty, make a new Alignment.
     */
    void getFreeAlignmentPointer(Alignment*& newAlignment) {
        if (!freeAlignments.empty()) {
            newAlignment = freeAlignments.back();
            freeAlignments.pop_back();
        } else {
            newAlignment = new Alignment();
        }
    }

    /**
     * receive alignment information from AlnSink3NSam::appendMate() and append it to alignment pool.
     */
    void append(Alignment *newAlignment) {

        newAlignment->extractFlagInfo();
        paired = newAlignment->paired;
        newAlignment->DNA = DNA;
        if (passThroughLines[newAlignment->pairSegment].empty()) {
            passThroughLines[newAlignment->pairSegment] = newAlignment->passThroughLine;
        }

        // check if the alignment is already exist. if exist, ignore it.
        if (!alignmentPositions.append(newAlignment)) {
            alignments.push_back(newAlignment);
            return;
        }

        // construct MD tag and check if the alignment has too many mismatch, if do, ignore it.
        if (newAlignment->repeat) {
            if (!newAlignment->constructRepeatMD(bitReference, alignmentPositions)) {
                alignmentPositions.badAligned();
                alignments.push_back(newAlignment);
                return;
            }
            nRepeatAlignment++; // for each repeat alignment, record it.
        } else {
            // check mismatch, update tags
            if (!newAlignment->constructMD(bitReference)) {
                alignmentPositions.badAligned();
                alignments.push_back(newAlignment);
                return;
            }
        }

        // update pair score or AS, for output using.
        // if the new alignment has lower paring score or AS than bestPairScore or bestAS, ignore it.
        if (paired) {
            if (!alignmentPositions.updatePairScore()) {
                alignments.push_back(newAlignment);
                return;
            }
            if (alignmentPositions.bestPairScore == maxPairScore && alignmentPositions.nBestPair > 1) {
                multipleAligned = true;
            }
        } else {
            if (!alignmentPositions.updateAS()) {
                alignments.push_back(newAlignment);
                return;
            }
            if (alignmentPositions.bestAS == 0 && alignmentPositions.nBestSingle > 1) {
                multipleAligned = true;
            }
        }
        alignments.push_back(newAlignment);
    }

    /**
     * if there is no alignment, output unAlignment result.
     * this function is important when hisat2 give mapped result, but it does not pass my filter (has too many mismatch).
     */
    void outputUnAlignmentRead(BTString& o) {
        if (paired) {
            for (int i = 0; i < 2; i++) {
                assert(!readName[i].empty());
                string flag = (i == 0) ? "77" : "141";
                o.append(readName[i].toZBuf());
                o.append("\t");
                o.append(flag.c_str());
                o.append("\t*\t0\t0\t*\t*\t0\t0\t");
                o.append(readSequence[i].toZBuf());
                o.append("\t");
                o.append(qualityScore[i].toZBuf());
                o.append("\tYT:Z:UP");
                o.append(passThroughLines[i].toZBuf());
                o.append('\n');
            }
        } else {
            assert(!readName[0].empty());
            o.append(readName[0].toZBuf());
            o.append("\t4\t*\t0\t0\t*\t*\t0\t0\t");
            o.append(readSequence[0].toZBuf());
            o.append("\t");
            o.append(qualityScore[0].toZBuf());
            o.append("\tYT:Z:UU");
            o.append(passThroughLines[0].toZBuf());
            o.append('\n');
        }
    }

    /**
     * report alignment statistics for single-end alignment
     */
    void reportStats_single(ReportingMetrics& met);


    /**
     * report alignment statistics for paired-end alignment
     */
    void reportStats_paired(ReportingMetrics& met);

    /**
     * output single-end alignment reuslts
     */
    void output_single(BTString& o,
                       ReportingMetrics& met) {

        reportStats_single(met);

        // output
        if (uniqueOutputOnly && (alignmentPositions.nBestSingle != 1 || multipleAligned)) {
            // do not output anything
        } else if (alignments.empty() || alignmentPositions.nBestSingle == 0) {
            // make a unalignment result and output it.
            outputUnAlignmentRead(o);
        } else {
            // output
            alignmentPositions.outputSingle(o);
        }
    }

    /**
     * output paired-end alignment reuslts
     */
    void output_paired(BTString& o,
                       ReportingMetrics& met) {

        reportStats_paired(met);

        if ((uniqueOutputOnly && (alignmentPositions.nBestPair != 1 || multipleAligned))) {
            // do not report anything
        } else if (alignments.empty() ||
                   alignmentPositions.nBestPair == 0 ||
                   alignmentPositions.bestPairScore == numeric_limits<int>::min()) {
            // make a unalignment result and output it.
            outputUnAlignmentRead(o);
        } else {
            // output
            alignmentPositions.outputPair(o);
        }
    }
    /**
     * output function will be redirected to output_single or output_paired
     */
    void output(ReportingMetrics& met,
                BTString& o) {

        if (paired) {
            output_paired(o, met);
        } else {
            output_single(o,met);
        }
        initialize();
    }
};

#endif //HISAT2_ALIGNMENT_3N_H
