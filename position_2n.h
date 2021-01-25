/*
 * Copyright 2020, Yun (Leo) Zhang <imzhangyun@gmail.com>
 *
 * This file is part of HISAT-2N.
 *
 * HISAT-2N is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT-2N is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT-2N.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HISAT2_POSITION_2N_H
#define HISAT2_POSITION_2N_H

#include "sstring.h"
#include "alignment_2n.h"

class Alignment;
class RepeatMappingPosition;

/**
 * the data structure to store existing mapping positions.
 */
class MappingPosition {
public:
    long long int* locations[2] = {NULL};
    BTString* chromosome;
    int AS = numeric_limits<int>::min();
    int pairScore; // score to decide which mapping position should be output.
    bool segmentExist[2] = {false}; // indicate whether we have the segment alignment information.
    bool badAlignment = false; // if the alignment result
    bool repeat = false; // whether the mapping position is belong to a expanded repeat alignment
    Alignment* alignments[2] = {NULL};
    RepeatMappingPosition* repeats[2] = {NULL};

    void initialize() {
        for (int i = 0; i < 2; i++) {
            locations[i] = NULL;
            segmentExist[i] = false;
            alignments[i] = NULL;
            repeats[i] = NULL;
        }
        chromosome = NULL;
        AS = numeric_limits<int>::min();
        pairScore = numeric_limits<int>::min();
        badAlignment = false;
        repeat = false;
    }

    MappingPosition() {

    }

    /**
     * constructor for non-repeat alignment, or original repeat alignment( with chromosomeName = rep*).
     */
    MappingPosition (Alignment* newAlignment);


    /**
     * constructor for expanded repeat alignment position.
     */
    MappingPosition (RepeatMappingPosition* repeat0, Alignment* newAlignment0, RepeatMappingPosition* repeat1, Alignment* newAlignment1);

    /**
     * return true if the MappingPosition has same information as input Alignment
     */
    bool operator==(Alignment* o);
};

/**
 * this is the data structure to store all MappingPosition
 */
class MappingPositions {
public:
    vector<MappingPosition> positions;
    int bestPairScore; // the best pairing score, for paired-end alignment output
    int nBestPair; // the number of pair have bestPairScore, should equal to NH.
    int bestAS; // the best AS score, for single-end alignment output.
    int nBestSingle; // the number of alignment have bestAS, should equal to NH.
    int index; // the index number on positions. should always point to the last or current MappingPosition.
    Alignment* oppositeAlignment; // the temporary pointer point to the opposite mate's Alignment. use in append function.
    bool concordantExist; // whether concordant alignment is exist. use for paired-end output statistics.

    void initialize() {
        positions.clear();
        bestPairScore = numeric_limits<int>::min();
        nBestPair = 0;
        bestAS = numeric_limits<int>::min();
        nBestSingle = 0;
        index = -1;
        oppositeAlignment = NULL;
        concordantExist = false;
    }

    MappingPositions() {
        initialize();
    };

    /**
     * return number of MappingPosition in positions.
     */
    int size() {
        return positions.size();
    }

    /**
     * recursively search the positions to find whether there is a Mapping position has target information.
     * if the opposite mate is exist, we will save it's Alignment address to oppositeAlignment.
     */
    bool findPosition (long long int* inputLocations[2], BTString& chromosome, int& pairSegment) {
        oppositeAlignment = NULL;
        for (int i = 0; i < positions.size(); i++) {
            if (positions[i].locations[1-pairSegment] == NULL ||
                *(positions[i].locations[1-pairSegment]) == *inputLocations[1-pairSegment]) {
                if (!positions[i].badAlignment) {
                    oppositeAlignment = positions[i].alignments[1-pairSegment];
                }
                if (*positions[i].locations[pairSegment] == *inputLocations[pairSegment] &&
                    (*positions[i].chromosome == chromosome)) {
                    index = i;
                    return positions[i].segmentExist[pairSegment];
                }
            }
        }
        return false;
    }

    /**
     * set current MappingPosition to a bad alignment.
     */
    void badAligned() {
        positions[index].badAlignment = true;
    }

    /**
     * return true if current positions is a bad alignment.
     */
    bool isBad() {
        return positions[index].badAlignment;
    }

    /**
     * return if both segment is exist for current MappingPosition
     */
    bool mateExist() {
        return positions[index].segmentExist[0] && positions[index].segmentExist[1];
    }

    /**
     * calculate the pairing score,
     * if one of mate is repeat, calculate the pairing score by knn and append the pair has best pairing score to positions.
     */
    bool updatePairScore();

    /**
     * calculate the pairing score for regular (non-repeat) alignment.
     */
    bool updatePairScore_regular();

    /**
     * calculate the pairing score for repeat alignment
     * append to positions if the new pair has better (or equal) pairing score.
     */
    bool updatePairScore_repeat();

    /**
     * redirect to updateAS_regular() or updateAS_repeat().
     */
    bool updateAS();

    /**
     * if AS is larger than bestAS, update bestAS.
     */
    bool updateAS_regular();

    /**
     * if AS in repeatPosition is larger than bestAS, add it to positions and update BestAS.
     */
    bool updateAS_repeat();

    /**
     *  return true if the is a MappingPosition has same information as input Alignment.
     *  first check the latest MappingPosition, if not same, search all MappingPositions.
     *  both mate and opposite mate position should be same to MappingPosition to return true.
     */
    bool positionExist(Alignment* newAlignment);

    /**
     *  return true if the is a MappingPosition has same information as position.
     *  this function is to check repeat mapping position.
     *  without checking it's mate, if the repeat mapping location is exist, return true.
     */
    bool positionExist (BTString& chromosome, long long int& location, int& segment) {
        for (int i = 0; i < positions.size(); i++) {
            if ((*(positions[i].locations[segment]) == location) &&
                (*(positions[i].chromosome) == chromosome)) {
                return true;
            }
        }
        return false;
    }

    /**
     * output paired-end alignment results.
     */
    void outputPair(BTString& o);

    /**
     * output single-end alignment results.
     */
    void outputSingle(BTString& o);

    /**
     * append new Alignment to positions.
     * return true if the new Alignment successfully append.
     * return false if the new Alignment is exist or it's mate is bad aligned.
     */
    bool append(Alignment* newAlignment);
};

/**
 * the data structure to store repeat information.
 */
class RepeatMappingPosition: public MappingPosition {
public:
    long long int repeatLocation;
    BTString MD;
    int XM;
    int NM;
    int YS;
    int Yf;
    char YZ;
    BTString refSequence;
    BTString repeatChromosome;
    bool outputted = false;
    int flagInfoIndex = -1;

    RepeatMappingPosition() {};

    /**
     * constructor for new repeat information.
     */
    RepeatMappingPosition (long long int& inputLocation,
                           BTString& inputChromosome,
                           BTString &inputRefSequence,
                           int &inputAS,
                           BTString &inputMD,
                           int &inputXM,
                           int &inputNM,
                           int &inputTC,
                           char &repeatYZ) {
        repeatLocation = inputLocation;
        repeatChromosome = inputChromosome;
        refSequence = inputRefSequence;
        AS = inputAS;
        MD = inputMD;
        XM = inputXM;
        NM = inputNM;
        Yf = inputTC;
        YZ = repeatYZ;
        pairScore = numeric_limits<int>::min();
        flagInfoIndex = -1;
    }

    /**
     * constructor for the repeat which has same reference sequence.
     * we save the index for pattern RepeatMappingPosition, because they should have same information except location and chromosome.
     */
    RepeatMappingPosition(long long int &inputLocation,
                          BTString &inputChromosome,
                          int& inputAS,
                          int& index) {
        repeatLocation = inputLocation;
        repeatChromosome = inputChromosome;
        AS = inputAS;
        flagInfoIndex = index;
    }
};

/**
 * this is the data structure to store all repeatMappingPosition after expansion.
 */
class RepeatMappingPositions {
public:
    vector<RepeatMappingPosition> positions;

    void initialize() {
        positions.clear();
    }

    /**
     * return number of MappingPosition in positions.
     */
    int size() {
        return positions.size();
    }

    /**
     * return true if reference sequence is exist, else, return false.
     */
    bool sequenceExist (BTString& refSequence, int &index) {

        for (int i = 0; i < positions.size(); i++) {
            if ((positions[i].flagInfoIndex == -1) && (refSequence == positions[i].refSequence)) {
                index = i;
                return true;
            }
        }
        return false;
    }

    /**
     * add repeat mapping information.
     */
    void append (long long int &location,
                 BTString &chromosome,
                 BTString &refSequence,
                 int &AS,
                 BTString &MD,
                 int &XM,
                 int &NM,
                 int &Yf,
                 char &repeatYZ) {
        positions.emplace_back(location, chromosome, refSequence, AS, MD, XM, NM, Yf, repeatYZ);
    }

    /**
     * add repeat mapping information.
     */
    void append(BTString &chromosome, long long int &location, int &index) {
        positions.emplace_back(location, chromosome, positions[index].AS, index);
    }
};

#endif //HISAT2_POSITION_2N_H

