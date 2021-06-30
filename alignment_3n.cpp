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

#include "alignment_3n.h"
#include "aln_sink.h"

/**
 * return true if two location is concordant.
 * return false, if there are not concordant or too far (>maxPairDistance).
 */
bool Alignment::isConcordant(long long int &location1, bool &forward1, long long int &location2, bool &forward2) {
    if (abs(location1-location2) > maxPairDistance) { return false; }
    if (location1 < location2) {
        if (forward1 && !forward2) { return true; }
    } else {
        if (!forward1 && forward2) { return true; }
    }
    return false;
}

/**
 * this is the basic function to calculate DNA pair score.
 * if the distance between 2 alignments is more than penaltyFreeDistance_DNA, we reduce the score by the distance/100.
 * if two alignment is concordant we add concordantScoreBounce to make sure to select the concordant pair as best pair.
 */
int Alignment::calculatePairScore_DNA (long long int &location0, int& AS0, bool& forward0, long long int &location1, int &AS1, bool &forward1, bool& concordant) {

    int score = ASPenalty*AS0 + ASPenalty*AS1;
    int distance = abs(location0 - location1);
    if (distance > maxPairDistance) { return numeric_limits<int>::min(); }
    if (distance > penaltyFreeDistance_DNA) { score -= distance/distancePenaltyFraction_DNA; }
    concordant = isConcordant(location0, forward0, location1, forward1);
    if (concordant) { score += concordantScoreBounce; }
    return score;
}

/**
 * this is the basic function to calculate RNA pair score.
 * if the distance between 2 alignments is more than penaltyFreeDistance_RNA, we reduce the score by the distance/1000.
 * if two alignment is concordant we add concordantScoreBounce to make sure to select the concordant pair as best pair.
 */
int Alignment::calculatePairScore_RNA (long long int &location0, int& XM0, bool& forward0, long long int &location1, int &XM1, bool &forward1, bool& concordant) {
    // this is the basic function to calculate pair score.
    // if the distance between 2 alignment is more than 100,000, we reduce the score by the distance/1000.
    // if two alignment is concordant we add 500,000 to make sure to select the concordant pair as best pair.
    int score = -ASPenalty*XM0 + -ASPenalty*XM1;
    int distance = abs(location0 - location1);
    if (distance > maxPairDistance) { return numeric_limits<int>::min(); }
    if (distance > penaltyFreeDistance_RNA) { score -= distance/distancePenaltyFraction_RNA; }
    concordant = isConcordant(location0, forward0, location1, forward1);
    if (concordant) { score += concordantScoreBounce; }
    return score;
}

/**
 * calculate the pairScore for a pair of alignment result. Output pair Score and number of pair.
 * Do not update their pairScore.
 */
int Alignment::calculatePairScore(Alignment *inputAlignment, int &nPair) {
    int pairScore = numeric_limits<int>::min();
    nPair = 0;
    if (pairSegment == inputAlignment->pairSegment){
        // when 2 alignment results are from same pair segment, output the lowest score and number of pair equal zero.
        pairScore = numeric_limits<int>::min();
    } else if (!mapped && !inputAlignment->mapped) {
        // both unmapped.
        pairScore = numeric_limits<int>::min()/2 - 1;
    } else if (!mapped || !inputAlignment->mapped) {
        // one of the segment unmapped.
        pairScore = numeric_limits<int>::min()/2;
        nPair = 1;
    } else if ((!repeat && !inputAlignment->repeat)){
        // both mapped and (both non-repeat or not expand repeat)
        bool concordant;
        if (DNA) {
            pairScore = calculatePairScore_DNA(location,
                                               AS,
                                               forward,
                                               inputAlignment->location,
                                               inputAlignment->AS,
                                               inputAlignment->forward,
                                               concordant);
        } else {
            pairScore = calculatePairScore_RNA(location,
                                               XM,
                                               forward,
                                               inputAlignment->location,
                                               inputAlignment->XM,
                                               inputAlignment->forward,
                                               concordant);
        }
        setConcordant(concordant);
        inputAlignment->setConcordant(concordant);
        nPair = 1;
    }
    return pairScore;
}

void Alignments::reportStats_single(ReportingMetrics& met) {

    int nAlignment = alignmentPositions.nBestSingle;
    if (nAlignment == 0) {
        met.nunp_0++;
    } else {
        met.nunp_uni++;
        if (nAlignment == 1) { met.nunp_uni1++; }
        else { met.nunp_uni2++; }
    }
}

void Alignments::reportStats_paired(ReportingMetrics& met) {
    if (!alignmentPositions.concordantExist) {
        met.nconcord_0++;
        if (alignmentPositions.nBestPair == 0) {
            met.nunp_0_0 += 2;
            return;
        }
        if (alignmentPositions.bestPairScore == numeric_limits<int>::min()/2) {
            // one mate is unmapped, one mate is mapped
            met.nunp_0_0++;
            met.nunp_0_uni++;
            if (alignmentPositions.nBestPair == 1) { met.nunp_0_uni1++; }
            else { met.nunp_0_uni2++; }
        } else { //both mate is mapped
            if (alignmentPositions.nBestPair == 1) {
                met.ndiscord++;
                return;
            }
            else {
                met.nunp_0_uni += 2;
                met.nunp_0_uni2 += 2;
            }
        }
    } else {
        assert(alignmentPositions.nBestPair > 0);
        met.nconcord_uni++;
        if (alignmentPositions.nBestPair == 1) { met.nconcord_uni1++; }
        else { met.nconcord_uni2++; }
    }
}

bool getConversion(char usrInputFrom, char usrInputTo, char& convertFrom, char& convertTo) {
    if ((usrInputFrom == 'A' && usrInputTo == 'T') ||
        (usrInputFrom == 'A' && usrInputTo == 'C') ||
        (usrInputFrom == 'C' && usrInputTo == 'G') ||
        (usrInputFrom == 'C' && usrInputTo == 'T')) {
        convertFrom = usrInputFrom;
        convertTo = usrInputTo;
        return false;
    }
    if ((usrInputFrom == 'A' && usrInputTo == 'G') ||
        (usrInputFrom == 'G' && usrInputTo == 'T')) {
        swap(usrInputFrom, usrInputTo);
        convertFrom = usrInputFrom;
        convertTo = usrInputTo;
        return true;
    }
    if ((usrInputFrom == 'C' && usrInputTo == 'A') ||
        (usrInputFrom == 'G' && usrInputTo == 'C') ||
        (usrInputFrom == 'T' && usrInputTo == 'A') ||
        (usrInputFrom == 'T' && usrInputTo == 'C')) {
        swap(usrInputFrom, usrInputTo);
        convertFrom = usrInputFrom;
        convertTo = usrInputTo;
        return false;
    }
    if ((usrInputFrom == 'G' && usrInputTo == 'A') ||
        (usrInputFrom == 'T' && usrInputTo == 'G')){
        convertFrom = usrInputFrom;
        convertTo = usrInputTo;
        return true;
    }
}

