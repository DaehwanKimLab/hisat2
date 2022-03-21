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
#include "position_3n.h"


/**
 * compare the MappingPosition to Alignment.
 * if they have same chromosome and location information, return true.
 */
bool MappingPosition::operator==(Alignment* o) {

    BTString* testChromosome;
    if (!o->repeat && o->pairToRepeat) {
        testChromosome = &o->pairToChromosome;
    } else {
        testChromosome = &o->chromosomeName;
    }

    if (locations[1] == NULL) {
        return (*locations[o->pairSegment] == o->location) &&
               (*chromosome == *testChromosome);
    } else {
        return (*locations[o->pairSegment] == o->location) &&
               (*locations[1-(o->pairSegment)] == o->pairToLocation) &&
               (*chromosome == *testChromosome);
    }
}

/**
 * constructor for non-repeat alignment, or original repeat alignment ( with chromosomeName = rep*).
 */
MappingPosition::MappingPosition(Alignment* newAlignment) {
    initialize();
    locations[newAlignment->pairSegment] = &newAlignment->location;
    locations[1-newAlignment->pairSegment] = &newAlignment->pairToLocation;
    segmentExist[newAlignment->pairSegment] = true;
    alignments[newAlignment->pairSegment] = newAlignment;
    if (!newAlignment->repeat && newAlignment->pairToRepeat) {
        chromosome = &newAlignment->pairToChromosome;
    } else {
        chromosome = &newAlignment->chromosomeName;
    }
    pairScore = numeric_limits<int>::min();
}

/**
 * constructor for expanded repeat alignment position.
 */
MappingPosition::MappingPosition (RepeatMappingPosition* repeat0, Alignment* newAlignment0, RepeatMappingPosition* repeat1=NULL, Alignment* newAlignment1=NULL) {
    initialize();
    locations[newAlignment0->pairSegment] = &repeat0->repeatLocation;
    chromosome = &repeat0->repeatChromosome;
    repeats[0] = repeat0;
    repeats[1] = repeat1;
    alignments[0] = newAlignment0;
    alignments[1] = newAlignment1;
    segmentExist[0] = true;
    if (alignments[1] != NULL) {
        locations[newAlignment1->pairSegment] = &repeat1->repeatLocation;
        segmentExist[1] = true;
    }
    AS = repeat0->AS;
    repeat = true;
}

/**
 *  return true if the is a MappingPosition has same information as input Alignment.
 *  first check the latest MappingPosition, if not same, search all MappingPositions.
 *  both mate and opposite mate position should be same to MappingPosition to return true.
 */
bool MappingPositions::positionExist (Alignment* newAlignment) {
    if (positions.empty()) {
        index = 0;
        return false;
    }

    if (positions[index] == newAlignment) {
        return positions[index].segmentExist[newAlignment->pairSegment];
    }

    int segment = newAlignment->pairSegment;
    long long int* targetLocations[2];
    targetLocations[segment] = &newAlignment->location;
    targetLocations[1-segment] = &newAlignment->pairToLocation;

    return findPosition (targetLocations,
                         (!newAlignment->repeat && newAlignment->pairToRepeat)?newAlignment->pairToChromosome:newAlignment->chromosomeName,
                         segment);
}

/**
 * append new Alignment to positions.
 * return true if the new Alignment successfully append.
 * return false if the new Alignment is exist or it's mate is bad aligned.
 */
bool MappingPositions::append(Alignment* newAlignment) {
    if (positionExist(newAlignment)) { // check if position is exist.
        return false;
    } else {
        int segment = newAlignment->pairSegment;
        if (!positions.empty() && positions[index] == newAlignment) {
            // check if current MappingPosition is same to new Alignment.
            positions[index].segmentExist[segment] = true;
            if (positions[index].badAlignment) {
                return false;
            }
            positions[index].alignments[segment] = newAlignment;
        } else {
            // add the new Alignment to positions.
            positions.emplace_back(newAlignment);
            index = positions.size()-1;
            if (oppositeAlignment != NULL) {
                positions[index].alignments[1-segment] = oppositeAlignment;
                positions[index].segmentExist[1-segment] = true;
            }
        }
        return true;
    }
}

/**
 * output paired-end alignment results.
 */
void MappingPositions::outputPair(BTString& o) {
    int outputCount = 0;
    bool primary = true; // for primary alignment flag.
    for (int i = 0; i < positions.size(); i++) {
        if (positions[i].pairScore == bestPairScore) {
            outputCount++;
            assert(positions[i].alignments[0] != NULL);
            assert(positions[i].alignments[1] != NULL);

            // change the NH tag
            positions[i].alignments[0]->updateNH(nBestPair);
            positions[i].alignments[1]->updateNH(nBestPair);

            // get concordant information and change the concordant flag.
            bool concordant;
            if (!positions[i].alignments[0]->mapped || !positions[i].alignments[1]->mapped) {
                concordant = false;
            } else {
                concordant = Alignment::isConcordant(*positions[i].locations[0],
                                                     positions[i].alignments[0]->forward,
                                                     positions[i].alignments[0]->readSequence.length(),
                                                     *positions[i].locations[1],
                                                     positions[i].alignments[1]->forward,
                                                     positions[i].alignments[1]->readSequence.length());
            }

            positions[i].alignments[0]->setConcordant(concordant);
            positions[i].alignments[1]->setConcordant(concordant);

            positions[i].alignments[0]->setMateMappingFlag(positions[i].alignments[1]->mapped ? positions[i].locations[1] : NULL);
            positions[i].alignments[1]->setMateMappingFlag(positions[i].alignments[0]->mapped ? positions[i].locations[0] : NULL);

            if (!positions[i].repeat) {
                // output regular alignment result

                // if both mate is outputted before, change the mate 1 status and output mate1.
                if (positions[i].alignments[0]->outputted && positions[i].alignments[1]->outputted) {
                    positions[i].alignments[1]->outputted = false;
                }
                // change YS tag.
                positions[i].alignments[0]->setYS(positions[i].alignments[1]);
                positions[i].alignments[1]->setYS(positions[i].alignments[0]);
                // output
                positions[i].alignments[0]->outputAlignment(o, NULL, positions[i].locations[1], primary);
                positions[i].alignments[1]->outputAlignment(o, NULL, positions[i].locations[0], primary);
            } else {
                //output repeat alignment result.

                // if both mate is outputted before, change the mate 1 status and output mate1.
                if (positions[i].repeats[0]->outputted && positions[i].repeats[1]->outputted) {
                    positions[i].repeats[1]->outputted = false;
                }
                // change YS tag.
                positions[i].alignments[0]->setYS(positions[i].repeats[1]);
                positions[i].alignments[1]->setYS(positions[i].repeats[0]);
                //output
                positions[i].alignments[0]->outputAlignment(o, positions[i].repeats[0], positions[i].locations[1], primary);
                positions[i].alignments[1]->outputAlignment(o, positions[i].repeats[1], positions[i].locations[0], primary);
            }
            primary = false; // after output the first pair for each read, set primary status to false;
        }
    }
    assert(outputCount == nBestPair);
}

/**
 * output single-end alignment results.
 */
void MappingPositions::outputSingle(BTString &o) {
    int outputCount = 0;
    bool primary = true; // for primary alignment flag.
    for (int i = 0; i < positions.size(); i++) {
        if (positions[i].AS == bestAS && !positions[i].badAlignment) {
            outputCount++;
            assert(positions[i].alignments[0] != NULL);
            // set NH tag
            positions[i].alignments[0]->updateNH(nBestSingle);

            if (!positions[i].repeat) { // output regular alignment result
                positions[i].alignments[0]->outputAlignment(o, NULL, NULL, primary);
            } else { // output repeat alignment result
                positions[i].alignments[0]->outputAlignment(o, positions[i].repeats[0], NULL, primary);
            }
            primary = false; // after output the first alignment for each read, set primary status to false;
        }
    }
    assert(outputCount == nBestSingle);
}

bool MappingPositions::updateAS_regular() {
    if (isBad()) { return false; }
    if (!positions[index].alignments[0]->mapped) { return true; }
    int AS = positions[index].alignments[0]->AS;
    if (AS > bestAS) {
        bestAS = AS;
        nBestSingle = 1;
    } else if (AS == bestAS) {
        nBestSingle++;
    } else {
        badAligned();
        return false;
    }
    positions[index].AS = AS;
    return true;
}

/**
 * if AS in repeatPosition is larger than bestAS, add it to positions and update BestAS.
 */
bool MappingPositions::updateAS_repeat() {
    if (isBad()) { return false; }
    Alignment* alignment = positions[index].alignments[0];
    RepeatMappingPosition* repeatPosition;
    badAligned(); // label this as bad alignment to avoid directly output.
    int AS;
    for (int i = 0; i < alignment->repeatPositions.size(); i++) {
        repeatPosition = &alignment->repeatPositions.positions[i];
        AS = (repeatPosition->flagInfoIndex == -1)?repeatPosition->AS : alignment->repeatPositions.positions[repeatPosition->flagInfoIndex].AS;
        if (AS >= bestAS) {
            positions.emplace_back(repeatPosition, alignment);
            if (AS > bestAS) {
                bestAS = AS;
                nBestSingle = 1;
            } else {
                nBestSingle++;
            }
        }
    }
    return true;
}

/**
 * redirect to updateAS_regular() or updateAS_repeat().
 */
bool MappingPositions::updateAS() {
    if (positions[index].alignments[0]->repeat) {
        return updateAS_repeat();
    } else {
        return updateAS_regular();
    }
}

/**
 * calculate the pairing score for regular (non-repeat) alignment.
 */
bool MappingPositions::updatePairScore_regular() {
    if (positions[index].alignments[0]->chromosomeName != positions[index].alignments[1]->chromosomeName) {
        badAligned();
        return false;
    }
    int nPair;
    int score;
    score = positions[index].alignments[0]->calculatePairScore(positions[index].alignments[1], nPair);
    if (score > bestPairScore) {
        bestPairScore = score;
        nBestPair = nPair;
        concordantExist = positions[index].alignments[0]->concordant;
    } else if (score == bestPairScore) {
        nBestPair += nPair;
    } else { // the newPair Score is less than bestPairScore, label it
        badAligned();
        return false;
    }
    positions[index].pairScore = score;
    return true;
}

/**
 * calculate the pairing score for repeat alignment
 * append to positions if the new pair has better (or equal) pairing score.
 */
bool MappingPositions::updatePairScore_repeat() {
    Alignment* alignments[2];
    alignments[0] = positions[index].alignments[0];
    alignments[1] = positions[index].alignments[1];
    if ((!alignments[0]->mapped || !alignments[1]->mapped) &&
        (bestPairScore >= (numeric_limits<int>::min()/2 - 1))) {
        badAligned();
        return false;
    }
    RepeatMappingPosition *repeatPosition0;
    RepeatMappingPosition *repeatPosition1;
    RepeatMappingPosition *repeatFlag0;
    RepeatMappingPosition *repeatFlag1;
    bool forward[2];
    forward[0] = alignments[0]->forward;
    forward[1] = alignments[1]->forward;
    bool DNA = alignments[0]->DNA;
    int score;
    bool concordant;
    for (int i = 0; i < alignments[0]->repeatPositions.size(); i++) {
        repeatPosition0 = &alignments[0]->repeatPositions.positions[i];
        repeatFlag0 = repeatPosition0->flagInfoIndex==-1 ? repeatPosition0 : &alignments[0]->repeatPositions.positions[repeatPosition0->flagInfoIndex];
        for (int j = 0; j < alignments[1]->repeatPositions.size(); j++) {
            repeatPosition1 = &alignments[1]->repeatPositions.positions[j];
            if (repeatPosition0->repeatChromosome == repeatPosition1->repeatChromosome) {
                repeatFlag1 = repeatPosition1->flagInfoIndex==-1 ? repeatPosition1 : &alignments[1]->repeatPositions.positions[repeatPosition1->flagInfoIndex];
                if (DNA) {
                    score = Alignment::calculatePairScore_DNA(repeatPosition0->repeatLocation,
                                                   repeatFlag0->AS,
                                                   forward[0],
                                                   alignments[0]->readSequence.length(),
                                                   repeatPosition1->repeatLocation,
                                                   repeatFlag1->AS,
                                                   forward[1],
                                                   alignments[1]->readSequence.length(),
                                                   concordant);
                } else {
                    score = Alignment::calculatePairScore_RNA(repeatPosition0->repeatLocation,
                                                   repeatFlag0->XM,
                                                   forward[0],
                                                   alignments[0]->readSequence.length(),
                                                   repeatPosition1->repeatLocation,
                                                   repeatFlag1->XM,
                                                   forward[1],
                                                   alignments[1]->readSequence.length(),
                                                   concordant);
                }
                if (score >= bestPairScore) {
                    positions.emplace_back(repeatPosition0, alignments[0], repeatPosition1, alignments[1]);
                    positions.back().pairScore = score;
                    if (score > bestPairScore) {
                        nBestPair = 1;
                        bestPairScore = score;
                        concordantExist = concordant;
                    } else {
                        nBestPair++;
                    }
                }
            }
        }
    }
    return true;
}

/**
 * calculate the pairing score,
 * if one of mate is repeat, calculate the pairing score by knn and append the pair has best pairing score to positions.
 */
bool MappingPositions::updatePairScore() {
    if (!mateExist()) { return true; }

    assert(positions[index].alignments[0] != NULL);
    assert(positions[index].alignments[1] != NULL);

    if (positions[index].alignments[0]->repeat || positions[index].alignments[1]->repeat) {
        return updatePairScore_repeat();
    } else {
        return updatePairScore_regular();
    }
}