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

#ifndef POSITION_3N_TABLE_H
#define POSITION_3N_TABLE_H

#include <string>
#include <vector>
#include <fstream>
#include <mutex>
#include <thread>
#include "alignment_3n_table.h"

using namespace std;

extern bool CG_only;
extern long long int loadingBlockSize;

/**
 * basic class to store reference position information
 */
class Position{
    mutex mutex_;
public:
    string chromosome; // reference chromosome name
    long long int location; // 1-based position
    char strand; // +(REF) or -(REF-RC)
    string convertedQualities; // each char is a mapping quality on this position for converted base.
    string unconvertedQualities; // each char is a mapping quality on this position for unconverted base.

    void initialize() {
        chromosome.clear();
        location = -1;
        strand = '?';
        convertedQualities.clear();
        unconvertedQualities.clear();
    }

    Position(){
        initialize();
    };

    /**
     * return true if there is mapping information in this reference position.
     */
    bool empty() {
        return convertedQualities.empty() && unconvertedQualities.empty();
    }

    /**
     * set the chromosome, location (position), and strand information.
     */

    void set (string& inputChr, long long int inputLoc) {
        chromosome = inputChr;
        location = inputLoc + 1;
    }

    void set(char inputStrand) {
        strand = inputStrand;
    }

    /**
     * append the SAM information into this position.
     */
    void appendBase (PosQuality& input) {
        mutex_.lock();
        if (input.converted) {
            convertedQualities += input.qual;
        } else {
            unconvertedQualities += input.qual;
        }
        mutex_.unlock();
    }
};

/**
 * store all reference position in this class.
 */
class Positions{
public:
    vector<Position*> refPositions; // the pool of all current reference position.
    string chromosome; // current reference chromosome name.
    long long int location; // current location (position) in reference chromosome.
    char lastBase = 'X'; // the last base of reference line. this is for CG_only mode.
    SafeQueue<string*> linePool; // pool to store unprocessed SAM line.
    SafeQueue<string*> freeLinePool; // pool to store free string pointer for SAM line.
    SafeQueue<Position*> freePositionPool; // pool to store free position pointer for reference position.
    SafeQueue<Position*> outputPositionPool; // pool to store the reference position which is loaded and ready to output.
    bool working;
    mutex mutex_;
    long long int refCoveredPosition; // this is the last position in reference chromosome we loaded in refPositions.
    ifstream refFile;
    vector<mutex*> workerLock; // one lock for one worker thread.
    int nThreads = 1;
    ChromosomeFilePositions chromosomePos; // store the chromosome name and it's streamPos. To quickly find new chromosome in file.

    Positions(string inputRefFileName, int inputNThreads) {
        working = true;
        nThreads = inputNThreads;
        for (int i = 0; i < nThreads; i++) {
            workerLock.push_back(new mutex);
        }
        refFile.open(inputRefFileName, ios_base::in);
    }

    ~Positions() {
        for (int i = 0; i < workerLock.size(); i++) {
            delete workerLock[i];
        }
        Position* pos;
        while(freePositionPool.popFront(pos)) {
            delete pos;
        }
    }

    /**
     * given the target Position output the corresponding position index in refPositions.
     */
    int getIndex(long long int &targetPos) {
        int firstPos = refPositions[0]->location;
        return targetPos - firstPos;
    }

    /**
     * given reference line (start with '>'), extract the chromosome information.
     * this is important when there is space in chromosome name. the SAM information only contain the first word.
     */
    string getChrName(string& inputLine) {
        size_t endPosition = inputLine.find(' ', 0);
        return inputLine.substr(1, endPosition-1);
    }

    /**
     * get a fasta line (not header), append the bases to positions.
     */
    void appendRefPosition(string& line) {
        Position* newPos;
        // check the base one by one
        char* b;
        for (int i = 0; i < line.size(); i++) {
            getFreePosition(newPos);
            newPos->set(chromosome, location+i);
            b = &line[i];
            if (CG_only) {
                if (i == 0 && lastBase == 'C' && *b == 'G') {
                    refPositions.back()->set('+');
                    newPos->set('-');
                }
            } else {
                if (*b == convertFrom) {
                    newPos->set('+');
                } else if (*b == convertFromComplement) {
                    newPos->set('-');
                }
            }
            refPositions.push_back(newPos);
            lastBase = *b;
        }
        location += line.size();
    }

    /**
     * if we can go through all the workerLock, that means no worker is appending new position.
     */
    void appendingFinished() {
        for (int i = 0; i < nThreads; i++) {
            workerLock[i]->lock();
            workerLock[i]->unlock();
        }
    }

    /**
     * the output function for output thread.
     */
    void outputFunction(string outputFileName) {
        ofstream tableFile;
        tableFile.open(outputFileName, ios_base::out);
        tableFile << "ref\tpos\tstrand\tconvertedBaseQualities\tconvertedBaseCount\tunconvertedBaseQualities\tunconvertedBaseCount\n";
        Position* pos;
        while (working) {
            if (outputPositionPool.popFront(pos)) {
                tableFile << pos->chromosome << '\t'
                          << to_string(pos->location) << '\t'
                          << pos->strand << '\t'
                          << pos->convertedQualities << '\t'
                          << to_string(pos->convertedQualities.size()) << '\t'
                          << pos->unconvertedQualities << '\t'
                          << to_string(pos->unconvertedQualities.size()) << '\n';
                returnPosition(pos);
            } else {
                this_thread::sleep_for (std::chrono::microseconds(1));
            }
        }
        tableFile.close();
    }

    /**
     * move the position which position smaller than refCoveredPosition - loadingBlockSize, output it.
     */
    void moveBlockToOutput() {
        if (refPositions.empty()) {
            return;
        }
        int index;
        for (index = 0; index < refPositions.size(); index++) {
            if (refPositions[index]->location < refCoveredPosition - loadingBlockSize) {
                if (refPositions[index]->empty() || refPositions[index]->strand == '?') {
                    returnPosition(refPositions[index]);
                } else {
                    outputPositionPool.push(refPositions[index]);
                }
            } else {
                break;
            }
        }
        if (index != 0) {
            refPositions.erase(refPositions.begin(), refPositions.begin()+index);
        }
    }

    /**
     * move all the refPosition into output pool.
     */
    void moveAllToOutput() {
        if (refPositions.empty()) {
            return;
        }
        for (int index = 0; index < refPositions.size(); index++) {
            if (refPositions[index]->empty() || refPositions[index]->strand == '?') {
                returnPosition(refPositions[index]);
            } else {
                outputPositionPool.push(refPositions[index]);
            }
        }
        refPositions.clear();
    }

    /**
     * initially load reference sequence for 2 million bp
     */
    void loadNewChromosome(string targetChromosome) {
        refFile.clear();
        // find the start position in file based on chromosome name.
        streampos startPos = chromosomePos.getStreamPos(targetChromosome);
        refFile.seekg(startPos, ios::beg);
        refCoveredPosition = 2 * loadingBlockSize;
        string line;
        bool load = false;
        while (refFile.good()) {
            getline(refFile, line);
            if (line.front() == '>') { // this line is chromosome name
                if (load) { // meet next chromosome, return it.
                    return;
                }
                chromosome = getChrName(line);
                if (chromosome == targetChromosome) {
                    load = true;
                } else {
                    streampos currentPos = refFile.tellg();
                    chromosomePos.append(chromosome, currentPos);
                }
                lastBase = 'X';
                location = 0;
            } else {
                if (!load) {
                    continue;
                }
                if (line.empty()) { continue; }
                // change all base to upper case
                for (int i = 0; i < line.size(); i++) {
                    line[i] = toupper(line[i]);
                }
                appendRefPosition(line);
                if (location >= refCoveredPosition) {
                    return;
                }
            }
        }
    }

    /**
     * load more Position (loadingBlockSize bp) to positions
     * if we meet next chromosome, return false. Else, return ture.
     */
    void loadMore() {
        refCoveredPosition += loadingBlockSize;
        string line;
        while (refFile.good()) {
            getline(refFile, line);
            if (line.front() == '>') { // meet next chromosome, return.
                return ;
            } else {
                if (line.empty()) { continue; }

                // change all base to upper case
                for (int i = 0; i < line.size(); i++) {
                    line[i] = toupper(line[i]);
                }

                appendRefPosition(line);
                if (location >= refCoveredPosition) {
                    return ;
                }
            }
        }
    }

    /**
     * add position information from Alignment into ref position.
     */
    void appendPositions(Alignment& newAlignment) {
        if (!newAlignment.mapped || newAlignment.bases.empty()) {
            return;
        }
        long long int startPos = newAlignment.location; // 1-based position
        // find the first reference position in pool.
        int index = getIndex(newAlignment.location);

        for (int i = 0; i < newAlignment.sequence.size(); i++) {
            PosQuality* b = &newAlignment.bases[i];
            if (b->remove) {
                continue;
            }

            Position* pos = refPositions[index+b->refPos];
            assert (pos->location == startPos + b->refPos);

            if (pos->strand == '?') {
                // this is for CG-only mode. read has a 'C' or 'G' but not 'CG'.
                continue;
            }
            pos->appendBase(newAlignment.bases[i]);
        }
    }

    /**
     * get a string pointer from freeLinePool, if freeLinePool is empty, make a new string pointer.
     */
    void getFreeStringPointer(string*& newLine) {
        if (freeLinePool.popFront(newLine)) {
            return;
        } else {
            newLine = new string();
        }
    }

    /**
     * get a Position pointer from freePositionPool, if freePositionPool is empty, make a new Position pointer.
     */
    void getFreePosition(Position*& newPosition) {
        if (freePositionPool.popFront(newPosition)) {
            return;
        } else {
            newPosition = new Position();
        }
    }

    /**
     * return the line to freeLinePool
     */
    void returnLine(string* line) {
        line->clear();
        freeLinePool.push(line);
    }

    /**
     * return the position to freePositionPool.
     */
    void returnPosition(Position* pos) {
        pos->initialize();
        freePositionPool.push(pos);
    }

    /**
     * this is the working function.
     * it take the SAM line from linePool, parse it.
     */
    void append(int threadID) {
        string* line;
        Alignment newAlignment;

        while (working) {
            workerLock[threadID]->lock();
            if(!linePool.popFront(line)) {
                workerLock[threadID]->unlock();
                this_thread::sleep_for (std::chrono::nanoseconds(1));
                continue;
            }
            while (refPositions.empty()) {
                this_thread::sleep_for (std::chrono::microseconds(1));
            }
            newAlignment.parse(line);
            returnLine(line);
            appendPositions(newAlignment);
            workerLock[threadID]->unlock();
        }
    }
};

#endif //POSITION_3N_TABLE_H
