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

#ifndef UTILITY_3N_TABLE_H
#define UTILITY_3N_TABLE_H

#include <mutex>
#include <queue>

using namespace std;

/**
 * return complement of input base.
 */
char asc2dnacomp[] = {
        /*   0 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*  16 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*  32 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,'-',  0,  0,
        /*  48 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /*  64 */ 0,'T','V','G','H',  0,  0,'C','D',  0,  0,'M',  0,'K','N',  0,
        /*    A   B   C   D           G   H           K       M   N */
        /*  80 */ 0,  0,'Y','S','A',  0,'B','W',  0,'R',  0,  0,  0,  0,  0,  0,
        /*        R   S   T       V   W       Y */
        /*  96 */ 0,'T','V','G','H',  0,  0,'C','D',  0,  0,'M',  0,'K','N',  0,
        /*   a   b   c   d           g   h           k       m   n */
        /* 112 */ 0,  0,'Y','S','A',  0,'B','W',  0,'R',  0,  0,  0,  0,  0,  0,
        /*        r   s   t       v   w       y */
        /* 128 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 144 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 160 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 176 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 192 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 208 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 224 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        /* 240 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

/**
 * the simple data structure to bind quality score and position (on reference) together.
 */
class PosQuality {
public:
    int readPos; // 0-based
    int refPos;  // 0-based
    char qual;
    bool converted;
    bool remove;

    PosQuality(int& inputPos) {
        readPos = inputPos;
        refPos = inputPos;
        remove = true;
    }

    void setQual (char& inputQual, bool inputConverted) {
        qual = inputQual;
        converted = inputConverted;
        remove = false;
    }
};

/**
 * the base class for string we need to search.
 */
class string_search {
public:
    int start;
    string s;
    int stringLen;

    void initialize() {
        start = 0;
        stringLen = 0;
        s.clear();
    }

    void loadString(string intputString) {
        s = intputString;
        stringLen = s.size();
        start = 0;
    }
};


/**
 * to store CIGAR string and search segments in it.
 */
class CIGAR : public string_search{
public:

    bool getNextSegment(int& len, char& symbol) {
        if (start == stringLen) {
            return false;
        }
        len = 0;
        int currentIndex = start;
        while (true) {
            if (isalpha(s[currentIndex])) {
                len = stoi(s.substr(start, currentIndex-start));
                symbol = s[currentIndex];
                start = currentIndex+1;
                return true;
            }
            currentIndex++;
        }
    }
};

/**
 * to store MD tag and search segments in it.
 */
class MD_tag : public string_search {
public:

    bool getNextSegment(string& seg) {
        if (start == stringLen) {
            return false;
        }
        seg.clear();
        int currentIndex = start;
        bool deletion = false;

        while (true) {
            if (start >= stringLen) {
                return !seg.empty();
            }
            if (seg.empty() && s[currentIndex] == '0') {
                currentIndex++;
                continue;
            }
            if (isalpha(s[currentIndex])) {
                if (seg.empty()) {
                    seg = s[currentIndex];
                    start = currentIndex+1;
                    return true;
                } else {
                    if (deletion) {
                        seg += s[currentIndex];
                        //currentIndex++;
                    } else {
                        start = currentIndex;
                        return true;
                    }
                }
            } else if (s[currentIndex] == '^') {
                if (seg.empty()) {
                    seg = s[currentIndex];
                    deletion = true;
                } else {
                    start = currentIndex;
                    return true;
                }
            } else { // number
                if (seg.empty()) {
                    seg = s[currentIndex];
                } else {
                    if (deletion || isalpha(seg.back())) {
                        start = currentIndex;
                        return true;
                    } else {
                        seg += s[currentIndex];
                    }
                }
            }
            currentIndex++;
        }
    }
};

/**
 * simple safe queue
 */
template <typename T>
class SafeQueue {
private:
    mutex mutex_;
    queue<T> queue_;
public:
    void pop() {
        mutex_.lock();
        queue_.pop();
        mutex_.unlock();
    }

    T front() {
        mutex_.lock();
        T value = queue_.front();
        mutex_.unlock();
        return value;
    }

    int size() {
        mutex_.lock();
        int s = queue_.size();
        mutex_.unlock();
        return s;
    }

    /**
     * return true if the queue is not empty and pop front and get value.
     * return false if the queue is empty.
     */
    bool popFront(T& value) {
        mutex_.lock();
        bool isEmpty = queue_.empty();
        if (!isEmpty) {
            value = queue_.front();
            queue_.pop();
        }
        mutex_.unlock();
        return !isEmpty;
    }

    void push(T value) {
        mutex_.lock();
        queue_.push(value);
        mutex_.unlock();
    }

    bool empty() {
        mutex_.lock();
        bool check = queue_.empty();
        mutex_.unlock();
        return check;
    }
};

/**
 * store one chromosome and it's stream position
 */
class ChromosomeFilePosition {
public:
    string chromosome;
    streampos linePos;
    ChromosomeFilePosition(string inputChromosome, streampos inputPos) {
        chromosome = inputChromosome;
        linePos = inputPos;
    }
};

/**
 * store all chromosome and it's stream position
 */
class ChromosomeFilePositions {
public:
    vector <ChromosomeFilePosition> pos;
    streampos largestPos;

    /**
     * input chromosome name, return it's streamPosition.
     * if the chromosome is not in pos, return 0 (first position).
     */
    streampos getStreamPos(string &chromosome) {
        if (pos.empty()) {
             return 0;
        }
        int index = searchChromosome(chromosome, 0, pos.size()-1);
        if (pos[index].chromosome == chromosome) {
            return pos[index].linePos;
        } else {
            return largestPos;
        }
    }

    /**
     * input the chromosome name and it's streamPos, if it is not in pos, add it.
     */
    void append (string &chromosome, streampos& linePos) {
        if (linePos > largestPos) {
            largestPos = linePos;
        }
        if (pos.empty() || chromosome > pos.back().chromosome) {
            pos.push_back(ChromosomeFilePosition(chromosome, linePos));
            return;
        }
        int index = searchChromosome(chromosome, 0, pos.size()-1);
        if (pos[index].chromosome != chromosome) {
            pos.insert(pos.begin()+index, ChromosomeFilePosition(chromosome, linePos));
        }
    }

    /**
     * make binary search on pos for target chromosome name
     */
    int searchChromosome(string &targetChromosome, int start, int end) {
        if (pos.empty() || targetChromosome > pos.back().chromosome) {
            return 0;
        }
        if (start <= end) {
            int middle = (start + end) / 2;
            if (pos[middle].chromosome == targetChromosome) {
                return middle;
            }
            if (pos[middle].chromosome > targetChromosome) {
                return searchChromosome(targetChromosome, start, middle-1);
            }
            return searchChromosome(targetChromosome, middle+1, end);
        }
        return start; // return the bigger one
    }
};
#endif //UTILITY_3N_TABLE_H
