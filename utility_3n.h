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

#ifndef HISAT2_UTILITY_3N_H
#define HISAT2_UTILITY_3N_H

#include <string>


using namespace std;

/**
 * this is the class to convert asc2dna, asc2dnacomp, and dnacomp matrix for hisat-3n conversion.
 * always save the conversion as convertFrom and convertTo.
 * this class is to convert matrix for hisat2-build with (--3N) or hisat2-repeat with (--3N).
 * hisat2 with (--base-change C,T) will not use this class.
 */
class ConvertMatrix3N {
    char convertFrom = 'A';
    char convertTo = 'A';
    string allBase = "ACGT";
    string allBaseLower = "acgt";

    /**
     * helper function to convert character to int presentation.
     */
    int charToInt(char inputChar) { return allBase.find(inputChar); }

    /**
     * return the complement nucleotide
     */
    int complement(char inputChar) { return allBase[3-charToInt(inputChar)]; }

    /**
     *  convert asc2dna, asc2dnacomp, and dnacomp matrix according to convertFrom and convertTo variable.
     *  always restore the matrix to original (4 letter) matrix.
     */
    void convertMatrix() ;
public:
    ConvertMatrix3N(){

    };

    /**
     * save the conversion information then convert the matrix.
     */
    void convert(char from, char to)  {
        convertFrom = from;
        convertTo = to;
        convertMatrix();
    }

    /**
     * change convertFrom and convertTO to it's complement nucleotide, then convert the matrix.
     */
    void inverseConversion() {
        convertFrom = complement(convertFrom);
        convertTo = complement(convertTo);
        convertMatrix();
    }

    /**
     * restore the asc2dna, asc2dnacomp, and dnacomp matrix to original (4 letter).
     * do not change the convertFrom and convertTO variable.
     */
    void restoreNormal() ;

    /**
     * convert the matrix according to the convertFrom and convertTo variable.
     */
    void restoreConversion() {
        convertMatrix();
    }
};

/**
 *  the simple data structure to store cigar information.
 */
class Cigar {
    int len;
    char label;
public:
    Cigar() { }

    Cigar(int inputLen, char inputLabel): len(inputLen), label(inputLabel) {
    }

    int& getLen() { return len; }

    char& getLabel() { return label; }
};

extern void getConversion(char usrInputFrom, char usrInputTo, char& convertFrom, char& convertTo);
extern bool fileExist (string name);


#endif //HISAT2_UTILITY_3N_H
