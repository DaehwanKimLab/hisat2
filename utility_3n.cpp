/*
 * Copyright 2021, Yun (Leo) Zhang <imzhangyun@gmail.com>
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

#include <unistd.h>
#include "utility_3n.h"
#include "alphabet.h"

void getConversion(char usrInputFrom, char usrInputTo, char& convertFrom, char& convertTo) {
    if ((usrInputFrom == 'A' && usrInputTo == 'T') ||
        (usrInputFrom == 'A' && usrInputTo == 'C') ||
        (usrInputFrom == 'C' && usrInputTo == 'G') ||
        (usrInputFrom == 'C' && usrInputTo == 'T') ||
        (usrInputFrom == 'G' && usrInputTo == 'A') ||
        (usrInputFrom == 'T' && usrInputTo == 'G')) {
        convertFrom = usrInputFrom;
        convertTo = usrInputTo;
        return;
    }
    if ((usrInputFrom == 'C' && usrInputTo == 'A') ||
        (usrInputFrom == 'G' && usrInputTo == 'C') ||
        (usrInputFrom == 'T' && usrInputTo == 'A') ||
        (usrInputFrom == 'T' && usrInputTo == 'C') ||
        (usrInputFrom == 'A' && usrInputTo == 'G') ||
        (usrInputFrom == 'G' && usrInputTo == 'T')) {
        convertFrom = usrInputTo;
        convertTo = usrInputFrom;
        return;
    }
    cerr << "Un-identified --base-change type: " << usrInputFrom << "," << usrInputTo << endl;
    throw 1;
}

bool fileExist (string name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}

void ConvertMatrix3N:: convertMatrix() {
    restoreNormal();
    for (int i = 0; i < 4; i++) {
        char base = allBase[i];
        char lowerBase = allBaseLower[i];
        if (convertFrom == base) {
            asc2dna[base] = charToInt(convertTo);
            asc2dna[lowerBase] = charToInt(convertTo);
        } else if (complement(convertFrom) == base) {
            asc2dnacomp[base] = convertTo;
            asc2dnacomp[lowerBase] = convertTo;
            dnacomp[i] = charToInt(convertTo);
        }
    }
}

void ConvertMatrix3N::restoreNormal() {
    for (int i = 0; i < 4; i++) {
        char base = allBase[i];
        char lowerBase = allBaseLower[i];
        asc2dna[base] = charToInt(base);
        asc2dna[lowerBase] = charToInt(base);
        asc2dnacomp[base] = complement(base);
        asc2dnacomp[lowerBase] = complement(base);
        dnacomp[i] = charToInt(complement(base));
    }
}

