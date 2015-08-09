/*
 * Copyright 2015, Daehwan Kim <infphilo@gmail.com>
 *
 * This file is part of HISAT 2.
 *
 * HISAT 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include "gfm.h"

using namespace std;

#ifdef BOWTIE_64BIT_INDEX

const std::string gfm_ext("ht2l");

#else

const std::string gfm_ext("ht2");

#endif  // BOWTIE_64BIT_INDEX

string gLastIOErrMsg;

/**
 * Try to find the Bowtie index specified by the user.  First try the
 * exact path given by the user.  Then try the user-provided string
 * appended onto the path of the "indexes" subdirectory below this
 * executable, then try the provided string appended onto
 * "$HISAT2_INDEXES/".
 */
string adjustEbwtBase(const string& cmdline,
                      const string& gfmFileBase,
                      bool verbose)
{
    string str = gfmFileBase;
    ifstream in;
    if(verbose) cout << "Trying " << str.c_str() << endl;
    in.open((str + ".1." + gfm_ext).c_str(), ios_base::in | ios::binary);
    if(!in.is_open()) {
        if(verbose) cout << "  didn't work" << endl;
        in.close();
        if(getenv("HISAT2_INDEXES") != NULL) {
            str = string(getenv("HISAT2_INDEXES")) + "/" + gfmFileBase;
            if(verbose) cout << "Trying " << str.c_str() << endl;
            in.open((str + ".1." + gfm_ext).c_str(), ios_base::in | ios::binary);
            if(!in.is_open()) {
                if(verbose) cout << "  didn't work" << endl;
                in.close();
            } else {
                if(verbose) cout << "  worked" << endl;
            }
        }
    }
    if(!in.is_open()) {
        cerr << "Could not locate a HISAT2 index corresponding to basename \"" << gfmFileBase.c_str() << "\"" << endl;
        throw 1;
    }
    return str;
}
