/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "bt2_idx.h"

using namespace std;

#ifdef HISAT_CLASS

const std::string gEbwt_ext("hc");

#else

#  ifdef BOWTIE_64BIT_INDEX

const std::string gEbwt_ext("bt2l");

#  else

const std::string gEbwt_ext("bt2");

#  endif  // BOWTIE_64BIT_INDEX

#endif

/**
 * Try to find the Bowtie index specified by the user.  First try the
 * exact path given by the user.  Then try the user-provided string
 * appended onto the path of the "indexes" subdirectory below this
 * executable, then try the provided string appended onto
 * "$BOWTIE2_INDEXES/".
 */
string adjustEbwtBase(const string& cmdline,
					  const string& ebwtFileBase,
					  bool verbose = false)
{
	string str = ebwtFileBase;
	ifstream in;
	if(verbose) cout << "Trying " << str.c_str() << endl;
	in.open((str + ".1." + gEbwt_ext).c_str(), ios_base::in | ios::binary);
	if(!in.is_open()) {
		if(verbose) cout << "  didn't work" << endl;
		in.close();
		if(getenv("BOWTIE2_INDEXES") != NULL) {
			str = string(getenv("BOWTIE2_INDEXES")) + "/" + ebwtFileBase;
			if(verbose) cout << "Trying " << str.c_str() << endl;
			in.open((str + ".1." + gEbwt_ext).c_str(), ios_base::in | ios::binary);
			if(!in.is_open()) {
				if(verbose) cout << "  didn't work" << endl;
				in.close();
			} else {
				if(verbose) cout << "  worked" << endl;
			}
		}
	}
	if(!in.is_open()) {
		cerr << "Could not locate a Bowtie index corresponding to basename \"" << ebwtFileBase.c_str() << "\"" << endl;
		throw 1;
	}
	return str;
}

string gLastIOErrMsg;
