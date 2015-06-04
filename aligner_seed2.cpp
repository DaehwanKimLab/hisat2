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

#include <limits>
#include <ctype.h>
#include "aligner_seed2.h"
#include "assert_helpers.h"
#include "gfm.h"

#ifdef ALIGNER_SEED2_MAIN

#include <string>
#include "sstring.h"

using namespace std;

/**
 * A way of feeding simply tests to the seed alignment infrastructure.
 */
int main(int argc, char **argv) {

    EList<string> strs;
    //                            GCTATATAGCGCGCTCGCATCATTTTGTGT
    strs.push_back(string("CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAACCA"
                          "NNNNNNNNNN"
                          "CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAACCA"));
    //                            GCTATATAGCGCGCTTGCATCATTTTGTGT
    //                                           ^
    bool packed = false;
    int color = 0;
	pair<GFM*, GFM*> gfms = GFM::fromStrings<SString<char> >(
		strs,
		packed,
		REF_READ_REVERSE,
		Ebwt::default_bigEndian,
		Ebwt::default_lineRate,
		Ebwt::default_offRate,
		Ebwt::default_ftabChars,
		".aligner_seed2.cpp.tmp",
		Ebwt::default_useBlockwise,
		Ebwt::default_bmax,
		Ebwt::default_bmaxMultSqrt,
		Ebwt::default_bmaxDivN,
		Ebwt::default_dcv,
		Ebwt::default_seed,
		false,  // verbose
		false,  // autoMem
		false); // sanity
    
    gfms.first->loadIntoMemory (-1, true, true, true, true, false);
    gfms.second->loadIntoMemory(1, true, true, true, true, false);
	
	int testnum = 0;

	// Query is longer than ftab and matches exactly twice
    for(int rc = 0; rc < 2; rc++) {
		for(int i = 0; i < 2; i++) {
			cerr << "Test " << (++testnum) << endl;
			cerr << "  Query with length greater than ftab" << endl;
			DescentMetrics mets;
			PerReadMetrics prm;
			DescentDriver dr;
			
			// Set up the read
			BTDnaString seq ("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
			BTString    qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
			if(rc) {
				seq.reverseComp();
				qual.reverse();
			}
			dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -30, 30);

			// Set up the DescentConfig
			DescentConfig conf;
			conf.cons.init(GFM::default_ftabChars, 1.0);
			conf.expol = DESC_EX_NONE;
			
			// Set up the search roots
			dr.addRoot(
				conf,   // DescentConfig
				(i == 0) ? 0 : (seq.length() - 1), // 5' offset into read of root
				(i == 0) ? true : false,           // left-to-right?
				rc == 0,   // forward?
				0.0f);   // root priority
			
			// Do the search
			Scoring sc = Scoring::base1();
			dr.go(sc, *gfms.first, *gfms.second, mets, prm);
			
			// Confirm that an exact-matching alignment was found
			assert_eq(1, dr.sink().nrange());
			assert_eq(2, dr.sink().nelt());
		}
	}
	
	// Query has length euqal to ftab and matches exactly twice
    for(int i = 0; i < 2; i++) {
		cerr << "Test " << (++testnum) << endl;
		cerr << "  Query with length equal to ftab" << endl;
        DescentMetrics mets;
		PerReadMetrics prm;
        DescentDriver dr;
        
        // Set up the read
        BTDnaString seq ("GCTATATAGC", true);
        BTString    qual("ABCDEFGHIa");
		dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -30, 30);
        
        // Set up the DescentConfig
        DescentConfig conf;
        conf.cons.init(GFM::default_ftabChars, 1.0);
        conf.expol = DESC_EX_NONE;
        
        // Set up the search roots
        dr.addRoot(
            conf,   // DescentConfig
            (i == 0) ? 0 : (seq.length() - 1), // 5' offset into read of root
            (i == 0) ? true : false,           // left-to-right?
            true,   // forward?
            0.0f);   // root priority
        
        // Do the search
        Scoring sc = Scoring::base1();
        dr.go(sc, *gfms.first, *gfms.second, mets, prm);
		
		// Confirm that an exact-matching alignment was found
		assert_eq(1, dr.sink().nrange());
		assert_eq(2, dr.sink().nelt());
    }

	// Query has length less than ftab length and matches exactly twice
    for(int i = 0; i < 2; i++) {
		cerr << "Test " << (++testnum) << endl;
		cerr << "  Query with length less than ftab" << endl;
        DescentMetrics mets;
		PerReadMetrics prm;
        DescentDriver dr;
        
        // Set up the read
        BTDnaString seq ("GCTATATAG", true);
        BTString    qual("ABCDEFGHI");
		dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -30, 30);
        
        // Set up the DescentConfig
        DescentConfig conf;
        conf.cons.init(GFM::default_ftabChars, 1.0);
        conf.expol = DESC_EX_NONE;
        
        // Set up the search roots
        dr.addRoot(
            conf,   // DescentConfig
            (i == 0) ? 0 : (seq.length() - 1), // 5' offset into read of root
            (i == 0) ? true : false,           // left-to-right?
            true,   // forward?
            0.0f);   // root priority
        
        // Do the search
        Scoring sc = Scoring::base1();
        dr.go(sc, *gfms.first, *gfms.second, mets, prm);
		
		// Confirm that an exact-matching alignment was found
		assert_eq(1, dr.sink().nrange());
		assert_eq(2, dr.sink().nelt());
    }
	
	// Search root is in the middle of the read, requiring a bounce
    for(int i = 0; i < 2; i++) {
		cerr << "Test " << (++testnum) << endl;
		cerr << "  Search root in middle of read" << endl;
        DescentMetrics mets;
		PerReadMetrics prm;
        DescentDriver dr;
        
        // Set up the read
		//                012345678901234567890123456789
        BTDnaString seq ("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
        BTString    qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
		TIndexOffU top, bot;
		top = bot = 0;
		bool ret = gfms.first->contains("GCGCTCGCATCATTTTGTGT", &top, &bot);
		cerr << ret << ", " << top << ", " << bot << endl;
		dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -30, 30);
        
        // Set up the DescentConfig
        DescentConfig conf;
        conf.cons.init(GFM::default_ftabChars, 1.0);
        conf.expol = DESC_EX_NONE;
        
        // Set up the search roots
        dr.addRoot(
            conf,   // DescentConfig
            (i == 0) ? 10 : (seq.length() - 1 - 10), // 5' offset into read of root
            (i == 0) ? true : false,                 // left-to-right?
            true,   // forward?
            0.0f);   // root priority
        
        // Do the search
        Scoring sc = Scoring::base1();
        dr.go(sc, *gfms.first, *gfms.second, mets, prm);
		
		// Confirm that an exact-matching alignment was found
		assert_eq(1, dr.sink().nrange());
		assert_eq(2, dr.sink().nelt());
    }

	delete gfms.first;
	delete gfms.second;
	
	strs.clear();
    strs.push_back(string("CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAACCA"
                          "NNNNNNNNNN"
                          "CATGTCAGCTATATAGCG"));
	gfms = GFM::fromStrings<SString<char> >(
		strs,
		packed,
		REF_READ_REVERSE,
		GFM::default_bigEndian,
		GFM::default_lineRate,
		GFM::default_offRate,
		GFM::default_ftabChars,
		".aligner_seed2.cpp.tmp",
		GFM::default_useBlockwise,
		GFM::default_bmax,
		GfM::default_bmaxMultSqrt,
		GFM::default_bmaxDivN,
		GFM::default_dcv,
		GFM::default_seed,
		false,  // verbose
		false,  // autoMem
		false); // sanity
    
    gfms.first->loadIntoMemory (-1, true, true, true, true, false);
    gfms.second->loadIntoMemory(1, true, true, true, true, false);
	
	// Query is longer than ftab and matches exactly once.  One search root for
	// forward read.
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for(int i = 0; i < 2; i++) {
			BTDnaString seq ("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
			BTString    qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
			for(size_t j = 0; j < seq.length(); j++) {
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query with length greater than ftab and matches exactly once" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				
				// Set up the read
				dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -30, 30);
				
				// Set up the DescentConfig
				DescentConfig conf;
				conf.cons.init(GFM::default_ftabChars, 1.0);
				conf.expol = DESC_EX_NONE;
				
				// Set up the search roots
				dr.addRoot(
					conf,   // DescentConfig
					j,      // 5' offset into read of root
					i == 0, // left-to-right?
					true,   // forward?
					0.0f);   // root priority
				
				// Do the search
				Scoring sc = Scoring::base1();
				dr.go(sc, *gfms.first, *gfms.second, mets, prm);
				
				// Confirm that an exact-matching alignment was found
				assert_eq(1, dr.sink().nrange());
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				assert_eq(1, dr.sink().nelt());
			}
		}
	}

	// Query is longer than ftab and its reverse complement matches exactly
	// once.  Search roots on forward and reverse-comp reads.
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for(int i = 0; i < 2; i++) {
			BTDnaString seq ("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
			BTString    qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
			for(size_t j = 0; j < seq.length(); j++) {
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query with length greater than ftab and reverse complement matches exactly once" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				
				// Set up the read
				dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -30, 30);
				
				// Set up the DescentConfig
				DescentConfig conf;
				conf.cons.init(GFM::default_ftabChars, 1.0);
				conf.expol = DESC_EX_NONE;
				
				// Set up the search roots
				dr.addRoot(
					conf,   // DescentConfig
					j,      // 5' offset into read of root
					i == 0, // left-to-right?
					true,   // forward?
					0.0f);   // root priority
				dr.addRoot(
					conf,   // DescentConfig
					j,      // 5' offset into read of root
					i == 0, // left-to-right?
					false,  // forward?
					1.0f);   // root priority
				
				// Do the search
				Scoring sc = Scoring::base1();
				dr.go(sc, *gfms.first, *gfms.second, mets, prm);
				
				// Confirm that an exact-matching alignment was found
				assert_eq(1, dr.sink().nrange());
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				assert_eq(1, dr.sink().nelt());
			}
		}
	}

	// Query is longer than ftab and matches exactly once with one mismatch
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for(int i = 0; i < 2; i++) {
			// Set up the read
			//    Ref: CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAACCA
			//                ||||||||||||||||||||||||||||||
			BTDnaString orig("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
			//                012345678901234567890123456789
			BTString    qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
			for(size_t k = 0; k < orig.length(); k++) {
				BTDnaString seq = orig;
				seq.set(seq[k] ^ 3, k);
				for(size_t j = 0; j < seq.length(); j++) {
					// Assume left-to-right
					size_t beg = j;
					size_t end = j + GFM::default_ftabChars;
					// Mismatch penalty is 3, so we have to skip starting
					// points that are within 2 from the mismatch
					if((i > 0 && j > 0) || j == seq.length()-1) {
						// Right-to-left
						if(beg < GFM::default_ftabChars) {
							beg = 0;
						} else {
							beg -= GFM::default_ftabChars;
						}
						end -= GFM::default_ftabChars;
					}
					size_t kk = k;
					//if(rc) {
					//	kk = seq.length() - k - 1;
					//}
					if(beg <= kk && end > kk) {
						continue;
					}
					if((j > kk) ? (j - kk <= 2) : (kk - j <= 2)) {
						continue;
					}
					cerr << "Test " << (++testnum) << endl;
					cerr << "  Query with length greater than ftab and matches exactly once with 1mm" << endl;
					DescentMetrics mets;
					PerReadMetrics prm;
					DescentDriver dr;
					
					dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -30, 30);
					
					// Set up the DescentConfig
					DescentConfig conf;
					// Changed 
					conf.cons.init(0, 1.0);
					conf.expol = DESC_EX_NONE;
					
					// Set up the search roots
					dr.addRoot(
						conf,    // DescentConfig
						j,       // 5' offset into read of root
						i == 0,  // left-to-right?
						true,    // forward?
						0.0f);    // root priority
					
					// Do the search
					Scoring sc = Scoring::base1();
					dr.go(sc, *gfms.first, *gfms.second, mets, prm);
					
					// Confirm that an exact-matching alignment was found
					assert_eq(1, dr.sink().nrange());
					assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
					assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
					cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
					assert_eq(1, dr.sink().nelt());
					last_topf = dr.sink()[0].topf;
					last_botf = dr.sink()[0].botf;
				}
			}
		}
    }

	// Query is longer than ftab and matches exactly once with one N mismatch
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for(int i = 0; i < 2; i++) {
			// Set up the read
			//    Ref: CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAACCA
			//                ||||||||||||||||||||||||||||||
			BTDnaString orig("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
			//                012345678901234567890123456789
			BTString    qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
			for(size_t k = 0; k < orig.length(); k++) {
				BTDnaString seq = orig;
				seq.set(4, k);
				for(size_t j = 0; j < seq.length(); j++) {
					// Assume left-to-right
					size_t beg = j;
					size_t end = j + GFM::default_ftabChars;
					// Mismatch penalty is 3, so we have to skip starting
					// points that are within 2 from the mismatch
					if((i > 0 && j > 0) || j == seq.length()-1) {
						// Right-to-left
						if(beg < GFM::default_ftabChars) {
							beg = 0;
						} else {
							beg -= GFM::default_ftabChars;
						}
						end -= GFM::default_ftabChars;
					}
					if(beg <= k && end > k) {
						continue;
					}
					if((j > k) ? (j - k <= 2) : (k - j <= 2)) {
						continue;
					}
					cerr << "Test " << (++testnum) << endl;
					cerr << "  Query with length greater than ftab and matches exactly once with 1mm" << endl;
					DescentMetrics mets;
					PerReadMetrics prm;
					DescentDriver dr;
					
					dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -30, 30);
					
					// Set up the DescentConfig
					DescentConfig conf;
					// Changed 
					conf.cons.init(0, 1.0);
					conf.expol = DESC_EX_NONE;
					
					// Set up the search roots
					dr.addRoot(
						conf,   // DescentConfig
						j,      // 5' offset into read of root
						i == 0, // left-to-right?
						true,   // forward?
						0.0f);   // root priority
					
					// Do the search
					Scoring sc = Scoring::base1();
					dr.go(sc, *gfms.first, *gfms.second, mets, prm);
					
					// Confirm that an exact-matching alignment was found
					assert_eq(1, dr.sink().nrange());
					assert_eq(sc.n(40), dr.sink()[0].pen);
					assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
					assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
					cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
					assert_eq(1, dr.sink().nelt());
					last_topf = dr.sink()[0].topf;
					last_botf = dr.sink()[0].botf;
				}
			}
		}
    }

	// Throw a bunch of queries with a bunch of Ns in and try to force an assert
	{
		RandomSource rnd(79);
		for(int i = 0; i < 2; i++) {
			// Set up the read
			//    Ref: CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAACCA
			//                ||||||||||||||||||||||||||||||
			BTDnaString orig("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
			//                012345678901234567890123456789
			BTString    qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
			if(i == 1) {
				orig.reverseComp();
				qual.reverse();
			}
			for(size_t trials = 0; trials < 100; trials++) {
				BTDnaString seq = orig;
				size_t ns = 10;
				for(size_t k = 0; k < ns; k++) {
					size_t pos = rnd.nextU32() % seq.length();
					seq.set(4, pos);
				}
				
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query with a bunch of Ns" << endl;
				
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				
				dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -30, 30);
				
				// Set up the DescentConfig
				DescentConfig conf;
				// Changed 
				conf.cons.init(GFM::default_ftabChars, 1.0);
				conf.expol = DESC_EX_NONE;
				
				// Set up the search roots
				for(size_t k = 0; k < ns; k++) {
					size_t j = rnd.nextU32() % seq.length();
					bool ltr = (rnd.nextU2() == 0) ? true : false;
					bool fw = (rnd.nextU2() == 0) ? true : false;
					dr.addRoot(
						conf,   // DescentConfig
						j,      // 5' offset into read of root
						ltr,    // left-to-right?
						fw,     // forward?
						0.0f);   // root priority
				}
				
				// Do the search
				Scoring sc = Scoring::base1();
				dr.go(sc, *gfms.first, *gfms.second, mets, prm);
			}
		}
    }

	// Query is longer than ftab and matches exactly once with one mismatch
	{
		RandomSource rnd(77);
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for(int i = 0; i < 2; i++) {
			// Set up the read
			//    Ref: CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAACCA
			//                ||||||||||||||||||||||||||||||
			BTDnaString orig("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
			//                012345678901234567890123456789
			BTString    qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
			//       revcomp: ACACAAAATGATGCGAGCGCGCTATATAGC
			//       revqual: cbaIHGFEDCBAihgfedcbaIHGFEDCBA
			bool fwi = (i == 0);
			if(!fwi) {
				orig.reverseComp();
			}
			for(size_t k = 0; k < orig.length(); k++) {
				BTDnaString seq = orig;
				seq.set(seq[k] ^ 3, k);
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query with length greater than ftab and matches exactly once with 1mm.  Many search roots." << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				
				dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -30, 30);
				
				// Set up the DescentConfig
				DescentConfig conf;
				// Changed 
				conf.cons.init(0, 1.0);
				conf.expol = DESC_EX_NONE;
				
				// Set up several random search roots
				bool onegood = false;
				for(size_t y = 0; y < 10; y++) {
					size_t j = rnd.nextU32() % seq.length();
					bool ltr = (rnd.nextU2() == 0) ? true : false;
					bool fw = (rnd.nextU2() == 0) ? true : false;
					dr.addRoot(
						conf,     // DescentConfig
						(TReadOff)j,        // 5' offset into read of root
						ltr,      // left-to-right?
						fw,       // forward?
						(float)((float)y * 1.0f)); // root priority
					// Assume left-to-right
					size_t beg = j;
					size_t end = j + GFM::default_ftabChars;
					// Mismatch penalty is 3, so we have to skip starting
					// points that are within 2 from the mismatch
					if(!ltr) {
						// Right-to-left
						if(beg < GFM::default_ftabChars) {
							beg = 0;
						} else {
							beg -= GFM::default_ftabChars;
						}
						end -= GFM::default_ftabChars;
					}
					bool good = true;
					if(fw != fwi) {
						good = false;
					}
					if(beg <= k && end > k) {
						good = false;
					}
					if((j > k) ? (j - k <= 2) : (k - j <= 2)) {
						good = false;
					}
					if(good) {
						onegood = true;
					}
				}
				if(!onegood) {
					continue;
				}
				
				// Do the search
				Scoring sc = Scoring::base1();
				dr.go(sc, *gfms.first, *gfms.second, mets, prm);
				
				// Confirm that an exact-matching alignment was found
				assert_eq(1, dr.sink().nrange());
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
				assert_eq(1, dr.sink().nelt());
				last_topf = dr.sink()[0].topf;
				last_botf = dr.sink()[0].botf;
			}
		}
    }

	// Query is longer than ftab and matches exactly once with one read gap
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for(int i = 0; i < 2; i++) {
		for(int k = 0; k < 2; k++) {
			// Set up the read
			//                GCTATATAGCGCGCCTGCATCATTTTGTGT
			//    Ref: CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAACCA
			//                |||||||||||||||///////////////
			BTDnaString seq ("GCTATATAGCGCGCTGCATCATTTTGTGT", true);
			//                01234567890123456789012345678
			//                87654321098765432109876543210
			BTString    qual("ABCDEFGHIabcdefghiABCDEFGHIab");
			if(k == 1) {
				seq.reverseComp();
				qual.reverse();
			}
			assert_eq(seq.length(), qual.length());
			// js iterate over offsets from 5' end for the search root
			for(size_t j = 0; j < seq.length(); j++) {
				// Assume left-to-right
				size_t beg = j;
				if(k == 1) {
					beg = seq.length() - beg - 1;
				}
				size_t end = beg + GFM::default_ftabChars;
				// Mismatch penalty is 3, so we have to skip starting
				// points that are within 2 from the mismatch
				if((i > 0 && j > 0) || j == seq.length()-1) {
					// Right-to-left
					if(beg < GFM::default_ftabChars) {
						beg = 0;
					} else {
						beg -= GFM::default_ftabChars;
					}
					end -= GFM::default_ftabChars;
				}
				assert_geq(end, beg);
				if(beg <= 15 && end >= 15) {
					continue;
				}
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query matches once with a read gap of length 1" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				
				Read q("test", seq.toZBuf(), qual.toZBuf());
				assert(q.repOk());
				dr.initRead(q, -30, 30);
				
				// Set up the DescentConfig
				DescentConfig conf;
				// Changed 
				conf.cons.init(0, 0.5);
				conf.expol = DESC_EX_NONE;
				
				// Set up the search roots
				dr.addRoot(
					conf,   // DescentConfig
					j,      // 5' offset into read of root
					i == 0, // left-to-right?
					k == 0, // forward?
					0.0f);  // root priority
				
				// Do the search
				Scoring sc = Scoring::base1();
				dr.go(sc, *gfms.first, *gfms.second, mets, prm);
				
				// Confirm that an exact-matching alignment was found
				assert_eq(1, dr.sink().nrange());
				assert_eq(sc.readGapOpen() + 0 * sc.readGapExtend(), dr.sink()[0].pen);
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
				assert_eq(1, dr.sink().nelt());
				last_topf = dr.sink()[0].topf;
				last_botf = dr.sink()[0].botf;
			}
		}}
    }

	// Query is longer than ftab and matches exactly once with one read gap of
	// length 3
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for(int i = 0; i < 2; i++) {
		for(int k = 0; k < 2; k++) {
			// Set up the read
			//                GCTATATAGCGCGCGCTCATCATTTTGTGT
			//    Ref: CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAACCA
			//                ||||||||||||||   |||||||||||||
			BTDnaString seq ("GCTATATAGCGCGC" "CATCATTTTGTGT", true);
			//                01234567890123   4567890123456
			//                65432109876543   2109876543210
			BTString    qual("ABCDEFGHIabcde" "fghiABCDEFGHI");
			if(k == 1) {
				seq.reverseComp();
				qual.reverse();
			}
			for(size_t j = 0; j < seq.length(); j++) {
				// Assume left-to-right
				size_t beg = j;
				if(k == 1) {
					beg = seq.length() - beg - 1;
				}
				size_t end = beg + GFM::default_ftabChars;
				// Mismatch penalty is 3, so we have to skip starting
				// points that are within 2 from the mismatch
				if((i > 0 && j > 0) || j == seq.length()-1) {
					// Right-to-left
					if(beg < GFM::default_ftabChars) {
						beg = 0;
					} else {
						beg -= GFM::default_ftabChars;
					}
					end -= GFM::default_ftabChars;
				}
				if(beg <= 14 && end >= 14) {
					continue;
				}
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query matches once with a read gap of length 3" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				
				dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -30, 30);
				
				// Set up the DescentConfig
				DescentConfig conf;
				// Changed
				conf.cons.init(0, 0.2);
				conf.expol = DESC_EX_NONE;
				
				// Set up the search roots
				dr.addRoot(
					conf,   // DescentConfig
					j,      // 5' offset into read of root
					i == 0, // left-to-right?
					k == 0, // forward?
					0.0f);  // root priority
				
				// Do the search
				Scoring sc = Scoring::base1();
				// Need to adjust the mismatch penalty up to avoid alignments
				// with lots of mismatches.
				sc.setMmPen(COST_MODEL_CONSTANT, 6, 6);
				dr.go(sc, *gfms.first, *gfms.second, mets, prm);
				
				// Confirm that an exact-matching alignment was found
				assert_eq(1, dr.sink().nrange());
				assert_eq(sc.readGapOpen() + 2 * sc.readGapExtend(), dr.sink()[0].pen);
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
				assert_eq(1, dr.sink().nelt());
				last_topf = dr.sink()[0].topf;
				last_botf = dr.sink()[0].botf;
			}
		}}
    }

	// Query is longer than ftab and matches exactly once with one reference gap
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for(int i = 0; i < 2; i++) {
			// Set up the read
			//    Ref: CATGTCAGCTATATAGCGCGC" "TCGCATCATTTTGTGTGTAAACCA
			//                ||||||||||||||   ||||||||||||||||
			BTDnaString seq ("GCTATATAGCGCGCA""TCGCATCATTTTGTGT", true);
			//                012345678901234  5678901234567890
			BTString    qual("ABCDEFGHIabcdef""ghiABCDEFGHIabcd");
			for(size_t j = 0; j < seq.length(); j++) {
				// Assume left-to-right
				size_t beg = j;
				size_t end = j + GFM::default_ftabChars;
				// Mismatch penalty is 3, so we have to skip starting
				// points that are within 2 from the mismatch
				if((i > 0 && j > 0) || j == seq.length()-1) {
					// Right-to-left
					if(beg < GFM::default_ftabChars) {
						beg = 0;
					} else {
						beg -= GFM::default_ftabChars;
					}
					end -= GFM::default_ftabChars;
				}
				if(beg <= 14 && end >= 14) {
					continue;
				}
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query matches once with a reference gap of length 1" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				
				dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -30, 30);
				
				// Set up the DescentConfig
				DescentConfig conf;
				// Changed 
				conf.cons.init(1, 0.5);
				conf.expol = DESC_EX_NONE;
				
				// Set up the search roots
				dr.addRoot(
					conf,   // DescentConfig
					j,      // 5' offset into read of root
					i == 0, // left-to-right?
					true,   // forward?
					0.0f);  // root priority
				
				// Do the search
				Scoring sc = Scoring::base1();
				// Need to adjust the mismatch penalty up to avoid alignments
				// with lots of mismatches.
				sc.setMmPen(COST_MODEL_CONSTANT, 6, 6);
				dr.go(sc, *gfms.first, *gfms.second, mets, prm);
				
				// Confirm that an exact-matching alignment was found
				assert_eq(1, dr.sink().nrange());
				assert_eq(sc.refGapOpen() + 0 * sc.refGapExtend(), dr.sink()[0].pen);
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
				assert_eq(1, dr.sink().nelt());
				last_topf = dr.sink()[0].topf;
				last_botf = dr.sink()[0].botf;
			}
		}
    }

	// Query is longer than ftab and matches exactly once with one reference gap
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for(int i = 0; i < 2; i++) {
			// Set up the read
			//    Ref: CATGTCAGCTATATAGCGCGC"   "TCGCATCATTTTGTGTGTAAACCA
			//                ||||||||||||||     ||||||||||||||||
			BTDnaString seq ("GCTATATAGCGCGCATG""TCGCATCATTTTGTGT", true);
			//                01234567890123456  7890123456789012
			BTString    qual("ABCDEFGHIabcdefgh""iABCDEFGHIabcdef");
			for(size_t j = 0; j < seq.length(); j++) {
				// Assume left-to-right
				size_t beg = j;
				size_t end = j + GFM::default_ftabChars;
				// Mismatch penalty is 3, so we have to skip starting
				// points that are within 2 from the mismatch
				if((i > 0 && j > 0) || j == seq.length()-1) {
					// Right-to-left
					if(beg < GFM::default_ftabChars) {
						beg = 0;
					} else {
						beg -= GFM::default_ftabChars;
					}
					end -= GFM::default_ftabChars;
				}
				if(beg <= 14 && end >= 14) {
					continue;
				}
				if(beg <= 15 && end >= 15) {
					continue;
				}
				if(beg <= 16 && end >= 16) {
					continue;
				}
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query matches once with a reference gap of length 1" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				
				dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -30, 30);
				
				// Set up the DescentConfig
				DescentConfig conf;
				// Changed 
				conf.cons.init(1, 0.25);
				conf.expol = DESC_EX_NONE;
				
				// Set up the search roots
				dr.addRoot(
					conf,   // DescentConfig
					j,      // 5' offset into read of root
					i == 0, // left-to-right?
					true,   // forward?
					0.0f);  // root priority
				
				// Do the search
				Scoring sc = Scoring::base1();
				// Need to adjust the mismatch penalty up to avoid alignments
				// with lots of mismatches.
				sc.setMmPen(COST_MODEL_CONSTANT, 6, 6);
				dr.go(sc, *gfms.first, *gfms.second, mets, prm);
				
				// Confirm that an exact-matching alignment was found
				assert_eq(1, dr.sink().nrange());
				assert_eq(sc.refGapOpen() + 2 * sc.refGapExtend(), dr.sink()[0].pen);
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
				assert_eq(1, dr.sink().nelt());
				last_topf = dr.sink()[0].topf;
				last_botf = dr.sink()[0].botf;
			}
		}
    }

	// Query is longer than ftab and matches exactly once with one read gap,
	// one ref gap, and one mismatch
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for(int i = 0; i < 2; i++) {
			// Set up the read
			//           Ref: CATGTCAGCT   ATATAGCGCGCT  CGCATCATTTTGTGTGTAAACCA
			//                ||||||||||   ||||||||||||   |||||| |||||||||||||
			BTDnaString seq ("CATGTCAGCT""GATATAGCGCGCT" "GCATCAATTTGTGTGTAAAC", true);
			//                0123456789  0123456789012   34567890123456789012
			BTString    qual("ABCDEFGHIa""bcdefghiACDEF" "GHIabcdefghijkABCDEF");
			for(size_t j = 0; j < seq.length(); j++) {
				// Assume left-to-right
				size_t beg = j;
				size_t end = j + GFM::default_ftabChars;
				// Mismatch penalty is 3, so we have to skip starting
				// points that are within 2 from the mismatch
				if((i > 0 && j > 0) || j == seq.length()-1) {
					// Right-to-left
					if(beg < GFM::default_ftabChars) {
						beg = 0;
					} else {
						beg -= GFM::default_ftabChars;
					}
					end -= GFM::default_ftabChars;
				}
				if(beg <= 10 && end >= 10) {
					continue;
				}
				if(beg <= 22 && end >= 22) {
					continue;
				}
				if(beg <= 30 && end >= 30) {
					continue;
				}
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query matches once with a read gap of length 1" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				
				dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -50, 50);
				
				// Set up the DescentConfig
				DescentConfig conf;
				// Changed 
				conf.cons.init(1, 0.5);
				conf.expol = DESC_EX_NONE;
				
				// Set up the search roots
				dr.addRoot(
					conf,   // DescentConfig
					j,      // 5' offset into read of root
					i == 0, // left-to-right?
					true,   // forward?
					0.0f);  // root priority
				
				// Do the search
				Scoring sc = Scoring::base1();
				dr.go(sc, *gfms.first, *gfms.second, mets, prm);
				
				// Confirm that an exact-matching alignment was found
				assert_eq(1, dr.sink().nrange());
				assert_eq(sc.readGapOpen() + sc.refGapOpen() + sc.mm((int)'d' - 33), dr.sink()[0].pen);
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
				assert_eq(1, dr.sink().nelt());
				last_topf = dr.sink()[0].topf;
				last_botf = dr.sink()[0].botf;
			}
		}
    }

	delete gfms.first;
	delete gfms.second;
	
	//  Ref CATGTCAGCT-ATATAGCGCGCTCGCATCATTTTGTGTGTAAAC
	//      |||||||||| |||||||||||| |||||| |||||||||||||
	//  Rd  CATGTCAGCTGATATAGCGCGCT-GCATCAATTTGTGTGTAAAC
	strs.clear();
    strs.push_back(string("CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAAC"
                          "NNNNNNNNNN"
                          "CATGTCAGCTGATATAGCGCGCTCGCATCATTTTGTGTGTAAAC" // same but without first ref gap
                          "N"
                          "CATGTCAGCTATATAGCGCGCTGCATCATTTTGTGTGTAAAC" // same but without first read gap
                          "N"
                          "CATGTCAGCTATATAGCGCGCTCGCATCAATTTGTGTGTAAAC" // same but without first mismatch
                          "N"
                          "CATGTCAGCTGATATAGCGCGCTGCATCAATTTGTGTGTAAAC" // Exact match for read
						  ));
	gfms = GFM::fromStrings<SString<char> >(
		strs,
		packed,
		REF_READ_REVERSE,
		GFM::default_bigEndian,
		GFM::default_lineRate,
		GFM::default_offRate,
		GFM::default_ftabChars,
		".aligner_seed2.cpp.tmp",
		GFM::default_useBlockwise,
		GFM::default_bmax,
		GFM::default_bmaxMultSqrt,
		GFM::default_bmaxDivN,
		GFM::default_dcv,
		GFM::default_seed,
		false,  // verbose
		false,  // autoMem
		false); // sanity
    
    gfms.first->loadIntoMemory (color, -1, true, true, true, true, false);
    gfms.second->loadIntoMemory(color,  1, true, true, true, true, false);

	// Query is longer than ftab and matches exactly once with one read gap,
	// one ref gap, and one mismatch
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for(int i = 0; i < 2; i++) {
			// Set up the read
			//           Ref: CATGTCAGCT   ATATAGCGCGCT  CGCATCATTTTGTGTGTAAACCA
			//                ||||||||||   ||||||||||||   |||||| |||||||||||||
			BTDnaString seq ("CATGTCAGCT""GATATAGCGCGCT" "GCATCAATTTGTGTGTAAAC", true);
			//                0123456789  0123456789012   34567890123456789012
			BTString    qual("ABCDEFGHIa""bcdefghiACDEF" "GHIabcdefghijkABCDEF");
			for(size_t j = 0; j < seq.length(); j++) {
				// Assume left-to-right
				size_t beg = j;
				size_t end = j + GFM::default_ftabChars;
				// Mismatch penalty is 3, so we have to skip starting
				// points that are within 2 from the mismatch
				if((i > 0 && j > 0) || j == seq.length()-1) {
					// Right-to-left
					if(beg < GFM::default_ftabChars) {
						beg = 0;
					} else {
						beg -= GFM::default_ftabChars;
					}
					end -= GFM::default_ftabChars;
				}
				if(beg <= 10 && end >= 10) {
					continue;
				}
				if(beg <= 22 && end >= 22) {
					continue;
				}
				if(beg <= 30 && end >= 30) {
					continue;
				}
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query matches once with a read gap of length 1" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				
				dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -50, 50);
				
				// Set up the DescentConfig
				DescentConfig conf;
				// Changed 
				conf.cons.init(1, 0.5);
				conf.expol = DESC_EX_NONE;
				
				// Set up the search roots
				dr.addRoot(
					conf,   // DescentConfig
					j,      // 5' offset into read of root
					i == 0, // left-to-right?
					true,   // forward?
					0.0f);  // root priority
				
				// Do the search
				Scoring sc = Scoring::base1();
				dr.go(sc, *gfms.first, *gfms.second, mets, prm);
				
				// Confirm that an exact-matching alignment was found
				assert_eq(5, dr.sink().nrange());
				assert_eq(0, dr.sink()[0].pen);
				assert_eq(min(sc.readGapOpen(), sc.refGapOpen()) + sc.mm((int)'d' - 33), dr.sink()[1].pen);
				assert_eq(max(sc.readGapOpen(), sc.refGapOpen()) + sc.mm((int)'d' - 33), dr.sink()[2].pen);
				assert_eq(sc.readGapOpen() + sc.refGapOpen(), dr.sink()[3].pen);
				assert_eq(sc.readGapOpen() + sc.refGapOpen() + sc.mm((int)'d' - 33), dr.sink()[4].pen);
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
				assert_eq(5, dr.sink().nelt());
				last_topf = dr.sink()[0].topf;
				last_botf = dr.sink()[0].botf;
			}
		}
    }

	// Query is longer than ftab and matches exactly once with one read gap,
	// one ref gap, one mismatch, and one N
	{
		size_t last_topf = std::numeric_limits<size_t>::max();
		size_t last_botf = std::numeric_limits<size_t>::max();
		for(int i = 0; i < 2; i++) {
			// Set up the read
			//           Ref: CATGTCAGCT   ATATAGCGCGCT  CGCATCATTTTGTGTGTAAACCA
			//                ||||||||||   ||||||||||||   |||||| |||||| ||||||
			BTDnaString seq ("CATGTCAGCT""GATATAGCGCGCT" "GCATCAATTTGTGNGTAAAC", true);
			//                0123456789  0123456789012   34567890123456789012
			BTString    qual("ABCDEFGHIa""bcdefghiACDEF" "GHIabcdefghijkABCDEF");
			for(size_t j = 0; j < seq.length(); j++) {
				// Assume left-to-right
				size_t beg = j;
				size_t end = j + GFM::default_ftabChars;
				// Mismatch penalty is 3, so we have to skip starting
				// points that are within 2 from the mismatch
				if((i > 0 && j > 0) || j == seq.length()-1) {
					// Right-to-left
					if(beg < GFM::default_ftabChars) {
						beg = 0;
					} else {
						beg -= GFM::default_ftabChars;
					}
					end -= GFM::default_ftabChars;
				}
				if(beg <= 10 && end >= 10) {
					continue;
				}
				if(beg <= 22 && end >= 22) {
					continue;
				}
				if(beg <= 30 && end >= 30) {
					continue;
				}
				if(beg <= 36 && end >= 36) {
					continue;
				}
				cerr << "Test " << (++testnum) << endl;
				cerr << "  Query matches with various patterns of gaps, mismatches and Ns" << endl;
				DescentMetrics mets;
				PerReadMetrics prm;
				DescentDriver dr;
				
				dr.initRead(Read("test", seq.toZBuf(), qual.toZBuf()), -50, 50);
				
				// Set up the DescentConfig
				DescentConfig conf;
				// Changed 
				conf.cons.init(1, 0.5);
				conf.expol = DESC_EX_NONE;
				
				// Set up the search roots
				dr.addRoot(
					conf,   // DescentConfig
					j,      // 5' offset into read of root
					i == 0, // left-to-right?
					true,   // forward?
					0.0f);  // root priority
				
				// Do the search
				Scoring sc = Scoring::base1();
				sc.setNPen(COST_MODEL_CONSTANT, 1);
				dr.go(sc, *gfms.first, *gfms.second, mets, prm);
				
				// Confirm that an exact-matching alignment was found
				assert_eq(5, dr.sink().nrange());
				assert_eq(sc.n(40), dr.sink()[0].pen);
				assert_eq(sc.n(40) + min(sc.readGapOpen(), sc.refGapOpen()) + sc.mm((int)'d' - 33), dr.sink()[1].pen);
				assert_eq(sc.n(40) + max(sc.readGapOpen(), sc.refGapOpen()) + sc.mm((int)'d' - 33), dr.sink()[2].pen);
				assert_eq(sc.n(40) + sc.readGapOpen() + sc.refGapOpen(), dr.sink()[3].pen);
				assert_eq(sc.n(40) + sc.readGapOpen() + sc.refGapOpen() + sc.mm((int)'d' - 33), dr.sink()[4].pen);
				assert(last_topf == std::numeric_limits<size_t>::max() || last_topf == dr.sink()[0].topf);
				assert(last_botf == std::numeric_limits<size_t>::max() || last_botf == dr.sink()[0].botf);
				cerr << dr.sink()[0].topf << ", " << dr.sink()[0].botf << endl;
				assert_eq(5, dr.sink().nelt());
				last_topf = dr.sink()[0].topf;
				last_botf = dr.sink()[0].botf;
			}
		}
    }

    delete gfms.first;
    delete gfms.second;
	
	cerr << "DONE" << endl;
}
#endif

