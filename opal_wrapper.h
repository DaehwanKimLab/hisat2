/*
 * Copyright 2018, Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>
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

#ifndef __OPAL_WRAPPER_H__
#define __OPAL_WRAPPER_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include "timer.h"
#include "aligner_sw.h"
#include "aligner_result.h"
#include "scoring.h"
#include "sstring.h"
#include "opal.h"

// wrapper
class OPAL {
public:
    OPAL() {
        match_ = 0;
        mismatch_ = -1;
        alphasize_ = 4;
        gap_open_ = 2;
        gap_extend_ = 1;

        score_matrix_ = NULL;

        buildScoreMatrix(); 
    }


    ~OPAL() {
        if (score_matrix_) {
            delete[] score_matrix_;
        }
    }

public:
    int *score_matrix_;

    int match_;
    int mismatch_;
    int gap_open_;
    int gap_extend_;
    int minsc_;

    size_t alphasize_;

    void setMinScore(int rpt_edit)
    {
        minsc_ = rpt_edit * mismatch_;
    }

    void buildScoreMatrix()
    {
        int *ptr = new int[alphasize_ * alphasize_];

        for(size_t i = 0; i < alphasize_; i++) {
            for(size_t j = 0; j < alphasize_; j++) {
                if(i == j) {
                    ptr[i * alphasize_ + j] = match_; // match
                } else {
                    ptr[i * alphasize_ + j] = mismatch_; // mismatch
                }
            }
        }
        score_matrix_ = ptr;
    }



    int alignStrings(const string &ref, const string &read, EList<Edit>& edits)
    {
        BTDnaString btread;
        BTDnaString btref;

        btread.install(read.c_str(), true);
        btref.install(ref.c_str(), true);

        unsigned char *query = (unsigned char *)btread.wbuf();
        size_t query_len = btread.length();

        size_t db_length = 1;

        unsigned char *db[1];
        db[0] = (unsigned char *)btref.wbuf();
        int db_seqs_len[1];
        db_seqs_len[0] = btref.length();


        OpalSearchResult *results[1];

        results[0] = new OpalSearchResult;
        opalInitSearchResult(results[0]);

        int result_code = opalSearchDatabase(
                query,
                query_len,
                db,
                db_length,
                db_seqs_len,
                gap_open_,
                gap_extend_,
                score_matrix_,
                alphasize_,
                results,
                OPAL_SEARCH_ALIGNMENT,
                OPAL_MODE_NW,
                OPAL_OVERFLOW_SIMPLE);

        //cerr << "OPAL " << result_code << endl;

        if (result_code == 0) {
            assert_gt(results[0]->alignmentLength, 0);

            //cerr << "score " << results[0]->score << endl;

            if (results[0]->score > minsc_) {

                // Show alignment
#if 0
                printAlignments(
                        cerr, 
                        ref, read,
                        (char *)results[0]->alignment,
                        results[0]->alignmentLength
                        );
#endif

                makeEdits(
                        (char *)results[0]->alignment,
                        results[0]->alignmentLength,
                        ref,
                        read,
                        edits);

                //Edit::print(cerr, edits); cerr << endl;
            }

            free(results[0]->alignment);

        }

        return 0;

    }


    void makeEdits(const char *alignments, const int alignment_length,
            const string& ref, const string& read,
            EList<Edit>& edits)
    {
        size_t read_pos = 0;
        size_t ref_pos = 0;

#if 0
        cerr << ref << endl;
        cerr << read << endl;
        for(size_t i = 0; i < alignment_length; i++) {
            cerr << (int)alignments[i];
        }
        cerr << endl;
#endif

        for(size_t i = 0; i < alignment_length; i++) {
            if(alignments[i] == OPAL_ALIGN_MATCH) {
                read_pos++;
                ref_pos++;
            } else if(alignments[i] == OPAL_ALIGN_MISMATCH) {
                edits.expand();
                edits.back().init(read_pos, 
                        ref[ref_pos],
                        read[read_pos],
                        EDIT_TYPE_MM);
                read_pos++;
                ref_pos++;
            } else if(alignments[i] == OPAL_ALIGN_DEL) {
                edits.expand();
                edits.back().init(read_pos,
                        '-',
                        read[read_pos],
                        EDIT_TYPE_REF_GAP);
                read_pos++;
            } else if(alignments[i] == OPAL_ALIGN_INS) {
                edits.expand();
                edits.back().init(read_pos,
                        ref[ref_pos],
                        '-',
                        EDIT_TYPE_READ_GAP);
                ref_pos++;
            }
        }

        assert_eq(read_pos, read.length());
        assert_eq(ref_pos, ref.length());
    }


    void printAlignments(ostream& os, const string& ref, const string& read,
            const char *alignment, size_t alignemt_len)
    {

        // print Ref
        size_t ref_pos = 0;
        size_t read_pos = 0;
        
        os << "REF : ";
        for(size_t i = 0; i < alignemt_len; i++) {
            if(alignment[i] == OPAL_ALIGN_DEL) {
                os << "-";
            } else {
                os << ref[ref_pos++];
            }
        }
        os << endl;
        os << "      ";
        for(size_t i = 0; i < alignemt_len; i++) {
            if(alignment[i] == OPAL_ALIGN_MATCH) {
                os << "|";
            } else {
                os << " ";
            }
        }
        os << endl;
        os << "READ: ";
        for(size_t i = 0; i < alignemt_len; i++) {
            if(alignment[i] == OPAL_ALIGN_INS) {
                os << "-";
            } else {
                os << read[read_pos++];
            }
        }
        os << endl;
    }

};

#endif /* __OPAL_WRAPPER_H__ */
