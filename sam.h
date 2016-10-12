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

#ifndef SAM_H_
#define SAM_H_

#include <string>
#include <sys/time.h>
#include "ds.h"
#include "read.h"
#include "util.h"
#include "aligner_result.h"
#include "scoring.h"
#include "alt.h"
#include "filebuf.h"

enum {
	// Comments use language from v1.4-r962 spec
	SAM_FLAG_PAIRED         = 1,   // templ. having mult. frag.s in sequencing
	SAM_FLAG_MAPPED_PAIRED  = 2,   // each frag properly aligned
	SAM_FLAG_UNMAPPED       = 4,   // fragment unmapped
	SAM_FLAG_MATE_UNMAPPED  = 8,   // next fragment in template unmapped
	SAM_FLAG_QUERY_STRAND   = 16,  // SEQ is reverse comp'ed from original
	SAM_FLAG_MATE_STRAND    = 32,  // next fragment SEQ reverse comp'ed
	SAM_FLAG_FIRST_IN_PAIR  = 64,  // first fragment in template
	SAM_FLAG_SECOND_IN_PAIR = 128, // last fragment in template
	SAM_FLAG_NOT_PRIMARY    = 256, // secondary alignment
	SAM_FLAG_FAILS_CHECKS   = 512, // not passing quality controls
	SAM_FLAG_DUPLICATE      = 1024 // PCR or optical duplicate
};

class AlnRes;
class AlnFlags;
class AlnSetSumm;

/**
 * Encapsulates all the various ways that a user may wish to customize SAM
 * output.
 */
template<typename index_t>
class SamConfig {

	typedef EList<std::string> StrList;
	typedef EList<size_t> LenList;

public:

	SamConfig(
              const StrList& refnames,  // reference sequence names
              const LenList& reflens,   // reference sequence lengths
              bool truncQname,          // truncate read name to 255?
              bool omitsec,             // omit secondary SEQ/QUAL
              bool noUnal,              // omit unaligned reads
              const std::string& pg_id, // id
              const std::string& pg_pn, // name
              const std::string& pg_vn, // version
              const std::string& pg_cl, // command-line
              const std::string& rgs,   // read groups string
              int rna_strandness,
              bool print_as,
              bool print_xs,
              bool print_xss,
              bool print_yn,
              bool print_xn,
              bool print_cs,
              bool print_cq,
              bool print_x0,
              bool print_x1,
              bool print_xm,
              bool print_xo,
              bool print_xg,
              bool print_nm,
              bool print_md,
              bool print_yf,
              bool print_yi,
              bool print_ym,
              bool print_yp,
              bool print_yt,
              bool print_ys,
              bool print_zs,
              bool print_xr,
              bool print_xt,
              bool print_xd,
              bool print_xu,
              bool print_ye, // streak of failed DPs at end
              bool print_yl, // longest streak of failed DPs
              bool print_yu, // index of last succeeded DP
              bool print_xp, // print seed hit information
              bool print_yr, // # redundant seed hits
              bool print_zb, // # Ftab lookups
              bool print_zr, // # redundant path checks
              bool print_zf, // # FM Index ops
              bool print_zm, // FM Index op string for best-first search
              bool print_zi, // # seed extend loop iters
              bool print_zp,
              bool print_zu,
              bool print_xs_a,
              bool print_nh) :
		truncQname_(truncQname),
		omitsec_(omitsec),
		noUnal_(noUnal),
		pg_id_(pg_id),
		pg_pn_(pg_pn),
		pg_vn_(pg_vn),
		pg_cl_(pg_cl),
		rgs_(rgs),
		refnames_(refnames),
		reflens_(reflens),
        rna_strandness_(rna_strandness),
		print_as_(print_as), // alignment score of best alignment
		print_xs_(print_xs), // alignment score of second-best alignment
		print_xss_(print_xss),
		print_yn_(print_yn), // minimum valid score and perfect score
		print_xn_(print_xn),
		print_cs_(print_cs),
		print_cq_(print_cq),
		print_x0_(print_x0),
		print_x1_(print_x1),
		print_xm_(print_xm),
		print_xo_(print_xo),
		print_xg_(print_xg),
		print_nm_(print_nm),
		print_md_(print_md),
		print_yf_(print_yf),
		print_yi_(print_yi),
		print_ym_(print_ym),
		print_yp_(print_yp),
		print_yt_(print_yt),
		print_ys_(print_ys),
		print_zs_(print_zs),
		print_xr_(print_xr),
		print_xt_(print_xt), // time elapsed in microseconds
		print_xd_(print_xd), // DP extend attempts
		print_xu_(print_xu), // ungapped extend attempts
		print_ye_(print_ye), // streak of failed DPs at end
		print_yl_(print_yl), // longest streak of failed DPs
		print_yu_(print_yu), // index of last succeeded DP
		print_xp_(print_xp), // print seed hit information
		print_yr_(print_yr), // index of last succeeded DP
		print_zb_(print_zb), // # Ftab lookups
		print_zr_(print_zr), // # redundant path checks
		print_zf_(print_zf), // # FM Index ops
		print_zm_(print_zm), // FM Index op string for best-first search
		print_zi_(print_zi), // # seed extend loop iters
		print_zp_(print_zp), // # seed extend loop iters
		print_zu_(print_zu), // # seed extend loop iters
        print_xs_a_(print_xs_a),
        print_nh_(print_nh)
	{
		assert_eq(refnames_.size(), reflens_.size());
	}

	/**
	 * Print a reference name in a way that doesn't violate SAM's character
	 * constraints. \*|[!-()+-<>-~][!-~]*
	 */
	void printRefName(
		BTString& o,
		const std::string& name)
		const;

	/**
	 * Print a :Z optional field where certain characters (whitespace, colon
	 * and percent) are escaped using % escapes.
	 */
	template<typename T>
	void printOptFieldEscapedZ(BTString& o, const T& s) const {
		size_t len = s.length();
		for(size_t i = 0; i < len; i++) {
			if(s[i] < 33 || s[i] > 126 || s[i] == ':' || s[i] == '%') {
				// percent-encode it
				o.append('%');
				int ms = s[i] >> 4;
				int ls = s[i] & 15;
				assert_range(0, 15, ms);
				assert_range(0, 15, ls);
				o.append("0123456789ABCDEF"[ms]);
				o.append("0123456789ABCDEF"[ls]);
			} else {
				o.append(s[i]);
			}
		}
	}

	/**
	 * Print a :Z optional field where newline characters are escaped using %
	 * escapes.
	 */
	template<typename T>
	void printOptFieldNewlineEscapedZ(BTString& o, const T& s) const {
		size_t len = s.length();
		for(size_t i = 0; i < len; i++) {
			if(s[i] == 10 || s[i] == 13 || s[i] == '%') {
				// percent-encode it
				o.append('%');
				int ms = s[i] >> 4;
				int ls = s[i] & 15;
				assert_range(0, 15, ms);
				assert_range(0, 15, ls);
				o.append("0123456789ABCDEF"[ms]);
				o.append("0123456789ABCDEF"[ls]);
			} else {
				o.append(s[i]);
			}
		}
	}
	
	/**
	 * Print a read name in a way that doesn't violate SAM's character
	 * constraints. [!-?A-~]{1,255} (i.e. [33, 63], [65, 126])
	 */
	template<typename TStr>
	void printReadName(
		BTString& o,
		const TStr& name,
		bool omitSlashMate)
		const
	{
		size_t namelen = name.length();
		if(omitSlashMate &&
		   namelen >= 2 &&
		   name[namelen-2] == '/' &&
		   (name[namelen-1] == '1' || name[namelen-1] == '2' || name[namelen-1] == '3'))
		{
			namelen -= 2;
		}
		if(truncQname_ && namelen > 255) {
			namelen = 255;
		}
		for(size_t i = 0; i < namelen; i++) {
			if(truncQname_ && isspace(name[i])) {
				return;
			}
			o.append(name[i]);
		}
	}

	/**
	 * Print a reference name given a reference index.
	 */
	void printRefNameFromIndex(
		BTString& o,
		size_t i)
		const;
	
	/**
	 * Print SAM header to given output buffer.
	 */
	void printHeader(
		BTString& o,
		const std::string& rgid,
		const std::string& rgs,
		bool printHd,
		bool printSq,
		bool printPg)
		const;

	/**
	 * Print the @HD header line to the given string.
	 */
	void printHdLine(BTString& o, const char *samver) const;

	/**
	 * Print the @SQ header lines to the given string.
	 */
	void printSqLines(BTString& o) const;

	/**
	 * Print the @PG header line to the given string.
	 */
	void printPgLine(BTString& o) const;

	/**
	 * Print the optional flags to the given string.
	 */
	void printAlignedOptFlags(
		BTString& o,                 // output buffer
		bool first,                  // first opt flag printed is first overall?
		const Read& rd,              // the read
		AlnRes& res,                 // individual alignment result
		StackedAln& staln,           // stacked alignment
		const AlnFlags& flags,       // alignment flags
		const AlnSetSumm& summ,      // summary of alignments for this read
		const SeedAlSumm& ssm,       // seed alignment summary
		const PerReadMetrics& prm,   // per-read metics
		const Scoring& sc,           // scoring scheme
		const char *mapqInp,         // inputs to MAPQ calculation
        const ALTDB<index_t>* altdb)
		const;

	/**
	 * Print the optional flags to the given string.
	 */
	void printEmptyOptFlags(
		BTString& o,               // output buffer
		bool first,                // first opt flag printed is first overall?
		const Read& rd,            // the read
		const AlnFlags& flags,     // alignment flags
		const AlnSetSumm& summ,    // summary of alignments for this read
		const SeedAlSumm& ssm,     // seed alignment summary
		const PerReadMetrics& prm, // per-read metrics
		const Scoring& sc)         // scoring scheme
		const;
	
	/**
	 * Return true iff we should try to obey the SAM spec's recommendations
	 * that:
	 *
	 * SEQ and QUAL of secondary alignments should be set to ‘*’ to reduce the
	 * file size.
	 */
	bool omitSecondarySeqQual() const {
		return omitsec_;
	}
	
	bool omitUnalignedReads() const {
		return noUnal_;
	}

protected:

	bool truncQname_;   // truncate QNAME to 255 chars?
	bool omitsec_;      // omit secondary 
	bool noUnal_;       // omit unaligned reads
	
	std::string pg_id_; // @PG ID: Program record identifier
	std::string pg_pn_; // @PG PN: Program name
	std::string pg_vn_; // @PG VN: Program version
	std::string pg_cl_; // @PG CL: Program command-line
	std::string rgs_;   // Read-group string to add to all records
	const StrList& refnames_; // reference sequence names
	const LenList& reflens_;  // reference sequence lengths
    
    int rna_strandness_;
	
	// Which alignment flags to print?

	// Following are printed by BWA-SW
	bool print_as_; // AS:i: Alignment score generated by aligner
	bool print_xs_; // XS:i: Suboptimal alignment score
	bool print_xss_;// Xs:i: Best invalid alignment score found
	bool print_yn_; // YN:i:, Yn:i: minimum valid score and perfect score
	bool print_xn_; // XN:i: Number of ambiguous bases in the referenece

	// Other optional flags
	bool print_cs_; // CS:Z: Color read sequence on the original strand
	bool print_cq_; // CQ:Z: Color read quality on the original strand

	// Following are printed by BWA
	bool print_x0_; // X0:i: Number of best hits
	bool print_x1_; // X1:i: Number of sub-optimal best hits
	bool print_xm_; // XM:i: Number of mismatches in the alignment
	bool print_xo_; // XO:i: Number of gap opens
	bool print_xg_; // XG:i: Number of gap extensions (incl. opens)
	bool print_nm_; // NM:i: Edit dist. to the ref, Ns count, clipping doesn't
	bool print_md_; // MD:Z: String for mms. [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*2

	// Following are Bowtie2-specific
	bool print_yf_; // YF:i: Read was filtered out?
	bool print_yi_; // YI:Z: Summary of inputs to MAPQ calculation
	bool print_ym_; // YM:i: Read was repetitive when aligned unpaired?
	bool print_yp_; // YP:i: Read was repetitive when aligned paired?
	bool print_yt_; // YT:Z: String representing alignment type
	bool print_ys_; // YS:i: Score of other mate
	bool print_zs_; // ZS:i: Pseudo-random seed
	
	bool print_xr_; // XR:Z: Original read string
	bool print_xt_; // XT:i: Time taken to align
	bool print_xd_; // XD:i: DP problems
	bool print_xu_; // XU:i: ungapped alignment
	bool print_ye_; // YE:i: streak of failed DPs at end
	bool print_yl_; // YL:i: longest streak of failed DPs
	bool print_yu_; // YU:i: index of last succeeded DP
	bool print_xp_; // XP:BI: seed hit information
	bool print_yr_; // YR:i: # redundant seed hits
	bool print_zb_; // ZB:i: # Ftab lookups
	bool print_zr_; // ZR:i: # redundant path checks
	bool print_zf_; // ZF:i: # FM Index ops
	bool print_zm_; // ZM:i: FM ops string for best-first search
	bool print_zi_; // ZI:i: # extend loop iters
	bool print_zp_; // ZP:i: Score of best/second-best paired-end alignment
	bool print_zu_; // ZU:i: Score of best/second-best unpaired alignment
    
    bool print_xs_a_; // XS:A:[+=] Sense/anti-sense strand splice sites correspond to
    bool print_nh_;   // NH:i: # alignments
};

/**
 * Print a reference name in a way that doesn't violate SAM's character
 * constraints. \*|[!-()+-<>-~][!-~]* (i.e. [33, 63], [65, 126])
 */
template<typename index_t>
void SamConfig<index_t>::printRefName(
                             BTString& o,
                             const std::string& name) const
{
    size_t namelen = name.length();
    for(size_t i = 0; i < namelen; i++) {
        if(isspace(name[i])) {
            return;
        }
        o.append(name[i]);
    }
}

/**
 * Print a reference name given a reference index.
 */
template<typename index_t>
void SamConfig<index_t>::printRefNameFromIndex(BTString& o, size_t i) const {
    printRefName(o, refnames_[i]);
}

/**
 * Print SAM header to given output buffer.
 */
template<typename index_t>
void SamConfig<index_t>::printHeader(
                            BTString& o,
                            const string& rgid,
                            const string& rgs,
                            bool printHd,
                            bool printSq,
                            bool printPg) const
{
    if(printHd) printHdLine(o, "1.0");
    if(printSq) printSqLines(o);
    if(!rgid.empty()) {
        o.append("@RG");
        o.append(rgid.c_str());
        o.append(rgs.c_str());
        o.append('\n');
    }
    if(printPg) printPgLine(o);
}

/**
 * Print the @HD header line to the given string.
 */
template<typename index_t>
void SamConfig<index_t>::printHdLine(BTString& o, const char *samver) const {
    o.append("@HD\tVN:");
    o.append(samver);
    o.append("\tSO:unsorted\n");
}

/**
 * Print the @SQ header lines to the given string.
 */
template<typename index_t>
void SamConfig<index_t>::printSqLines(BTString& o) const {
    char buf[1024];
    for(size_t i = 0; i < refnames_.size(); i++) {
        o.append("@SQ\tSN:");
        printRefName(o, refnames_[i]);
        o.append("\tLN:");
        itoa10<size_t>(reflens_[i], buf);
        o.append(buf);
        o.append('\n');
    }
}

/**
 * Print the @PG header line to the given string.
 */
template<typename index_t>
void SamConfig<index_t>::printPgLine(BTString& o) const {
    o.append("@PG\tID:");
    o.append(pg_id_.c_str());
    o.append("\tPN:");
    o.append(pg_pn_.c_str());
    o.append("\tVN:");
    o.append(pg_vn_.c_str());
    o.append("\tCL:\"");
    o.append(pg_cl_.c_str());
    o.append('"');
    o.append('\n');
}

#define WRITE_SEP() { \
if(!first) o.append('\t'); \
first = false; \
}

/**
 * Print the optional flags to the given string.
 */
template<typename index_t>
void SamConfig<index_t>::printAlignedOptFlags(
                                     BTString& o,               // output buffer
                                     bool first,                // first opt flag printed is first overall?
                                     const Read& rd,            // the read
                                     AlnRes& res,               // individual alignment result
                                     StackedAln& staln,         // stacked alignment buffer
                                     const AlnFlags& flags,     // alignment flags
                                     const AlnSetSumm& summ,    // summary of alignments for this read
                                     const SeedAlSumm& ssm,     // seed alignment summary
                                     const PerReadMetrics& prm, // per-read metrics
                                     const Scoring& sc,         // scoring scheme
                                     const char *mapqInp,       // inputs to MAPQ calculation
                                     const ALTDB<index_t>* altdb)
const
{
    char buf[1024];
    if(print_as_) {
        // AS:i: Alignment score generated by aligner
        itoa10<TAlScore>(res.score().score(), buf);
        WRITE_SEP();
        o.append("AS:i:");
        o.append(buf);
    }
    
    // Do not output suboptimal alignment score, which conflicts with Cufflinks and StringTie
    if(print_xs_) {
        // XS:i: Suboptimal alignment score
        // Use ZS:i: to avoid conflict with XS:A:
        AlnScore sco = summ.secbestMate(rd.mate < 2);
        if(sco.valid()) {
            itoa10<TAlScore>(sco.score(), buf);
            WRITE_SEP();
            o.append("ZS:i:");
            o.append(buf);
        }
    }
    if(print_xn_) {
        // XN:i: Number of ambiguous bases in the referenece
        itoa10<size_t>(res.refNs(), buf);
        WRITE_SEP();
        o.append("XN:i:");
        o.append(buf);
    }
    if(print_x0_) {
        // X0:i: Number of best hits
    }
    if(print_x1_) {
        // X1:i: Number of sub-optimal best hits
    }
    size_t num_mm = 0;
    size_t num_go = 0;
    size_t num_gx = 0;
    for(size_t i = 0; i < res.ned().size(); i++) {
        if(res.ned()[i].isMismatch()) {
            if(res.ned()[i].snpID >= altdb->alts().size()) {
                num_mm++;
            }
        } else if(res.ned()[i].isReadGap()) {
            if(res.ned()[i].snpID >= altdb->alts().size()) {
                num_go++;
                num_gx++;
            }
            while(i < res.ned().size()-1 &&
                  res.ned()[i+1].pos == res.ned()[i].pos &&
                  res.ned()[i+1].isReadGap())
            {
                i++;
                if(res.ned()[i].snpID >= altdb->alts().size()) {
                    num_gx++;
                }
            }
        } else if(res.ned()[i].isRefGap()) {
            if(res.ned()[i].snpID >= altdb->alts().size()) {
                num_go++;
                num_gx++;
            }
            while(i < res.ned().size()-1 &&
                  res.ned()[i+1].pos == res.ned()[i].pos+1 &&
                  res.ned()[i+1].isRefGap())
            {
                i++;
                if(res.ned()[i].snpID >= altdb->alts().size()) {
                    num_gx++;
                }
            }
        }
    }
    if(print_xm_) {
        // XM:i: Number of mismatches in the alignment
        itoa10<size_t>(num_mm, buf);
        WRITE_SEP();
        o.append("XM:i:");
        o.append(buf);
    }
    if(print_xo_) {
        // XO:i: Number of gap opens
        itoa10<size_t>(num_go, buf);
        WRITE_SEP();
        o.append("XO:i:");
        o.append(buf);
    }
    if(print_xg_) {
        // XG:i: Number of gap extensions (incl. opens)
        itoa10<size_t>(num_gx, buf);
        WRITE_SEP();
        o.append("XG:i:");
        o.append(buf);
    }
    if(print_nm_) {
        // NM:i: Edit dist. to the ref, Ns count, clipping doesn't
        size_t NM = 0;
        for(size_t i = 0; i < res.ned().size(); i++) {
            if(res.ned()[i].type != EDIT_TYPE_SPL) {
                if(res.ned()[i].snpID >= altdb->alts().size()) {
                    NM++;
                }
            }
        }
        itoa10<size_t>(NM, buf);
        WRITE_SEP();
        o.append("NM:i:");
        o.append(buf);
    }
    if(print_md_) {
        // MD:Z: String for mms. [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*2
        WRITE_SEP();
        o.append("MD:Z:");
        staln.buildMdz();
        staln.writeMdz(
                       &o,        // output buffer
                       NULL);     // no char buffer
    }
    if(print_ys_ && summ.paired()) {
        // YS:i: Alignment score of opposite mate
        assert(res.oscore().valid());
        itoa10<TAlScore>(res.oscore().score(), buf);
        WRITE_SEP();
        o.append("YS:i:");
        o.append(buf);
    }
    if(print_yn_) {
        // YN:i: Minimum valid score for this mate
        TAlScore mn = sc.scoreMin.f<TAlScore>(rd.length());
        itoa10<TAlScore>(mn, buf);
        WRITE_SEP();
        o.append("YN:i:");
        o.append(buf);
        // Yn:i: Perfect score for this mate
        TAlScore pe = sc.perfectScore(rd.length());
        itoa10<TAlScore>(pe, buf);
        WRITE_SEP();
        o.append("Yn:i:");
        o.append(buf);
    }
    if(print_xss_) {
        // Xs:i: Best invalid alignment score of this mate
        bool one = true;
        if(flags.partOfPair() && !flags.readMate1()) {
            one = false;
        }
        TAlScore bst = one ? prm.bestLtMinscMate1 : prm.bestLtMinscMate2;
        if(bst > std::numeric_limits<TAlScore>::min()) {
            itoa10<TAlScore>(bst, buf);
            WRITE_SEP();
            o.append("Xs:i:");
            o.append(buf);
        }
        if(flags.partOfPair()) {
            // Ys:i: Best invalid alignment score of opposite mate
            bst = one ? prm.bestLtMinscMate2 : prm.bestLtMinscMate1;
            if(bst > std::numeric_limits<TAlScore>::min()) {
                itoa10<TAlScore>(bst, buf);
                WRITE_SEP();
                o.append("Ys:i:");
                o.append(buf);
            }
        }
    }
    if(print_zs_) {
        // ZS:i: Pseudo-random seed for read
        itoa10<uint32_t>(rd.seed, buf);
        WRITE_SEP();
        o.append("ZS:i:");
        o.append(buf);
    }
    if(print_yt_) {
        // YT:Z: String representing alignment type
        WRITE_SEP();
        flags.printYT(o);
    }
    if(print_yp_ && flags.partOfPair() && flags.canMax()) {
        // YP:i: Read was repetitive when aligned paired?
        WRITE_SEP();
        flags.printYP(o);
    }
    if(print_ym_ && flags.canMax() && (flags.isMixedMode() || !flags.partOfPair())) {
        // YM:i: Read was repetitive when aligned unpaired?
        WRITE_SEP();
        flags.printYM(o);
    }
    if(print_yf_ && flags.filtered()) {
        // YF:i: Read was filtered?
        first = flags.printYF(o, first) && first;
    }
    if(print_yi_) {
        // Print MAPQ calibration info
        if(mapqInp[0] != '\0') {
            // YI:i: Suboptimal alignment score
            WRITE_SEP();
            o.append("YI:Z:");
            o.append(mapqInp);
        }
    }
    if(flags.partOfPair() && print_zp_) {
        // ZP:i: Score of best concordant paired-end alignment
        WRITE_SEP();
        o.append("ZP:Z:");
        if(summ.bestPaired().valid()) {
            itoa10<TAlScore>(summ.bestPaired().score(), buf);
            o.append(buf);
        } else {
            o.append("NA");
        }
        // Zp:i: Second-best concordant paired-end alignment score
        WRITE_SEP();
        o.append("Zp:Z:");
        if(summ.secbestPaired().valid()) {
            itoa10<TAlScore>(summ.secbestPaired().score(), buf);
            o.append(buf);
        } else {
            o.append("NA");
        }
    }
    if(print_zu_) {
        // ZU:i: Score of best unpaired alignment
        AlnScore best    = (rd.mate <= 1 ? summ.best1()    : summ.best2());
        AlnScore secbest = (rd.mate <= 1 ? summ.secbest1() : summ.secbest2());
        WRITE_SEP();
        o.append("ZU:i:");
        if(best.valid()) {
            itoa10<TAlScore>(best.score(), buf);
            o.append(buf);
        } else {
            o.append("NA");
        }
        // Zu:i: Score of second-best unpaired alignment
        WRITE_SEP();
        o.append("Zu:i:");
        if(secbest.valid()) {
            itoa10<TAlScore>(secbest.score(), buf);
            o.append(buf);
        } else {
            o.append("NA");
        }
    }
    if(!rgs_.empty()) {
        WRITE_SEP();
        o.append(rgs_.c_str());
    }
    if(print_xt_) {
        // XT:i: Timing
        WRITE_SEP();
        struct timeval  tv_end;
        struct timezone tz_end;
        gettimeofday(&tv_end, &tz_end);
        size_t total_usecs =
        (tv_end.tv_sec  - prm.tv_beg.tv_sec) * 1000000 +
        (tv_end.tv_usec - prm.tv_beg.tv_usec);
        itoa10<size_t>(total_usecs, buf);
        o.append("XT:i:");
        o.append(buf);
    }
    if(print_xd_) {
        // XD:i: Extend DPs
        WRITE_SEP();
        itoa10<uint64_t>(prm.nExDps, buf);
        o.append("XD:i:");
        o.append(buf);
        // Xd:i: Mate DPs
        WRITE_SEP();
        itoa10<uint64_t>(prm.nMateDps, buf);
        o.append("Xd:i:");
        o.append(buf);
    }
    if(print_xu_) {
        // XU:i: Extend ungapped tries
        WRITE_SEP();
        itoa10<uint64_t>(prm.nExUgs, buf);
        o.append("XU:i:");
        o.append(buf);
        // Xu:i: Mate ungapped tries
        WRITE_SEP();
        itoa10<uint64_t>(prm.nMateUgs, buf);
        o.append("Xu:i:");
        o.append(buf);
    }
    if(print_ye_) {
        // YE:i: Streak of failed DPs at end
        WRITE_SEP();
        itoa10<uint64_t>(prm.nDpFail, buf);
        o.append("YE:i:");
        o.append(buf);
        // Ye:i: Streak of failed ungaps at end
        WRITE_SEP();
        itoa10<uint64_t>(prm.nUgFail, buf);
        o.append("Ye:i:");
        o.append(buf);
    }
    if(print_yl_) {
        // YL:i: Longest streak of failed DPs
        WRITE_SEP();
        itoa10<uint64_t>(prm.nDpFailStreak, buf);
        o.append("YL:i:");
        o.append(buf);
        // Yl:i: Longest streak of failed ungaps
        WRITE_SEP();
        itoa10<uint64_t>(prm.nUgFailStreak, buf);
        o.append("Yl:i:");
        o.append(buf);
    }
    if(print_yu_) {
        // YU:i: Index of last succesful DP
        WRITE_SEP();
        itoa10<uint64_t>(prm.nDpLastSucc, buf);
        o.append("YU:i:");
        o.append(buf);
        // Yu:i: Index of last succesful DP
        WRITE_SEP();
        itoa10<uint64_t>(prm.nUgLastSucc, buf);
        o.append("Yu:i:");
        o.append(buf);
    }
    if(print_xp_) {
        // XP:Z: String describing seed hits
        WRITE_SEP();
        o.append("XP:B:I,");
        itoa10<uint64_t>(prm.nSeedElts, buf);
        o.append(buf);
        o.append(',');
        itoa10<uint64_t>(prm.nSeedEltsFw, buf);
        o.append(buf);
        o.append(',');
        itoa10<uint64_t>(prm.nSeedEltsRc, buf);
        o.append(buf);
        o.append(',');
        itoa10<uint64_t>(prm.seedMean, buf);
        o.append(buf);
        o.append(',');
        itoa10<uint64_t>(prm.seedMedian, buf);
        o.append(buf);
    }
    if(print_yr_) {
        // YR:i: Redundant seed hits
        WRITE_SEP();
        itoa10<uint64_t>(prm.nRedundants, buf);
        o.append("YR:i:");
        o.append(buf);
    }
    if(print_zb_) {
        // ZB:i: Ftab ops for seed alignment
        WRITE_SEP();
        itoa10<uint64_t>(prm.nFtabs, buf);
        o.append("ZB:i:");
        o.append(buf);
    }
    if(print_zr_) {
        // ZR:Z: Redundant path skips in seed alignment
        WRITE_SEP();
        o.append("ZR:Z:");
        itoa10<uint64_t>(prm.nRedSkip, buf); o.append(buf);
        o.append(',');
        itoa10<uint64_t>(prm.nRedFail, buf); o.append(buf);
        o.append(',');
        itoa10<uint64_t>(prm.nRedIns, buf); o.append(buf);
    }
    if(print_zf_) {
        // ZF:i: FM Index ops for seed alignment
        WRITE_SEP();
        itoa10<uint64_t>(prm.nSdFmops, buf);
        o.append("ZF:i:");
        o.append(buf);
        // Zf:i: FM Index ops for offset resolution
        WRITE_SEP();
        itoa10<uint64_t>(prm.nExFmops, buf);
        o.append("Zf:i:");
        o.append(buf);
    }
    if(print_zm_) {
        // ZM:Z: Print FM index op string for best-first search
        WRITE_SEP();
        o.append("ZM:Z:");
        prm.fmString.print(o, buf);
    }
    if(print_zi_) {
        // ZI:i: Seed extend loop iterations
        WRITE_SEP();
        itoa10<uint64_t>(prm.nExIters, buf);
        o.append("ZI:i:");
        o.append(buf);
    }
    if(print_xs_a_) {
        if(rna_strandness_ == RNA_STRANDNESS_UNKNOWN) {
            uint8_t whichsense = res.spliced_whichsense_transcript();
            if(whichsense != SPL_UNKNOWN) {
                WRITE_SEP();
                o.append("XS:A:");
                if(whichsense == SPL_FW || whichsense == SPL_SEMI_FW) {
                    o.append('+');
                } else {
                    assert(whichsense == SPL_RC || whichsense == SPL_SEMI_RC);
                    o.append('-');
                }
            }
        } else {
            WRITE_SEP();
            o.append("XS:A:");
            char strandness = '+';
            if(res.readMate1()) {
                if(res.orient()) {
                    if(rna_strandness_ == RNA_STRANDNESS_R || rna_strandness_ == RNA_STRANDNESS_RF) {
                        strandness = '-';
                    }
                } else {
                    if(rna_strandness_ == RNA_STRANDNESS_F || rna_strandness_ == RNA_STRANDNESS_FR) {
                        strandness = '-';
                    }
                }
            } else {
                assert(res.readMate2());
                assert(rna_strandness_ == RNA_STRANDNESS_FR || rna_strandness_ == RNA_STRANDNESS_RF);
                if(res.orient()) {
                    if(rna_strandness_ == RNA_STRANDNESS_FR) {
                        strandness = '-';
                    }
                } else {
                    if(rna_strandness_ == RNA_STRANDNESS_RF) {
                        strandness = '-';
                    }
                }
            }
            o.append(strandness);
        }
    }
    if(print_nh_) {
        if(flags.alignedPaired()) {
            WRITE_SEP();
            itoa10<uint64_t>(summ.numAlnsPaired(), buf);
            o.append("NH:i:");
            o.append(buf);
        } else if(flags.alignedUnpaired() || flags.alignedUnpairedMate()) {
            WRITE_SEP();
            itoa10<uint64_t>((flags.alignedUnpaired() || flags.readMate1()) ?
                             summ.numAlns1() : summ.numAlns2(), buf);
            o.append("NH:i:");
            o.append(buf);
        }
    }
    
    bool snp_first = true;
    index_t prev_snp_idx = INDEX_MAX;
    size_t len_trimmed = rd.length() - res.trimmed5p(true) - res.trimmed3p(true);
    if(!res.fw()) {
        Edit::invertPoss(const_cast<EList<Edit>&>(res.ned()), len_trimmed, false);
    }
    for(size_t i = 0; i < res.ned().size(); i++) {
        if(res.ned()[i].snpID < altdb->alts().size()) {
            index_t snp_idx = res.ned()[i].snpID;
            assert_lt(snp_idx, altdb->alts().size());
            const ALT<index_t>& snp = altdb->alts()[snp_idx];
            const string& snpID = altdb->altnames()[snp_idx];
            if(snp_idx == prev_snp_idx) continue;
            if(snp_first) {
                WRITE_SEP();
                o.append("Zs:Z:");
            }
            if(!snp_first) o.append(",");
            uint64_t pos = res.ned()[i].pos;
            size_t j = i;
            while(j > 0) {
                if(res.ned()[j-1].snpID < altdb->alts().size()) {
                    const ALT<index_t>& snp2 = altdb->alts()[res.ned()[j-1].snpID];
                    if(snp2.type == ALT_SNP_SGL) {
                        pos -= (res.ned()[j-1].pos + 1);
                    } else if(snp2.type == ALT_SNP_DEL) {
                        pos -= res.ned()[j-1].pos;
                    } else if(snp2.type == ALT_SNP_INS) {
                        pos -= (res.ned()[j-1].pos + snp.len);
                    }
                    break;
                }
                j--;
            }
            itoa10<uint64_t>(pos, buf);
            o.append(buf);
            o.append("|");
            if(snp.type == ALT_SNP_SGL) {
                o.append("S");
            } else if(snp.type == ALT_SNP_DEL) {
                o.append("D");
            } else {
                assert_eq(snp.type, ALT_SNP_INS);
                o.append("I");
            }
            o.append("|");
            o.append(snpID.c_str());
            
            if(snp_first) snp_first = false;
            prev_snp_idx = snp_idx;
        }
    }
    if(!res.fw()) {
        Edit::invertPoss(const_cast<EList<Edit>&>(res.ned()), len_trimmed, false);
    }
    
    if(print_xr_) {
        // Original read string
        o.append("\n");
        printOptFieldNewlineEscapedZ(o, rd.readOrigBuf);
    }
}

/**
 * Print the optional flags to the given string.
 */
template<typename index_t>
void SamConfig<index_t>::printEmptyOptFlags(
                                   BTString& o,               // output buffer
                                   bool first,                // first opt flag printed is first overall?
                                   const Read& rd,            // read
                                   const AlnFlags& flags,     // alignment flags
                                   const AlnSetSumm& summ,    // summary of alignments for this read
                                   const SeedAlSumm& ssm,     // seed alignment summary
                                   const PerReadMetrics& prm, // per-read metrics
                                   const Scoring& sc)         // scoring scheme
const
{
    char buf[1024];
    if(print_yn_) {
        // YN:i: Minimum valid score for this mate
        TAlScore mn = sc.scoreMin.f<TAlScore>(rd.length());
        itoa10<TAlScore>(mn, buf);
        WRITE_SEP();
        o.append("YN:i:");
        o.append(buf);
        // Yn:i: Perfect score for this mate
        TAlScore pe = sc.perfectScore(rd.length());
        itoa10<TAlScore>(pe, buf);
        WRITE_SEP();
        o.append("Yn:i:");
        o.append(buf);
    }
    if(print_zs_) {
        // ZS:i: Pseudo-random seed for read
        itoa10<uint32_t>(rd.seed, buf);
        WRITE_SEP();
        o.append("ZS:i:");
        o.append(buf);
    }
    if(print_yt_) {
        // YT:Z: String representing alignment type
        WRITE_SEP();
        flags.printYT(o);
    }
    if(print_yp_ && flags.partOfPair() && flags.canMax()) {
        // YP:i: Read was repetitive when aligned paired?
        WRITE_SEP();
        flags.printYP(o);
    }
    if(print_ym_ && flags.canMax() && (flags.isMixedMode() || !flags.partOfPair())) {
        // YM:i: Read was repetitive when aligned unpaired?
        WRITE_SEP();
        flags.printYM(o);
    }
    if(print_yf_ && flags.filtered()) {
        // YM:i: Read was repetitive when aligned unpaired?
        first = flags.printYF(o, first) && first;
    }
    if(!rgs_.empty()) {
        WRITE_SEP();
        o.append(rgs_.c_str());
    }
    if(print_xt_) {
        // XT:i: Timing
        WRITE_SEP();
        struct timeval  tv_end;
        struct timezone tz_end;
        gettimeofday(&tv_end, &tz_end);
        size_t total_usecs =
        (tv_end.tv_sec  - prm.tv_beg.tv_sec) * 1000000 +
        (tv_end.tv_usec - prm.tv_beg.tv_usec);
        itoa10<size_t>(total_usecs, buf);
        o.append("XT:i:");
        o.append(buf);
    }
    if(print_xd_) {
        // XD:i: Extend DPs
        WRITE_SEP();
        itoa10<uint64_t>(prm.nExDps, buf);
        o.append("XD:i:");
        o.append(buf);
        // Xd:i: Mate DPs
        WRITE_SEP();
        itoa10<uint64_t>(prm.nMateDps, buf);
        o.append("Xd:i:");
        o.append(buf);
    }
    if(print_xu_) {
        // XU:i: Extend ungapped tries
        WRITE_SEP();
        itoa10<uint64_t>(prm.nExUgs, buf);
        o.append("XU:i:");
        o.append(buf);
        // Xu:i: Mate ungapped tries
        WRITE_SEP();
        itoa10<uint64_t>(prm.nMateUgs, buf);
        o.append("Xu:i:");
        o.append(buf);
    }
    if(print_ye_) {
        // YE:i: Streak of failed DPs at end
        WRITE_SEP();
        itoa10<uint64_t>(prm.nDpFail, buf);
        o.append("YE:i:");
        o.append(buf);
        // Ye:i: Streak of failed ungaps at end
        WRITE_SEP();
        itoa10<uint64_t>(prm.nUgFail, buf);
        o.append("Ye:i:");
        o.append(buf);
    }
    if(print_yl_) {
        // YL:i: Longest streak of failed DPs
        WRITE_SEP();
        itoa10<uint64_t>(prm.nDpFailStreak, buf);
        o.append("YL:i:");
        o.append(buf);
        // Yl:i: Longest streak of failed ungaps
        WRITE_SEP();
        itoa10<uint64_t>(prm.nUgFailStreak, buf);
        o.append("Yl:i:");
        o.append(buf);
    }
    if(print_yu_) {
        // YU:i: Index of last succesful DP
        WRITE_SEP();
        itoa10<uint64_t>(prm.nDpLastSucc, buf);
        o.append("YU:i:");
        o.append(buf);
        // Yu:i: Index of last succesful DP
        WRITE_SEP();
        itoa10<uint64_t>(prm.nUgLastSucc, buf);
        o.append("Yu:i:");
        o.append(buf);
    }
    if(print_xp_) {
        // XP:Z: String describing seed hits
        WRITE_SEP();
        o.append("XP:B:I,");
        itoa10<uint64_t>(prm.nSeedElts, buf);
        o.append(buf);
        o.append(',');
        itoa10<uint64_t>(prm.nSeedEltsFw, buf);
        o.append(buf);
        o.append(',');
        itoa10<uint64_t>(prm.nSeedEltsRc, buf);
        o.append(buf);
        o.append(',');
        itoa10<uint64_t>(prm.seedMean, buf);
        o.append(buf);
        o.append(',');
        itoa10<uint64_t>(prm.seedMedian, buf);
        o.append(buf);
    }
    if(print_yr_) {
        // YR:i: Redundant seed hits
        WRITE_SEP();
        itoa10<uint64_t>(prm.nRedundants, buf);
        o.append("YR:i:");
        o.append(buf);
    }
    if(print_zb_) {
        // ZB:i: Ftab ops for seed alignment
        WRITE_SEP();
        itoa10<uint64_t>(prm.nFtabs, buf);
        o.append("ZB:i:");
        o.append(buf);
    }
    if(print_zr_) {
        // ZR:Z: Redundant path skips in seed alignment
        WRITE_SEP();
        o.append("ZR:Z:");
        itoa10<uint64_t>(prm.nRedSkip, buf); o.append(buf);
        o.append(',');
        itoa10<uint64_t>(prm.nRedFail, buf); o.append(buf);
        o.append(',');
        itoa10<uint64_t>(prm.nRedIns, buf); o.append(buf);
    }
    if(print_zf_) {
        // ZF:i: FM Index ops for seed alignment
        WRITE_SEP();
        itoa10<uint64_t>(prm.nSdFmops, buf);
        o.append("ZF:i:");
        o.append(buf);
        // Zf:i: FM Index ops for offset resolution
        WRITE_SEP();
        itoa10<uint64_t>(prm.nExFmops, buf);
        o.append("Zf:i:");
        o.append(buf);
    }
    if(print_zm_) {
        // ZM:Z: Print FM index op string for best-first search
        WRITE_SEP();
        o.append("ZM:Z:");
        prm.fmString.print(o, buf);
    }
    if(print_zi_) {
        // ZI:i: Seed extend loop iterations
        WRITE_SEP();
        itoa10<uint64_t>(prm.nExIters, buf);
        o.append("ZI:i:");
        o.append(buf);
    }
    if(print_xr_) {
        // Original read string
        o.append("\n");
        printOptFieldNewlineEscapedZ(o, rd.readOrigBuf);
    }
}

#endif /* SAM_H_ */
