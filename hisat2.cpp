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

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <stdexcept>
#include <getopt.h>
#include <math.h>
#include <utility>
#include <limits>
#include "alphabet.h"
#include "assert_helpers.h"
#include "endian_swap.h"
#include "hgfm.h"
#include "rfm.h"
#include "formats.h"
#include "sequence_io.h"
#include "tokenize.h"
#include "aln_sink.h"
#include "pat.h"
#include "threading.h"
#include "ds.h"
#include "aligner_metrics.h"
#include "sam.h"
#include "aligner_seed.h"
#include "splice_site.h"
#include "spliced_aligner.h"
#include "aligner_seed_policy.h"
#include "aligner_sw.h"
#include "aligner_sw_driver.h"
#include "aligner_cache.h"
#include "util.h"
#include "pe.h"
#include "tp.h"
#include "gp.h"
#include "simple_func.h"
#include "presets.h"
#include "opts.h"
#include "outq.h"
#include "repeat_kmer.h"

using namespace std;

MemoryTally gMemTally;

static EList<string> mates1;  // mated reads (first mate)
static EList<string> mates2;  // mated reads (second mate)
static EList<string> mates12; // mated reads (1st/2nd interleaved in 1 file)
static string adjIdxBase;
bool gColor;              // colorspace (not supported)
int gVerbose;             // be talkative
static bool startVerbose; // be talkative at startup
int gQuiet;               // print nothing but the alignments
static int sanityCheck;   // enable expensive sanity checks
static int format;        // default read format is FASTQ
static string origString; // reference text, or filename(s)
static int seed;          // srandom() seed
static int timing;        // whether to report basic timing data
static int metricsIval;   // interval between alignment metrics messages (0 = no messages)
static string metricsFile;// output file to put alignment metrics in
static bool metricsStderr;// output file to put alignment metrics in
static bool metricsPerRead; // report a metrics tuple for every read
static bool allHits;      // for multihits, report just one
static bool showVersion;  // just print version and quit?
static int ipause;        // pause before maching?
static uint32_t qUpto;    // max # of queries to read
int gTrim5;               // amount to trim from 5' end
int gTrim3;               // amount to trim from 3' end
static int offRate;       // keep default offRate
static bool solexaQuals;  // quality strings are solexa quals, not phred, and subtract 64 (not 33)
static bool phred64Quals; // quality chars are phred, but must subtract 64 (not 33)
static bool integerQuals; // quality strings are space-separated strings of integers, not ASCII
static int nthreads;      // number of pthreads operating concurrently
static int outType;       // style of output
static bool noRefNames;   // true -> print reference indexes; not names
static uint32_t khits;    // number of hits per read; >1 is much slower
static uint32_t mhits;    // don't report any hits if there are > mhits
static int partitionSz;   // output a partitioning key in first field
static bool useSpinlock;  // false -> don't use of spinlocks even if they're #defines
static bool fileParallel; // separate threads read separate input files in parallel
static bool useShmem;     // use shared memory to hold the index
static bool useMm;        // use memory-mapped files to hold the index
static bool mmSweep;      // sweep through memory-mapped files immediately after mapping
int gMinInsert;           // minimum insert size
int gMaxInsert;           // maximum insert size
bool gMate1fw;            // -1 mate aligns in fw orientation on fw strand
bool gMate2fw;            // -2 mate aligns in rc orientation on fw strand
bool gFlippedMatesOK;     // allow mates to be in wrong order
bool gDovetailMatesOK;    // allow one mate to extend off the end of the other
bool gContainMatesOK;     // allow one mate to contain the other in PE alignment
bool gOlapMatesOK;        // allow mates to overlap in PE alignment
bool gExpandToFrag;       // incr max frag length to =larger mate len if necessary
bool gReportDiscordant;   // find and report discordant paired-end alignments
bool gReportMixed;        // find and report unpaired alignments for paired reads
static uint32_t cacheLimit;      // ranges w/ size > limit will be cached
static uint32_t cacheSize;       // # words per range cache
static uint32_t skipReads;       // # reads/read pairs to skip
bool gNofw; // don't align fw orientation of read
bool gNorc; // don't align rc orientation of read
static uint32_t fastaContLen;
static uint32_t fastaContFreq;
static bool hadoopOut; // print Hadoop status and summary messages
static bool fuzzy;
static bool fullRef;
static bool samTruncQname; // whether to truncate QNAME to 255 chars
static bool samOmitSecSeqQual; // omit SEQ/QUAL for 2ndary alignments?
static bool samNoUnal; // don't print records for unaligned reads
static bool samNoHead; // don't print any header lines in SAM output
static bool samNoSQ;   // don't print @SQ header lines
static bool sam_print_as;
static bool sam_print_xs;  // XS:i
static bool sam_print_xss; // Xs:i and Ys:i
static bool sam_print_yn;  // YN:i and Yn:i
static bool sam_print_xn;
static bool sam_print_cs;
static bool sam_print_cq;
static bool sam_print_x0;
static bool sam_print_x1;
static bool sam_print_xm;
static bool sam_print_xo;
static bool sam_print_xg;
static bool sam_print_nm;
static bool sam_print_md;
static bool sam_print_yf;
static bool sam_print_yi;
static bool sam_print_ym;
static bool sam_print_yp;
static bool sam_print_yt;
static bool sam_print_ys;
static bool sam_print_zs;
static bool sam_print_xr;
static bool sam_print_xt;
static bool sam_print_xd;
static bool sam_print_xu;
static bool sam_print_yl;
static bool sam_print_ye;
static bool sam_print_yu;
static bool sam_print_xp;
static bool sam_print_yr;
static bool sam_print_zb;
static bool sam_print_zr;
static bool sam_print_zf;
static bool sam_print_zm;
static bool sam_print_zi;
static bool sam_print_zp;
static bool sam_print_zu;
static bool sam_print_xs_a;
static bool sam_print_nh;
static bool bwaSwLike;
static float bwaSwLikeC;
static float bwaSwLikeT;
static bool qcFilter;
static bool sortByScore;      // prioritize alignments to report by score?
bool gReportOverhangs;        // false -> filter out alignments that fall off the end of a reference sequence
static string rgid;           // ID: setting for @RG header line
static string rgs;            // SAM outputs for @RG header line
static string rgs_optflag;    // SAM optional flag to add corresponding to @RG ID
static bool msample;          // whether to report a random alignment when maxed-out via -m/-M
int      gGapBarrier;         // # diags on top/bot only to be entered diagonally
static EList<string> qualities;
static EList<string> qualities1;
static EList<string> qualities2;
static string polstr;         // temporary holder for policy string
static bool  msNoCache;       // true -> disable local cache
static int   bonusMatchType;  // how to reward matches
static int   bonusMatch;      // constant reward if bonusMatchType=constant
static int   penMmcType;      // how to penalize mismatches
static int   penMmcMax;       // max mm penalty
static int   penMmcMin;       // min mm penalty
static int   penScMax;       // max sc penalty
static int   penScMin;       // min sc penalty
static int   penNType;        // how to penalize Ns in the read
static int   penN;            // constant if N pelanty is a constant
static bool  penNCatPair;     // concatenate mates before N filtering?
static bool  localAlign;      // do local alignment in DP steps
static bool  noisyHpolymer;   // set to true if gap penalties should be reduced to be consistent with a sequencer that under- and overcalls homopolymers
static int   penRdGapConst;   // constant cost of extending a gap in the read
static int   penRfGapConst;   // constant cost of extending a gap in the reference
static int   penRdGapLinear;  // coeff of linear term for cost of gap extension in read
static int   penRfGapLinear;  // coeff of linear term for cost of gap extension in ref
static SimpleFunc scoreMin;   // minimum valid score as function of read len
static SimpleFunc nCeil;      // max # Ns allowed as function of read len
static SimpleFunc msIval;     // interval between seeds as function of read len
static double descConsExp;    // how to adjust score minimum as we descent further into index-assisted alignment
static size_t descentLanding; // don't place a search root if it's within this many positions of end
static SimpleFunc descentTotSz;    // maximum space a DescentDriver can use in bytes
static SimpleFunc descentTotFmops; // maximum # FM ops a DescentDriver can perform
static int    multiseedMms;   // mismatches permitted in a multiseed seed
static int    multiseedLen;   // length of multiseed seeds
static size_t multiseedOff;   // offset to begin extracting seeds
static uint32_t seedCacheLocalMB;   // # MB to use for non-shared seed alignment cacheing
static uint32_t seedCacheCurrentMB; // # MB to use for current-read seed hit cacheing
static uint32_t exactCacheCurrentMB; // # MB to use for current-read seed hit cacheing
static size_t maxhalf;        // max width on one side of DP table
static bool seedSumm;         // print summary information about seed hits, not alignments
static bool doUngapped;       // do ungapped alignment
static size_t maxIters;       // stop after this many extend loop iterations
static size_t maxUg;          // stop after this many ungap extends
static size_t maxDp;          // stop after this many DPs
static size_t maxItersIncr;   // amt to add to maxIters for each -k > 1
static size_t maxEeStreak;    // stop after this many end-to-end fails in a row
static size_t maxUgStreak;    // stop after this many ungap fails in a row
static size_t maxDpStreak;    // stop after this many dp fails in a row
static size_t maxStreakIncr;  // amt to add to streak for each -k > 1
static size_t maxMateStreak;  // stop seed range after this many mate-find fails
static bool doExtend;         // extend seed hits
static bool enable8;          // use 8-bit SSE where possible?
static size_t cminlen;        // longer reads use checkpointing
static size_t cpow2;          // checkpoint interval log2
static bool doTri;            // do triangular mini-fills?
static string defaultPreset;  // default preset; applied immediately
static bool ignoreQuals;      // all mms incur same penalty, regardless of qual
static string wrapper;        // type of wrapper script, so we can print correct usage
static EList<string> queries; // list of query files
static string outfile;        // write SAM output to this file
static int mapqv;             // MAPQ calculation version
static int tighten;           // -M tighten mode (0=none, 1=best, 2=secbest+1)
static bool doExactUpFront;   // do exact search up front if seeds seem good enough
static bool do1mmUpFront;     // do 1mm search up front if seeds seem good enough
static size_t do1mmMinLen;    // length below which we disable 1mm e2e search
static int seedBoostThresh;   // if average non-zero position has more than this many elements
static size_t maxSeeds;       // maximum number of seeds allowed
static size_t nSeedRounds;    // # seed rounds
static bool reorder;          // true -> reorder SAM recs in -p mode
static float sampleFrac;      // only align random fraction of input reads
static bool arbitraryRandom;  // pseudo-randoms no longer a function of read properties
static bool bowtie2p5;
static bool useTempSpliceSite;
static int penCanSplice;
static int penNoncanSplice;
static int penConflictSplice;
static SimpleFunc penCanIntronLen;
static SimpleFunc penNoncanIntronLen;
static size_t minIntronLen;
static size_t maxIntronLen;
static string knownSpliceSiteInfile;  //
static string novelSpliceSiteInfile;  //
static string novelSpliceSiteOutfile; //
static bool secondary;
static bool no_spliced_alignment;
static int rna_strandness; //
static bool splicesite_db_only; //

static bool anchorStop;
static bool pseudogeneStop;
static bool tranMapOnly; // transcriptome mapping only
static bool tranAssm;    // alignments selected for downstream transcript assembly such as StringTie and Cufflinks
static string tranAssm_program;
static bool avoid_pseudogene;

#ifdef USE_SRA
static EList<string> sra_accs;
#endif

static string bt2index;      // read Bowtie 2 index from files with this prefix
static EList<pair<int, string> > extra_opts;
static size_t extra_opts_cur;

static EList<uint64_t> thread_rids;
static MUTEX_T         thread_rids_mutex;
static uint64_t        thread_rids_mindist;

static bool rmChrName;  // remove "chr" from reference names (e.g., chr18 to 18)
static bool addChrName; // add "chr" to reference names (e.g., 18 to chr18)

static size_t max_alts_tried;
static bool use_haplotype;
static bool enable_codis;

static bool templateLenAdjustment;
static string alignSumFile; // write alignment summary stat. to this file
static bool newAlignSummary;

static int bowtie2_dp; // Bowtie2's dynamic programming alignment (0: no dynamic programming, 1: conditional dynamic programming, and 2: uncoditional dynamic programming)
static bool fast;           // --fast
static bool sensitive;      // --sensitive
static bool very_sensitive; // --very-sensitive

static bool repeat;
static bool use_repeat_index;
static EList<size_t> readLens;


#define DMAX std::numeric_limits<double>::max()

static void resetOptions() {
	mates1.clear();
	mates2.clear();
	mates12.clear();
	adjIdxBase	            = "";
	gColor                  = false;
	gVerbose                = 0;
	startVerbose			= 0;
	gQuiet					= false;
	sanityCheck				= 0;  // enable expensive sanity checks
	format					= FASTQ; // default read format is FASTQ
	origString				= ""; // reference text, or filename(s)
	seed					= 0; // srandom() seed
	timing					= 0; // whether to report basic timing data
	metricsIval				= 1; // interval between alignment metrics messages (0 = no messages)
	metricsFile             = ""; // output file to put alignment metrics in
	metricsStderr           = false; // print metrics to stderr (in addition to --metrics-file if it's specified
	metricsPerRead          = false; // report a metrics tuple for every read?
	allHits					= false; // for multihits, report just one
	showVersion				= false; // just print version and quit?
	ipause					= 0; // pause before maching?
	qUpto					= 0xffffffff; // max # of queries to read
	gTrim5					= 0; // amount to trim from 5' end
	gTrim3					= 0; // amount to trim from 3' end
	offRate					= -1; // keep default offRate
	solexaQuals				= false; // quality strings are solexa quals, not phred, and subtract 64 (not 33)
	phred64Quals			= false; // quality chars are phred, but must subtract 64 (not 33)
	integerQuals			= false; // quality strings are space-separated strings of integers, not ASCII
	nthreads				= 1;     // number of pthreads operating concurrently
	outType					= OUTPUT_SAM;  // style of output
	noRefNames				= false; // true -> print reference indexes; not names
	khits					= 10;    // number of hits per read; >1 is much slower
	mhits					= 0;     // stop after finding this many alignments+1
	partitionSz				= 0;     // output a partitioning key in first field
	useSpinlock				= true;  // false -> don't use of spinlocks even if they're #defines
	fileParallel			= false; // separate threads read separate input files in parallel
	useShmem				= false; // use shared memory to hold the index
	useMm					= false; // use memory-mapped files to hold the index
	mmSweep					= false; // sweep through memory-mapped files immediately after mapping
	gMinInsert				= 0;     // minimum insert size
	gMaxInsert				= 1000;   // maximum insert size
	gMate1fw				= true;  // -1 mate aligns in fw orientation on fw strand
	gMate2fw				= false; // -2 mate aligns in rc orientation on fw strand
	gFlippedMatesOK         = false; // allow mates to be in wrong order
	gDovetailMatesOK        = false; // allow one mate to extend off the end of the other
	gContainMatesOK         = true;  // allow one mate to contain the other in PE alignment
	gOlapMatesOK            = true;  // allow mates to overlap in PE alignment
	gExpandToFrag           = true;  // incr max frag length to =larger mate len if necessary
	gReportDiscordant       = true;  // find and report discordant paired-end alignments
	gReportMixed            = true;  // find and report unpaired alignments for paired reads

	cacheLimit				= 5;     // ranges w/ size > limit will be cached
	cacheSize				= 0;     // # words per range cache
	skipReads				= 0;     // # reads/read pairs to skip
	gNofw					= false; // don't align fw orientation of read
	gNorc					= false; // don't align rc orientation of read
	fastaContLen			= 0;
	fastaContFreq			= 0;
	hadoopOut				= false; // print Hadoop status and summary messages
	fuzzy					= false; // reads will have alternate basecalls w/ qualities
	fullRef					= false; // print entire reference name instead of just up to 1st space
	samTruncQname           = true;  // whether to truncate QNAME to 255 chars
	samOmitSecSeqQual       = false; // omit SEQ/QUAL for 2ndary alignments?
	samNoUnal               = false; // omit SAM records for unaligned reads
	samNoHead				= false; // don't print any header lines in SAM output
	samNoSQ					= false; // don't print @SQ header lines
	sam_print_as            = true;
	sam_print_xs            = true;
	sam_print_xss           = false; // Xs:i and Ys:i
	sam_print_yn            = false; // YN:i and Yn:i
	sam_print_xn            = true;
	sam_print_cs            = false;
	sam_print_cq            = false;
	sam_print_x0            = true;
	sam_print_x1            = true;
	sam_print_xm            = true;
	sam_print_xo            = true;
	sam_print_xg            = true;
	sam_print_nm            = true;
	sam_print_md            = true;
	sam_print_yf            = true;
	sam_print_yi            = false;
	sam_print_ym            = false;
	sam_print_yp            = false;
	sam_print_yt            = true;
	sam_print_ys            = true;
	sam_print_zs            = false;
	sam_print_xr            = false;
	sam_print_xt            = false;
	sam_print_xd            = false;
	sam_print_xu            = false;
	sam_print_yl            = false;
	sam_print_ye            = false;
	sam_print_yu            = false;
	sam_print_xp            = false;
	sam_print_yr            = false;
	sam_print_zb            = false;
	sam_print_zr            = false;
	sam_print_zf            = false;
	sam_print_zm            = false;
	sam_print_zi            = false;
	sam_print_zp            = false;
	sam_print_zu            = false;
    sam_print_xs_a          = true;
    sam_print_nh            = true;
	bwaSwLike               = false;
	bwaSwLikeC              = 5.5f;
	bwaSwLikeT              = 20.0f;
	qcFilter                = false; // don't believe upstream qc by default
	sortByScore             = true;  // prioritize alignments to report by score?
	rgid					= "";    // SAM outputs for @RG header line
	rgs						= "";    // SAM outputs for @RG header line
	rgs_optflag				= "";    // SAM optional flag to add corresponding to @RG ID
	msample				    = true;
	gGapBarrier				= 4;     // disallow gaps within this many chars of either end of alignment
	qualities.clear();
	qualities1.clear();
	qualities2.clear();
	polstr.clear();
	msNoCache       = true; // true -> disable local cache
	bonusMatchType  = DEFAULT_MATCH_BONUS_TYPE;
	bonusMatch      = DEFAULT_MATCH_BONUS;
	penMmcType      = DEFAULT_MM_PENALTY_TYPE;
	penMmcMax       = DEFAULT_MM_PENALTY_MAX;
	penMmcMin       = DEFAULT_MM_PENALTY_MIN;
    penScMax        = DEFAULT_SC_PENALTY_MAX;
    penScMin        = DEFAULT_SC_PENALTY_MIN;
	penNType        = DEFAULT_N_PENALTY_TYPE;
	penN            = DEFAULT_N_PENALTY;
	penNCatPair     = DEFAULT_N_CAT_PAIR; // concatenate mates before N filtering?
	localAlign      = false;     // do local alignment in DP steps
	noisyHpolymer   = false;
	penRdGapConst   = DEFAULT_READ_GAP_CONST;
	penRfGapConst   = DEFAULT_REF_GAP_CONST;
	penRdGapLinear  = DEFAULT_READ_GAP_LINEAR;
	penRfGapLinear  = DEFAULT_REF_GAP_LINEAR;
    scoreMin.init  (SIMPLE_FUNC_LINEAR, 0.0f, -0.2f);
    // scoreMin.init  (SIMPLE_FUNC_CONST, -18, 0);
	nCeil.init     (SIMPLE_FUNC_LINEAR, 0.0f, DMAX, 2.0f, 0.1f);
	msIval.init    (SIMPLE_FUNC_LINEAR, 1.0f, DMAX, DEFAULT_IVAL_B, DEFAULT_IVAL_A);
	descConsExp     = 2.0;
	descentLanding  = 20;
	descentTotSz.init(SIMPLE_FUNC_LINEAR, 1024.0, DMAX, 0.0, 1024.0);
	descentTotFmops.init(SIMPLE_FUNC_LINEAR, 100.0, DMAX, 0.0, 10.0);
	multiseedMms    = DEFAULT_SEEDMMS;
	multiseedLen    = DEFAULT_SEEDLEN;
	multiseedOff    = 0;
	seedCacheLocalMB   = 32; // # MB to use for non-shared seed alignment cacheing
	seedCacheCurrentMB = 20; // # MB to use for current-read seed hit cacheing
	exactCacheCurrentMB = 20; // # MB to use for current-read seed hit cacheing
	maxhalf            = 15; // max width on one side of DP table
	seedSumm           = false; // print summary information about seed hits, not alignments
	doUngapped         = true;  // do ungapped alignment
	maxIters           = 400;   // max iterations of extend loop
	maxUg              = 300;   // stop after this many ungap extends
	maxDp              = 300;   // stop after this many dp extends
	maxItersIncr       = 20;    // amt to add to maxIters for each -k > 1
	maxEeStreak        = 15;    // stop after this many end-to-end fails in a row
	maxUgStreak        = 15;    // stop after this many ungap fails in a row
	maxDpStreak        = 15;    // stop after this many dp fails in a row
	maxStreakIncr      = 10;    // amt to add to streak for each -k > 1
	maxMateStreak      = 10;    // in PE: abort seed range after N mate-find fails
	doExtend           = true;  // do seed extensions
	enable8            = true;  // use 8-bit SSE where possible?
	cminlen            = 2000;  // longer reads use checkpointing
	cpow2              = 4;     // checkpoint interval log2
	doTri              = false; // do triangular mini-fills?
	defaultPreset      = "sensitive%LOCAL%"; // default preset; applied immediately
	extra_opts.clear();
	extra_opts_cur = 0;
	bt2index.clear();        // read Bowtie 2 index from files with this prefix
	ignoreQuals = false;     // all mms incur same penalty, regardless of qual
	wrapper.clear();         // type of wrapper script, so we can print correct usage
	queries.clear();         // list of query files
	outfile.clear();         // write SAM output to this file
	mapqv = 2;               // MAPQ calculation version
	tighten = 3;             // -M tightening mode
	doExactUpFront = true;   // do exact search up front if seeds seem good enough
	do1mmUpFront = true;     // do 1mm search up front if seeds seem good enough
	seedBoostThresh = 300;   // if average non-zero position has more than this many elements
	nSeedRounds = 2;         // # rounds of seed searches to do for repetitive reads
    maxSeeds = 0;            // maximum number of seeds allowed
	do1mmMinLen = 60;        // length below which we disable 1mm search
	reorder = false;         // reorder SAM records with -p > 1
	sampleFrac = 1.1f;       // align all reads
	arbitraryRandom = false; // let pseudo-random seeds be a function of read properties
	bowtie2p5 = false;
    useTempSpliceSite = true;
    penCanSplice = 0;
    penNoncanSplice = 12;
    penConflictSplice = 1000000;
    penCanIntronLen.init(SIMPLE_FUNC_LOG, -8, 1);
    penNoncanIntronLen.init(SIMPLE_FUNC_LOG, -8, 1);
    minIntronLen = 20;
    maxIntronLen = 500000;
    knownSpliceSiteInfile = "";
    novelSpliceSiteInfile = "";
    novelSpliceSiteOutfile = "";
    secondary = false;       // allow secondary alignments
    no_spliced_alignment = false;
    rna_strandness = RNA_STRANDNESS_UNKNOWN;
    splicesite_db_only = false;
    anchorStop = true;
    pseudogeneStop = true;
    tranMapOnly = false;
    tranAssm = false;
    tranAssm_program = "";
    avoid_pseudogene = false;
    
#ifdef USE_SRA
    sra_accs.clear();
#endif
    
    rmChrName = false;
    addChrName = false;
    
    max_alts_tried = 16;
    use_haplotype = false;
    enable_codis = false;
    
    templateLenAdjustment = true;
    alignSumFile = "";
    newAlignSummary = false;
    
    bowtie2_dp = 0; // disable Bowtie2's dynamic programming alignment
    fast = false;
    sensitive = false;
    very_sensitive = false;
    
    repeat = false; // true iff alignments to repeat sequences are directly reported.
    use_repeat_index = true;
    readLens.clear();
}

static const char *short_options = "fF:qbzhcu:rv:s:aP:t3:5:w:p:k:M:1:2:I:X:CQ:N:i:L:U:x:S:g:O:D:R:";

static struct option long_options[] = {
	{(char*)"verbose",      no_argument,       0,            ARG_VERBOSE},
	{(char*)"startverbose", no_argument,       0,            ARG_STARTVERBOSE},
	{(char*)"quiet",        no_argument,       0,            ARG_QUIET},
	{(char*)"sanity",       no_argument,       0,            ARG_SANITY},
	{(char*)"pause",        no_argument,       &ipause,      1},
	{(char*)"orig",         required_argument, 0,            ARG_ORIG},
	{(char*)"all",          no_argument,       0,            'a'},
	{(char*)"solexa-quals", no_argument,       0,            ARG_SOLEXA_QUALS},
	{(char*)"integer-quals",no_argument,       0,            ARG_INTEGER_QUALS},
	{(char*)"int-quals",    no_argument,       0,            ARG_INTEGER_QUALS},
	{(char*)"metrics",      required_argument, 0,            ARG_METRIC_IVAL},
	{(char*)"metrics-file", required_argument, 0,            ARG_METRIC_FILE},
	{(char*)"metrics-stderr",no_argument,      0,            ARG_METRIC_STDERR},
	{(char*)"metrics-per-read", no_argument,   0,            ARG_METRIC_PER_READ},
	{(char*)"met-read",     no_argument,       0,            ARG_METRIC_PER_READ},
	{(char*)"met",          required_argument, 0,            ARG_METRIC_IVAL},
	{(char*)"met-file",     required_argument, 0,            ARG_METRIC_FILE},
	{(char*)"met-stderr",   no_argument,       0,            ARG_METRIC_STDERR},
	{(char*)"time",         no_argument,       0,            't'},
	{(char*)"trim3",        required_argument, 0,            '3'},
	{(char*)"trim5",        required_argument, 0,            '5'},
	{(char*)"seed",         required_argument, 0,            ARG_SEED},
	{(char*)"qupto",        required_argument, 0,            'u'},
	{(char*)"upto",         required_argument, 0,            'u'},
	{(char*)"version",      no_argument,       0,            ARG_VERSION},
	{(char*)"filepar",      no_argument,       0,            ARG_FILEPAR},
	{(char*)"help",         no_argument,       0,            'h'},
	{(char*)"threads",      required_argument, 0,            'p'},
	{(char*)"khits",        required_argument, 0,            'k'},
	{(char*)"minins",       required_argument, 0,            'I'},
	{(char*)"maxins",       required_argument, 0,            'X'},
	{(char*)"quals",        required_argument, 0,            'Q'},
	{(char*)"Q1",           required_argument, 0,            ARG_QUALS1},
	{(char*)"Q2",           required_argument, 0,            ARG_QUALS2},
	{(char*)"refidx",       no_argument,       0,            ARG_REFIDX},
	{(char*)"partition",    required_argument, 0,            ARG_PARTITION},
	{(char*)"ff",           no_argument,       0,            ARG_FF},
	{(char*)"fr",           no_argument,       0,            ARG_FR},
	{(char*)"rf",           no_argument,       0,            ARG_RF},
	{(char*)"cachelim",     required_argument, 0,            ARG_CACHE_LIM},
	{(char*)"cachesz",      required_argument, 0,            ARG_CACHE_SZ},
	{(char*)"nofw",         no_argument,       0,            ARG_NO_FW},
	{(char*)"norc",         no_argument,       0,            ARG_NO_RC},
	{(char*)"skip",         required_argument, 0,            's'},
	{(char*)"12",           required_argument, 0,            ARG_ONETWO},
	{(char*)"tab5",         required_argument, 0,            ARG_TAB5},
	{(char*)"tab6",         required_argument, 0,            ARG_TAB6},
	{(char*)"phred33-quals", no_argument,      0,            ARG_PHRED33},
	{(char*)"phred64-quals", no_argument,      0,            ARG_PHRED64},
	{(char*)"phred33",       no_argument,      0,            ARG_PHRED33},
	{(char*)"phred64",      no_argument,       0,            ARG_PHRED64},
	{(char*)"solexa1.3-quals", no_argument,    0,            ARG_PHRED64},
	{(char*)"mm",           no_argument,       0,            ARG_MM},
	{(char*)"shmem",        no_argument,       0,            ARG_SHMEM},
	{(char*)"mmsweep",      no_argument,       0,            ARG_MMSWEEP},
	{(char*)"hadoopout",    no_argument,       0,            ARG_HADOOPOUT},
	{(char*)"fuzzy",        no_argument,       0,            ARG_FUZZY},
	{(char*)"fullref",      no_argument,       0,            ARG_FULLREF},
	{(char*)"usage",        no_argument,       0,            ARG_USAGE},
	{(char*)"sam-no-qname-trunc", no_argument, 0,            ARG_SAM_NO_QNAME_TRUNC},
	{(char*)"sam-omit-sec-seq", no_argument,   0,            ARG_SAM_OMIT_SEC_SEQ},
	{(char*)"omit-sec-seq", no_argument,       0,            ARG_SAM_OMIT_SEC_SEQ},
	{(char*)"sam-no-head",  no_argument,       0,            ARG_SAM_NOHEAD},
	{(char*)"sam-nohead",   no_argument,       0,            ARG_SAM_NOHEAD},
	{(char*)"sam-noHD",     no_argument,       0,            ARG_SAM_NOHEAD},
	{(char*)"sam-no-hd",    no_argument,       0,            ARG_SAM_NOHEAD},
	{(char*)"sam-nosq",     no_argument,       0,            ARG_SAM_NOSQ},
	{(char*)"sam-no-sq",    no_argument,       0,            ARG_SAM_NOSQ},
	{(char*)"sam-noSQ",     no_argument,       0,            ARG_SAM_NOSQ},
	{(char*)"no-head",      no_argument,       0,            ARG_SAM_NOHEAD},
	{(char*)"no-hd",        no_argument,       0,            ARG_SAM_NOHEAD},
	{(char*)"no-sq",        no_argument,       0,            ARG_SAM_NOSQ},
	{(char*)"no-HD",        no_argument,       0,            ARG_SAM_NOHEAD},
	{(char*)"no-SQ",        no_argument,       0,            ARG_SAM_NOSQ},
	{(char*)"no-unal",      no_argument,       0,            ARG_SAM_NO_UNAL},
	{(char*)"color",        no_argument,       0,            'C'},
	{(char*)"sam-RG",       required_argument, 0,            ARG_SAM_RG},
	{(char*)"sam-rg",       required_argument, 0,            ARG_SAM_RG},
	{(char*)"sam-rg-id",    required_argument, 0,            ARG_SAM_RGID},
	{(char*)"RG",           required_argument, 0,            ARG_SAM_RG},
	{(char*)"rg",           required_argument, 0,            ARG_SAM_RG},
	{(char*)"rg-id",        required_argument, 0,            ARG_SAM_RGID},
	{(char*)"snpphred",     required_argument, 0,            ARG_SNPPHRED},
	{(char*)"snpfrac",      required_argument, 0,            ARG_SNPFRAC},
	{(char*)"gbar",         required_argument, 0,            ARG_GAP_BAR},
	{(char*)"qseq",         no_argument,       0,            ARG_QSEQ},
	{(char*)"policy",       required_argument, 0,            ARG_ALIGN_POLICY},
	{(char*)"preset",       required_argument, 0,            'P'},
	{(char*)"seed-summ",    no_argument,       0,            ARG_SEED_SUMM},
	{(char*)"seed-summary", no_argument,       0,            ARG_SEED_SUMM},
	{(char*)"overhang",     no_argument,       0,            ARG_OVERHANG},
	{(char*)"no-cache",     no_argument,       0,            ARG_NO_CACHE},
	{(char*)"cache",        no_argument,       0,            ARG_USE_CACHE},
	{(char*)"454",          no_argument,       0,            ARG_NOISY_HPOLY},
	{(char*)"ion-torrent",  no_argument,       0,            ARG_NOISY_HPOLY},
	{(char*)"no-mixed",     no_argument,       0,            ARG_NO_MIXED},
	{(char*)"no-discordant",no_argument,       0,            ARG_NO_DISCORDANT},
	// {(char*)"local",        no_argument,       0,            ARG_LOCAL},
	{(char*)"end-to-end",   no_argument,       0,            ARG_END_TO_END},
	{(char*)"ungapped",     no_argument,       0,            ARG_UNGAPPED},
	{(char*)"no-ungapped",  no_argument,       0,            ARG_UNGAPPED_NO},
	{(char*)"sse8",         no_argument,       0,            ARG_SSE8},
	{(char*)"no-sse8",      no_argument,       0,            ARG_SSE8_NO},
	{(char*)"scan-narrowed",no_argument,       0,            ARG_SCAN_NARROWED},
	{(char*)"qc-filter",    no_argument,       0,            ARG_QC_FILTER},
	{(char*)"bwa-sw-like",  no_argument,       0,            ARG_BWA_SW_LIKE},
	{(char*)"multiseed",        required_argument, 0,        ARG_MULTISEED_IVAL},
	{(char*)"ma",               required_argument, 0,        ARG_SCORE_MA},
	{(char*)"mp",               required_argument, 0,        ARG_SCORE_MMP},
    {(char*)"sp",               required_argument, 0,        ARG_SCORE_SCP},
    {(char*)"no-softclip",      no_argument,       0,        ARG_NO_SOFTCLIP},
	{(char*)"np",               required_argument, 0,        ARG_SCORE_NP},
	{(char*)"rdg",              required_argument, 0,        ARG_SCORE_RDG},
	{(char*)"rfg",              required_argument, 0,        ARG_SCORE_RFG},
	{(char*)"score-min",        required_argument, 0,        ARG_SCORE_MIN},
	{(char*)"min-score",        required_argument, 0,        ARG_SCORE_MIN},
	{(char*)"n-ceil",           required_argument, 0,        ARG_N_CEIL},
	{(char*)"dpad",             required_argument, 0,        ARG_DPAD},
	{(char*)"mapq-print-inputs",no_argument,       0,        ARG_SAM_PRINT_YI},
	{(char*)"very-fast",        no_argument,       0,        ARG_PRESET_VERY_FAST},
	{(char*)"fast",             no_argument,       0,        ARG_PRESET_FAST},
	{(char*)"sensitive",        no_argument,       0,        ARG_PRESET_SENSITIVE},
	{(char*)"very-sensitive",   no_argument,       0,        ARG_PRESET_VERY_SENSITIVE},
	// {(char*)"very-fast-local",      no_argument,   0,        ARG_PRESET_VERY_FAST_LOCAL},
	// {(char*)"fast-local",           no_argument,   0,        ARG_PRESET_FAST_LOCAL},
	// {(char*)"sensitive-local",      no_argument,   0,        ARG_PRESET_SENSITIVE_LOCAL},
	// {(char*)"very-sensitive-local", no_argument,   0,        ARG_PRESET_VERY_SENSITIVE_LOCAL},
	{(char*)"no-score-priority",no_argument,       0,        ARG_NO_SCORE_PRIORITY},
	{(char*)"seedlen",          required_argument, 0,        'L'},
	{(char*)"seedmms",          required_argument, 0,        'N'},
	{(char*)"seedival",         required_argument, 0,        'i'},
	{(char*)"ignore-quals",     no_argument,       0,        ARG_IGNORE_QUALS},
	{(char*)"index",            required_argument, 0,        'x'},
	{(char*)"arg-desc",         no_argument,       0,        ARG_DESC},
	{(char*)"wrapper",          required_argument, 0,        ARG_WRAPPER},
	{(char*)"unpaired",         required_argument, 0,        'U'},
	{(char*)"output",           required_argument, 0,        'S'},
	{(char*)"mapq-v",           required_argument, 0,        ARG_MAPQ_V},
	{(char*)"dovetail",         no_argument,       0,        ARG_DOVETAIL},
	{(char*)"no-dovetail",      no_argument,       0,        ARG_NO_DOVETAIL},
	{(char*)"contain",          no_argument,       0,        ARG_CONTAIN},
	{(char*)"no-contain",       no_argument,       0,        ARG_NO_CONTAIN},
	{(char*)"overlap",          no_argument,       0,        ARG_OVERLAP},
	{(char*)"no-overlap",       no_argument,       0,        ARG_NO_OVERLAP},
	{(char*)"tighten",          required_argument, 0,        ARG_TIGHTEN},
	{(char*)"exact-upfront",    no_argument,       0,        ARG_EXACT_UPFRONT},
	{(char*)"1mm-upfront",      no_argument,       0,        ARG_1MM_UPFRONT},
	{(char*)"no-exact-upfront", no_argument,       0,        ARG_EXACT_UPFRONT_NO},
	{(char*)"no-1mm-upfront",   no_argument,       0,        ARG_1MM_UPFRONT_NO},
	{(char*)"1mm-minlen",       required_argument, 0,        ARG_1MM_MINLEN},
	{(char*)"seed-off",         required_argument, 0,        'O'},
	{(char*)"seed-boost",       required_argument, 0,        ARG_SEED_BOOST_THRESH},
    {(char*)"max-seeds",        required_argument, 0,        ARG_MAX_SEEDS},
	{(char*)"read-times",       no_argument,       0,        ARG_READ_TIMES},
	{(char*)"show-rand-seed",   no_argument,       0,        ARG_SHOW_RAND_SEED},
	{(char*)"dp-fail-streak",   required_argument, 0,        ARG_DP_FAIL_STREAK_THRESH},
	{(char*)"ee-fail-streak",   required_argument, 0,        ARG_EE_FAIL_STREAK_THRESH},
	{(char*)"ug-fail-streak",   required_argument, 0,        ARG_UG_FAIL_STREAK_THRESH},
	{(char*)"fail-streak",      required_argument, 0,        'D'},
	{(char*)"dp-fails",         required_argument, 0,        ARG_DP_FAIL_THRESH},
	{(char*)"ug-fails",         required_argument, 0,        ARG_UG_FAIL_THRESH},
	{(char*)"extends",          required_argument, 0,        ARG_EXTEND_ITERS},
	{(char*)"no-extend",        no_argument,       0,        ARG_NO_EXTEND},
	{(char*)"mapq-extra",       no_argument,       0,        ARG_MAPQ_EX},
	{(char*)"seed-rounds",      required_argument, 0,        'R'},
	{(char*)"reorder",          no_argument,       0,        ARG_REORDER},
	{(char*)"passthrough",      no_argument,       0,        ARG_READ_PASSTHRU},
	{(char*)"sample",           required_argument, 0,        ARG_SAMPLE},
	{(char*)"cp-min",           required_argument, 0,        ARG_CP_MIN},
	{(char*)"cp-ival",          required_argument, 0,        ARG_CP_IVAL},
	{(char*)"tri",              no_argument,       0,        ARG_TRI},
	{(char*)"nondeterministic", no_argument,       0,        ARG_NON_DETERMINISTIC},
	{(char*)"non-deterministic", no_argument,      0,        ARG_NON_DETERMINISTIC},
	// {(char*)"local-seed-cache-sz", required_argument, 0,     ARG_LOCAL_SEED_CACHE_SZ},
	{(char*)"seed-cache-sz",       required_argument, 0,     ARG_CURRENT_SEED_CACHE_SZ},
	{(char*)"no-unal",          no_argument,       0,        ARG_SAM_NO_UNAL},
	{(char*)"test-25",          no_argument,       0,        ARG_TEST_25},
	// TODO: following should be a function of read length?
	{(char*)"desc-kb",          required_argument, 0,        ARG_DESC_KB},
	{(char*)"desc-landing",     required_argument, 0,        ARG_DESC_LANDING},
	{(char*)"desc-exp",         required_argument, 0,        ARG_DESC_EXP},
	{(char*)"desc-fmops",       required_argument, 0,        ARG_DESC_FMOPS},
    {(char*)"no-temp-splicesite",  no_argument,    0,        ARG_NO_TEMPSPLICESITE},
    {(char*)"pen-cansplice",  required_argument,   0,        ARG_PEN_CANSPLICE},
    {(char*)"pen-noncansplice",  required_argument, 0,       ARG_PEN_NONCANSPLICE},
    {(char*)"pen-conflictsplice",  required_argument, 0,     ARG_PEN_CONFLICTSPLICE},
    {(char*)"pen-intronlen",  required_argument,   0,        ARG_PEN_CANINTRONLEN},
    {(char*)"pen-canintronlen",  required_argument, 0,       ARG_PEN_CANINTRONLEN},
    {(char*)"pen-noncanintronlen",  required_argument, 0,    ARG_PEN_NONCANINTRONLEN},
    {(char*)"min-intronlen",  required_argument,   0,        ARG_MIN_INTRONLEN},
    {(char*)"max-intronlen",  required_argument,   0,        ARG_MAX_INTRONLEN},
    {(char*)"known-splicesite-infile",       required_argument, 0,        ARG_KNOWN_SPLICESITE_INFILE},
    {(char*)"novel-splicesite-infile",       required_argument, 0,        ARG_NOVEL_SPLICESITE_INFILE},
    {(char*)"novel-splicesite-outfile",      required_argument, 0,        ARG_NOVEL_SPLICESITE_OUTFILE},
    {(char*)"secondary",        no_argument,       0,        ARG_SECONDARY},
    {(char*)"no-spliced-alignment",   no_argument, 0,        ARG_NO_SPLICED_ALIGNMENT},
    {(char*)"rna-strandness",   required_argument, 0,        ARG_RNA_STRANDNESS},
    {(char*)"splicesite-db-only",   no_argument,   0,        ARG_SPLICESITE_DB_ONLY},
    {(char*)"no-anchorstop",   no_argument,        0,        ARG_NO_ANCHORSTOP},
    {(char*)"transcriptome-mapping-only", no_argument, 0,    ARG_TRANSCRIPTOME_MAPPING_ONLY},
    {(char*)"tmo",             no_argument,        0,        ARG_TRANSCRIPTOME_MAPPING_ONLY},
    {(char*)"downstream-transcriptome-assembly",   no_argument, 0,        ARG_TRANSCRIPTOME_ASSEMBLY},
    {(char*)"dta",             no_argument,        0,        ARG_TRANSCRIPTOME_ASSEMBLY},
    {(char*)"dta-cufflinks",   no_argument,        0,        ARG_TRANSCRIPTOME_ASSEMBLY_CUFFLINKS},
    {(char*)"avoid-pseudogene",no_argument,        0,        ARG_AVOID_PSEUDOGENE},
    {(char*)"no-templatelen-adjustment",    no_argument,        0,        ARG_NO_TEMPLATELEN_ADJUSTMENT},
#ifdef USE_SRA
    {(char*)"sra-acc",         required_argument,  0,        ARG_SRA_ACC},
#endif
    {(char*)"remove-chrname",  no_argument,        0,        ARG_REMOVE_CHRNAME},
    {(char*)"add-chrname",     no_argument,        0,        ARG_ADD_CHRNAME},
    {(char*)"max-altstried",   required_argument,  0,        ARG_MAX_ALTSTRIED},
    {(char*)"haplotype",       no_argument,        0,        ARG_HAPLOTYPE},
    {(char*)"enable-codis",    no_argument,        0,        ARG_CODIS},
    {(char*)"summary-file",    required_argument,  0,        ARG_SUMMARY_FILE},
    {(char*)"new-summary",     no_argument,        0,        ARG_NEW_SUMMARY},
    {(char*)"enable-dp",       no_argument,        0,        ARG_DP},
    {(char*)"bowtie2-dp",      required_argument,  0,        ARG_DP},
    {(char*)"repeat",          no_argument,        0,        ARG_REPEAT},
    {(char*)"no-repeat-index", no_argument,        0,        ARG_NO_REPEAT_INDEX},
    {(char*)"read-lengths",    required_argument,  0,        ARG_READ_LENGTHS},
	{(char*)0, 0, 0, 0} // terminator
};

/**
 * Print out a concise description of what options are taken and whether they
 * take an argument.
 */
static void printArgDesc(ostream& out) {
	// struct option {
	//   const char *name;
	//   int has_arg;
	//   int *flag;
	//   int val;
	// };
	size_t i = 0;
	while(long_options[i].name != 0) {
		out << long_options[i].name << "\t"
		    << (long_options[i].has_arg == no_argument ? 0 : 1)
		    << endl;
		i++;
	}
	size_t solen = strlen(short_options);
	for(i = 0; i < solen; i++) {
		// Has an option?  Does if next char is :
		if(i == solen-1) {
			assert_neq(':', short_options[i]);
			cout << (char)short_options[i] << "\t" << 0 << endl;
		} else {
			if(short_options[i+1] == ':') {
				// Option with argument
				cout << (char)short_options[i] << "\t" << 1 << endl;
				i++; // skip the ':'
			} else {
				// Option with no argument
				cout << (char)short_options[i] << "\t" << 0 << endl;
			}
		}
	}
}

/**
 * Print a summary usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "HISAT2 version " << string(HISAT2_VERSION).c_str() << " by Daehwan Kim (infphilo@gmail.com, www.ccb.jhu.edu/people/infphilo)" << endl;
	string tool_name = "hisat2-align";
	if(wrapper == "basic-0") {
		tool_name = "hisat2";
	}
	out << "Usage: " << endl
#ifdef USE_SRA
	    << "  " << tool_name.c_str() << " [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]" << endl
#else
    << "  " << tool_name.c_str() << " [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]" << endl
#endif
	    << endl
		<<     "  <ht2-idx>  Index filename prefix (minus trailing .X." << gfm_ext << ")." << endl
	    <<     "  <m1>       Files with #1 mates, paired with files in <m2>." << endl;
	if(wrapper == "basic-0") {
		out << "             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)." << endl;
	}
	out <<     "  <m2>       Files with #2 mates, paired with files in <m1>." << endl;
	if(wrapper == "basic-0") {
		out << "             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)." << endl;
	}
	out <<     "  <r>        Files with unpaired reads." << endl;
	if(wrapper == "basic-0") {
		out << "             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)." << endl;
	}
#ifdef USE_SRA
    out <<     "  <SRA accession number>        Comma-separated list of SRA accession numbers, e.g. --sra-acc SRR353653,SRR353654." << endl;
#endif
	out <<     "  <sam>      File for SAM output (default: stdout)" << endl
	    << endl
	    << "  <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be" << endl
		<< "  specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'." << endl
		// Wrapper script should write <bam> line next
		<< endl
	    << "Options (defaults in parentheses):" << endl
		<< endl
	    << " Input:" << endl
	    << "  -q                 query input files are FASTQ .fq/.fastq (default)" << endl
	    << "  --qseq             query input files are in Illumina's qseq format" << endl
	    << "  -f                 query input files are (multi-)FASTA .fa/.mfa" << endl
	    << "  -r                 query input files are raw one-sequence-per-line" << endl
	    << "  -c                 <m1>, <m2>, <r> are sequences themselves, not files" << endl
	    << "  -s/--skip <int>    skip the first <int> reads/pairs in the input (none)" << endl
	    << "  -u/--upto <int>    stop after first <int> reads/pairs (no limit)" << endl
	    << "  -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)" << endl
	    << "  -3/--trim3 <int>   trim <int> bases from 3'/right end of reads (0)" << endl
	    << "  --phred33          qualities are Phred+33 (default)" << endl
	    << "  --phred64          qualities are Phred+64" << endl
	    << "  --int-quals        qualities encoded as space-delimited integers" << endl
#ifdef USE_SRA
        << "  --sra-acc          SRA accession ID" << endl
#endif
		<< endl

	    << " Presets:                 Same as:" << endl
		// << "  For --end-to-end:" << endl
		// << "   --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50" << endl
		// << "   --fast                 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50" << endl
		// << "   --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)" << endl
		// << "   --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50" << endl
        << "   --fast                 --no-repeat-index" << endl
        << "   --sensitive            --bowtie2-dp 1 -k 30 --score-min L,0,-0.5" << endl
        << "   --very-sensitive       --bowtie2-dp 2 -k 50 --score-min L,0,-1" << endl
		<< endl
	    << " Alignment:" << endl
		//<< "  -N <int>           max # mismatches in seed alignment; can be 0 or 1 (0)" << endl
		//<< "  -L <int>           length of seed substrings; must be >3, <32 (22)" << endl
		//<< "  -i <func>          interval between seed substrings w/r/t read len (S,1,1.15)" << endl
        << "  --bowtie2-dp <int> use Bowtie2's dynamic programming alignment algorithm (0) - 0: no dynamic programming, 1: conditional dynamic programming, and 2: unconditional dynamic programming (slowest)" << endl
		<< "  --n-ceil <func>    func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)" << endl
		//<< "  --dpad <int>       include <int> extra ref chars on sides of DP table (15)" << endl
		//<< "  --gbar <int>       disallow gaps within <int> nucs of read extremes (4)" << endl
		<< "  --ignore-quals     treat all quality values as 30 on Phred scale (off)" << endl
	    << "  --nofw             do not align forward (original) version of read (off)" << endl
	    << "  --norc             do not align reverse-complement version of read (off)" << endl
        << "  --no-repeat-index  do not use repeat index" << endl
		<< endl
        << " Spliced Alignment:" << endl
        << "  --pen-cansplice <int>              penalty for a canonical splice site (0)" << endl
        << "  --pen-noncansplice <int>           penalty for a non-canonical splice site (12)" << endl
        // << "  --pen-conflictsplice <int>         penalty for conflicting splice sites (1000000)" << endl
        << "  --pen-canintronlen <func>          penalty for long introns (G,-8,1) with canonical splice sites" << endl
        << "  --pen-noncanintronlen <func>       penalty for long introns (G,-8,1) with noncanonical splice sites" << endl
        << "  --min-intronlen <int>              minimum intron length (20)" << endl
        << "  --max-intronlen <int>              maximum intron length (500000)" << endl
        << "  --known-splicesite-infile <path>   provide a list of known splice sites" << endl
        << "  --novel-splicesite-outfile <path>  report a list of splice sites" << endl
        << "  --novel-splicesite-infile <path>   provide a list of novel splice sites" << endl
        << "  --no-temp-splicesite               disable the use of splice sites found" << endl
        << "  --no-spliced-alignment             disable spliced alignment" << endl
        << "  --rna-strandness <string>          specify strand-specific information (unstranded)" << endl
        << "  --tmo                              reports only those alignments within known transcriptome" << endl
        << "  --dta                              reports alignments tailored for transcript assemblers" << endl
        << "  --dta-cufflinks                    reports alignments tailored specifically for cufflinks" << endl
        << "  --avoid-pseudogene                 tries to avoid aligning reads to pseudogenes (experimental option)" << endl
        << "  --no-templatelen-adjustment        disables template length adjustment for RNA-seq reads" << endl
        << endl
		<< " Scoring:" << endl
		//<< "  --ma <int>         match bonus (0 for --end-to-end, 2 for --local) " << endl
		<< "  --mp <int>,<int>   max and min penalties for mismatch; lower qual = lower penalty <6,2>" << endl
        << "  --sp <int>,<int>   max and min penalties for soft-clipping; lower qual = lower penalty <2,1>" << endl
        << "  --no-softclip      no soft-clipping" << endl
		<< "  --np <int>         penalty for non-A/C/G/Ts in read/ref (1)" << endl
		<< "  --rdg <int>,<int>  read gap open, extend penalties (5,3)" << endl
		<< "  --rfg <int>,<int>  reference gap open, extend penalties (5,3)" << endl
		<< "  --score-min <func> min acceptable alignment score w/r/t read length" << endl
		<< "                     (L,0.0,-0.2)" << endl
		<< endl
	    << " Reporting:" << endl
	    << "  -k <int>           It searches for at most <int> distinct, primary alignments for each read. Primary alignments mean " << endl
        << "                     alignments whose alignment score is equal to or higher than any other alignments. The search terminates " << endl
        << "                     when it cannot find more distinct valid alignments, or when it finds <int>, whichever happens first. " << endl
        << "                     The alignment score for a paired-end alignment equals the sum of the alignment scores of " << endl 
        << "                     the individual mates. Each reported read or pair alignment beyond the first has the SAM ‘secondary’ bit " << endl
        << "                     (which equals 256) set in its FLAGS field. For reads that have more than <int> distinct, " << endl
        << "                     valid alignments, hisat2 does not guarantee that the <int> alignments reported are the best possible " << endl 
        << "                     in terms of alignment score. Default: 5 (linear index) or 10 (graph index)." << endl
        << "                     Note: HISAT2 is not designed with large values for -k in mind, and when aligning reads to long, " << endl
        << "                     repetitive genomes, large -k could make alignment much slower." << endl
        << "  --max-seeds <int>  HISAT2, like other aligners, uses seed-and-extend approaches. HISAT2 tries to extend seeds to " << endl
        << "                     full-length alignments. In HISAT2, --max-seeds is used to control the maximum number of seeds that " << endl
        << "                     will be extended. For DNA-read alignment (--no-spliced-alignment), HISAT2 extends up to these many seeds" << endl
        << "                     and skips the rest of the seeds. For RNA-read alignment, HISAT2 skips extending seeds and reports " << endl
        << "                     no alignments if the number of seeds is larger than the number specified with the option, " << endl
        << "                     to be compatible with previous versions of HISAT2. Large values for --max-seeds may improve alignment " << endl
        << "                     sensitivity, but HISAT2 is not designed with large values for --max-seeds in mind, and when aligning " << endl
        << "                     reads to long, repetitive genomes, large --max-seeds could make alignment much slower. " << endl
        << "                     The default value is the maximum of 5 and the value that comes with -k times 2." << endl
        << "  -a/--all           HISAT2 reports all alignments it can find. Using the option is equivalent to using both --max-seeds " << endl
        << "                     and -k with the maximum value that a 64-bit signed integer can represent (9,223,372,036,854,775,807)." << endl 
        << "  --repeat           report alignments to repeat sequences directly" << endl
		<< endl
	    //<< " Effort:" << endl
	    //<< "  -D <int>           give up extending after <int> failed extends in a row (15)" << endl
	    //<< "  -R <int>           for reads w/ repetitive seeds, try <int> sets of seeds (2)" << endl
		//<< endl
		<< " Paired-end:" << endl
	    << "  -I/--minins <int>  minimum fragment length (0), only valid with --no-spliced-alignment" << endl
	    << "  -X/--maxins <int>  maximum fragment length (500), only valid with --no-spliced-alignment" << endl
	    << "  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)" << endl
		<< "  --no-mixed         suppress unpaired alignments for paired reads" << endl
		<< "  --no-discordant    suppress discordant alignments for paired reads" << endl
		<< endl
	    << " Output:" << endl;
	//if(wrapper == "basic-0") {
	//	out << "  --bam              output directly to BAM (by piping through 'samtools view')" << endl;
	//}
	out << "  -t/--time          print wall-clock time taken by search phases" << endl;
	if(wrapper == "basic-0") {
	out << "  --un <path>           write unpaired reads that didn't align to <path>" << endl
	    << "  --al <path>           write unpaired reads that aligned at least once to <path>" << endl
	    << "  --un-conc <path>      write pairs that didn't align concordantly to <path>" << endl
	    << "  --al-conc <path>      write pairs that aligned concordantly at least once to <path>" << endl
	    << "  (Note: for --un, --al, --un-conc, or --al-conc, add '-gz' to the option name, e.g." << endl
		<< "  --un-gz <path>, to gzip compress output, or add '-bz2' to bzip2 compress output.)" << endl;
	}
    out << "  --summary-file <path> print alignment summary to this file." << endl
        << "  --new-summary         print alignment summary in a new style, which is more machine-friendly." << endl
        << "  --quiet               print nothing to stderr except serious errors" << endl
	//  << "  --refidx              refer to ref. seqs by 0-based index rather than name" << endl
		<< "  --met-file <path>     send metrics to file at <path> (off)" << endl
		<< "  --met-stderr          send metrics to stderr (off)" << endl
		<< "  --met <int>           report internal counters & metrics every <int> secs (1)" << endl
	// Following is supported in the wrapper instead
	//  << "  --no-unal             suppress SAM records for unaligned reads" << endl
	    << "  --no-head             suppress header lines, i.e. lines starting with @" << endl
	    << "  --no-sq               suppress @SQ header lines" << endl
	    << "  --rg-id <text>        set read group id, reflected in @RG line and RG:Z: opt field" << endl
	    << "  --rg <text>           add <text> (\"lab:value\") to @RG line of SAM header." << endl
	    << "                        Note: @RG line only printed when --rg-id is set." << endl
	    << "  --omit-sec-seq        put '*' in SEQ and QUAL fields for secondary alignments." << endl
		<< endl
	    << " Performance:" << endl
	    << "  -o/--offrate <int> override offrate of index; must be >= index's offrate" << endl
	    << "  -p/--threads <int> number of alignment threads to launch (1)" << endl
	    << "  --reorder          force SAM output order to match order of input reads" << endl
#ifdef BOWTIE_MM
	    << "  --mm               use memory-mapped I/O for index; many 'hisat2's can share" << endl
#endif
#ifdef BOWTIE_SHARED_MEM
		//<< "  --shmem            use shared mem for index; many 'hisat2's can share" << endl
#endif
		<< endl
	    << " Other:" << endl;
	if(wrapper == "basic-0") {
    out << "  --temp-directory   set the directory for holding temporary files (/tmp)" << endl;
  }
		out << "  --qc-filter        filter out reads that are bad according to QSEQ filter" << endl
	    << "  --seed <int>       seed for random number generator (0)" << endl
	    << "  --non-deterministic seed rand. gen. arbitrarily instead of using read attributes" << endl
        << "  --remove-chrname   remove 'chr' from reference names in alignment" << endl
        << "  --add-chrname      add 'chr' to reference names in alignment " << endl
	//  << "  --verbose          verbose output for debugging" << endl
	    << "  --version          print version information and quit" << endl
	    << "  -h/--help          print this usage message" << endl
	    ;
	if(wrapper.empty()) {
		cerr << endl
		     << "*** Warning ***" << endl
			 << "'hisat2-align' was run directly.  It is recommended that you run the wrapper script 'hisat2' instead." << endl
			 << endl;
	}
}

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static int parseInt(int lower, int upper, const char *errmsg, const char *arg) {
	long l;
	char *endPtr= NULL;
	l = strtol(arg, &endPtr, 10);
	if (endPtr != NULL) {
		if (l < lower || l > upper) {
			cerr << errmsg << endl;
			printUsage(cerr);
			throw 1;
		}
		return (int32_t)l;
	}
	cerr << errmsg << endl;
	printUsage(cerr);
	throw 1;
	return -1;
}

/**
 * Upper is maximum int by default.
 */
static int parseInt(int lower, const char *errmsg, const char *arg) {
	return parseInt(lower, std::numeric_limits<int>::max(), errmsg, arg);
}

/**
 * Parse a T string 'str'.
 */
template<typename T>
T parse(const char *s) {
	T tmp;
	stringstream ss(s);
	ss >> tmp;
	return tmp;
}

/**
 * Parse a pair of Ts from a string, 'str', delimited with 'delim'.
 */
template<typename T>
pair<T, T> parsePair(const char *str, char delim) {
	string s(str);
	EList<string> ss;
	tokenize(s, delim, ss);
	pair<T, T> ret;
	ret.first = parse<T>(ss[0].c_str());
	ret.second = parse<T>(ss[1].c_str());
	return ret;
}

/**
 * Parse a pair of Ts from a string, 'str', delimited with 'delim'.
 */
template<typename T>
void parseTuple(const char *str, char delim, EList<T>& ret) {
	string s(str);
	EList<string> ss;
	tokenize(s, delim, ss);
	for(size_t i = 0; i < ss.size(); i++) {
		ret.push_back(parse<T>(ss[i].c_str()));
	}
}

static string applyPreset(const string& sorig, Presets& presets) {
	string s = sorig;
	size_t found = s.find("%LOCAL%");
	if(found != string::npos) {
		s.replace(found, strlen("%LOCAL%"), localAlign ? "-local" : "");
	}
	if(gVerbose) {
		cerr << "Applying preset: '" << s.c_str() << "' using preset menu '"
			 << presets.name() << "'" << endl;
	}
	string pol;
	presets.apply(s, pol, extra_opts);
	return pol;
}

static bool saw_M;
static bool saw_a;
static bool saw_k;
static EList<string> presetList;

/**
 * TODO: Argument parsing is very, very flawed.  The biggest problem is that
 * there are two separate worlds of arguments, the ones set via polstr, and
 * the ones set directly in variables.  This makes for nasty interactions,
 * e.g., with the -M option being resolved at an awkward time relative to
 * the -k and -a options.
 */
static void parseOption(int next_option, const char *arg) {
	switch (next_option) {
		case ARG_TEST_25: bowtie2p5 = true; break;
		case ARG_DESC_KB: descentTotSz = SimpleFunc::parse(arg, 0.0, 1024.0, 1024.0, DMAX); break;
		case ARG_DESC_FMOPS: descentTotFmops = SimpleFunc::parse(arg, 0.0, 10.0, 100.0, DMAX); break;
		case ARG_DESC_LANDING: descentLanding = parse<int>(arg); break;
		case ARG_DESC_EXP: {
			descConsExp = parse<double>(arg);
			if(descConsExp < 0.0) {
				cerr << "Error: --desc-exp must be greater than or equal to 0" << endl;
				throw 1;
			}
			break;
		}
		case '1': tokenize(arg, ",", mates1); break;
		case '2': tokenize(arg, ",", mates2); break;
		case ARG_ONETWO: tokenize(arg, ",", mates12); format = TAB_MATE5; break;
		case ARG_TAB5:   tokenize(arg, ",", mates12); format = TAB_MATE5; break;
		case ARG_TAB6:   tokenize(arg, ",", mates12); format = TAB_MATE6; break;
		case 'f': format = FASTA; break;
		case 'F': {
			format = FASTA_CONT;
			pair<uint32_t, uint32_t> p = parsePair<uint32_t>(arg, ',');
			fastaContLen = p.first;
			fastaContFreq = p.second;
			break;
		}
		case ARG_BWA_SW_LIKE: {
			bwaSwLikeC = 5.5f;
			bwaSwLikeT = 30;
			bwaSwLike = true;
			localAlign = true;
			// -a INT   Score of a match [1]
			// -b INT   Mismatch penalty [3]
			// -q INT   Gap open penalty [5]
			// -r INT   Gap extension penalty. The penalty for a contiguous
			//          gap of size k is q+k*r. [2] 
			polstr += ";MA=1;MMP=C3;RDG=5,2;RFG=5,2";
			break;
		}
		case 'q': format = FASTQ; break;
		case 'r': format = RAW; break;
		case 'c': format = CMDLINE; break;
		case ARG_QSEQ: format = QSEQ; break;
		case 'C': {
			cerr << "Error: -C specified but Bowtie 2 does not support colorspace input." << endl;
			throw 1;
			break;
		}
		case 'I':
			gMinInsert = parseInt(0, "-I arg must be positive", arg);
			break;
		case 'X':
			gMaxInsert = parseInt(1, "-X arg must be at least 1", arg);
			break;
		case ARG_NO_DISCORDANT: gReportDiscordant = false; break;
		case ARG_NO_MIXED: gReportMixed = false; break;
		case 's':
			skipReads = (uint32_t)parseInt(0, "-s arg must be positive", arg);
			break;
		case ARG_FF: gMate1fw = true;  gMate2fw = true;  break;
		case ARG_RF: gMate1fw = false; gMate2fw = true;  break;
		case ARG_FR: gMate1fw = true;  gMate2fw = false; break;
		case ARG_SHMEM: useShmem = true; break;
		case ARG_SEED_SUMM: seedSumm = true; break;
		case ARG_MM: {
#ifdef BOWTIE_MM
			useMm = true;
			break;
#else
			cerr << "Memory-mapped I/O mode is disabled because bowtie was not compiled with" << endl
				 << "BOWTIE_MM defined.  Memory-mapped I/O is not supported under Windows.  If you" << endl
				 << "would like to use memory-mapped I/O on a platform that supports it, please" << endl
				 << "refrain from specifying BOWTIE_MM=0 when compiling Bowtie." << endl;
			throw 1;
#endif
		}
		case ARG_MMSWEEP: mmSweep = true; break;
		case ARG_HADOOPOUT: hadoopOut = true; break;
		case ARG_SOLEXA_QUALS: solexaQuals = true; break;
		case ARG_INTEGER_QUALS: integerQuals = true; break;
		case ARG_PHRED64: phred64Quals = true; break;
		case ARG_PHRED33: solexaQuals = false; phred64Quals = false; break;
		case ARG_OVERHANG: gReportOverhangs = true; break;
		case ARG_NO_CACHE: msNoCache = true; break;
		case ARG_USE_CACHE: msNoCache = false; break;
		case ARG_LOCAL_SEED_CACHE_SZ:
			seedCacheLocalMB = (uint32_t)parseInt(1, "--local-seed-cache-sz arg must be at least 1", arg);
			break;
		case ARG_CURRENT_SEED_CACHE_SZ:
			seedCacheCurrentMB = (uint32_t)parseInt(1, "--seed-cache-sz arg must be at least 1", arg);
			break;
		case ARG_REFIDX: noRefNames = true; break;
		case ARG_FUZZY: fuzzy = true; break;
		case ARG_FULLREF: fullRef = true; break;
		case ARG_GAP_BAR:
			gGapBarrier = parseInt(1, "--gbar must be no less than 1", arg);
			break;
		case ARG_SEED:
			seed = parseInt(0, "--seed arg must be at least 0", arg);
			break;
		case ARG_NON_DETERMINISTIC:
			arbitraryRandom = true;
			break;
		case 'u':
			qUpto = (uint32_t)parseInt(1, "-u/--qupto arg must be at least 1", arg);
			break;
		case 'Q':
			tokenize(arg, ",", qualities);
			integerQuals = true;
			break;
		case ARG_QUALS1:
			tokenize(arg, ",", qualities1);
			integerQuals = true;
			break;
		case ARG_QUALS2:
			tokenize(arg, ",", qualities2);
			integerQuals = true;
			break;
		case ARG_CACHE_LIM:
			cacheLimit = (uint32_t)parseInt(1, "--cachelim arg must be at least 1", arg);
			break;
		case ARG_CACHE_SZ:
			cacheSize = (uint32_t)parseInt(1, "--cachesz arg must be at least 1", arg);
			cacheSize *= (1024 * 1024); // convert from MB to B
			break;
		case ARG_WRAPPER: wrapper = arg; break;
		case 'p':
			nthreads = parseInt(1, "-p/--threads arg must be at least 1", arg);
			break;
		case ARG_FILEPAR:
			fileParallel = true;
			break;
		case '3': gTrim3 = parseInt(0, "-3/--trim3 arg must be at least 0", arg); break;
		case '5': gTrim5 = parseInt(0, "-5/--trim5 arg must be at least 0", arg); break;
		case 'h': printUsage(cout); throw 0; break;
		case ARG_USAGE: printUsage(cout); throw 0; break;
		//
		// NOTE that unlike in Bowtie 1, -M, -a and -k are mutually
		// exclusive here.
		//
		case 'M': {
			msample = true;
			mhits = parse<uint32_t>(arg);
			if(saw_a || saw_k) {
				cerr << "Warning: -M, -k and -a are mutually exclusive. "
					 << "-M will override" << endl;
				khits = 1;
			}
			assert_eq(1, khits);
			saw_M = true;
			cerr << "Warning: -M is deprecated.  Use -D and -R to adjust " <<
			        "effort instead." << endl;
			break;
		}
		case ARG_EXTEND_ITERS: {
			maxIters = parse<size_t>(arg);
			break;
		}
		case ARG_NO_EXTEND: {
			doExtend = false;
			break;
		}
		case 'R': { polstr += ";ROUNDS="; polstr += arg; break; }
		case 'D': { polstr += ";DPS=";    polstr += arg; break; }
		case ARG_DP_MATE_STREAK_THRESH: {
			maxMateStreak = parse<size_t>(arg);
			break;
		}
		case ARG_DP_FAIL_STREAK_THRESH: {
			maxDpStreak = parse<size_t>(arg);
			break;
		}
		case ARG_EE_FAIL_STREAK_THRESH: {
			maxEeStreak = parse<size_t>(arg);
			break;
		}
		case ARG_UG_FAIL_STREAK_THRESH: {
			maxUgStreak = parse<size_t>(arg);
			break;
		}
		case ARG_DP_FAIL_THRESH: {
			maxDp = parse<size_t>(arg);
			break;
		}
		case ARG_UG_FAIL_THRESH: {
			maxUg = parse<size_t>(arg);
			break;
		}
        case ARG_MAX_SEEDS: {
            maxSeeds = parse<size_t>(arg);
            break;
        }
		case ARG_SEED_BOOST_THRESH: {
			seedBoostThresh = parse<int>(arg);
			break;
		}
		case 'a': {
			msample = false;
			allHits = true;
			mhits = 0; // disable -M
			if(saw_M || saw_k) {
				cerr << "Warning: -M, -k and -a are mutually exclusive. "
					 << "-a will override" << endl;
			}
			saw_a = true;
			break;
		}
		case 'k': {
			msample = false;
			khits = (uint32_t)parseInt(1, "-k arg must be at least 1", arg);
			mhits = 0; // disable -M
			if(saw_M || saw_a) {
				cerr << "Warning: -M, -k and -a are mutually exclusive. "
					 << "-k will override" << endl;
			}
			saw_k = true;
			break;
		}
		case ARG_VERBOSE: gVerbose = 1; break;
		case ARG_STARTVERBOSE: startVerbose = true; break;
		case ARG_QUIET: gQuiet = true; break;
		case ARG_SANITY: sanityCheck = true; break;
		case 't': timing = true; break;
		case ARG_METRIC_IVAL: {
			metricsIval = parseInt(1, "--metrics arg must be at least 1", arg);
			break;
		}
		case ARG_METRIC_FILE: metricsFile = arg; break;
		case ARG_METRIC_STDERR: metricsStderr = true; break;
		case ARG_METRIC_PER_READ: metricsPerRead = true; break;
		case ARG_NO_FW: gNofw = true; break;
		case ARG_NO_RC: gNorc = true; break;
		case ARG_SAM_NO_QNAME_TRUNC: samTruncQname = false; break;
		case ARG_SAM_OMIT_SEC_SEQ: samOmitSecSeqQual = true; break;
		case ARG_SAM_NO_UNAL: samNoUnal = true; break;
		case ARG_SAM_NOHEAD: samNoHead = true; break;
		case ARG_SAM_NOSQ: samNoSQ = true; break;
		case ARG_SAM_PRINT_YI: sam_print_yi = true; break;
		case ARG_REORDER: reorder = true; break;
		case ARG_MAPQ_EX: {
			sam_print_zp = true;
			sam_print_zu = true;
			sam_print_xp = true;
			sam_print_xss = true;
			sam_print_yn = true;
			break;
		}
		case ARG_SHOW_RAND_SEED: {
			sam_print_zs = true;
			break;
		}
		case ARG_SAMPLE:
			sampleFrac = parse<float>(arg);
			break;
		case ARG_CP_MIN:
			cminlen = parse<size_t>(arg);
			break;
		case ARG_CP_IVAL:
			cpow2 = parse<size_t>(arg);
			break;
		case ARG_TRI:
			doTri = true;
			break;
		case ARG_READ_PASSTHRU: {
			sam_print_xr = true;
			break;
		}
		case ARG_READ_TIMES: {
			sam_print_xt = true;
			sam_print_xd = true;
			sam_print_xu = true;
			sam_print_yl = true;
			sam_print_ye = true;
			sam_print_yu = true;
			sam_print_yr = true;
			sam_print_zb = true;
			sam_print_zr = true;
			sam_print_zf = true;
			sam_print_zm = true;
			sam_print_zi = true;
			break;
		}
		case ARG_SAM_RG: {
			string argstr = arg;
			if(argstr.substr(0, 3) == "ID:") {
				rgid = "\t";
				rgid += argstr;
				rgs_optflag = "RG:Z:" + argstr.substr(3);
			} else {
				rgs += '\t';
				rgs += argstr;
			}
			break;
		}
		case ARG_SAM_RGID: {
			string argstr = arg;
			rgid = "\t";
			rgid = "\tID:" + argstr;
			rgs_optflag = "RG:Z:" + argstr;
			break;
		}
		case ARG_PARTITION: partitionSz = parse<int>(arg); break;
		case ARG_DPAD:
			maxhalf = parseInt(0, "--dpad must be no less than 0", arg);
			break;
		case ARG_ORIG:
			if(arg == NULL || strlen(arg) == 0) {
				cerr << "--orig arg must be followed by a string" << endl;
				printUsage(cerr);
				throw 1;
			}
			origString = arg;
			break;
		case ARG_LOCAL: localAlign = true; break;
		case ARG_END_TO_END: localAlign = false; break;
		case ARG_SSE8: enable8 = true; break;
		case ARG_SSE8_NO: enable8 = false; break;
		case ARG_UNGAPPED: doUngapped = true; break;
		case ARG_UNGAPPED_NO: doUngapped = false; break;
		// case ARG_NO_DOVETAIL: gDovetailMatesOK = false; break;
		// case ARG_NO_CONTAIN:  gContainMatesOK  = false; break;
		// case ARG_NO_OVERLAP:  gOlapMatesOK     = false; break;
		// case ARG_DOVETAIL:    gDovetailMatesOK = true;  break;
		// case ARG_CONTAIN:     gContainMatesOK  = true;  break;
		// case ARG_OVERLAP:     gOlapMatesOK     = true;  break;
		case ARG_QC_FILTER: qcFilter = true; break;
		case ARG_NO_SCORE_PRIORITY: sortByScore = false; break;
		case ARG_IGNORE_QUALS: ignoreQuals = true; break;
		case ARG_MAPQ_V: mapqv = parse<int>(arg); break;
		case ARG_TIGHTEN: tighten = parse<int>(arg); break;
		case ARG_EXACT_UPFRONT:    doExactUpFront = true; break;
		case ARG_1MM_UPFRONT:      do1mmUpFront   = true; break;
		case ARG_EXACT_UPFRONT_NO: doExactUpFront = false; break;
		case ARG_1MM_UPFRONT_NO:   do1mmUpFront   = false; break;
		case ARG_1MM_MINLEN:       do1mmMinLen = parse<size_t>(arg); break;
		case ARG_NOISY_HPOLY: noisyHpolymer = true; break;
		case 'x': bt2index = arg; break;
		case ARG_PRESET_VERY_FAST_LOCAL: localAlign = true;
		case ARG_PRESET_VERY_FAST: {
			presetList.push_back("very-fast%LOCAL%"); break;
		}
		case ARG_PRESET_FAST_LOCAL: localAlign = true;
		case ARG_PRESET_FAST: {
            fast = true;
			presetList.push_back("fast%LOCAL%"); break;
		}
		case ARG_PRESET_SENSITIVE_LOCAL: localAlign = true;
		case ARG_PRESET_SENSITIVE: {
            sensitive = true;
			presetList.push_back("sensitive%LOCAL%"); break;
		}
		case ARG_PRESET_VERY_SENSITIVE_LOCAL: localAlign = true;
		case ARG_PRESET_VERY_SENSITIVE: {
            very_sensitive = true;
			presetList.push_back("very-sensitive%LOCAL%"); break;
		}
		case 'P': { presetList.push_back(arg); break; }
		case ARG_ALIGN_POLICY: {
			if(strlen(arg) > 0) {
				polstr += ";"; polstr += arg;
			}
			break;
		}
		case 'N': { polstr += ";SEED="; polstr += arg; break; }
		case 'L': {
			int64_t len = parse<size_t>(arg);
			if(len < 0) {
				cerr << "Error: -L argument must be >= 0; was " << arg << endl;
				throw 1;
			}
			if(len > 32) {
				cerr << "Error: -L argument must be <= 32; was" << arg << endl;
				throw 1;
			}
			polstr += ";SEEDLEN="; polstr += arg; break;
		}
		case 'O':
			multiseedOff = parse<size_t>(arg);
			break;
		case 'i': {
			EList<string> args;
			tokenize(arg, ",", args);
			if(args.size() > 3 || args.size() == 0) {
				cerr << "Error: expected 3 or fewer comma-separated "
					 << "arguments to -i option, got "
					 << args.size() << endl;
				throw 1;
			}
			// Interval-settings arguments
			polstr += (";IVAL=" + args[0]); // Function type
			if(args.size() > 1) {
				polstr += ("," + args[1]);  // Constant term
			}
			if(args.size() > 2) {
				polstr += ("," + args[2]);  // Coefficient
			}
			break;
		}
		case ARG_MULTISEED_IVAL: {
			polstr += ";";
			// Split argument by comma
			EList<string> args;
			tokenize(arg, ",", args);
			if(args.size() > 5 || args.size() == 0) {
				cerr << "Error: expected 5 or fewer comma-separated "
					 << "arguments to --multiseed option, got "
					 << args.size() << endl;
				throw 1;
			}
			// Seed mm and length arguments
			polstr += "SEED=";
			polstr += (args[0]); // # mismatches
			if(args.size() >  1) polstr += ("," + args[ 1]); // length
			if(args.size() >  2) polstr += (";IVAL=" + args[2]); // Func type
			if(args.size() >  3) polstr += ("," + args[ 3]); // Constant term
			if(args.size() >  4) polstr += ("," + args[ 4]); // Coefficient
			break;
		}
		case ARG_N_CEIL: {
			// Split argument by comma
			EList<string> args;
			tokenize(arg, ",", args);
			if(args.size() > 3) {
				cerr << "Error: expected 3 or fewer comma-separated "
					 << "arguments to --n-ceil option, got "
					 << args.size() << endl;
				throw 1;
			}
			if(args.size() == 0) {
				cerr << "Error: expected at least one argument to --n-ceil option" << endl;
				throw 1;
			}
			polstr += ";NCEIL=";
			if(args.size() == 3) {
				polstr += (args[0] + "," + args[1] + "," + args[2]);
			} else {
                if(args.size() == 1) {
                    polstr += ("C," + args[0]);
                } else {
					polstr += (args[0] + "," + args[1]);
				}
			}
			break;
		}
		case ARG_SCORE_MA:  polstr += ";MA=";    polstr += arg; break;
		case ARG_SCORE_MMP: {
			EList<string> args;
			tokenize(arg, ",", args);
			if(args.size() > 2 || args.size() == 0) {
				cerr << "Error: expected 1 or 2 comma-separated "
					 << "arguments to --mp option, got " << args.size() << endl;
				throw 1;
			}
			if(args.size() >= 1) {
				polstr += ";MMP=Q,";
				polstr += args[0];
				if(args.size() >= 2) {
					polstr += ",";
					polstr += args[1];
				}
			}
			break;
		}
        case ARG_SCORE_SCP: {
            EList<string> args;
            tokenize(arg, ",", args);
            if(args.size() > 2 || args.size() == 0) {
                cerr << "Error: expected 1 or 2 comma-separated "
                << "arguments to --sp option, got " << args.size() << endl;
                throw 1;
            }
            if(args.size() >= 1) {
                polstr += ";SCP=Q,";
                polstr += args[0];
                if(args.size() >= 2) {
                    polstr += ",";
                    polstr += args[1];
                }
            }
            break;
        }
        case ARG_NO_SOFTCLIP: {
            ostringstream convert;
            convert << std::numeric_limits<int>::max();
            polstr += ";SCP=Q,";
            polstr += convert.str();
            polstr += ",";
            polstr += convert.str();
            break;
        }
		case ARG_SCORE_NP:  polstr += ";NP=C";   polstr += arg; break;
		case ARG_SCORE_RDG: polstr += ";RDG=";   polstr += arg; break;
		case ARG_SCORE_RFG: polstr += ";RFG=";   polstr += arg; break;
		case ARG_SCORE_MIN: {
			polstr += ";";
			EList<string> args;
			tokenize(arg, ",", args);
			if(args.size() > 3 && args.size() == 0) {
				cerr << "Error: expected 3 or fewer comma-separated "
					 << "arguments to --n-ceil option, got "
					 << args.size() << endl;
				throw 1;
			}
			polstr += ("MIN=" + args[0]);
			if(args.size() > 1) {
				polstr += ("," + args[1]);
			}
			if(args.size() > 2) {
				polstr += ("," + args[2]);
			}
			break;
		}
		case ARG_DESC: printArgDesc(cout); throw 0;
		case 'S': outfile = arg; break;
		case 'U': {
			EList<string> args;
			tokenize(arg, ",", args);
			for(size_t i = 0; i < args.size(); i++) {
				queries.push_back(args[i]);
			}
			break;
		}
		case ARG_VERSION: showVersion = 1; break;
        case ARG_NO_TEMPSPLICESITE: useTempSpliceSite = false; break;
        case ARG_PEN_CANSPLICE: {
            penCanSplice = parseInt(0, "--pen-cansplice arg must be at least 0", arg);
            break;
        }
        case ARG_PEN_NONCANSPLICE: {
            penNoncanSplice = parseInt(0, "--pen-noncansplice arg must be at least 0", arg);
            break;
        }
        case ARG_PEN_CONFLICTSPLICE: {
            penConflictSplice = parseInt(0, "--pen-conflictsplice arg must be at least 0", arg);
            break;
        }
        case ARG_PEN_CANINTRONLEN: {
			polstr += ";";
			EList<string> args;
			tokenize(arg, ",", args);
			if(args.size() > 3 && args.size() == 0) {
				cerr << "Error: expected 3 or fewer comma-separated "
                << "arguments to --n-ceil option, got "
                << args.size() << endl;
				throw 1;
			}
			polstr += ("CANINTRONLEN=" + args[0]);
			if(args.size() > 1) {
				polstr += ("," + args[1]);
			}
			if(args.size() > 2) {
				polstr += ("," + args[2]);
			}
			break;
		}
        case ARG_PEN_NONCANINTRONLEN: {
            polstr += ";";
            EList<string> args;
            tokenize(arg, ",", args);
            if(args.size() > 3 && args.size() == 0) {
                cerr << "Error: expected 3 or fewer comma-separated "
                << "arguments to --n-ceil option, got "
                << args.size() << endl;
                throw 1;
            }
            polstr += ("NONCANINTRONLEN=" + args[0]);
            if(args.size() > 1) {
                polstr += ("," + args[1]);
            }
            if(args.size() > 2) {
                polstr += ("," + args[2]);
            }
            break;
        }
        case ARG_MIN_INTRONLEN: {
            minIntronLen = parseInt(20, "--min-intronlen arg must be at least 20", arg);
            break;
        }
        case ARG_MAX_INTRONLEN: {
            maxIntronLen = parseInt(20, "--max-intronlen arg must be at least 20", arg);
            break;
        }
        case ARG_KNOWN_SPLICESITE_INFILE: knownSpliceSiteInfile = arg; break;
        case ARG_NOVEL_SPLICESITE_INFILE: novelSpliceSiteInfile = arg; break;
        case ARG_NOVEL_SPLICESITE_OUTFILE: novelSpliceSiteOutfile = arg; break;
        case ARG_SECONDARY: secondary = true; break;
        case ARG_NO_SPLICED_ALIGNMENT: no_spliced_alignment = true; break;
        case ARG_RNA_STRANDNESS: {
            string strandness = arg;
            if(strandness == "F")       rna_strandness = RNA_STRANDNESS_F;
            else if(strandness == "R")  rna_strandness = RNA_STRANDNESS_R;
            else if(strandness == "FR") rna_strandness = RNA_STRANDNESS_FR;
            else if(strandness == "RF") rna_strandness = RNA_STRANDNESS_RF;
            else {
                cerr << "Error: should be one of F, R, FR, or RF " << endl;
				throw 1;
            }
            break;
        }
        case ARG_SPLICESITE_DB_ONLY: {
            splicesite_db_only = true;
            break;
        }
        case ARG_NO_ANCHORSTOP: {
            anchorStop = false;
            break;
        }
        case ARG_TRANSCRIPTOME_MAPPING_ONLY: {
            tranMapOnly = true;
            break;
        }
        case ARG_TRANSCRIPTOME_ASSEMBLY: {
            tranAssm = true;
            break;
        }
        case ARG_TRANSCRIPTOME_ASSEMBLY_CUFFLINKS: {
            tranAssm = true;
            tranAssm_program = "cufflinks";
            break;
        }
        case ARG_AVOID_PSEUDOGENE: {
            avoid_pseudogene = true;
            break;
        }
#ifdef USE_SRA
        case ARG_SRA_ACC: {
            tokenize(arg, ",", sra_accs); format = SRA_FASTA;
            break;
        }
#endif
        case ARG_REMOVE_CHRNAME: {
            rmChrName = true;
            break;
        }
        case ARG_ADD_CHRNAME: {
            addChrName = true;
            break;
        }
        case ARG_MAX_ALTSTRIED: {
            max_alts_tried = parseInt(8, "--max-altstried arg must be at least 8", arg);
            break;
        }
        case ARG_HAPLOTYPE: {
            use_haplotype = true;
            break;
        }
        case ARG_CODIS: {
            enable_codis = true;
            break;
        }
        case ARG_NO_TEMPLATELEN_ADJUSTMENT: {
            templateLenAdjustment = false;
            break;
        }
        case ARG_SUMMARY_FILE: {
            alignSumFile = arg;
            break;
        }
        case ARG_NEW_SUMMARY: {
            newAlignSummary = true;
            break;
        }
        case ARG_DP: {
            bowtie2_dp = parseInt(0, "--bowtie2-dp arg must be 0, 1, or 2", arg);
            break;
        }
        case ARG_REPEAT: {
            repeat = true;
            break;
        }
        case ARG_NO_REPEAT_INDEX: {
            use_repeat_index = false;
            break;
        }
        case ARG_READ_LENGTHS: {
            EList<string> str_readLens;
            tokenize(arg, ",", str_readLens);
            for(size_t i = 0; i < str_readLens.size(); i++) {
                int readLen = parseInt(0, "--read-lengths arg must be at least 0", str_readLens[i].c_str());
                readLens.push_back(readLen);
            }
            readLens.sort();
            break;
        }
		default:
			printUsage(cerr);
			throw 1;
	}
}

/**
 * Read command-line arguments
 */
static void parseOptions(int argc, const char **argv) {
	int option_index = 0;
	int next_option;
	saw_M = false;
	saw_a = false;
	saw_k = false;
	presetList.clear();
	if(startVerbose) { cerr << "Parsing options: "; logTime(cerr, true); }
	while(true) {
		next_option = getopt_long(
			argc, const_cast<char**>(argv),
			short_options, long_options, &option_index);
		const char * arg = optarg;
		if(next_option == EOF) {
			if(extra_opts_cur < extra_opts.size()) {
				next_option = extra_opts[extra_opts_cur].first;
				arg = extra_opts[extra_opts_cur].second.c_str();
				extra_opts_cur++;
			} else {
				break;
			}
		}
		parseOption(next_option, arg);
	}
	// Now parse all the presets.  Might want to pick which presets version to
	// use according to other parameters.
	auto_ptr<Presets> presets(new PresetsV0());
	// Apply default preset
	if(!defaultPreset.empty()) {
		polstr = applyPreset(defaultPreset, *presets.get()) + polstr;
	}
	// Apply specified presets
	for(size_t i = 0; i < presetList.size(); i++) {
		polstr += applyPreset(presetList[i], *presets.get());
	}
	for(size_t i = 0; i < extra_opts.size(); i++) {
		next_option = extra_opts[extra_opts_cur].first;
		const char *arg = extra_opts[extra_opts_cur].second.c_str();
		parseOption(next_option, arg);
	}
	// Remove initial semicolons
	while(!polstr.empty() && polstr[0] == ';') {
		polstr = polstr.substr(1);
	}
	if(gVerbose) {
		cerr << "Final policy string: '" << polstr.c_str() << "'" << endl;
	}
    
    size_t failStreakTmp = 0;
	SeedAlignmentPolicy::parseString(
                                     polstr,
                                     localAlign,
                                     noisyHpolymer,
                                     ignoreQuals,
                                     bonusMatchType,
                                     bonusMatch,
                                     penMmcType,
                                     penMmcMax,
                                     penMmcMin,
                                     penScMax,
                                     penScMin,
                                     penNType,
                                     penN,
                                     penRdGapConst,
                                     penRfGapConst,
                                     penRdGapLinear,
                                     penRfGapLinear,
                                     scoreMin,
                                     nCeil,
                                     penNCatPair,
                                     multiseedMms,
                                     multiseedLen,
                                     msIval,
                                     failStreakTmp,
                                     nSeedRounds,
                                     &penCanIntronLen,
                                     &penNoncanIntronLen);
	if(failStreakTmp > 0) {
		maxEeStreak = failStreakTmp;
		maxUgStreak = failStreakTmp;
		maxDpStreak = failStreakTmp;
	}
	if(saw_a || saw_k || true) {
		msample = false;
		mhits = 0;
	} else {
		assert_gt(mhits, 0);
		msample = true;
	}
    
    if(fast) {
        use_repeat_index = false;
    } else if(sensitive) {
        if(bowtie2_dp == 0) {
            bowtie2_dp = 1;
        }
        
        if(khits < 10) {
            khits = 10;
            saw_k = true;
        }
        scoreMin.init(SIMPLE_FUNC_LINEAR, 0.0f, -0.5f);
    } else if(very_sensitive) {
        bowtie2_dp = 2;
        if(khits < 30) {
            khits = 30;
            saw_k = true;
        }
        scoreMin.init(SIMPLE_FUNC_LINEAR, 0.0f, -1.0f);
    }
    
	if(mates1.size() != mates2.size()) {
		cerr << "Error: " << mates1.size() << " mate files/sequences were specified with -1, but " << mates2.size() << endl
		     << "mate files/sequences were specified with -2.  The same number of mate files/" << endl
		     << "sequences must be specified with -1 and -2." << endl;
		throw 1;
	}
	if(qualities.size() && format != FASTA) {
		cerr << "Error: one or more quality files were specified with -Q but -f was not" << endl
		     << "enabled.  -Q works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities1.size() && format != FASTA) {
		cerr << "Error: one or more quality files were specified with --Q1 but -f was not" << endl
		     << "enabled.  --Q1 works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities2.size() && format != FASTA) {
		cerr << "Error: one or more quality files were specified with --Q2 but -f was not" << endl
		     << "enabled.  --Q2 works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities1.size() > 0 && mates1.size() != qualities1.size()) {
		cerr << "Error: " << mates1.size() << " mate files/sequences were specified with -1, but " << qualities1.size() << endl
		     << "quality files were specified with --Q1.  The same number of mate and quality" << endl
		     << "files must sequences must be specified with -1 and --Q1." << endl;
		throw 1;
	}
	if(qualities2.size() > 0 && mates2.size() != qualities2.size()) {
		cerr << "Error: " << mates2.size() << " mate files/sequences were specified with -2, but " << qualities2.size() << endl
		     << "quality files were specified with --Q2.  The same number of mate and quality" << endl
		     << "files must sequences must be specified with -2 and --Q2." << endl;
		throw 1;
	}
	if(!rgs.empty() && rgid.empty()) {
		cerr << "Warning: --rg was specified without --rg-id also "
		     << "being specified.  @RG line is not printed unless --rg-id "
			 << "is specified." << endl;
	}
	// Check for duplicate mate input files
	if(format != CMDLINE) {
		for(size_t i = 0; i < mates1.size(); i++) {
			for(size_t j = 0; j < mates2.size(); j++) {
				if(mates1[i] == mates2[j] && !gQuiet) {
					cerr << "Warning: Same mate file \"" << mates1[i].c_str() << "\" appears as argument to both -1 and -2" << endl;
				}
			}
		}
	}
	// If both -s and -u are used, we need to adjust qUpto accordingly
	// since it uses rdid to know if we've reached the -u limit (and
	// rdids are all shifted up by skipReads characters)
	if(qUpto + skipReads > qUpto) {
		qUpto += skipReads;
	}
	if(useShmem && useMm && !gQuiet) {
		cerr << "Warning: --shmem overrides --mm..." << endl;
		useMm = false;
	}
	if(gGapBarrier < 1) {
		cerr << "Warning: --gbar was set less than 1 (=" << gGapBarrier
		     << "); setting to 1 instead" << endl;
		gGapBarrier = 1;
	}
	if(multiseedMms >= multiseedLen) {
		assert_gt(multiseedLen, 0);
		cerr << "Warning: seed mismatches (" << multiseedMms
		     << ") is less than seed length (" << multiseedLen
			 << "); setting mismatches to " << (multiseedMms-1)
			 << " instead" << endl;
		multiseedMms = multiseedLen-1;
	}
	sam_print_zm = sam_print_zm && bowtie2p5;
#ifndef NDEBUG
	if(!gQuiet) {
		cerr << "Warning: Running in debug mode.  Please use debug mode only "
			 << "for diagnosing errors, and not for typical use of HISAT2."
			 << endl;
	}
#endif
}

static const char *argv0 = NULL;

/// Create a PatternSourcePerThread for the current thread according
/// to the global params and return a pointer to it
static PatternSourcePerThreadFactory*
createPatsrcFactory(PairedPatternSource& _patsrc, int tid) {
	PatternSourcePerThreadFactory *patsrcFact;
	patsrcFact = new WrappedPatternSourcePerThreadFactory(_patsrc);
	assert(patsrcFact != NULL);
	return patsrcFact;
}

#define PTHREAD_ATTRS (PTHREAD_CREATE_JOINABLE | PTHREAD_CREATE_DETACHED)

typedef TIndexOffU index_t;
typedef uint16_t local_index_t;
static PairedPatternSource*              multiseed_patsrc;
static HGFM<index_t>*                    multiseed_gfm;
static RFM<index_t>*                     multiseed_rgfm;
static Scoring*                          multiseed_sc;
static BitPairReference*                 multiseed_refs;
static BitPairReference*                 multiseed_rrefs;
static AlnSink<index_t>*                 multiseed_msink;
static OutFileBuf*                       multiseed_metricsOfb;
static SpliceSiteDB*                     ssdb;
static ALTDB<index_t>*                   altdb;
static RepeatDB<index_t>*                repeatdb;
static ALTDB<index_t>*                   raltdb;
static TranscriptomePolicy*              multiseed_tpol;
static GraphPolicy*                      gpol;

/**
 * Metrics for measuring the work done by the outer read alignment
 * loop.
 */
struct OuterLoopMetrics {

	OuterLoopMetrics() {
	    reset();
	}

	/**
	 * Set all counters to 0.
	 */
	void reset() {
		reads = bases = srreads = srbases =
		freads = fbases = ureads = ubases = 0;
	}

	/**
	 * Sum the counters in m in with the conters in this object.  This
	 * is the only safe way to update an OuterLoopMetrics that's shared
	 * by multiple threads.
	 */
	void merge(
		const OuterLoopMetrics& m,
		bool getLock = false)
	{
		ThreadSafe ts(&mutex_m, getLock);
		reads += m.reads;
		bases += m.bases;
		srreads += m.srreads;
		srbases += m.srbases;
		freads += m.freads;
		fbases += m.fbases;
		ureads += m.ureads;
		ubases += m.ubases;
	}

	uint64_t reads;   // total reads
	uint64_t bases;   // total bases
	uint64_t srreads; // same-read reads
	uint64_t srbases; // same-read bases
	uint64_t freads;  // filtered reads
	uint64_t fbases;  // filtered bases
	uint64_t ureads;  // unfiltered reads
	uint64_t ubases;  // unfiltered bases
	MUTEX_T mutex_m;
};

/**
 * Collection of all relevant performance metrics when aligning in
 * multiseed mode.
 */
struct PerfMetrics {

	PerfMetrics() : first(true) { reset(); }

	/**
	 * Set all counters to 0.
	 */
	void reset() {
		olm.reset();
		sdm.reset();
		wlm.reset();
		swmSeed.reset();
		swmMate.reset();
		rpm.reset();
		dpSse8Seed.reset();   // 8-bit SSE seed extensions
		dpSse8Mate.reset();   // 8-bit SSE mate finds
		dpSse16Seed.reset();  // 16-bit SSE seed extensions
		dpSse16Mate.reset();  // 16-bit SSE mate finds
		nbtfiltst = 0;
		nbtfiltsc = 0;
		nbtfiltdo = 0;
		
		olmu.reset();
		sdmu.reset();
		wlmu.reset();
		swmuSeed.reset();
		swmuMate.reset();
		rpmu.reset();
		dpSse8uSeed.reset();  // 8-bit SSE seed extensions
		dpSse8uMate.reset();  // 8-bit SSE mate finds
		dpSse16uSeed.reset(); // 16-bit SSE seed extensions
		dpSse16uMate.reset(); // 16-bit SSE mate finds
		nbtfiltst_u = 0;
		nbtfiltsc_u = 0;
		nbtfiltdo_u = 0;
        
        him.reset();
	}

	/**
	 * Merge a set of specific metrics into this object.
	 */
	void merge(
		const OuterLoopMetrics *ol,
		const SeedSearchMetrics *sd,
		const WalkMetrics *wl,
		const SwMetrics *swSeed,
		const SwMetrics *swMate,
		const ReportingMetrics *rm,
		const SSEMetrics *dpSse8Ex,
		const SSEMetrics *dpSse8Ma,
		const SSEMetrics *dpSse16Ex,
		const SSEMetrics *dpSse16Ma,
		uint64_t nbtfiltst_,
		uint64_t nbtfiltsc_,
		uint64_t nbtfiltdo_,
        const HIMetrics *hi,
		bool getLock)
	{
		ThreadSafe ts(&mutex_m, getLock);
		if(ol != NULL) {
			olmu.merge(*ol, false);
		}
		if(sd != NULL) {
			sdmu.merge(*sd, false);
		}
		if(wl != NULL) {
			wlmu.merge(*wl, false);
		}
		if(swSeed != NULL) {
			swmuSeed.merge(*swSeed, false);
		}
		if(swMate != NULL) {
			swmuMate.merge(*swMate, false);
		}
		if(rm != NULL) {
			rpmu.merge(*rm, false);
		}
		if(dpSse8Ex != NULL) {
			dpSse8uSeed.merge(*dpSse8Ex, false);
		}
		if(dpSse8Ma != NULL) {
			dpSse8uMate.merge(*dpSse8Ma, false);
		}
		if(dpSse16Ex != NULL) {
			dpSse16uSeed.merge(*dpSse16Ex, false);
		}
		if(dpSse16Ma != NULL) {
			dpSse16uMate.merge(*dpSse16Ma, false);
		}
		nbtfiltst_u += nbtfiltst_;
		nbtfiltsc_u += nbtfiltsc_;
		nbtfiltdo_u += nbtfiltdo_;
        if(hi != NULL) {
            him.merge(*hi, false);
        }
	}

	/**
	 * Reports a matrix of results, incl. column labels, to an OutFileBuf.
	 * Optionally also sends results to stderr (unbuffered).  Can optionally
	 * print a per-read record with the read name at the beginning.
	 */
	void reportInterval(
		OutFileBuf* o,        // file to send output to
		bool metricsStderr,   // additionally output to stderr?
		bool total,           // true -> report total, otherwise incremental
		bool sync,            //  synchronize output
		const BTString *name) // non-NULL name pointer if is per-read record
	{
		ThreadSafe ts(&mutex_m, sync);
		ostringstream stderrSs;
		time_t curtime = time(0);
		char buf[1024];
		if(first) {
			const char *str =
				/*  1 */ "Time"           "\t"
				/*  2 */ "Read"           "\t"
				/*  3 */ "Base"           "\t"
				/*  4 */ "SameRead"       "\t"
				/*  5 */ "SameReadBase"   "\t"
				/*  6 */ "UnfilteredRead" "\t"
				/*  7 */ "UnfilteredBase" "\t"
				
				/*  8 */ "Paired"         "\t"
				/*  9 */ "Unpaired"       "\t"
				/* 10 */ "AlConUni"       "\t"
				/* 11 */ "AlConRep"       "\t"
				/* 12 */ "AlConFail"      "\t"
				/* 13 */ "AlDis"          "\t"
				/* 14 */ "AlConFailUni"   "\t"
				/* 15 */ "AlConFailRep"   "\t"
				/* 16 */ "AlConFailFail"  "\t"
				/* 17 */ "AlConRepUni"    "\t"
				/* 18 */ "AlConRepRep"    "\t"
				/* 19 */ "AlConRepFail"   "\t"
				/* 20 */ "AlUnpUni"       "\t"
				/* 21 */ "AlUnpRep"       "\t"
				/* 22 */ "AlUnpFail"      "\t"
				
				/* 23 */ "SeedSearch"     "\t"
				/* 24 */ "IntraSCacheHit" "\t"
				/* 25 */ "InterSCacheHit" "\t"
				/* 26 */ "OutOfMemory"    "\t"
				/* 27 */ "AlBWOp"         "\t"
				/* 28 */ "AlBWBranch"     "\t"
				/* 29 */ "ResBWOp"        "\t"
				/* 30 */ "ResBWBranch"    "\t"
				/* 31 */ "ResResolve"     "\t"
				/* 34 */ "ResReport"      "\t"
				/* 35 */ "RedundantSHit"  "\t"

				/* 36 */ "BestMinEdit0"   "\t"
				/* 37 */ "BestMinEdit1"   "\t"
				/* 38 */ "BestMinEdit2"   "\t"

				/* 39 */ "ExactAttempts"  "\t"
				/* 40 */ "ExactSucc"      "\t"
				/* 41 */ "ExactRanges"    "\t"
				/* 42 */ "ExactRows"      "\t"
				/* 43 */ "ExactOOMs"      "\t"

				/* 44 */ "1mmAttempts"    "\t"
				/* 45 */ "1mmSucc"        "\t"
				/* 46 */ "1mmRanges"      "\t"
				/* 47 */ "1mmRows"        "\t"
				/* 48 */ "1mmOOMs"        "\t"

				/* 49 */ "UngappedSucc"   "\t"
				/* 50 */ "UngappedFail"   "\t"
				/* 51 */ "UngappedNoDec"  "\t"

				/* 52 */ "DPExLt10Gaps"   "\t"
				/* 53 */ "DPExLt5Gaps"    "\t"
				/* 54 */ "DPExLt3Gaps"    "\t"

				/* 55 */ "DPMateLt10Gaps" "\t"
				/* 56 */ "DPMateLt5Gaps"  "\t"
				/* 57 */ "DPMateLt3Gaps"  "\t"

				/* 58 */ "DP16ExDps"      "\t"
				/* 59 */ "DP16ExDpSat"    "\t"
				/* 60 */ "DP16ExDpFail"   "\t"
				/* 61 */ "DP16ExDpSucc"   "\t"
				/* 62 */ "DP16ExCol"      "\t"
				/* 63 */ "DP16ExCell"     "\t"
				/* 64 */ "DP16ExInner"    "\t"
				/* 65 */ "DP16ExFixup"    "\t"
				/* 66 */ "DP16ExGathSol"  "\t"
				/* 67 */ "DP16ExBt"       "\t"
				/* 68 */ "DP16ExBtFail"   "\t"
				/* 69 */ "DP16ExBtSucc"   "\t"
				/* 70 */ "DP16ExBtCell"   "\t"
				/* 71 */ "DP16ExCoreRej"  "\t"
				/* 72 */ "DP16ExNRej"     "\t"

				/* 73 */ "DP8ExDps"       "\t"
				/* 74 */ "DP8ExDpSat"     "\t"
				/* 75 */ "DP8ExDpFail"    "\t"
				/* 76 */ "DP8ExDpSucc"    "\t"
				/* 77 */ "DP8ExCol"       "\t"
				/* 78 */ "DP8ExCell"      "\t"
				/* 79 */ "DP8ExInner"     "\t"
				/* 80 */ "DP8ExFixup"     "\t"
				/* 81 */ "DP8ExGathSol"   "\t"
				/* 82 */ "DP8ExBt"        "\t"
				/* 83 */ "DP8ExBtFail"    "\t"
				/* 84 */ "DP8ExBtSucc"    "\t"
				/* 85 */ "DP8ExBtCell"    "\t"
				/* 86 */ "DP8ExCoreRej"   "\t"
				/* 87 */ "DP8ExNRej"      "\t"

				/* 88 */ "DP16MateDps"     "\t"
				/* 89 */ "DP16MateDpSat"   "\t"
				/* 90 */ "DP16MateDpFail"  "\t"
				/* 91 */ "DP16MateDpSucc"  "\t"
				/* 92 */ "DP16MateCol"     "\t"
				/* 93 */ "DP16MateCell"    "\t"
				/* 94 */ "DP16MateInner"   "\t"
				/* 95 */ "DP16MateFixup"   "\t"
				/* 96 */ "DP16MateGathSol" "\t"
				/* 97 */ "DP16MateBt"      "\t"
				/* 98 */ "DP16MateBtFail"  "\t"
				/* 99 */ "DP16MateBtSucc"  "\t"
				/* 100 */ "DP16MateBtCell"  "\t"
				/* 101 */ "DP16MateCoreRej" "\t"
				/* 102 */ "DP16MateNRej"    "\t"

				/* 103 */ "DP8MateDps"     "\t"
				/* 104 */ "DP8MateDpSat"   "\t"
				/* 105 */ "DP8MateDpFail"  "\t"
				/* 106 */ "DP8MateDpSucc"  "\t"
				/* 107 */ "DP8MateCol"     "\t"
				/* 108 */ "DP8MateCell"    "\t"
				/* 109 */ "DP8MateInner"   "\t"
				/* 110 */ "DP8MateFixup"   "\t"
				/* 111 */ "DP8MateGathSol" "\t"
				/* 112 */ "DP8MateBt"      "\t"
				/* 113 */ "DP8MateBtFail"  "\t"
				/* 114 */ "DP8MateBtSucc"  "\t"
				/* 115 */ "DP8MateBtCell"  "\t"
				/* 116 */ "DP8MateCoreRej" "\t"
				/* 117 */ "DP8MateNRej"    "\t"

				/* 118 */ "DPBtFiltStart"  "\t"
				/* 119 */ "DPBtFiltScore"  "\t"
				/* 120 */ "DpBtFiltDom"    "\t"

				/* 121 */ "MemPeak"        "\t"
				/* 122 */ "UncatMemPeak"   "\t" // 0
				/* 123 */ "EbwtMemPeak"    "\t" // EBWT_CAT
				/* 124 */ "CacheMemPeak"   "\t" // CA_CAT
				/* 125 */ "ResolveMemPeak" "\t" // GW_CAT
				/* 126 */ "AlignMemPeak"   "\t" // AL_CAT
				/* 127 */ "DPMemPeak"      "\t" // DP_CAT
				/* 128 */ "MiscMemPeak"    "\t" // MISC_CAT
				/* 129 */ "DebugMemPeak"   "\t" // DEBUG_CAT
            
                /* 130 */ "LocalSearch"         "\t"
                /* 131 */ "AnchorSearch"        "\t"
                /* 132 */ "LocalIndexSearch"    "\t"
                /* 133 */ "LocalExtSearch"      "\t"
                /* 134 */ "LocalSearchRecur"    "\t"
                /* 135 */ "GlobalGenomeCoords"  "\t"
                /* 136 */ "LocalGenomeCoords"   "\t"
            
            
				"\n";
			
			if(name != NULL) {
				if(o != NULL) o->writeChars("Name\t");
				if(metricsStderr) stderrSs << "Name\t";
			}
			
			if(o != NULL) o->writeChars(str);
			if(metricsStderr) stderrSs << str;
			first = false;
		}
		
		if(total) mergeIncrementals();
		
		// 0. Read name, if needed
		if(name != NULL) {
			if(o != NULL) {
				o->writeChars(name->toZBuf());
				o->write('\t');
			}
			if(metricsStderr) {
				stderrSs << (*name) << '\t';
			}
		}
			
		// 1. Current time in secs
		itoa10<time_t>(curtime, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		
		const OuterLoopMetrics& ol = total ? olm : olmu;
		
		// 2. Reads
		itoa10<uint64_t>(ol.reads, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 3. Bases
		itoa10<uint64_t>(ol.bases, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 4. Same-read reads
		itoa10<uint64_t>(ol.srreads, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 5. Same-read bases
		itoa10<uint64_t>(ol.srbases, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 6. Unfiltered reads
		itoa10<uint64_t>(ol.ureads, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 7. Unfiltered bases
		itoa10<uint64_t>(ol.ubases, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }

		const ReportingMetrics& rp = total ? rpm : rpmu;

		// 8. Paired reads
		itoa10<uint64_t>(rp.npaired, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 9. Unpaired reads
		itoa10<uint64_t>(rp.nunpaired, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 10. Pairs with unique concordant alignments
		itoa10<uint64_t>(rp.nconcord_uni, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 11. Pairs with repetitive concordant alignments
		itoa10<uint64_t>(rp.nconcord_rep, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 12. Pairs with 0 concordant alignments
		itoa10<uint64_t>(rp.nconcord_0, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 13. Pairs with 1 discordant alignment
		itoa10<uint64_t>(rp.ndiscord, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 14. Mates from unaligned pairs that align uniquely
		itoa10<uint64_t>(rp.nunp_0_uni, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 15. Mates from unaligned pairs that align repetitively
		itoa10<uint64_t>(rp.nunp_0_rep, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 16. Mates from unaligned pairs that fail to align
		itoa10<uint64_t>(rp.nunp_0_0, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 17. Mates from repetitive pairs that align uniquely
		itoa10<uint64_t>(rp.nunp_rep_uni, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 18. Mates from repetitive pairs that align repetitively
		itoa10<uint64_t>(rp.nunp_rep_rep, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 19. Mates from repetitive pairs that fail to align
		itoa10<uint64_t>(rp.nunp_rep_0, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 20. Unpaired reads that align uniquely
		itoa10<uint64_t>(rp.nunp_uni, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 21. Unpaired reads that align repetitively
		itoa10<uint64_t>(rp.nunp_rep, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 22. Unpaired reads that fail to align
		itoa10<uint64_t>(rp.nunp_0, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }

		const SeedSearchMetrics& sd = total ? sdm : sdmu;
		
		// 23. Seed searches
		itoa10<uint64_t>(sd.seedsearch, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 24. Hits in 'current' cache
		itoa10<uint64_t>(sd.intrahit, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 25. Hits in 'local' cache
		itoa10<uint64_t>(sd.interhit, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 26. Out of memory
		itoa10<uint64_t>(sd.ooms, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 27. Burrows-Wheeler ops in aligner
		itoa10<uint64_t>(sd.bwops, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 28. Burrows-Wheeler branches (edits) in aligner
		itoa10<uint64_t>(sd.bweds, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		
		const WalkMetrics& wl = total ? wlm : wlmu;
		
		// 29. Burrows-Wheeler ops in resolver
		itoa10<uint64_t>(wl.bwops, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 30. Burrows-Wheeler branches in resolver
		itoa10<uint64_t>(wl.branches, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 31. Burrows-Wheeler offset resolutions
		itoa10<uint64_t>(wl.resolves, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 34. Offset reports
		itoa10<uint64_t>(wl.reports, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		
		// 35. Redundant seed hit
		itoa10<uint64_t>(total ? swmSeed.rshit : swmuSeed.rshit, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }

		// 36. # times the best (out of fw/rc) minimum # edits was 0
		itoa10<uint64_t>(total ? sdm.bestmin0 : sdmu.bestmin0, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 37. # times the best (out of fw/rc) minimum # edits was 1
		itoa10<uint64_t>(total ? sdm.bestmin1 : sdmu.bestmin1, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 38. # times the best (out of fw/rc) minimum # edits was 2
		itoa10<uint64_t>(total ? sdm.bestmin2 : sdmu.bestmin2, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		
		// 39. Exact aligner attempts
		itoa10<uint64_t>(total ? swmSeed.exatts : swmuSeed.exatts, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 40. Exact aligner successes
		itoa10<uint64_t>(total ? swmSeed.exsucc : swmuSeed.exsucc, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 41. Exact aligner ranges
		itoa10<uint64_t>(total ? swmSeed.exranges : swmuSeed.exranges, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 42. Exact aligner rows
		itoa10<uint64_t>(total ? swmSeed.exrows : swmuSeed.exrows, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 43. Exact aligner OOMs
		itoa10<uint64_t>(total ? swmSeed.exooms : swmuSeed.exooms, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }

		// 44. 1mm aligner attempts
		itoa10<uint64_t>(total ? swmSeed.mm1atts : swmuSeed.mm1atts, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 45. 1mm aligner successes
		itoa10<uint64_t>(total ? swmSeed.mm1succ : swmuSeed.mm1succ, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 46. 1mm aligner ranges
		itoa10<uint64_t>(total ? swmSeed.mm1ranges : swmuSeed.mm1ranges, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 47. 1mm aligner rows
		itoa10<uint64_t>(total ? swmSeed.mm1rows : swmuSeed.mm1rows, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 48. 1mm aligner OOMs
		itoa10<uint64_t>(total ? swmSeed.mm1ooms : swmuSeed.mm1ooms, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }

		// 49 Ungapped aligner success
		itoa10<uint64_t>(total ? swmSeed.ungapsucc : swmuSeed.ungapsucc, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 50. Ungapped aligner fail
		itoa10<uint64_t>(total ? swmSeed.ungapfail : swmuSeed.ungapfail, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 51. Ungapped aligner no decision
		itoa10<uint64_t>(total ? swmSeed.ungapnodec : swmuSeed.ungapnodec, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }

		// 52. # seed-extend DPs with < 10 gaps
		itoa10<uint64_t>(total ? swmSeed.sws10 : swmuSeed.sws10, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 53. # seed-extend DPs with < 5 gaps
		itoa10<uint64_t>(total ? swmSeed.sws5 : swmuSeed.sws5, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 54. # seed-extend DPs with < 3 gaps
		itoa10<uint64_t>(total ? swmSeed.sws3 : swmuSeed.sws3, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }

		// 55. # seed-extend DPs with < 10 gaps
		itoa10<uint64_t>(total ? swmMate.sws10 : swmuMate.sws10, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 56. # seed-extend DPs with < 5 gaps
		itoa10<uint64_t>(total ? swmMate.sws5 : swmuMate.sws5, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 57. # seed-extend DPs with < 3 gaps
		itoa10<uint64_t>(total ? swmMate.sws3 : swmuMate.sws3, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		
		const SSEMetrics& dpSse16s = total ? dpSse16Seed : dpSse16uSeed;
		
		// 58. 16-bit SSE seed-extend DPs tried
		itoa10<uint64_t>(dpSse16s.dp, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 59. 16-bit SSE seed-extend DPs saturated
		itoa10<uint64_t>(dpSse16s.dpsat, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 60. 16-bit SSE seed-extend DPs failed
		itoa10<uint64_t>(dpSse16s.dpfail, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 61. 16-bit SSE seed-extend DPs succeeded
		itoa10<uint64_t>(dpSse16s.dpsucc, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 62. 16-bit SSE seed-extend DP columns completed
		itoa10<uint64_t>(dpSse16s.col, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 63. 16-bit SSE seed-extend DP cells completed
		itoa10<uint64_t>(dpSse16s.cell, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 64. 16-bit SSE seed-extend DP inner loop iters completed
		itoa10<uint64_t>(dpSse16s.inner, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 65. 16-bit SSE seed-extend DP fixup loop iters completed
		itoa10<uint64_t>(dpSse16s.fixup, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 66. 16-bit SSE seed-extend DP gather, cells with potential solutions
		itoa10<uint64_t>(dpSse16s.gathsol, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 67. 16-bit SSE seed-extend DP backtrace attempts
		itoa10<uint64_t>(dpSse16s.bt, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 68. 16-bit SSE seed-extend DP failed backtrace attempts
		itoa10<uint64_t>(dpSse16s.btfail, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 69. 16-bit SSE seed-extend DP succesful backtrace attempts
		itoa10<uint64_t>(dpSse16s.btsucc, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 70. 16-bit SSE seed-extend DP backtrace cells
		itoa10<uint64_t>(dpSse16s.btcell, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 71. 16-bit SSE seed-extend DP core-diag rejections
		itoa10<uint64_t>(dpSse16s.corerej, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 72. 16-bit SSE seed-extend DP N rejections
		itoa10<uint64_t>(dpSse16s.nrej, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		
		const SSEMetrics& dpSse8s = total ? dpSse8Seed : dpSse8uSeed;
		
		// 73. 8-bit SSE seed-extend DPs tried
		itoa10<uint64_t>(dpSse8s.dp, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 74. 8-bit SSE seed-extend DPs saturated
		itoa10<uint64_t>(dpSse8s.dpsat, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 75. 8-bit SSE seed-extend DPs failed
		itoa10<uint64_t>(dpSse8s.dpfail, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 76. 8-bit SSE seed-extend DPs succeeded
		itoa10<uint64_t>(dpSse8s.dpsucc, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 77. 8-bit SSE seed-extend DP columns completed
		itoa10<uint64_t>(dpSse8s.col, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 78. 8-bit SSE seed-extend DP cells completed
		itoa10<uint64_t>(dpSse8s.cell, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 79. 8-bit SSE seed-extend DP inner loop iters completed
		itoa10<uint64_t>(dpSse8s.inner, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 80. 8-bit SSE seed-extend DP fixup loop iters completed
		itoa10<uint64_t>(dpSse8s.fixup, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 81. 16-bit SSE seed-extend DP gather, cells with potential solutions
		itoa10<uint64_t>(dpSse8s.gathsol, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 82. 16-bit SSE seed-extend DP backtrace attempts
		itoa10<uint64_t>(dpSse8s.bt, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 83. 16-bit SSE seed-extend DP failed backtrace attempts
		itoa10<uint64_t>(dpSse8s.btfail, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 84. 16-bit SSE seed-extend DP succesful backtrace attempts
		itoa10<uint64_t>(dpSse8s.btsucc, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 85. 16-bit SSE seed-extend DP backtrace cells
		itoa10<uint64_t>(dpSse8s.btcell, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 86. 16-bit SSE seed-extend DP core-diag rejections
		itoa10<uint64_t>(dpSse8s.corerej, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 87. 16-bit SSE seed-extend DP N rejections
		itoa10<uint64_t>(dpSse8s.nrej, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		
		const SSEMetrics& dpSse16m = total ? dpSse16Mate : dpSse16uMate;
		
		// 88. 16-bit SSE mate-finding DPs tried
		itoa10<uint64_t>(dpSse16m.dp, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 89. 16-bit SSE mate-finding DPs saturated
		itoa10<uint64_t>(dpSse16m.dpsat, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 90. 16-bit SSE mate-finding DPs failed
		itoa10<uint64_t>(dpSse16m.dpfail, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 91. 16-bit SSE mate-finding DPs succeeded
		itoa10<uint64_t>(dpSse16m.dpsucc, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 92. 16-bit SSE mate-finding DP columns completed
		itoa10<uint64_t>(dpSse16m.col, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 93. 16-bit SSE mate-finding DP cells completed
		itoa10<uint64_t>(dpSse16m.cell, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 94. 16-bit SSE mate-finding DP inner loop iters completed
		itoa10<uint64_t>(dpSse16m.inner, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 95. 16-bit SSE mate-finding DP fixup loop iters completed
		itoa10<uint64_t>(dpSse16m.fixup, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 96. 16-bit SSE mate-finding DP gather, cells with potential solutions
		itoa10<uint64_t>(dpSse16m.gathsol, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 97. 16-bit SSE mate-finding DP backtrace attempts
		itoa10<uint64_t>(dpSse16m.bt, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 98. 16-bit SSE mate-finding DP failed backtrace attempts
		itoa10<uint64_t>(dpSse16m.btfail, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 99. 16-bit SSE mate-finding DP succesful backtrace attempts
		itoa10<uint64_t>(dpSse16m.btsucc, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 100. 16-bit SSE mate-finding DP backtrace cells
		itoa10<uint64_t>(dpSse16m.btcell, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 101. 16-bit SSE mate-finding DP core-diag rejections
		itoa10<uint64_t>(dpSse16m.corerej, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 102. 16-bit SSE mate-finding DP N rejections
		itoa10<uint64_t>(dpSse16m.nrej, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		
		const SSEMetrics& dpSse8m = total ? dpSse8Mate : dpSse8uMate;
		
		// 103. 8-bit SSE mate-finding DPs tried
		itoa10<uint64_t>(dpSse8m.dp, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 104. 8-bit SSE mate-finding DPs saturated
		itoa10<uint64_t>(dpSse8m.dpsat, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 105. 8-bit SSE mate-finding DPs failed
		itoa10<uint64_t>(dpSse8m.dpfail, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 106. 8-bit SSE mate-finding DPs succeeded
		itoa10<uint64_t>(dpSse8m.dpsucc, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 107. 8-bit SSE mate-finding DP columns completed
		itoa10<uint64_t>(dpSse8m.col, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 108. 8-bit SSE mate-finding DP cells completed
		itoa10<uint64_t>(dpSse8m.cell, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 109. 8-bit SSE mate-finding DP inner loop iters completed
		itoa10<uint64_t>(dpSse8m.inner, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 110. 8-bit SSE mate-finding DP fixup loop iters completed
		itoa10<uint64_t>(dpSse8m.fixup, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 111. 16-bit SSE mate-finding DP gather, cells with potential solutions
		itoa10<uint64_t>(dpSse8m.gathsol, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 112. 16-bit SSE mate-finding DP backtrace attempts
		itoa10<uint64_t>(dpSse8m.bt, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 113. 16-bit SSE mate-finding DP failed backtrace attempts
		itoa10<uint64_t>(dpSse8m.btfail, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 114. 16-bit SSE mate-finding DP succesful backtrace attempts
		itoa10<uint64_t>(dpSse8m.btsucc, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 115. 16-bit SSE mate-finding DP backtrace cells
		itoa10<uint64_t>(dpSse8m.btcell, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 116. 16-bit SSE mate-finding DP core rejections
		itoa10<uint64_t>(dpSse8m.corerej, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 117. 16-bit SSE mate-finding N rejections
		itoa10<uint64_t>(dpSse8m.nrej, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		
		// 118. Backtrace candidates filtered due to starting cell
		itoa10<uint64_t>(total ? nbtfiltst : nbtfiltst_u, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 119. Backtrace candidates filtered due to low score
		itoa10<uint64_t>(total ? nbtfiltsc : nbtfiltsc_u, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 120. Backtrace candidates filtered due to domination
		itoa10<uint64_t>(total ? nbtfiltdo : nbtfiltdo_u, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		
		// 121. Overall memory peak
		itoa10<size_t>(gMemTally.peak() >> 20, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 122. Uncategorized memory peak
		itoa10<size_t>(gMemTally.peak(0) >> 20, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 123. Ebwt memory peak
		itoa10<size_t>(gMemTally.peak(EBWT_CAT) >> 20, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 124. Cache memory peak
		itoa10<size_t>(gMemTally.peak(CA_CAT) >> 20, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 125. Resolver memory peak
		itoa10<size_t>(gMemTally.peak(GW_CAT) >> 20, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 126. Seed aligner memory peak
		itoa10<size_t>(gMemTally.peak(AL_CAT) >> 20, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 127. Dynamic programming aligner memory peak
		itoa10<size_t>(gMemTally.peak(DP_CAT) >> 20, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 128. Miscellaneous memory peak
		itoa10<size_t>(gMemTally.peak(MISC_CAT) >> 20, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 129. Debug memory peak
		itoa10<size_t>(gMemTally.peak(DEBUG_CAT) >> 20, buf);
        if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
        
        // 130
        itoa10<size_t>(him.localatts, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
        // 131
        itoa10<size_t>(him.anchoratts, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
        // 132
        itoa10<size_t>(him.localindexatts, buf);
		if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
        // 133
        itoa10<size_t>(him.localextatts, buf);
        if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
        // 134
        itoa10<size_t>(him.localsearchrecur, buf);
        if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
        // 135
        itoa10<size_t>(him.globalgenomecoords, buf);
        if(metricsStderr) stderrSs << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
        // 136
        itoa10<size_t>(him.localgenomecoords, buf);
        if(metricsStderr) stderrSs << buf;
		if(o != NULL) { o->writeChars(buf); }

		if(o != NULL) { o->write('\n'); }
		if(metricsStderr) cerr << stderrSs.str().c_str() << endl;
		if(!total) mergeIncrementals();
	}
	
	void mergeIncrementals() {
		olm.merge(olmu, false);
		sdm.merge(sdmu, false);
		wlm.merge(wlmu, false);
		swmSeed.merge(swmuSeed, false);
		swmMate.merge(swmuMate, false);
		dpSse8Seed.merge(dpSse8uSeed, false);
		dpSse8Mate.merge(dpSse8uMate, false);
		dpSse16Seed.merge(dpSse16uSeed, false);
		dpSse16Mate.merge(dpSse16uMate, false);
		nbtfiltst_u += nbtfiltst;
		nbtfiltsc_u += nbtfiltsc;
		nbtfiltdo_u += nbtfiltdo;

		olmu.reset();
		sdmu.reset();
		wlmu.reset();
		swmuSeed.reset();
		swmuMate.reset();
		rpmu.reset();
		dpSse8uSeed.reset();
		dpSse8uMate.reset();
		dpSse16uSeed.reset();
		dpSse16uMate.reset();
		nbtfiltst_u = 0;
		nbtfiltsc_u = 0;
		nbtfiltdo_u = 0;
	}

	// Total over the whole job
	OuterLoopMetrics  olm;   // overall metrics
	SeedSearchMetrics sdm;   // metrics related to seed alignment
	WalkMetrics       wlm;   // metrics related to walking left (i.e. resolving reference offsets)
	SwMetrics         swmSeed;  // metrics related to DP seed-extend alignment
	SwMetrics         swmMate;  // metrics related to DP mate-finding alignment
	ReportingMetrics  rpm;   // metrics related to reporting
	SSEMetrics        dpSse8Seed;  // 8-bit SSE seed extensions
	SSEMetrics        dpSse8Mate;    // 8-bit SSE mate finds
	SSEMetrics        dpSse16Seed; // 16-bit SSE seed extensions
	SSEMetrics        dpSse16Mate;   // 16-bit SSE mate finds
	uint64_t          nbtfiltst;
	uint64_t          nbtfiltsc;
	uint64_t          nbtfiltdo;

	// Just since the last update
	OuterLoopMetrics  olmu;  // overall metrics
	SeedSearchMetrics sdmu;  // metrics related to seed alignment
	WalkMetrics       wlmu;  // metrics related to walking left (i.e. resolving reference offsets)
	SwMetrics         swmuSeed; // metrics related to DP seed-extend alignment
	SwMetrics         swmuMate; // metrics related to DP mate-finding alignment
	ReportingMetrics  rpmu;  // metrics related to reporting
	SSEMetrics        dpSse8uSeed;  // 8-bit SSE seed extensions
	SSEMetrics        dpSse8uMate;  // 8-bit SSE mate finds
	SSEMetrics        dpSse16uSeed; // 16-bit SSE seed extensions
	SSEMetrics        dpSse16uMate; // 16-bit SSE mate finds
	uint64_t          nbtfiltst_u;
	uint64_t          nbtfiltsc_u;
	uint64_t          nbtfiltdo_u;
    
    //
    HIMetrics         him;

	MUTEX_T           mutex_m;  // lock for when one ob
	bool              first; // yet to print first line?
	time_t            lastElapsed; // used in reportInterval to measure time since last call
};

static PerfMetrics metrics;

// Cyclic rotations
#define ROTL(n, x) (((x) << (n)) | ((x) >> (32-n)))
#define ROTR(n, x) (((x) >> (n)) | ((x) << (32-n)))

static inline void printMmsSkipMsg(
	const PatternSourcePerThread& ps,
	bool paired,
	bool mate1,
	int seedmms)
{
	ostringstream os;
	if(paired) {
		os << "Warning: skipping mate #" << (mate1 ? '1' : '2')
		   << " of read '" << (mate1 ? ps.bufa().name : ps.bufb().name)
		   << "' because length (" << (mate1 ? ps.bufa().patFw.length() : ps.bufb().patFw.length())
		   << ") <= # seed mismatches (" << seedmms << ")" << endl;
	} else {
		os << "Warning: skipping read '" << (mate1 ? ps.bufa().name : ps.bufb().name)
		   << "' because length (" << (mate1 ? ps.bufa().patFw.length() : ps.bufb().patFw.length())
		   << ") <= # seed mismatches (" << seedmms << ")" << endl;
	}
	cerr << os.str().c_str();
}

static inline void printLenSkipMsg(
	const PatternSourcePerThread& ps,
	bool paired,
	bool mate1)
{
	ostringstream os;
	if(paired) {
		os << "Warning: skipping mate #" << (mate1 ? '1' : '2')
		   << " of read '" << (mate1 ? ps.bufa().name : ps.bufb().name)
		   << "' because it was < 2 characters long" << endl;
	} else {
		os << "Warning: skipping read '" << (mate1 ? ps.bufa().name : ps.bufb().name)
		   << "' because it was < 2 characters long" << endl;
	}
	cerr << os.str().c_str();
}

static inline void printLocalScoreMsg(
	const PatternSourcePerThread& ps,
	bool paired,
	bool mate1)
{
	ostringstream os;
	if(paired) {
		os << "Warning: minimum score function gave negative number in "
		   << "--local mode for mate #" << (mate1 ? '1' : '2')
		   << " of read '" << (mate1 ? ps.bufa().name : ps.bufb().name)
		   << "; setting to 0 instead" << endl;
	} else {
		os << "Warning: minimum score function gave negative number in "
		   << "--local mode for read '" << (mate1 ? ps.bufa().name : ps.bufb().name)
		   << "; setting to 0 instead" << endl;
	}
	cerr << os.str().c_str();
}

static inline void printEEScoreMsg(
	const PatternSourcePerThread& ps,
	bool paired,
	bool mate1)
{
	ostringstream os;
	if(paired) {
		os << "Warning: minimum score function gave positive number in "
		   << "--end-to-end mode for mate #" << (mate1 ? '1' : '2')
		   << " of read '" << (mate1 ? ps.bufa().name : ps.bufb().name)
		   << "; setting to 0 instead" << endl;
	} else {
		os << "Warning: minimum score function gave positive number in "
		   << "--end-to-end mode for read '" << (mate1 ? ps.bufa().name : ps.bufb().name)
		   << "; setting to 0 instead" << endl;
	}
	cerr << os.str().c_str();
}


#define MERGE_METRICS(met, sync) { \
	msink.mergeMetrics(rpm); \
	met.merge( \
		&olm, \
		&sdm, \
		&wlm, \
		&swmSeed, \
		&swmMate, \
		&rpm, \
		&sseU8ExtendMet, \
		&sseU8MateMet, \
		&sseI16ExtendMet, \
		&sseI16MateMet, \
		nbtfiltst, \
		nbtfiltsc, \
		nbtfiltdo, \
        &him, \
		sync); \
	olm.reset(); \
	sdm.reset(); \
	wlm.reset(); \
	swmSeed.reset(); \
	swmMate.reset(); \
	rpm.reset(); \
	sseU8ExtendMet.reset(); \
	sseU8MateMet.reset(); \
	sseI16ExtendMet.reset(); \
	sseI16MateMet.reset(); \
    him.reset(); \
}

#define MERGE_SW(x) { \
	x.merge( \
		sseU8ExtendMet, \
		sseU8MateMet, \
		sseI16ExtendMet, \
		sseI16MateMet, \
		nbtfiltst, \
		nbtfiltsc, \
		nbtfiltdo); \
	x.resetCounters(); \
}

/**
 * Called once per thread.  Sets up per-thread pointers to the shared global
 * data structures, creates per-thread structures, then enters the alignment
 * loop.  The general flow of the alignment loop is:
 *
 * - If it's been a while and we're the master thread, report some alignment
 *   metrics
 * - Get the next read/pair
 * - Check if this read/pair is identical to the previous
 *   + If identical, check whether we can skip any or all alignment stages.  If
 *     we can skip all stages, report the result immediately and move to next
 *     read/pair
 *   + If not identical, continue
 * -
 */
static void multiseedSearchWorker_hisat2(void *vp) {
	int tid = *((int*)vp);
	assert(multiseed_gfm != NULL);
	assert(multiseedMms == 0);
	PairedPatternSource&             patsrc   = *multiseed_patsrc;
	const HGFM<index_t>&             gfm      = *multiseed_gfm;
    const RFM<index_t>*              rgfm     = multiseed_rgfm;
	const Scoring&                   sc       = *multiseed_sc;
    const BitPairReference&          ref      = *multiseed_refs;
    const BitPairReference*          rref     = multiseed_rrefs;
	AlnSink<index_t>&                msink    = *multiseed_msink;
	OutFileBuf*                      metricsOfb = multiseed_metricsOfb;
    
	// Sinks: these are so that we can print tables encoding counts for
	// events of interest on a per-read, per-seed, per-join, or per-SW
	// level.  These in turn can be used to diagnose performance
	// problems, or generally characterize performance.
	
	//const BitPairReference& refs   = *multiseed_refs;
	auto_ptr<PatternSourcePerThreadFactory> patsrcFact(createPatsrcFactory(patsrc, tid));
	auto_ptr<PatternSourcePerThread> ps(patsrcFact->create());
	
    // Instantiate an object for holding reporting-related parameters.
    if(maxSeeds == 0) {
        maxSeeds = max<size_t>(5, khits * 2);
    }
    ReportingParams rp(
                       (allHits ? std::numeric_limits<THitInt>::max() : khits), // -k
                       (allHits ? std::numeric_limits<THitInt>::max() : maxSeeds), // --max-seeds
                       mhits,             // -m/-M
                       0,                 // penalty gap (not used now)
                       msample,           // true -> -M was specified, otherwise assume -m
                       gReportDiscordant, // report discordang paired-end alignments?
                       gReportMixed,      // report unpaired alignments for paired reads?
                       secondary,
                       localAlign,
                       bowtie2_dp,
                       sensitive | very_sensitive,
                       repeat);
    
	// Instantiate a mapping quality calculator
	auto_ptr<Mapq> bmapq(new_mapq(mapqv, scoreMin, sc));
	
	// Make a per-thread wrapper for the global MHitSink object.
	AlnSinkWrap<index_t> msinkwrap(
                                   msink,         // global sink
                                   rp,            // reporting parameters
                                   *bmapq.get(),  // MAPQ calculator
                                   (size_t)tid,   // thread id
                                   secondary,     // secondary alignments
                                   no_spliced_alignment ? NULL : ssdb,
                                   thread_rids_mindist);
    
    SplicedAligner<index_t, local_index_t> splicedAligner(
                                                          gfm,
                                                          anchorStop,
                                                          thread_rids_mindist);
	SwAligner sw;
	OuterLoopMetrics olm;
	SeedSearchMetrics sdm;
	WalkMetrics wlm;
	SwMetrics swmSeed, swmMate;
	ReportingMetrics rpm;
	RandomSource rnd, rndArb;
	SSEMetrics sseU8ExtendMet;
	SSEMetrics sseU8MateMet;
	SSEMetrics sseI16ExtendMet;
	SSEMetrics sseI16MateMet;
   	DescentMetrics descm;
	uint64_t nbtfiltst = 0; // TODO: find a new home for these
	uint64_t nbtfiltsc = 0; // TODO: find a new home for these
	uint64_t nbtfiltdo = 0; // TODO: find a new home for these
    HIMetrics him;
    
	ASSERT_ONLY(BTDnaString tmp);
    
	int pepolFlag;
	if(gMate1fw && gMate2fw) {
		pepolFlag = PE_POLICY_FF;
	} else if(gMate1fw && !gMate2fw) {
		pepolFlag = PE_POLICY_FR;
	} else if(!gMate1fw && gMate2fw) {
		pepolFlag = PE_POLICY_RF;
	} else {
		pepolFlag = PE_POLICY_RR;
	}
	assert_geq(gMaxInsert, gMinInsert);
	assert_geq(gMinInsert, 0);
	PairedEndPolicy pepol(
                          pepolFlag,
                          gMaxInsert,
                          gMinInsert,
                          localAlign,
                          gFlippedMatesOK,
                          gDovetailMatesOK,
                          gContainMatesOK,
                          gOlapMatesOK,
                          gExpandToFrag);
    
  	PerfMetrics metricsPt; // per-thread metrics object; for read-level metrics
	BTString nametmp;
	
	PerReadMetrics prm;
    
	// Used by thread with threadid == 1 to measure time elapsed
	time_t iTime = time(0);
    
	// Keep track of whether last search was exhaustive for mates 1 and 2
	bool exhaustive[2] = { false, false };
	// Keep track of whether mates 1/2 were filtered out last time through
	bool filt[2]    = { true, true };
	// Keep track of whether mates 1/2 were filtered out due Ns last time
	bool nfilt[2]   = { true, true };
	// Keep track of whether mates 1/2 were filtered out due to not having
	// enough characters to rise about the score threshold.
	bool scfilt[2]  = { true, true };
	// Keep track of whether mates 1/2 were filtered out due to not having
	// more characters than the number of mismatches permitted in a seed.
	bool lenfilt[2] = { true, true };
	// Keep track of whether mates 1/2 were filtered out by upstream qc
	bool qcfilt[2]  = { true, true };
    
	rndArb.init((uint32_t)time(0));
	int mergei = 0;
	int mergeival = 16;
	while(true) {
		bool success = false, done = false, paired = false;
		ps->nextReadPair(success, done, paired, outType != OUTPUT_SAM);
		if(!success && done) {
			break;
		} else if(!success) {
			continue;
		}
		TReadId rdid = ps->rdid();
        if(nthreads > 1 && useTempSpliceSite) {
            assert_gt(tid, 0);
            assert_leq(tid, thread_rids.size());
            assert(thread_rids[tid - 1] == 0 || rdid > thread_rids[tid - 1]);
            thread_rids[tid - 1] = (rdid > 0 ? rdid - 1 : 0);
            while(true) {
                uint64_t min_rdid = thread_rids[0];
                {
                    for(size_t i = 1; i < thread_rids.size(); i++) {
                        if(thread_rids[i] < min_rdid) {
                            min_rdid = thread_rids[i];
                        }
                    }
                }
                
                if(min_rdid + thread_rids_mindist < rdid) {
#if defined(_TTHREAD_WIN32_)
                    Sleep(0);
#elif defined(_TTHREAD_POSIX_)
                    sched_yield();
#endif
                } else break;
            }
        }
        
		bool sample = true;
		if(arbitraryRandom) {
			ps->bufa().seed = rndArb.nextU32();
			ps->bufb().seed = rndArb.nextU32();
		}
		if(sampleFrac < 1.0f) {
			rnd.init(ROTL(ps->bufa().seed, 2));
			sample = rnd.nextFloat() < sampleFrac;
		}
		if(rdid >= skipReads && rdid < qUpto && sample) {
			// Align this read/pair
			bool retry = true;
			//
			// Check if there is metrics reporting for us to do.
			//
			if(metricsIval > 0 &&
			   (metricsOfb != NULL || metricsStderr) &&
			   !metricsPerRead &&
			   ++mergei == mergeival)
			{
				// Do a periodic merge.  Update global metrics, in a
				// synchronized manner if needed.
				MERGE_METRICS(metrics, nthreads > 1);
				mergei = 0;
				// Check if a progress message should be printed
				if(tid == 0) {
					// Only thread 1 prints progress messages
					time_t curTime = time(0);
					if(curTime - iTime >= metricsIval) {
						metrics.reportInterval(metricsOfb, metricsStderr, false, true, NULL);
						iTime = curTime;
					}
				}
			}
			prm.reset(); // per-read metrics
			prm.doFmString = false;
			if(sam_print_xt) {
				gettimeofday(&prm.tv_beg, &prm.tz_beg);
			}
			// Try to align this read
			while(retry) {
				retry = false;
				assert_eq(ps->bufa().color, false);
				olm.reads++;
				bool pair = paired;
				const size_t rdlen1 = ps->bufa().length();
				const size_t rdlen2 = pair ? ps->bufb().length() : 0;
				olm.bases += (rdlen1 + rdlen2);
				msinkwrap.nextRead(
                                   &ps->bufa(),
                                   pair ? &ps->bufb() : NULL,
                                   rdid,
                                   sc.qualitiesMatter());
				assert(msinkwrap.inited());
				size_t rdlens[2] = { rdlen1, rdlen2 };
				// Calculate the minimum valid score threshold for the read
				TAlScore minsc[2], maxpen[2];
				maxpen[0] = maxpen[1] = 0;
				minsc[0] = minsc[1] = std::numeric_limits<TAlScore>::max();
				if(bwaSwLike) {
					// From BWA-SW manual: "Given an l-long query, the
					// threshold for a hit to be retained is
					// a*max{T,c*log(l)}."  We try to recreate that here.
					float a = (float)sc.match(30);
					float T = bwaSwLikeT, c = bwaSwLikeC;
					minsc[0] = (TAlScore)max<float>(a*T, a*c*log(rdlens[0]));
					if(paired) {
						minsc[1] = (TAlScore)max<float>(a*T, a*c*log(rdlens[1]));
					}
				} else {
					minsc[0] = scoreMin.f<TAlScore>(rdlens[0]);
					if(paired) minsc[1] = scoreMin.f<TAlScore>(rdlens[1]);
					if(localAlign) {
						if(minsc[0] < 0) {
							if(!gQuiet) printLocalScoreMsg(*ps, paired, true);
							minsc[0] = 0;
						}
						if(paired && minsc[1] < 0) {
							if(!gQuiet) printLocalScoreMsg(*ps, paired, false);
							minsc[1] = 0;
						}
					} else {
						if(minsc[0] > 0) {
							if(!gQuiet) printEEScoreMsg(*ps, paired, true);
							minsc[0] = 0;
						}
						if(paired && minsc[1] > 0) {
							if(!gQuiet) printEEScoreMsg(*ps, paired, false);
							minsc[1] = 0;
						}
					}
				}
                
				// N filter; does the read have too many Ns?
				size_t readns[2] = {0, 0};
				sc.nFilterPair(
                               &ps->bufa().patFw,
                               pair ? &ps->bufb().patFw : NULL,
                               readns[0],
                               readns[1],
                               nfilt[0],
                               nfilt[1]);
				// Score filter; does the read enough character to rise above
				// the score threshold?
				scfilt[0] = sc.scoreFilter(minsc[0], rdlens[0]);
				scfilt[1] = sc.scoreFilter(minsc[1], rdlens[1]);
				lenfilt[0] = lenfilt[1] = true;
				if(rdlens[0] <= (size_t)multiseedMms || rdlens[0] < 2) {
					if(!gQuiet) printMmsSkipMsg(*ps, paired, true, multiseedMms);
					lenfilt[0] = false;
				}
				if((rdlens[1] <= (size_t)multiseedMms || rdlens[1] < 2) && paired) {
					if(!gQuiet) printMmsSkipMsg(*ps, paired, false, multiseedMms);
					lenfilt[1] = false;
				}
				if(rdlens[0] < 2) {
					if(!gQuiet) printLenSkipMsg(*ps, paired, true);
					lenfilt[0] = false;
				}
				if(rdlens[1] < 2 && paired) {
					if(!gQuiet) printLenSkipMsg(*ps, paired, false);
					lenfilt[1] = false;
				}
				qcfilt[0] = qcfilt[1] = true;
				if(qcFilter) {
					qcfilt[0] = (ps->bufa().filter != '0');
					qcfilt[1] = (ps->bufb().filter != '0');
				}
				filt[0] = (nfilt[0] && scfilt[0] && lenfilt[0] && qcfilt[0]);
				filt[1] = (nfilt[1] && scfilt[1] && lenfilt[1] && qcfilt[1]);
				prm.nFilt += (filt[0] ? 0 : 1) + (filt[1] ? 0 : 1);
				Read* rds[2] = { &ps->bufa(), &ps->bufb() };
				// For each mate...
				assert(msinkwrap.empty());
				//size_t minedfw[2] = { 0, 0 };
				//size_t minedrc[2] = { 0, 0 };
				// Calcualte nofw / no rc
				bool nofw[2] = { false, false };
				bool norc[2] = { false, false };
				nofw[0] = paired ? (gMate1fw ? gNofw : gNorc) : gNofw;
				norc[0] = paired ? (gMate1fw ? gNorc : gNofw) : gNorc;
				nofw[1] = paired ? (gMate2fw ? gNofw : gNorc) : gNofw;
				norc[1] = paired ? (gMate2fw ? gNorc : gNofw) : gNorc;
				// Calculate nceil
				int nceil[2] = { 0, 0 };
				nceil[0] = nCeil.f<int>((double)rdlens[0]);
				nceil[0] = min(nceil[0], (int)rdlens[0]);
				if(paired) {
					nceil[1] = nCeil.f<int>((double)rdlens[1]);
					nceil[1] = min(nceil[1], (int)rdlens[1]);
				}
				exhaustive[0] = exhaustive[1] = false;
				//size_t matemap[2] = { 0, 1 };
				bool pairPostFilt = filt[0] && filt[1];
				if(pairPostFilt) {
					rnd.init(ps->bufa().seed ^ ps->bufb().seed);
				} else {
					rnd.init(ps->bufa().seed);
				}
				// Calculate interval length for both mates
				int interval[2] = { 0, 0 };
				for(size_t mate = 0; mate < (pair ? 2:1); mate++) {
					interval[mate] = msIval.f<int>((double)rdlens[mate]);
					if(filt[0] && filt[1]) {
						// Boost interval length by 20% for paired-end reads
						interval[mate] = (int)(interval[mate] * 1.2 + 0.5);
					}
					interval[mate] = max(interval[mate], 1);
				}
				// Calculate streak length
				size_t streak[2]    = { maxDpStreak,   maxDpStreak };
				size_t mtStreak[2]  = { maxMateStreak, maxMateStreak };
				size_t mxDp[2]      = { maxDp,         maxDp       };
				size_t mxUg[2]      = { maxUg,         maxUg       };
				size_t mxIter[2]    = { maxIters,      maxIters    };
				if(allHits) {
					streak[0]   = streak[1]   = std::numeric_limits<size_t>::max();
					mtStreak[0] = mtStreak[1] = std::numeric_limits<size_t>::max();
					mxDp[0]     = mxDp[1]     = std::numeric_limits<size_t>::max();
					mxUg[0]     = mxUg[1]     = std::numeric_limits<size_t>::max();
					mxIter[0]   = mxIter[1]   = std::numeric_limits<size_t>::max();
				} else if(khits > 1) {
					for(size_t mate = 0; mate < 2; mate++) {
						streak[mate]   += (khits-1) * maxStreakIncr;
						mtStreak[mate] += (khits-1) * maxStreakIncr;
						mxDp[mate]     += (khits-1) * maxItersIncr;
						mxUg[mate]     += (khits-1) * maxItersIncr;
						mxIter[mate]   += (khits-1) * maxItersIncr;
					}
				}
				if(filt[0] && filt[1]) {
					streak[0] = (size_t)ceil((double)streak[0] / 2.0);
					streak[1] = (size_t)ceil((double)streak[1] / 2.0);
					assert_gt(streak[1], 0);
				}
				assert_gt(streak[0], 0);
				// Calculate # seed rounds for each mate
				size_t nrounds[2] = { nSeedRounds, nSeedRounds };
				if(filt[0] && filt[1]) {
					nrounds[0] = (size_t)ceil((double)nrounds[0] / 2.0);
					nrounds[1] = (size_t)ceil((double)nrounds[1] / 2.0);
					assert_gt(nrounds[1], 0);
				}
				assert_gt(nrounds[0], 0);
				// Increment counters according to what got filtered
				for(size_t mate = 0; mate < (pair ? 2:1); mate++) {
					if(!filt[mate]) {
						// Mate was rejected by N filter
						olm.freads++;               // reads filtered out
						olm.fbases += rdlens[mate]; // bases filtered out
					} else {
						//shs[mate].clear();
						//shs[mate].nextRead(mate == 0 ? ps->bufa() : ps->bufb());
						//assert(shs[mate].empty());
						olm.ureads++;               // reads passing filter
						olm.ubases += rdlens[mate]; // bases passing filter
					}
				}
				//size_t eePeEeltLimit = std::numeric_limits<size_t>::max();
				// Whether we're done with mate1 / mate2
                bool done[2] = { !filt[0], !filt[1] };
				// size_t nelt[2] = {0, 0};
                if(filt[0] && filt[1]) {
                    splicedAligner.initReads(rds, nofw, norc, minsc, maxpen);
                } else if(filt[0]) {
                    splicedAligner.initRead(rds[0], nofw[0], norc[0], minsc[0], maxpen[0], false);
                } else if(filt[1]) {
                    splicedAligner.initRead(rds[1], nofw[1], norc[1], minsc[1], maxpen[1], true);
                }
                if(filt[0] || filt[1]) {
                    int ret = splicedAligner.go(
                                                sc,
                                                pepol,
                                                *multiseed_tpol,
                                                *gpol,
                                                gfm,
                                                rgfm,
                                                *altdb,
                                                *repeatdb,
                                                *raltdb,
                                                ref,
                                                rref,
                                                sw,
                                                *ssdb,
                                                wlm,
                                                prm,
                                                swmSeed,
                                                him,
                                                rnd,
                                                msinkwrap);
                    MERGE_SW(sw);
                    // daehwan
                    size_t mate = 0;
                    
                    assert_gt(ret, 0);
                    // Clear out the exact hits so that we don't try to
                    // extend them again later!
                    if(ret == EXTEND_EXHAUSTED_CANDIDATES) {
                        // Not done yet
                    } else if(ret == EXTEND_POLICY_FULFILLED) {
                        // Policy is satisfied for this mate at least
                        if(msinkwrap.state().doneWithMate(mate == 0)) {
                            done[mate] = true;
                        }
                        if(msinkwrap.state().doneWithMate(mate == 1)) {
                            done[mate^1] = true;
                        }
                    } else if(ret == EXTEND_PERFECT_SCORE) {
                        // We exhausted this mode at least
                        done[mate] = true;
                    } else if(ret == EXTEND_EXCEEDED_HARD_LIMIT) {
                        // We exceeded a per-read limit
                        done[mate] = true;
                    } else if(ret == EXTEND_EXCEEDED_SOFT_LIMIT) {
                        // Not done yet
                    } else {
                        //
                        cerr << "Bad return value: " << ret << endl;
                        throw 1;
                    }
                    if(!done[mate]) {
                        TAlScore perfectScore = sc.perfectScore(rdlens[mate]);
                        if(!done[mate] && minsc[mate] == perfectScore) {
                            done[mate] = true;
                        }
                    }
                }

                for(size_t i = 0; i < 2; i++) {
                    assert_leq(prm.nExIters, mxIter[i]);
                    assert_leq(prm.nExDps,   mxDp[i]);
                    assert_leq(prm.nMateDps, mxDp[i]);
                    assert_leq(prm.nExUgs,   mxUg[i]);
                    assert_leq(prm.nMateUgs, mxUg[i]);
                    assert_leq(prm.nDpFail,  streak[i]);
                    assert_leq(prm.nUgFail,  streak[i]);
                    assert_leq(prm.nEeFail,  streak[i]);
                }
                
				// Commit and report paired-end/unpaired alignments
				msinkwrap.finishRead(
                                     NULL,
                                     NULL,
                                     exhaustive[0],        // exhausted seed hits for mate 1?
                                     exhaustive[1],        // exhausted seed hits for mate 2?
                                     nfilt[0],
                                     nfilt[1],
                                     scfilt[0],
                                     scfilt[1],
                                     lenfilt[0],
                                     lenfilt[1],
                                     qcfilt[0],
                                     qcfilt[1],
                                     sortByScore,          // prioritize by alignment score
                                     rnd,                  // pseudo-random generator
                                     rpm,                  // reporting metrics
                                     prm,                  // per-read metrics
                                     sc,                   // scoring scheme
                                     !seedSumm,            // suppress seed summaries?
                                     seedSumm,             // suppress alignments?
                                     templateLenAdjustment);
				assert(!retry || msinkwrap.empty());
			} // while(retry)
		} // if(rdid >= skipReads && rdid < qUpto)
		else if(rdid >= qUpto) {
			break;
		}
		if(metricsPerRead) {
			MERGE_METRICS(metricsPt, nthreads > 1);
			nametmp = ps->bufa().name;
			metricsPt.reportInterval(
                                     metricsOfb, metricsStderr, true, true, &nametmp);
			metricsPt.reset();
		}
	} // while(true)
	
	// One last metrics merge
	MERGE_METRICS(metrics, nthreads > 1);
    
	return;
}

/**
 * Called once per alignment job.  Sets up global pointers to the
 * shared global data structures, creates per-thread structures, then
 * enters the search loop.
 */
static void multiseedSearch(
                            Scoring& sc,
                            TranscriptomePolicy& tpol,
                            GraphPolicy& gp,
                            PairedPatternSource& patsrc,  // pattern source
                            AlnSink<index_t>& msink,      // hit sink
                            HGFM<index_t>& gfm,           // index of original text
                            RFM<index_t>* rgfm,           // index of repeat sequences
                            BitPairReference* refs,       // base reference
                            BitPairReference* rrefs,      // repeat reference
                            OutFileBuf *metricsOfb)
{
    multiseed_patsrc       = &patsrc;
	multiseed_msink        = &msink;
	multiseed_gfm          = &gfm;
    multiseed_rgfm         = rgfm;
	multiseed_sc           = &sc;
    multiseed_tpol         = &tpol;
    gpol                   = &gp;
	multiseed_metricsOfb   = metricsOfb;
	multiseed_refs         = refs;
    multiseed_rrefs        = rrefs;
	AutoArray<tthread::thread*> threads(nthreads);
	AutoArray<int> tids(nthreads);	
	// Start the metrics thread
	{
		Timer _t(cerr, "Multiseed full-index search: ", timing);
        
        thread_rids.resize(nthreads);
        thread_rids.fill(0);
        thread_rids_mindist = (nthreads == 1 || !useTempSpliceSite ? 0 : 1000 * nthreads);        
		for(int i = 0; i < nthreads; i++) {
			// Thread IDs start at 1
			tids[i] = i+1;
            threads[i] = new tthread::thread(multiseedSearchWorker_hisat2, (void*)&tids[i]);
		}

        for (int i = 0; i < nthreads; i++)
            threads[i]->join();

	}
	if(!metricsPerRead && (metricsOfb != NULL || metricsStderr)) {
		metrics.reportInterval(metricsOfb, metricsStderr, true, false, NULL);
	}
}

static string argstr;

extern void initializeCntLut();
extern void initializeCntBit();

template<typename TStr>
static void driver(
	const char * type,
	const string& bt2indexBase,
	const string& outfile)
{
	if(gVerbose || startVerbose)  {
		cerr << "Entered driver(): "; logTime(cerr, true);
	}
    
    initializeCntLut();
    initializeCntBit();
    
	// Vector of the reference sequences; used for sanity-checking
	EList<SString<char> > names, os;
	EList<size_t> nameLens, seqLens;
	// Read reference sequences from the command-line or from a FASTA file
	if(!origString.empty()) {
		// Read fasta file(s)
		EList<string> origFiles;
		tokenize(origString, ",", origFiles);
		parseFastas(origFiles, names, nameLens, os, seqLens);
	}
	PatternParams pp(
		format,        // file format
		fileParallel,  // true -> wrap files with separate PairedPatternSources
		seed,          // pseudo-random seed
		useSpinlock,   // use spin locks instead of pthreads
		solexaQuals,   // true -> qualities are on solexa64 scale
		phred64Quals,  // true -> qualities are on phred64 scale
		integerQuals,  // true -> qualities are space-separated numbers
		fuzzy,         // true -> try to parse fuzzy fastq
		fastaContLen,  // length of sampled reads for FastaContinuous...
		fastaContFreq, // frequency of sampled reads for FastaContinuous...
		skipReads      // skip the first 'skip' patterns
	);
	if(gVerbose || startVerbose) {
		cerr << "Creating PatternSource: "; logTime(cerr, true);
	}
	PairedPatternSource *patsrc = PairedPatternSource::setupPatternSources(
		queries,     // singles, from argv
		mates1,      // mate1's, from -1 arg
		mates2,      // mate2's, from -2 arg
		mates12,     // both mates on each line, from --12 arg
#ifdef USE_SRA
        sra_accs,    // SRA accessions
#endif
		qualities,   // qualities associated with singles
		qualities1,  // qualities associated with m1
		qualities2,  // qualities associated with m2
		pp,          // read read-in parameters
        nthreads,
		gVerbose || startVerbose); // be talkative
	// Open hit output file
	if(gVerbose || startVerbose) {
		cerr << "Opening hit output file: "; logTime(cerr, true);
	}
	OutFileBuf *fout;
	if(!outfile.empty()) {
		fout = new OutFileBuf(outfile.c_str(), false);
	} else {
		fout = new OutFileBuf();
	}
	// Initialize GFM object and read in header
	if(gVerbose || startVerbose) {
		cerr << "About to initialize fw GFM: "; logTime(cerr, true);
	}
    altdb = new ALTDB<index_t>();
    repeatdb = new RepeatDB<index_t>();
    raltdb = new ALTDB<index_t>();
	adjIdxBase = adjustEbwtBase(argv0, bt2indexBase, gVerbose);
	HGFM<index_t, local_index_t> gfm(
                                     adjIdxBase,
                                     altdb,
                                     NULL,
                                     NULL,
                                     -1,       // fw index
                                     true,     // index is for the forward direction
                                     /* overriding: */ offRate,
                                     0, // amount to add to index offrate or <= 0 to do nothing
                                     useMm,    // whether to use memory-mapped files
                                     useShmem, // whether to use shared memory
                                     mmSweep,  // sweep memory-mapped files
                                     !noRefNames, // load names?
                                     true,        // load SA sample?
                                     true,        // load ftab?
                                     true,        // load rstarts?
                                     !no_spliced_alignment, // load splice sites?
                                     gVerbose, // whether to be talkative
                                     startVerbose, // talkative during initialization
                                     false /*passMemExc*/,
                                     sanityCheck,
                                     use_haplotype); //use haplotypes?
	if(sanityCheck && !os.empty()) {
		// Sanity check number of patterns and pattern lengths in GFM
		// against original strings
		assert_eq(os.size(), gfm.nPat());
		for(size_t i = 0; i < os.size(); i++) {
			assert_eq(os[i].length(), gfm.plen()[i]);
		}
	}
	// Sanity-check the restored version of the GFM
	if(sanityCheck && !os.empty()) {
		gfm.loadIntoMemory(
			-1, // fw index
			true, // load SA sample
			true, // load ftab
			true, // load rstarts
			!noRefNames,
			startVerbose);
		gfm.checkOrigs(os, false);
		gfm.evictFromMemory();
	}
    {
        // Load the other half of the index into memory
        assert(!gfm.isInMemory());
        Timer _t(cerr, "Time loading forward index: ", timing);
        gfm.loadIntoMemory(
                           -1, // not the reverse index
                           true,         // load SA samp? (yes, need forward index's SA samp)
                           true,         // load ftab (in forward index)
                           true,         // load rstarts (in forward index)
                           !noRefNames,  // load names?
                           startVerbose);
    }
    RFM<index_t>* rgfm = NULL;
    string rep_adjIdxBase = adjIdxBase + ".rep";
    bool rep_index_exists = false;
    {
        std::ifstream infile((rep_adjIdxBase + ".1." + gfm_ext.c_str()).c_str());
        rep_index_exists = infile.good();
    }
    if(rep_index_exists && use_repeat_index) {
        rgfm = new RFM<index_t>(
                                rep_adjIdxBase,
                                raltdb,
                                repeatdb,
                                &readLens,
                                -1,       // fw index
                                true,     // index is for the forward direction
                                /* overriding: */ offRate,
                                0, // amount to add to index offrate or <= 0 to do nothing
                                useMm,    // whether to use memory-mapped files
                                useShmem, // whether to use shared memory
                                mmSweep,  // sweep memory-mapped files
                                !noRefNames, // load names?
                                true,        // load SA sample?
                                true,        // load ftab?
                                true,        // load rstarts?
                                !no_spliced_alignment, // load splice sites?
                                gVerbose, // whether to be talkative
                                startVerbose, // talkative during initialization
                                false /*passMemExc*/,
                                sanityCheck,
                                false); //use haplotypes?

        // CP to do
#if 0
        if(sanityCheck && !os.empty()) {
            // Sanity check number of patterns and pattern lengths in GFM
            // against original strings
            assert_eq(os.size(), gfm.nPat());
            for(size_t i = 0; i < os.size(); i++) {
                assert_eq(os[i].length(), rgfm->plen()[i]);
            }
        }
        // Sanity-check the restored version of the GFM
        if(sanityCheck && !os.empty()) {
            rgfm->loadIntoMemory(
                               -1, // fw index
                               true, // load SA sample
                               true, // load ftab
                               true, // load rstarts
                               !noRefNames,
                               startVerbose);
            rgfm->checkOrigs(os, false);
            rgfm->evictFromMemory();
        }
#endif
        {
            // Load the other half of the index into memory
            assert(!rgfm->isInMemory());
            Timer _t(cerr, "Time loading forward index: ", timing);
            rgfm->loadIntoMemory(
                                 -1, // not the reverse index
                                 true,         // load SA samp? (yes, need forward index's SA samp)
                                 true,         // load ftab (in forward index)
                                 true,         // load rstarts (in forward index)
                                 !noRefNames,  // load names?
                                 startVerbose);
            
            repeatdb->construct(gfm.rstarts(), gfm.nFrag());
        }
    }
    
    if(!saw_k) {
        if(gfm.gh().linearFM()) khits = 5;
        else                    khits = 10;
    }
	OutputQueue oq(
		*fout,                   // out file buffer
		reorder && nthreads > 1, // whether to reorder when there's >1 thread
		nthreads,                // # threads
		nthreads > 1,            // whether to be thread-safe
		skipReads);              // first read will have this rdid
	{
		Timer _t(cerr, "Time searching: ", timing);
		// Set up penalities
		if(bonusMatch > 0 && !localAlign) {
			cerr << "Warning: Match bonus always = 0 in --end-to-end mode; ignoring user setting" << endl;
			bonusMatch = 0;
		}
        if(tranAssm) {
            penNoncanIntronLen.init(SIMPLE_FUNC_LOG, -8, 2);
        }
		Scoring sc(
                   bonusMatch,     // constant reward for match
                   penMmcType,     // how to penalize mismatches
                   penMmcMax,      // max mm penalty
                   penMmcMin,      // min mm penalty
                   penScMax,       // max sc penalty
                   penScMin,       // min sc penalty
                   scoreMin,       // min score as function of read len
                   nCeil,          // max # Ns as function of read len
                   penNType,       // how to penalize Ns in the read
                   penN,           // constant if N pelanty is a constant
                   penNCatPair,    // whether to concat mates before N filtering
                   penRdGapConst,  // constant coeff for read gap cost
                   penRfGapConst,  // constant coeff for ref gap cost
                   penRdGapLinear, // linear coeff for read gap cost
                   penRfGapLinear, // linear coeff for ref gap cost
                   gGapBarrier,    // # rows at top/bot only entered diagonally
                   penCanSplice,   // canonical splicing penalty
                   penNoncanSplice,// non-canonical splicing penalty
                   penConflictSplice, // conflicting splice site penalty
                   &penCanIntronLen,      // penalty as to intron length
                   &penNoncanIntronLen);  // penalty as to intron length
        
		EList<size_t> reflens;
		for(size_t i = 0; i < gfm.nPat(); i++) {
			reflens.push_back(gfm.plen()[i]);
		}
		EList<string> refnames;
		readEbwtRefnames<index_t>(adjIdxBase, refnames);
        EList<size_t> replens;
        EList<string> repnames;
        if(rep_index_exists && use_repeat_index) {
            rgfm->getReferenceNames(repnames);
            rgfm->getReferenceLens(replens);
        }
        if(rmChrName && addChrName) {
            cerr << "Error: --remove-chrname and --add-chrname cannot be used at the same time" << endl;
            throw 1;
        }
        if(rmChrName) {
            for(size_t i = 0; i < refnames.size(); i++) {
                string& refname = refnames[i];
                if(refname.find("chr") == 0) {
                    refname = refname.substr(3);
                }
            }
        } else if(addChrName) {
            for(size_t i = 0; i < refnames.size(); i++) {
                string& refname = refnames[i];
                if(refname.find("chr") != 0) {
                    refname = string("chr") + refname;
                }
            }
        }

        EList<size_t> empty_replens;
        EList<string> empty_repnames;
		SamConfig<index_t> samc(
			refnames,               // reference sequence names
			reflens,                // reference sequence lengths
            repeat ? repnames : empty_repnames, // repeat sequence names
            repeat ? replens : empty_replens,   // repeat sequence lengths
			samTruncQname,          // whether to truncate QNAME to 255 chars
			samOmitSecSeqQual,      // omit SEQ/QUAL for 2ndary alignments?
			samNoUnal,              // omit unaligned-read records?
			string("hisat2"),       // program id
			string("hisat2"),       // program name
			string(HISAT2_VERSION), // program version
			argstr,                 // command-line
			rgs_optflag,            // read-group string
            rna_strandness,
			sam_print_as,
			sam_print_xs,
			sam_print_xss,
			sam_print_yn,
			sam_print_xn,
			sam_print_cs,
			sam_print_cq,
			sam_print_x0,
			sam_print_x1,
			sam_print_xm,
			sam_print_xo,
			sam_print_xg,
			sam_print_nm,
			sam_print_md,
			sam_print_yf,
			sam_print_yi,
			sam_print_ym,
			sam_print_yp,
			sam_print_yt,
			sam_print_ys,
			sam_print_zs,
			sam_print_xr,
			sam_print_xt,
			sam_print_xd,
			sam_print_xu,
			sam_print_yl,
			sam_print_ye,
			sam_print_yu,
			sam_print_xp,
			sam_print_yr,
			sam_print_zb,
			sam_print_zr,
			sam_print_zf,
			sam_print_zm,
			sam_print_zi,
			sam_print_zp,
			sam_print_zu,
            sam_print_xs_a,
            sam_print_nh);
		// Set up hit sink; if sanityCheck && !os.empty() is true,
		// then instruct the sink to "retain" hits in a vector in
		// memory so that we can easily sanity check them later on
		AlnSink<index_t> *mssink = NULL;
        Timer *_tRef = new Timer(cerr, "Time loading reference: ", timing);
        auto_ptr<BitPairReference> refs(
                                        new BitPairReference(
                                                             adjIdxBase,
                                                             NULL,
                                                             false,
                                                             sanityCheck,
                                                             NULL,
                                                             NULL,
                                                             false,
                                                             useMm,
                                                             useShmem,
                                                             mmSweep,
                                                             gVerbose,
                                                             startVerbose)
                                        );
        delete _tRef;
        if(!refs->loaded()) throw 1;
        
        BitPairReference* rrefs = NULL;
        if(rep_index_exists && use_repeat_index) {
            const EList<uint8_t>& included = rgfm->getRepeatIncluded();
            rrefs = new BitPairReference(
                                         rep_adjIdxBase,
                                         &included,
                                         false,
                                         sanityCheck,
                                         NULL,
                                         NULL,
                                         false,
                                         useMm,
                                         useShmem,
                                         mmSweep,
                                         gVerbose,
                                         startVerbose);
            if(!rrefs->loaded()) throw 1;
        }
        
        bool xsOnly = (tranAssm_program == "cufflinks");
        TranscriptomePolicy tpol(minIntronLen,
                                 maxIntronLen,
                                 tranAssm ? 15 : 7,
                                 tranAssm ? 20 : 14,
                                 no_spliced_alignment,
                                 tranMapOnly,
                                 tranAssm,
                                 xsOnly,
                                 avoid_pseudogene);
        
        GraphPolicy gpol(max_alts_tried,
                         use_haplotype,
                         altdb->haplotypes().size() > 0 && use_haplotype,
                         enable_codis);
        
        init_junction_prob();
        bool write = novelSpliceSiteOutfile != "" || useTempSpliceSite;
        bool read = knownSpliceSiteInfile != "" || novelSpliceSiteInfile != "" || useTempSpliceSite || altdb->hasSpliceSites();
        ssdb = new SpliceSiteDB(
                                *(refs.get()),
                                refnames,
                                nthreads > 1, // thread-safe
                                write, // write?
                                read);  // read?
        ssdb->read(gfm, altdb->alts());
        if(knownSpliceSiteInfile != "") {
            ifstream ssdb_file(knownSpliceSiteInfile.c_str(), ios::in);
            if(ssdb_file.is_open()) {
                ssdb->read(ssdb_file,
                           true); // known splice sites
                ssdb_file.close();
            }
        }
        if(novelSpliceSiteInfile != "") {
            ifstream ssdb_file(novelSpliceSiteInfile.c_str(), ios::in);
            if(ssdb_file.is_open()) {
                ssdb->read(ssdb_file,
                           false); // novel splice sites
                ssdb_file.close();
            }
        }
		switch(outType) {
			case OUTPUT_SAM: {
				mssink = new AlnSinkSam<index_t>(
                                                 oq,           // output queue
                                                 samc,         // settings & routines for SAM output
                                                 refnames,     // reference names
                                                 repnames,     // repeat names
                                                 gQuiet,       // don't print alignment summary at end
                                                 altdb,
                                                 ssdb);
				if(!samNoHead) {
					bool printHd = true, printSq = true;
					BTString buf;
					samc.printHeader(buf, rgid, rgs, printHd, !samNoSQ, printSq);
					fout->writeString(buf);
				}
				break;
			}
			default:
				cerr << "Invalid output type: " << outType << endl;
				throw 1;
		}
		if(gVerbose || startVerbose) {
			cerr << "Dispatching to search driver: "; logTime(cerr, true);
		}
		// Set up global constraint
		OutFileBuf *metricsOfb = NULL;
		if(!metricsFile.empty() && metricsIval > 0) {
			metricsOfb = new OutFileBuf(metricsFile);
		}
		// Do the search for all input reads
		assert(patsrc != NULL);
		assert(mssink != NULL);
		multiseedSearch(
                        sc,      // scoring scheme
                        tpol,
                        gpol,
                        *patsrc, // pattern source
                        *mssink, // hit sink
                        gfm,     // BWT
                        rgfm,
                        refs.get(),
                        rrefs,
                        metricsOfb);
		// Evict any loaded indexes from memory
		if(gfm.isInMemory()) {
			gfm.evictFromMemory();
		}
		if(!gQuiet && !seedSumm) {
			size_t repThresh = mhits;
			if(repThresh == 0) {
				repThresh = std::numeric_limits<size_t>::max();
			}
			mssink->finish(cerr,
                           repThresh,
                           gReportDiscordant,
                           gReportMixed,
                           newAlignSummary,
                           hadoopOut);
            if(alignSumFile != "") {
                ofstream sumfile(alignSumFile.c_str(), ios::out);
                if(sumfile.is_open()) {
                    mssink->finish(sumfile,
                                   repThresh,
                                   gReportDiscordant,
                                   gReportMixed,
                                   newAlignSummary,
                                   false); // hadoopOut
                    sumfile.close();
                }
            }
		}
        if(ssdb != NULL) {
            if(novelSpliceSiteOutfile != "") {
                ofstream ssdb_file(novelSpliceSiteOutfile.c_str(), ios::out);
                if(ssdb_file.is_open()) {
                    ssdb->print(ssdb_file);
                    ssdb_file.close();
                }
            }
        }
		oq.flush(true);
		assert_eq(oq.numStarted(), oq.numFinished());
		assert_eq(oq.numStarted(), oq.numFlushed());
		delete patsrc;
		delete mssink;
        delete altdb;
        delete repeatdb;
        delete raltdb;
        delete ssdb;
		delete metricsOfb;
        delete rgfm;
        delete rrefs;
		if(fout != NULL) {
			delete fout;
		}
	}
}

// C++ name mangling is disabled for the bowtie() function to make it
// easier to use Bowtie as a library.
extern "C" {

/**
 * Main bowtie entry function.  Parses argc/argv style command-line
 * options, sets global configuration variables, and calls the driver()
 * function.
 */
int hisat2(int argc, const char **argv) {
	try {
		// Reset all global state, including getopt state
		opterr = optind = 1;
		resetOptions();
		for(int i = 0; i < argc; i++) {
			argstr += argv[i];
			if(i < argc-1) argstr += " ";
		}
		if(startVerbose) { cerr << "Entered main(): "; logTime(cerr, true); }
		parseOptions(argc, argv);
		argv0 = argv[0];
		if(showVersion) {
			cout << argv0 << " version " << HISAT2_VERSION << endl;
			if(sizeof(void*) == 4) {
				cout << "32-bit" << endl;
			} else if(sizeof(void*) == 8) {
				cout << "64-bit" << endl;
			} else {
				cout << "Neither 32- nor 64-bit: sizeof(void*) = " << sizeof(void*) << endl;
			}
			cout << "Built on " << BUILD_HOST << endl;
			cout << BUILD_TIME << endl;
			cout << "Compiler: " << COMPILER_VERSION << endl;
			cout << "Options: " << COMPILER_OPTIONS << endl;
			cout << "Sizeof {int, long, long long, void*, size_t, off_t}: {"
				 << sizeof(int)
				 << ", " << sizeof(long) << ", " << sizeof(long long)
				 << ", " << sizeof(void *) << ", " << sizeof(size_t)
				 << ", " << sizeof(off_t) << "}" << endl;
			return 0;
		}
		{
			Timer _t(cerr, "Overall time: ", timing);
			if(startVerbose) {
				cerr << "Parsing index and read arguments: "; logTime(cerr, true);
			}

			// Get index basename (but only if it wasn't specified via --index)
			if(bt2index.empty()) {
				if(optind >= argc) {
					cerr << "No index, query, or output file specified!" << endl;
					printUsage(cerr);
					return 1;
				}
				bt2index = argv[optind++];
			}

			// Get query filename
			bool got_reads = !queries.empty() || !mates1.empty() || !mates12.empty();
#ifdef USE_SRA
            got_reads = got_reads || !sra_accs.empty();
#endif
            if(minIntronLen > maxIntronLen) {
                cerr << "--min-intronlen(" << minIntronLen << ") should not be greater than --max-intronlen("
                     << maxIntronLen << ")" << endl;
                printUsage(cerr);
                return 1;
            }
			if(optind >= argc) {
				if(!got_reads) {
					printUsage(cerr);
					cerr << "***" << endl
#ifdef USE_SRA
					     << "Error: Must specify at least one read input with -U/-1/-2/--sra-acc" << endl;
#else
                    << "Error: Must specify at least one read input with -U/-1/-2" << endl;

#endif
					return 1;
				}
			} else if(!got_reads) {
				// Tokenize the list of query files
				tokenize(argv[optind++], ",", queries);
				if(queries.empty()) {
					cerr << "Tokenized query file list was empty!" << endl;
					printUsage(cerr);
					return 1;
				}
			}

			// Get output filename
			if(optind < argc && outfile.empty()) {
				outfile = argv[optind++];
				cerr << "Warning: Output file '" << outfile.c_str()
				     << "' was specified without -S.  This will not work in "
					 << "future HISAT 2 versions.  Please use -S instead."
					 << endl;
			}

			// Extra parametesr?
			if(optind < argc) {
				cerr << "Extra parameter(s) specified: ";
				for(int i = optind; i < argc; i++) {
					cerr << "\"" << argv[i] << "\"";
					if(i < argc-1) cerr << ", ";
				}
				cerr << endl;
				if(mates1.size() > 0) {
					cerr << "Note that if <mates> files are specified using -1/-2, a <singles> file cannot" << endl
						 << "also be specified.  Please run HISAT2 separately for mates and singles." << endl;
				}
				throw 1;
			}

			// Optionally summarize
			if(gVerbose) {
				cout << "Input bt2 file: \"" << bt2index.c_str() << "\"" << endl;
				cout << "Query inputs (DNA, " << file_format_names[format].c_str() << "):" << endl;
				for(size_t i = 0; i < queries.size(); i++) {
					cout << "  " << queries[i].c_str() << endl;
				}
				cout << "Quality inputs:" << endl;
				for(size_t i = 0; i < qualities.size(); i++) {
					cout << "  " << qualities[i].c_str() << endl;
				}
				cout << "Output file: \"" << outfile.c_str() << "\"" << endl;
				cout << "Local endianness: " << (currentlyBigEndian()? "big":"little") << endl;
				cout << "Sanity checking: " << (sanityCheck? "enabled":"disabled") << endl;
			#ifdef NDEBUG
				cout << "Assertions: disabled" << endl;
			#else
				cout << "Assertions: enabled" << endl;
			#endif
			}
			if(ipause) {
				cout << "Press key to continue..." << endl;
				getchar();
			}
			driver<SString<char> >("DNA", bt2index, outfile);
		}
		return 0;
	} catch(std::exception& e) {
		cerr << "Error: Encountered exception: '" << e.what() << "'" << endl;
		cerr << "Command: ";
		for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
		cerr << endl;
		return 1;
	} catch(int e) {
		if(e != 0) {
			cerr << "Error: Encountered internal HISAT2 exception (#" << e << ")" << endl;
			cerr << "Command: ";
			for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
			cerr << endl;
		}
		return e;
	}
} // bowtie()
} // extern "C"
