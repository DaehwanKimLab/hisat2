#!/usr/bin/env python

#
# Copyright 2017, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT-genotype.
#
# HISAT-genotype is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT-genotype is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT-genotype.  If not, see <http://www.gnu.org/licenses/>.
#


import sys, os, subprocess, re, resource
import inspect, random
import math
from datetime import datetime, date, time
from argparse import ArgumentParser, FileType
import hisatgenotype_typing_common as typing_common


"""
Align reads, and sort the alignments into a BAM file
"""
def align_reads(base_fname,
                read_fnames,
                fastq,
                threads,
                verbose):
    aligner_cmd = ["hisat2",
                   "--no-unal",
                   "-p", str(threads),
                   "--no-spliced-alignment",
                   "--max-altstried", "64"]
    # aligner_cmd += ["--mm"]
    aligner_cmd += ["-x", "%s" % base_fname]

    assert len(read_fnames) > 0
    if not fastq:
        aligner_cmd += ["-f"]
    single = len(read_fnames) == 1
    if single:
        aligner_cmd += ["-U", read_fnames[0]]
    else:
        aligner_cmd += ["-1", read_fnames[0],
                        "-2", read_fnames[1]]

    out_base_fname = read_fnames[0].split('/')[-1].split('.')[0]

    print >> sys.stderr, "%s Aligning %s to %s ..." % (str(datetime.now()), ' '.join(read_fnames), base_fname)
    if verbose:
        print >> sys.stderr, "\t%s" % (' '.join(aligner_cmd))

    align_proc = subprocess.Popen(aligner_cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=open("/dev/null", 'w'))

    unsorted_bam_fname = "%s_unsorted.bam" % out_base_fname
    sambam_cmd = ["samtools",
                  "view",
                  "-bS",
                  "-"]
    sambam_proc = subprocess.Popen(sambam_cmd,
                                   stdin=align_proc.stdout,
                                   stdout=open(unsorted_bam_fname, 'w'))
    sambam_proc.communicate()

    # Increase the maximum number of files that can be opened
    resource.setrlimit(resource.RLIMIT_NOFILE, (10000, 10240))
    
    print >> sys.stderr, "%s Sorting %s ..." % (str(datetime.now()), unsorted_bam_fname)
    bam_fname = "%s.bam" % out_base_fname
    bamsort_cmd = ["samtools",
                   "sort",
                   "--threads", str(threads),
                   "-m", "1536M",
                   unsorted_bam_fname,
                   "-o", bam_fname]    
    if verbose:
        print >> sys.stderr, "\t%s" % ' '.join(bamsort_cmd)
    bamsort_proc = subprocess.call(bamsort_cmd)
    os.remove(unsorted_bam_fname)

    index_bam(bam_fname,
              verbose)
    
    return bam_fname


"""
"""
def index_bam(bam_fname,
              verbose):
    print >> sys.stderr, "%s Indexing %s ..." % (str(datetime.now()), bam_fname)
    bamindex_cmd = ["samtools",
                    "index",
                    bam_fname]
    if verbose:
        print >> sys.stderr, "\t%s" % ' '.join(bamindex_cmd)
    bamindex_proc = subprocess.call(bamindex_cmd)


"""
"""
def extract_reads(bam_fname,
                  chr,
                  left,
                  right,
                  read_base_fname, # sample => sample.1.fq.gz and sample.2.fq.gz
                  paired,
                  fastq,
                  verbose):
    out_read_dname = "hisatgenotype_out"
    if not os.path.exists(out_read_dname):
        os.mkdir(out_read_dname)
        
    read_fnames = []
    if paired:
        read_fnames = [out_read_dname + "/" + read_base_fname + ".1.fq.gz",
                       out_read_dname + "/" + read_base_fname + ".2.fq.gz"]
    else:
        read_fnames = [out_read_dname + "/" + read_base_fname + ".fq.gz"]

    if paired:
        gzip1_proc = subprocess.Popen(["gzip"],
                                      stdin=subprocess.PIPE,
                                      stdout=open(read_fnames[0], 'w'),
                                      stderr=open("/dev/null", 'w'))

        gzip2_proc = subprocess.Popen(["gzip"],
                                      stdin=subprocess.PIPE,
                                      stdout=open(read_fnames[1], 'w'),
                                      stderr=open("/dev/null", 'w'))
    else:
        gzip1_proc = subprocess.Popen(["gzip"],
                                      stdin=subprocess.PIPE,
                                      stdout=open(read_fnames[0], 'w'),
                                      stderr=open("/dev/null", 'w'))

    def write_read(gzip_proc, read_name, seq, qual):
        if fastq:
            gzip_proc.stdin.write("@%s\n" % read_name)
            gzip_proc.stdin.write("%s\n" % seq)
            gzip_proc.stdin.write("+\n")
            gzip_proc.stdin.write("%s\n" % qual)
        else:
            gzip_proc.stdin.write(">%s\n" % prev_read_name)
            gzip_proc.stdin.write("%s\n" % seq)                    

    bamview_cmd = ["samtools", "view", bam_fname, "%s:%d-%d" % (chr, left+1, right+1)]
    if verbose:
        print >> sys.stderr, "\t%s" % ' '.join(bamview_cmd)
    bamview_proc = subprocess.Popen(bamview_cmd,
                                    stdout=subprocess.PIPE,
                                    stderr=open("/dev/null", 'w'))

    sort_read_cmd = ["sort", "-k", "1,1", "-s"] # -s for stable sorting
    alignview_proc = subprocess.Popen(sort_read_cmd,
                                      stdin=bamview_proc.stdout,
                                      stdout=subprocess.PIPE,
                                      stderr=open("/dev/null", 'w'))

    prev_read_name, extract_read, read1, read2 = "", False, [], []
    for line in alignview_proc.stdout:
        if line.startswith('@'):
            continue
        line = line.strip()
        cols = line.split()
        read_name, flag, chr, pos, mapQ, cigar, _, _, _, read, qual = cols[:11]
        flag, pos = int(flag), int(pos)
        strand = '-' if flag & 0x10 else '+'                   
        AS, NH = "", ""
        for i in range(11, len(cols)):
            col = cols[i]
            if col.startswith("AS"):
                AS = int(col[5:])
            elif col.startswith("NH"):
                NH = int(col[5:])

        # DK - check this out
        simulation = True
        if (not simulation and read_name != prev_read_name) or \
           (simulation and read_name.split('|')[0] != prev_read_name.split('|')[0]):
            if extract_read:
                if paired:
                    if len(read1) == 2 and len(read2) == 2:
                        write_read(gzip1_proc, prev_read_name, read1[0], read1[1])
                        write_read(gzip2_proc, prev_read_name, read2[0], read2[1])
                else:                    
                    write_read(gzip1_proc, prev_read_name, read1[0], read1[1])
            prev_read_name, extract_read, read1, read2 = read_name, False, [], []

        if NH == 1:
            extract_read = True

        if flag & 0x40 or not paired: # left read
            if not read1:
                if flag & 0x10: # reverse complement
                    read1 = [typing_common.reverse_complement(read), qual[::-1]]
                else:
                    read1 = [read, qual]
        else:
            assert flag & 0x80 # right read
            if flag & 0x10: # reverse complement
                read2 = [typing_common.reverse_complement(read), qual[::-1]]
            else:
                read2 = [read, qual]

    if extract_read:
        if paired:
            if len(read1) == 2 and len(read2) == 2:
                write_read(gzip1_proc, prev_read_name, read1[0], read1[1])
                write_read(gzip2_proc, prev_read_name, read2[0], read2[1])
        else:                    
            write_read(gzip1_proc, prev_read_name, read1[0], read1[1])

    gzip1_proc.stdin.close()
    if paired:
        gzip2_proc.stdin.close()

    return read_fnames


"""
"""
def perform_genotyping(base_fname,
                       database,
                       locus_list,
                       read_fnames,
                       fastq,
                       num_editdist,
                       assembly,
                       local_database,
                       threads,
                       verbose):
    genotype_cmd = ["hisatgenotype_locus.py"]
    if not local_database:
        genotype_cmd += ["--genotype-genome", base_fname]
    genotype_cmd += ["--base", database]
    if len(locus_list) > 0:
        genotype_cmd += ["--locus-list", ','.join(locus_list)]
    genotype_cmd += ["-p", str(threads),
                     "--num-editdist", str(num_editdist)]
    if not fastq:
        genotype_cmd += ["-f"]

    if len(read_fnames) == 2: # paired
        genotype_cmd += ["-1", read_fnames[0],
                         "-2", read_fnames[1]]
    elif len(read_fnames) == 1:
        genotype_cmd += ["-U", read_fnames[0]] 
    else:
        assert len(read_fnames) == 0

    if assembly:
        genotype_cmd += ["--assembly"]

    if verbose:
        print >> sys.stderr, "\t%s" % ' '.join(genotype_cmd)
    genotype_proc = subprocess.Popen(genotype_cmd)
    genotype_proc.communicate()
        

"""
"""
def genotype(base_fname,
             target_region_list,
             fastq,
             read_fnames,
             alignment_fname,
             threads,
             num_editdist,
             assembly,
             local_database,
             verbose,
             debug):
    # variants, backbone sequence, and other sequeces
    genotype_fnames = ["%s.fa" % base_fname,
                       "%s.locus" % base_fname,
                       "%s.snp" % base_fname,
                       "%s.index.snp" % base_fname,
                       "%s.haplotype" % base_fname,
                       "%s.link" % base_fname,
                       "%s.coord" % base_fname,
                       "%s.clnsig" % base_fname]
    # hisat2 graph index files
    genotype_fnames += ["%s.%d.ht2" % (base_fname, i+1) for i in range(8)]
    if not typing_common.check_files(genotype_fnames):
        print >> sys.stderr, "Error: some of the following files are missing!"
        for fname in genotype_fnames:
            print >> sys.stderr, "\t%s" % fname
        sys.exit(1)

    # Read region alleles (names and sequences)
    regions, region_loci = {}, {}
    for line in open("%s.locus" % base_fname):
        family, allele_name, chr, left, right = line.strip().split()[:5]
        family = family.lower()
        if len(target_region_list) > 0 and \
           family not in target_region_list:
            continue
        
        locus_name = allele_name.split('*')[0]
        if family in target_region_list and \
           len(target_region_list[family]) > 0 and \
           locus_name not in target_region_list[family]:
            continue
        
        left, right = int(left), int(right)
        if family not in region_loci:
            region_loci[family] = []
        region_loci[family].append([locus_name, allele_name, chr, left, right])

    if len(region_loci) <= 0:
        print >> sys.stderr, "Warning: no region exists!"
        sys.exit(1)

    # Align reads, and sort the alignments into a BAM file
    if len(read_fnames) > 0:
        alignment_fname = align_reads(base_fname,
                                      read_fnames,
                                      fastq,
                                      threads,
                                      verbose)
    assert alignment_fname != "" and os.path.exists(alignment_fname)
    if not os.path.exists(alignment_fname + ".bai"):
        index_bam(alignment_fname,
                  verbose)
    assert os.path.exists(alignment_fname + ".bai")

    # Extract reads and perform genotyping
    for family, loci in region_loci.items():
        print >> sys.stderr, "Analyzing %s ..." % family.upper()
        for locus_name, allele_name, chr, left, right in loci:
            out_read_fname = "%s.%s" % (family, locus_name)
            if verbose:
                print >> sys.stderr, "\tExtracting reads beloning to %s-%s ..." % \
                    (family, locus_name)

            extracted_read_fnames = extract_reads(alignment_fname,
                                                  chr,
                                                  left,
                                                  right,
                                                  out_read_fname,
                                                  len(read_fnames) != 1, # paired?
                                                  fastq,
                                                  verbose)

            perform_genotyping(base_fname,
                               family,
                               [locus_name],
                               extracted_read_fnames,
                               fastq,
                               num_editdist,
                               assembly,
                               local_database,
                               threads,
                               verbose)
        print >> sys.stderr

    
                
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='HISAT-genotype')
    parser.add_argument("--base", "--base-name",
                        dest="base_fname",
                        type=str,
                        default="genotype_genome",
                        help="base filename for genotype genome")
    parser.add_argument("--region-list",
                        dest="region_list",
                        type=str,
                        default="",
                        help="A comma-separated list of regions (default: empty)")
    parser.add_argument("-f", "--fasta",
                        dest='fastq',
                        action='store_false',
                        help='FASTA file')    
    parser.add_argument("-U",
                        dest="read_fname_U",
                        type=str,
                        default="",
                        help="filename for single-end reads")
    parser.add_argument("-1",
                        dest="read_fname_1",
                        type=str,
                        default="",
                        help="filename for paired-end reads")
    parser.add_argument("-2",
                        dest="read_fname_2",
                        type=str,
                        default="",
                        help="filename for paired-end reads")
    parser.add_argument("--alignment-file",
                        dest="alignment_fname",
                        type=str,
                        default="",
                        help="Sorted BAM alignment file name")
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    parser.add_argument("--num-editdist",
                        dest="num_editdist",
                        type=int,
                        default=2,
                        help="Maximum number of mismatches per read alignment to be considered (default: 2)")
    parser.add_argument('--assembly',
                        dest='assembly',
                        action='store_true',
                        help='Perform assembly')
    parser.add_argument('--local-database',
                        dest='local_database',
                        action='store_true',
                        help='Use local database')    
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')
    parser.add_argument("--debug",
                        dest="debug",
                        type=str,
                        default="",
                        help="e.g., test_id:10,read_id:10000,basic_test")

    args = parser.parse_args()
    region_list = {}
    if args.region_list != "":
        for region in args.region_list.split(','):
            region = region.split('.')
            if len(region) < 1 or len(region) > 2:
                print >> sys.stderr, "Error: --region-list is incorrectly formatted."
                sys.exit(1)
                
            family = region[0].lower()
            if len(region) == 2:
                locus_name = region[1].upper()
            if family not in region_list:
                region_list[family] = set()
            if len(region) == 2:
                region_list[family].add(locus_name)

    read_fnames = []
    if args.alignment_fname != "":
        if not os.path.exists(args.alignment_fname):
            print >> sys.stderr, "Error: %s does not exist." % args.alignment_fname
    elif args.read_fname_U != "":
        read_fnames = [args.read_fname_U]
    else:
        if args.read_fname_1 == "" or args.read_fname_2 == "":
            print >> sys.stderr, "Error: please specify read file names correctly: -U or -1 and -2"
            sys.exit(1)
        read_fnames = [args.read_fname_1, args.read_fname_2]

    debug = {}
    if args.debug != "":
        for item in args.debug.split(','):
            if ':' in item:
                key, value = item.split(':')
                debug[key] = value
            else:
                debug[item] = 1

    genotype(args.base_fname,
             region_list,
             args.fastq,
             read_fnames,
             args.alignment_fname,
             args.threads,
             args.num_editdist,
             args.assembly,
             args.local_database,
             args.verbose,
             debug)


