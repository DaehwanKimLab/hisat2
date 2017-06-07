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


import sys, os, subprocess, re
import inspect
import random
import glob
from argparse import ArgumentParser, FileType
import hisatgenotype_typing_common as typing_common


"""
"""
def parallel_work(pids, 
                  work, 
                  fq_fname_base, 
                  fq_fname, 
                  fq_fname2, 
                  ranges,
                  simulation,
                  verbose):
    child = -1
    for i in range(len(pids)):
        if pids[i] == 0:
            child = i
            break

    while child == -1:
        status = os.waitpid(0, 0)
        for i in range(len(pids)):
            if status[0] == pids[i]:
                child = i
                pids[i] = 0
                break

    child_id = os.fork()
    if child_id == 0:
        work(fq_fname_base, 
             fq_fname, 
             fq_fname2, 
             ranges,
             simulation,
             verbose)
        os._exit(os.EX_OK)
    else:
        # print >> sys.stderr, '\t\t>> thread %d: %d' % (child, child_id)
        pids[child] = child_id

        
"""
"""
def wait_pids(pids):
    for pid in pids:
        if pid > 0:
            os.waitpid(pid, 0)
            

"""
"""
def extract_reads(base_fname,
                  database_list,
                  read_dir,
                  out_dir,
                  suffix,
                  read_fname,
                  fastq,
                  paired,
                  simulation,
                  threads,
                  max_sample,
                  job_range,
                  verbose):
    genotype_fnames = ["%s.fa" % base_fname,
                       "%s.locus" % base_fname,
                       "%s.snp" % base_fname,
                       "%s.haplotype" % base_fname,
                       "%s.link" % base_fname,
                       "%s.coord" % base_fname,
                       "%s.clnsig" % base_fname]
    # hisat2 graph index files
    genotype_fnames += ["%s.%d.ht2" % (base_fname, i+1) for i in range(8)]
    if not typing_common.check_files(genotype_fnames):        
        print >> sys.stderr, "Error: %s related files do not exist as follows:" % base_fname
        for fname in genotype_fnames:
            print >> sys.stderr, "\t%s" % fname
        sys.exit(1)

    filter_region = len(database_list) > 0
    ranges = []
    regions, region_loci = {}, {}
    for line in open("%s.locus" % base_fname):
        family, allele_name, chr, left, right = line.strip().split()[:5]
        if filter_region and family.lower() not in database_list:
            continue
        region_name = "%s-%s" % (family, allele_name.split('*')[0])
        assert region_name not in regions
        regions[region_name] = allele_name
        left, right = int(left), int(right)
        """
        exons = []
        for exon in exon_str.split(','):
            exon_left, exon_right = exon.split('-')
            exons.append([int(exon_left), int(exon_right)])
        """
        if chr not in region_loci:
            region_loci[chr] = {}
        region_loci[chr][region_name] = [allele_name, chr, left, right]
        database_list.add(family.lower())

    if out_dir != "" and not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Extract reads
    if len(read_fname) > 0:
        if paired:
            fq_fnames = [read_fname[0]]
            fq_fnames2 = [read_fname[1]]
        else:
            fq_fnames = read_fname
    else:
        if paired:
            fq_fnames = glob.glob("%s/*.1.%s" % (read_dir, suffix))
        else:
            fq_fnames = glob.glob("%s/*.%s" % (read_dir, suffix))
    count = 0
    pids = [0 for i in range(threads)]
    for file_i in range(len(fq_fnames)):
        if file_i >= max_sample:
            break
        fq_fname = fq_fnames[file_i]
        if job_range[1] > 1:
            if job_range[0] != (file_i % job_range[1]):
                continue

        fq_fname_base = fq_fname.split('/')[-1]
        one_suffix = ".1." + suffix
        if fq_fname_base.find(one_suffix) != -1:
            fq_fname_base = fq_fname_base[:fq_fname_base.find(one_suffix)]
        else:
            fq_fname_base = fq_fname_base.split('.')[0]
            
        if paired:
            if read_dir == "":
                fq_fname2 = fq_fnames2[file_i]
            else:
                fq_fname2 = "%s/%s.2.%s" % (read_dir, fq_fname_base, suffix)
            if not os.path.exists(fq_fname2):
                print >> sys.stderr, "%s does not exist." % fq_fname2
                continue
        else:
            fq_fname2 = ""

        if paired:
            if out_dir != "":
                if os.path.exists("%s/%s.extracted.1.fq.gz" % (out_dir, fq_fname_base)):
                    continue
        else:
            if out_dir != "":
                if os.path.exists("%s/%s.extracted.fq.gz" % (out_dir, fq_fname_base)):
                    continue
        count += 1

        print >> sys.stderr, "\t%d: Extracting reads from %s" % (count, fq_fname_base)
        def work(fq_fname_base,
                 fq_fname, 
                 fq_fname2, 
                 ranges,
                 simulation,
                 verbose):
            aligner_cmd = ["hisat2"]
            if not fastq:
                aligner_cmd += ["-f"]
            aligner_cmd += ["-x", base_fname]
            aligner_cmd += ["--no-spliced-alignment",
                            "--max-altstried", "64"]
            if paired:
                aligner_cmd += ["-1", fq_fname,
                                "-2", fq_fname2]
            else:
                aligner_cmd += ["-U", fq_fname]
            if verbose:
                print >> sys.stderr, "\t\trunning", ' '.join(aligner_cmd)
            align_proc = subprocess.Popen(aligner_cmd,
                                          stdout=subprocess.PIPE,
                                          stderr=open("/dev/null", 'w'))

            gzip_dic = {}
            out_dir_slash = out_dir
            if out_dir != "":
                out_dir_slash += "/"
            for database in database_list:
                if paired:
                    # LP6005041-DNA_A01.extracted.1.fq.gz
                    gzip1_proc = subprocess.Popen(["gzip"],
                                                  stdin=subprocess.PIPE,
                                                  stdout=open("%s%s.%s.extracted.1.fq.gz" % (out_dir_slash, fq_fname_base, database), 'w'),
                                                  stderr=open("/dev/null", 'w'))

                    # LP6005041-DNA_A01.extracted.2.fq.gz
                    gzip2_proc = subprocess.Popen(["gzip"],
                                                  stdin=subprocess.PIPE,
                                                  stdout=open("%s%s.%s.extracted.2.fq.gz" % (out_dir_slash, fq_fname_base, database), 'w'),
                                                  stderr=open("/dev/null", 'w'))
                else:
                    # LP6005041-DNA_A01.extracted.fq.gz
                    gzip1_proc = subprocess.Popen(["gzip"],
                                                  stdin=subprocess.PIPE,
                                                  stdout=open("%s%s.%s.extracted.fq.gz" % (out_dir_slash, fq_fname_base, database), 'w'),
                                                  stderr=open("/dev/null", 'w'))
                gzip_dic[database] = [gzip1_proc, gzip2_proc if paired else None]

            def write_read(gzip_proc, read_name, seq, qual):
                if fastq:
                    gzip_proc.stdin.write("@%s\n" % read_name)
                    gzip_proc.stdin.write("%s\n" % seq)
                    gzip_proc.stdin.write("+\n")
                    gzip_proc.stdin.write("%s\n" % qual)
                else:
                    gzip_proc.stdin.write(">%s\n" % prev_read_name)
                    gzip_proc.stdin.write("%s\n" % seq)                    

            prev_read_name, extract_read, read1, read2 = "", False, [], []
            for line in align_proc.stdout:
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

                if (not simulation and read_name != prev_read_name) or \
                   (simulation and read_name.split('|')[0] != prev_read_name.split('|')[0]):
                    if extract_read:
                        write_read(gzip_dic[region][0], prev_read_name, read1[0], read1[1])
                        if paired:
                            write_read(gzip_dic[region][1], prev_read_name, read2[0], read2[1])
                    prev_read_name, extract_read, read1, read2 = read_name, False, [], []

                if flag & 0x4 == 0 and NH == 1 and chr in region_loci:                    
                    for region, loci in region_loci[chr].items():
                        region = region.split('-')[0].lower()
                        _, _, loci_left, loci_right = loci
                        if pos >= loci_left and pos < loci_right:
                            extract_read = True
                            break

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
                write_read(gzip_dic[region][0], prev_read_name, read1[0], read1[1])
                if paired:
                    write_read(gzip_dic[region][1], prev_read_name, read2[0], read2[1])

            for gzip1_proc, gzip2_proc in gzip_dic.values():
                gzip1_proc.stdin.close()
                if paired:
                    gzip2_proc.stdin.close()                        

        if threads <= 1:
            work(fq_fname_base, 
                 fq_fname, 
                 fq_fname2,
                 ranges,
                 simulation,
                 verbose)
        else:
            parallel_work(pids, 
                          work, 
                          fq_fname_base, 
                          fq_fname, 
                          fq_fname2, 
                          ranges,
                          simulation,
                          verbose)

    if threads > 1:
        wait_pids(pids)


"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract reads')
    parser.add_argument("--base", "--base-fname",
                        dest="base_fname",
                        type=str,
                        default="genotype_genome",
                        help="base filename for genotype genome")
    parser.add_argument("--read-dir",
                        dest="read_dir",
                        type=str,
                        default="",
                        help="Directory name for read files")
    parser.add_argument("--out-dir",
                        dest="out_dir",
                        type=str,
                        default="",
                        help="Directory name for extracted read files")
    parser.add_argument("--suffix",
                        dest="suffix",
                        type=str,
                        default="fq.gz",
                        help="Read file suffix (Default: fq.gz)")
    parser.add_argument('-f', '--fasta',
                        dest='fastq',
                        action='store_false',
                        help='FASTA format')
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
    parser.add_argument("--database-list",
                        dest="database_list",
                        type=str,
                        default="",
                        help="A comma-separated list of database (default: empty)")
    parser.add_argument('--simulation',
                        dest='simulation',
                        action='store_true',
                        help='Simulated reads (Default: False)')    
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    parser.add_argument("--max-sample",
                        dest="max_sample",
                        type=int,
                        default=sys.maxint,
                        help="Number of samples to be extracted (default: sys.maxint)")
    parser.add_argument("--job-range",
                        dest="job_range",
                        type=str,
                        default="0,1",
                        help="two numbers (e.g. 1,3)")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()

    database_list = set()
    if args.database_list != "":
        for region in args.database_list.split(','):
            database_list.add(region)
    if args.read_fname_U != "":
        args.read_fname = [args.read_fname_U]
    elif args.read_fname_1 != "" or args.read_fname_2 != "":
        if args.read_fname_1 == "" or args.read_fname_2 == "":
            print >> sys.stderr, "Error: please specify both -1 and -2."
            sys.exit(1)
        args.read_fname = [args.read_fname_1, args.read_fname_2]
    else:
        args.read_fname = []
    if len(args.read_fname) == 0:
        if args.read_dir == "" or not os.path.exists(args.read_dir):
            print >> sys.stderr, "Error: please specify --read-dir with an existing directory."
            sys.exit(1)
        if args.out_dir == "":
            print >> sys.stderr, "Error: please specify --out-dir with a directory name."
            sys.exit(1)
    job_range = []
    for num in args.job_range.split(','):
        job_range.append(int(num))
        
    extract_reads(args.base_fname,
                  database_list,
                  args.read_dir,
                  args.out_dir,
                  args.suffix,
                  args.read_fname,
                  args.fastq,
                  False if args.read_fname_U != "" else True,
                  args.simulation,
                  args.threads,
                  args.max_sample,
                  job_range,
                  args.verbose)

