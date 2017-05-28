#!/usr/bin/env python

#
# Copyright 2015, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
#


import sys, os, subprocess, re, threading
import inspect
import random
import glob
from argparse import ArgumentParser, FileType
import hisatgenotype_typing_common as typing_common


"""
"""
class myThread(threading.Thread):
    def __init__(self,
                 lock, 
                 paths,
                 reference_type,
                 region_list,
                 num_editdist,
                 max_sample,
                 assembly,
                 out_dir,
                 verbose):
        threading.Thread.__init__(self)
        self.lock = lock
        self.paths = paths
        self.reference_type = reference_type
        self.region_list = region_list
        self.num_editdist = num_editdist
        self.max_sample = max_sample
        self.assembly = assembly
        self.out_dir = out_dir
        self.verbose = verbose

    def run(self):
        global work_idx
        while True:
            self.lock.acquire()
            my_work_idx = work_idx
            work_idx += 1
            self.lock.release()
            if my_work_idx >= len(self.paths) or \
               my_work_idx >= self.max_sample:
                return
            worker(self.lock,
                   self.paths[my_work_idx],
                   self.reference_type,
                   self.region_list,
                   self.num_editdist,
                   self.assembly,
                   self.out_dir,
                   self.verbose)

            
"""
"""
work_idx = 0
def worker(lock,
           path,
           reference_type,
           region_list,
           num_editdist,
           assembly,
           out_dir,
           verbose):
    fq_name = path.split('/')[-1]
    read_dir = '/'.join(path.split('/')[:-1])
    genome = fq_name.split('.')[0]
    if not fq_name.endswith("extracted.1.fq.gz"):
        return
    read_basename = fq_name[:fq_name.find("extracted.1.fq.gz")]
    read_fname_1, read_fname_2 = "%s/%sextracted.1.fq.gz" % \
                                 (read_dir, read_basename), "%s/%sextracted.2.fq.gz" % (read_dir, read_basename)

    if not os.path.exists(read_fname_1) or not os.path.exists(read_fname_2):
        return
    lock.acquire()
    print >> sys.stderr, genome
    lock.release()

    for family, loci in region_list.items():
        test_hla_cmd = ["hisatgenotype_locus.py",
                        "--base", family]
        if len(loci) > 0:
            test_hla_cmd += ["--locus", ','.join(loci)]
        test_hla_cmd += ["--num-editdist", str(num_editdist)]
        test_hla_cmd += ["-1", read_fname_1, "-2", read_fname_2]
        test_hla_cmd += ["--assembly-base"]
        if out_dir != "":
            test_hla_cmd += ["%s/%s" % (out_dir, genome)]
        else:
            test_hla_cmd += [genome]
        if assembly:
            test_hla_cmd += ["--assembly"]

        if verbose:
            lock.acquire()
            print >> sys.stderr, ' '.join(test_hla_cmd)
            lock.release()

        proc = subprocess.Popen(test_hla_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        test_alleles = set()
        output_list = []
        for line in proc.stdout:
            line = line.strip()
            if line.find("abundance") == -1:
                continue

            rank, _, allele, _, abundance = line.split()        
            output_list.append([allele, abundance[:-2]])

    lock.acquire()
    for output in output_list:
        allele, abundance = output
        print >> sys.stdout, "%s\t%s\t%s" % (genome, allele, abundance)
    sys.stdout.flush()
    lock.release()


"""
"""
def genotyping(read_dir,
               reference_type,
               region_list,
               num_editdist,
               nthreads,
               max_sample,
               assembly,
               out_dir,
               verbose):
    for database_name in region_list:
        # Extract variants, backbone sequence, and other sequeces
        typing_common.extract_database_if_not_exists(database_name,
                                                     [])            # locus_list
        # Build HISAT2's graph index
        typing_common.build_index_if_not_exists(database_name,
                                                "hisat2",
                                                "graph",
                                                1,            # threads
                                                verbose)
    
    if not os.path.exists(read_dir):
        print >> sys.stderr, "Error: %s does not exist." % read_dir
        sys.exit(1)

    if out_dir != "" and not os.path.exists(out_dir):
        os.mkdir(out_dir)        

    # fastq files
    fq_fnames = glob.glob("%s/*.extracted.1.fq.gz" % read_dir)

    lock = threading.Lock()
    threads = []
    for t in range(nthreads):
        thread = myThread(lock,
                          fq_fnames,
                          reference_type,
                          region_list,
                          num_editdist,
                          max_sample,
                          assembly,
                          out_dir,
                          verbose)
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()


"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='genotyping on many samples')
    parser.add_argument("--reference-type",
                        dest="reference_type",
                        type=str,
                        default="gene",
                        help="Reference type: gene, chromosome, and genome (default: gene)")
    parser.add_argument("--region-list",
                        dest="region_list",
                        type=str,
                        default="",
                        help="A comma-separated list of regions (default: empty)")
    parser.add_argument('--read-dir',
                        dest="read_dir",
                        type=str,
                        default="",
                        help='read directory (e.g. read_input)')
    parser.add_argument("--num-editdist",
                        dest="num_editdist",
                        type=int,
                        default=2,
                        help="Maximum number of mismatches per read alignment to be considered (default: 2)")
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    parser.add_argument('--assembly',
                        dest='assembly',
                        action='store_true',
                        help='Perform assembly')
    parser.add_argument("--max-sample",
                        dest="max_sample",
                        type=int,
                        default=sys.maxint,
                        help="Number of samples to be analyzed (default: sys.maxint)")
    parser.add_argument("--out-dir",
                        dest="out_dir",
                        type=str,
                        default="",
                        help='Output directory (default: (empty))')
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()

    if args.read_dir == "":
        print >> sys.stderr, "Error: please specify --read-dir."
        sys.exit(1)

    if not args.reference_type in ["gene", "chromosome", "genome"]:
        print >> sys.stderr, "Error: --reference-type (%s) must be one of gene, chromosome, and genome." % (args.reference_type)
        sys.exit(1)

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

    genotyping(args.read_dir,
               args.reference_type,
               region_list,
               args.num_editdist,
               args.threads,
               args.max_sample,
               args.assembly,
               args.out_dir,
               args.verbose)

