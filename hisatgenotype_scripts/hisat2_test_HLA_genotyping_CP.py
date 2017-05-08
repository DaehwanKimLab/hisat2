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

class myThread(threading.Thread):
    def __init__(self,
                 ex_path,
                 lock, 
                 paths,
                 reference_type,
                 hla_list,
                 partial,
                 aligners,
                 num_editdist,
                 max_sample,
                 assembly,
                 sample_frac,
                 out_dir,
                 verbose,
                 exclude_allele_list):
        threading.Thread.__init__(self)
        self.ex_path = ex_path
        self.lock = lock
        self.paths = paths
        self.reference_type = reference_type
        self.hla_list = hla_list
        self.partial = partial
        self.aligners = aligners
        self.num_editdist = num_editdist
        self.max_sample = max_sample
        self.assembly = assembly
        self.sample_frac = sample_frac
        self.out_dir = out_dir
        self.verbose = verbose
        self.exclude_allele_list = exclude_allele_list

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
            worker(self.ex_path,
                   self.lock,
                   self.paths[my_work_idx],
                   self.reference_type,
                   self.hla_list,
                   self.partial,
                   self.aligners,
                   self.num_editdist,
                   self.assembly,
                   self.sample_frac,
                   self.out_dir,
                   self.verbose,
                   self.exclude_allele_list)

work_idx = 0
def worker(ex_path,
           lock,
           path,
           reference_type,
           hla_list,
           partial,
           aligners,
           num_editdist,
           assembly,
           sample_frac,
           out_dir,
           verbose,
           exclude_allele_list):
    fq_name = path.split('/')[-1]
    read_dir = '/'.join(path.split('/')[:-1])
    genome = fq_name.split('.')[0]
    read_fname_1, read_fname_2 = "%s/%s.extracted.1.fq.gz" % (read_dir, genome), "%s/%s.extracted.2.fq.gz" % (read_dir, genome)

    if not os.path.exists(read_fname_1) or not os.path.exists(read_fname_2):
        return
    lock.acquire()
    print >> sys.stderr, genome
    lock.release()
    cmd_aligners = ['.'.join(aligners[i]) for i in range(len(aligners))]
    test_hla_cmd = ["hisatgenotype_locus.py",,
                    "--base", "hla",
                    "--locus", ','.join(hla_list),
                    "--aligner-list", ','.join(cmd_aligners),
                    "--reads", "%s,%s" % (read_fname_1, read_fname_2),
                    "--num-editdist", str(num_editdist)]
    test_hla_cmd += ["--assembly-base"]
    if out_dir != "":
        test_hla_cmd += ["%s/%s" % (out_dir, genome)]
    else:
        test_hla_cmd += [genome]

    if not assembly:
        test_hla_cmd += ["--no-assembly"]

    if len(exclude_allele_list) > 0:
        test_hla_cmd += ["--exclude-allele-list", ','.join(exclude_allele_list)]

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
def test_HLA_genotyping(read_dir,
                        reference_type,
                        hla_list,
                        partial,
                        aligners,
                        num_mismatch,
                        nthreads,
                        max_sample,
                        assembly,
                        sample_frac,
                        out_dir,
                        verbose,
                        exclude_allele_list):
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(test_HLA_genotyping))
    ex_path = os.path.dirname(curr_script)

    if not os.path.exists(read_dir):
        print >> sys.stderr, "Error: %s data is needed." % read_dir
        sys.exit(1)

    if out_dir != "" and not os.path.exists(out_dir):
        os.mkdir(out_dir)        

    # fastq files
    fq_fnames = glob.glob("%s/*.extracted.1.fq.gz" % read_dir)

    lock = threading.Lock()
    threads = []
    for t in range(nthreads):
        thread = myThread(ex_path,
                          lock,
                          fq_fnames,
                          reference_type,
                          hla_list,
                          partial,
                          aligners,
                          num_mismatch,
                          max_sample,
                          assembly,
                          sample_frac,
                          out_dir,
                          verbose,
                          exclude_allele_list)
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()


"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='test HLA genotyping for Genomes')
    parser.add_argument('read_dir',
                        nargs='?',
                        type=str,
                        help='read directory (e.g. CP)')
    parser.add_argument("--reference-type",
                        dest="reference_type",
                        type=str,
                        default="gene",
                        help="Reference type: gene, chromosome, and genome (default: gene)")
    parser.add_argument("--hla-list",
                        dest="hla_list",
                        type=str,
                        default="A,B,C,DQA1,DQB1,DRB1",
                        help="A comma-separated list of HLA genes (default: A,B,C,DQA1,DQB1,DRB1)")
    parser.add_argument('--no-partial',
                        dest='partial',
                        action='store_false',
                        help='Exclude partial alleles (e.g. A_nuc.fasta)')
    parser.add_argument("--aligner-list",
                        dest="aligners",
                        type=str,
                        default="hisat2.graph",
                        help="A comma-separated list of aligners (default: hisat2.graph)")
    parser.add_argument("--num-editdist",
                        dest="num_editdist",
                        type=int,
                        default=0,
                        help="Maximum number of mismatches per read alignment to be considered (default: 0)")
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
    parser.add_argument("--sample",
                        dest="sample_frac",
                        type=float,
                        default=1.1,
                        help="Reads sample rate (default: 1.1)")
    parser.add_argument("--out-dir",
                        dest="out_dir",
                        type=str,
                        default="",
                        help='Output directory (default: (empty))')
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')
    parser.add_argument("--exclude-allele-list",
                        dest="exclude_allele_list",
                        type=str,
                        default="",
                        help="A comma-separated list of allleles to be excluded")


    args = parser.parse_args()

    if not args.reference_type in ["gene", "chromosome", "genome"]:
        print >> sys.stderr, "Error: --reference-type (%s) must be one of gene, chromosome, and genome." % (args.reference_type)
        sys.exit(1)
    args.hla_list = args.hla_list.split(',')
    if args.aligners == "":
        print >> sys.stderr, "Error: --aligners must be non-empty."
        sys.exit(1)    
    args.aligners = args.aligners.split(',')
    for i in range(len(args.aligners)):
        args.aligners[i] = args.aligners[i].split('.')
    if args.exclude_allele_list != "":
        args.exclude_allele_list = args.exclude_allele_list.split(',')
    else:
        args.exclude_allele_list = []        

    test_HLA_genotyping(args.read_dir,
                        args.reference_type,
                        args.hla_list,
                        args.partial,
                        args.aligners,
                        args.num_editdist,
                        args.threads,
                        args.max_sample,
                        args.assembly,
                        args.sample_frac,
                        args.out_dir,
                        args.verbose,
                        args.exclude_allele_list)

