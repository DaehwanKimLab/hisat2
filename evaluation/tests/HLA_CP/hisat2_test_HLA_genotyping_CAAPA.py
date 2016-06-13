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


import sys, os, subprocess, re
import inspect
import random
import glob
from argparse import ArgumentParser, FileType

"""
"""
def test_HLA_genotyping(reference_type,
                        hla_list,
                        aligners,
                        exclude_allele_list,
                        num_mismatch,
                        sample_frac,
                        verbose):
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(test_HLA_genotyping))
    ex_path = os.path.dirname(curr_script)

    if not os.path.exists("CP"):
        print >> sys.stderr, "Error: CP data is needed."
        sys.exit(1)

    # fastq files
    fq_fnames = glob.glob("CP/*.1.gz")
    for fq_name in fq_fnames:
        fq_name = fq_name.split('/')[-1]
        genome = fq_name.split('.')[0]
        read_fname_1, read_fname_2 = "CP/%s.extracted.fq.1.gz" % genome, "CP/%s.extracted.fq.2.gz" % genome
        if not os.path.exists(read_fname_1) or not os.path.exists(read_fname_2):
            continue
        print >> sys.stderr, genome
        cmd_aligners = ['.'.join(aligners[i]) for i in range(len(aligners))]
        test_hla_script = os.path.join(ex_path, "hisat2_test_HLA_genotyping.py")
        test_hla_cmd = [test_hla_script,
                        "--reference-type", reference_type,
                        "--hla-list", ','.join(hla_list),
                        "--aligner-list", ','.join(cmd_aligners),
                        "--reads", "%s,%s" % (read_fname_1, read_fname_2),
                        # "--exclude-allele-list", ','.join(exclude_allele_list),
                        "--num-mismatch", str(num_mismatch)]
        
        if verbose:
            print >> sys.stderr, ' '.join(test_hla_cmd)
            
        proc = subprocess.Popen(test_hla_cmd, stdout=subprocess.PIPE, stderr=open("/dev/null", 'w'))
        test_alleles = set()
        for line in proc.stdout:
            print "\t\t", line,


"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='test HLA genotyping for Genomes')
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
    parser.add_argument("--aligner-list",
                        dest="aligners",
                        type=str,
                        default="hisat2.graph",
                        help="A comma-separated list of aligners (default: hisat2.graph)")
    parser.add_argument("--exclude-allele-list",
                        dest="exclude_allele_list",
                        type=str,
                        default="",
                        help="A comma-separated list of allleles to be excluded")
    parser.add_argument("--num-mismatch",
                        dest="num_mismatch",
                        type=int,
                        default=0,
                        help="Maximum number of mismatches per read alignment to be considered (default: 0)")
    parser.add_argument("--sample",
                        dest="sample_frac",
                        type=float,
                        default=1.1,
                        help="Reads sample rate (default: 1.1)")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

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
    args.exclude_allele_list = args.exclude_allele_list.split(',')

    test_HLA_genotyping(args.reference_type,
                        args.hla_list,
                        args.aligners,
                        args.exclude_allele_list,
                        args.num_mismatch,
                        args.sample_frac,
                        args.verbose)
