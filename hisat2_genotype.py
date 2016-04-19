#!/usr/bin/env python

#
# Copyright 2016, Daehwan Kim <infphilo@gmail.com>
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
import inspect, random
import math
from argparse import ArgumentParser, FileType


"""
"""
def genotype(reference_type,
             hla_list,
             partial,
             aligners,
             read_fname,
             alignment_fname,
             threads,
             simulate_interval,
             enable_coverage,
             best_alleles,
             exclude_allele_list,
             num_mismatch,
             verbose,
             daehwan_debug):
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(genotype))
    ex_path = os.path.dirname(curr_script)

        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='HISAT genotyping')
    parser.add_argument("--reference-type",
                        dest="reference_type",
                        type=str,
                        default="gene",
                        help="Reference type: gene, chromosome, and genome (default: gene)")
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    parser.add_argument("--simulate-interval",
                        dest="simulate_interval",
                        type=int,
                        default=1,
                        help="Reads simulated at every these base pairs (default: 1)")
    parser.add_argument("--num-mismatch",
                        dest="num_mismatch",
                        type=int,
                        default=0,
                        help="Maximum number of mismatches per read alignment to be considered (default: 0)")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')
    parser.add_argument("--daehwan-debug",
                        dest="daehwan_debug",
                        type=str,
                        default="",
                        help="e.g., test_id:10,read_id:10000,basic_test")

    args = parser.parse_args()
    if not args.reference_type in ["gene", "chromosome", "genome"]:
        print >> sys.stderr, "Error: --reference-type (%s) must be one of gene, chromosome, and genome." % (args.reference_type)
        sys.exit(1)
    daehwan_debug = {}
    if args.daehwan_debug != "":
        for item in args.daehwan_debug.split(','):
            if ':' in item:
                key, value = item.split(':')
                daehwan_debug[key] = value
            else:
                daehwan_debug[item] = 1

    random.seed(1)
    genotype(args.reference_type,
             args.threads,
             args.simulate_interval,
             args.num_mismatch,
             args.verbose,
             daehwan_debug)
