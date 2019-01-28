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
from argparse import ArgumentParser, FileType
from hisatgenotype_typing_process import extract_reads
import hisatgenotype_args as arguments

"""
This is the Wrapper script that runs the processing code found in hisatgenotype_modules/hisatgenotype_typing_process
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract reads')

    # Add Arguments
    arguments.args_databases(parser, 
                             True) # Add option to change genotype_genome name
    arguments.args_input_output(parser)
    arguments.args_aligner_inputs(parser)
    arguments.args_set_aligner(parser,
                               False) # No missmatch option
    arguments.args_single_end(parser)
    arguments.args_extract_reads(parser)
    arguments.args_common(parser)

    args = parser.parse_args()

    if args.locus_list:
        print "--locus-list option not implemented in this script yet, skipping"
    if not args.graph_genome:
        print "--linear-index not implemented yet, skipping"

    if not args.genotype_genome:
        args.genotype_genome = 'genotype_genome'

    database_list = []
    if args.base_fname != "":
        for region in args.base_fname.split(','):
            database_list.append(region)
            
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

    if args.aligner not in ["hisat2", "bowtie2"]:
        print >> sys.stderr, "Error: --aligner should be either hisat2 or bowtie2."
        sys.exit(1)        
    block_size = 20000000 if args.extract_whole else 0
    
    if args.read_fname_U != "":
        paired = False
    else:
        paired = args.paired
        if (args.read_fname_1 != "" and args.read_fname_2 != "") and not paired:
            print >> sys.stderr, "Error: Don't set --single-end when using -1 and -2 options."
            exit(1)

    _ = extract_reads(args.genotype_genome,
                  database_list, # base_fname
                  args.read_dir,
                  args.out_dir,
                  args.suffix,
                  args.read_fname,
                  args.fastq,
                  paired,
                  args.simulation,
                  args.threads,
                  args.threads_aprocess,
                  args.max_sample,
                  job_range,
                  args.aligner,
                  block_size,
                  args.verbose)

