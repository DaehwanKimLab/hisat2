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

import os, sys, subprocess, re
import multiprocessing
from argparse import ArgumentParser, FileType
from hisatgenotype_typing_process import extract_vars
import hisatgenotype_typing_common as typing_common
import hisatgenotype_args as arguments
        
"""
This is the Wrapper script that runs the processing code found in hisatgenotype_modules/hisatgenotype_typing_process
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Extract variants from multiple sequence alignments")

    # Add Arguments
    arguments.args_databases(parser)
    arguments.args_var_gaps(parser)
    arguments.args_extract_vars(parser)
    arguments.args_no_partial(parser)
    arguments.args_common(parser)

    args = parser.parse_args()
    if args.locus_list == "":
        locus_list = []
    else:
        locus_list = args.locus_list.split(',')
    if args.inter_gap > args.intra_gap:
        print >> sys.stderr, "Error: --inter-gap (%d) must be smaller than --intra-gap (%d)" % (args.inter_gap, args.intra_gap)
        sys.exit(1)
    
    # Clone hisatgenotype database from git
    if not os.path.exists("hisatgenotype_db"):
        typing_common.clone_hisatgenotype_database()
    
    if not args.base_fname:
        base_fname = []
        for database in os.listdir("hisatgenotype_db"):
            if database in ['.git', 'README.md']:
                continue
            base_fname.append(database.lower())
    else:
        base_fname = args.base_fname.split(',')

    def init(l):
        global lock
        lock = l
    l = multiprocessing.Lock()

    pool = multiprocessing.Pool(int(args.threads), initializer=init, initargs=(l,))
    for base in base_fname:
        if base.find('/') != -1:
            elems = base.split('/')
            base = elems[-1]
            base_dname = '/'.join(elems[:-1])
        else:
            base_dname = ""        

        pool.apply_async(extract_vars, 
                         args=(base,
                               base_dname,
                               locus_list,
                               args.inter_gap,
                               args.intra_gap,
                               args.whole_haplotype,
                               args.min_var_freq,
                               args.ext_seq_len,
                               args.leftshift,
                               args.partial,
                               args.verbose))
    pool.close()
    pool.join()

