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


import os, sys, subprocess, re
import inspect
from argparse import ArgumentParser, FileType



"""
"""
def build_genotype_genome(reference,
                          base_fname,                          
                          gap,
                          split,
                          verbose):    
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(build_genotype_genome))
    ex_path = os.path.dirname(curr_script)

    # Clone a git repository, IMGTHLA
    if not os.path.exists("IMGTHLA"):
        os.system("git clone https://github.com/jrob119/IMGTHLA.git")

    def check_files(fnames):
        for fname in fnames:
            if not os.path.exists(fname):
                return False
        return True

    # Download HISAT2 index
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]
    if not check_files(HISAT2_fnames):
        os.system("wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz; tar xvzf grch38.tar.gz; rm grch38.tar.gz")
        hisat2_inspect = os.path.join(ex_path, "hisat2-inspect")
        os.system("%s grch38/genome > genome.fa" % hisat2_inspect)
        os.system("samtools faidx genome.fa")

    # Check if the pre-existing files (hla*) are compatible with the current parameter setting
    if os.path.exists("hla.ref"):
        left = 0
        HLA_genes = set()
        for line in open("hla.ref"):
            HLA_name, chr, left, _, exon_str = line.strip().split()
            HLA_gene = HLA_name.split('*')[0]
            HLA_genes.add(HLA_gene)
            left = int(left)
        delete_hla_files = False
        if reference_type == "gene":
            if left > 0:
                delete_hla_files = True
        elif reference_type == "chromosome":
            if left == 0:
                delete_hla_files = True
        else:
            assert reference_type == "genome"
        if not set(hla_list).issubset(HLA_genes):
            delete_hla_files = True
        if delete_hla_files:
            os.system("rm hla*")
    
    # Extract HLA variants, backbone sequence, and other sequeces
    HLA_fnames = ["hla_backbone.fa",
                  "hla_sequences.fa",
                  "hla.ref",
                  "hla.snp",
                  "hla.haplotype",
                  "hla.link"]

    if not check_files(HLA_fnames):
        extract_hla_script = os.path.join(ex_path, "hisat2_extract_HLA_vars.py")
        extract_cmd = [extract_hla_script,
                       "--reference-type", "gene",
                       "--hla-list", "A,B,C,DRB1,DQB1,DRB1"]
        # if partial:
        #    extract_cmd += ["--partial"]
        extract_cmd += ["--gap", "30",
                        "--split", "50"]
        if verbose:
            print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
        proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()
        if not check_files(HLA_fnames):
            print >> sys.stderr, "Error: extract_HLA_vars failed!"
            sys.exit(1)

    # Write genotype genome
            
    # (1) genotype_genome.fa

    # (2) genotype_genome.snp

    # (3) genotype_genome.haplotype

    # (4) genotype_genome.link

        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Extract HLA variants from HLA multiple sequence alignments")
    parser.add_argument("reference",
                        nargs='?',
                        type=str,
                        help="Reference genome")
    parser.add_argument("base_fname",
                        nargs='?',
                        type=str,
                        help="base filename for genotype genome")
    parser.add_argument("-g", "--gap",
                        dest="gap",
                        type=int,
                        default=30,
                        help="Maximum distance for variants to be in the same haplotype")
    parser.add_argument("-s", "--split",
                        dest="split",
                        type=int,
                        default=50,
                        help="Break a haplotype into several haplotypes")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="also print some statistics to stderr")

    args = parser.parse_args()
    if not args.reference or not args.base_fname:
        parser.print_help()
        sys.exit(1)
    if args.gap > args.split:
        print >> sys.stderr, "Error: -g/--gap (%d) must be smaller than -s/--split (%d)" % (args.gap, args.split)
        sys.exit(1)
    build_genotype_genome(args.reference,
                          args.base_fname,
                          args.gap,
                          args.split,
                          args.verbose)
