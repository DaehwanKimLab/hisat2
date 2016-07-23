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
import inspect
import random
import glob
from argparse import ArgumentParser, FileType


"""
"""
def parallel_work(pids, work, cmd):
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
        work(cmd)
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
def test_HLA_genotyping(reference_type,
                        hla_list,
                        partial,
                        threads,
                        verbose):
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(test_HLA_genotyping))
    ex_path = "../../.."

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
                       "--reference-type", reference_type,
                       "--hla-list", ','.join(hla_list),
                       "--DRB1-REF"]
        if partial:
            extract_cmd += ["--partial"]
        extract_cmd += ["--gap", "30",
                        "--split", "50"]
        print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
        proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()
        if not check_files(HLA_fnames):
            print >> sys.stderr, "Error: extract_HLA_vars failed!"
            sys.exit(1)

    # Build HISAT2 graph indexes based on the above information
    HLA_hisat2_graph_index_fnames = ["hla.graph.%d.ht2" % (i+1) for i in range(8)]
    if not check_files(HLA_hisat2_graph_index_fnames):
        hisat2_build = os.path.join(ex_path, "hisat2-build")
        build_cmd = [hisat2_build,
                     # "-p", str(threads),
                     "--snp", "hla.snp",
                     "--haplotype", "hla.haplotype",
                     "hla_backbone.fa",
                     "hla.graph"]
        print >> sys.stderr, "\tRunning:", ' '.join(build_cmd)
        proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()        
        if not check_files(HLA_hisat2_graph_index_fnames):
            print >> sys.stderr, "Error: indexing HLA failed!  Perhaps, you may have forgotten to build hisat2 executables?"
            sys.exit(1)

    # Extract reads
    fq_fnames = glob.glob("*_1.fq.gz")
    count = 0
    pids = [0 for i in range(threads)]
    for fq_fname in fq_fnames:
        fq_fname_base = fq_fname.split('.')[0][:-2]
        fq_fname2 = "%s_2.fq.gz" % fq_fname_base
        if not os.path.exists(fq_fname2):
            print >> sys.stderr, "%s does not exist." % fq_fname2
            continue
        if os.path.exists("%s.extracted.fq.1.gz" % fq_fname_base):
            continue        
        count += 1

        print >> sys.stderr, "\t%d: Extracting reads from %s" % (count, fq_fname_base)
        hisat2 = os.path.join(ex_path, "hisat2")
        aligner_cmd = [hisat2,
                       "--al-conc-disc-gz", "%s.extracted.fq.gz" % fq_fname_base,
                       "--mm",
                       "-x", "hla.graph",
                       "-1", fq_fname,
                       "-2", fq_fname2]
        print >> sys.stderr, "\t\trunning", ' '.join(aligner_cmd)            
            
        def work(aligner_cmd):
            align_proc = subprocess.Popen(aligner_cmd,
                                          stdout=open("/dev/null", 'w'),
                                          stderr=open("/dev/null", 'w'))
            align_proc.communicate()

        if threads <= 1:
            work(aligner_cmd)
        else:
            parallel_work(pids, work, aligner_cmd)

    if threads > 1:
        wait_pids(pids)


"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='test HLA genotyping for Platinum Genomes')
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
    parser.add_argument('--partial',
                        dest='partial',
                        action='store_true',
                        help='Include partial alleles (e.g. A_nuc.fasta)')
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()

    if not args.reference_type in ["gene", "chromosome", "genome"]:
        print >> sys.stderr, "Error: --reference-type (%s) must be one of gene, chromosome, and genome." % (args.reference_type)
        sys.exit(1)
    args.hla_list = args.hla_list.split(',')
    test_HLA_genotyping(args.reference_type,
                        args.hla_list,
                        args.partial,
                        args.threads,
                        args.verbose)
