#!/usr/bin/env python

import sys, os
from argparse import ArgumentParser, FileType

def get_data(small = False):
    data_root = "http://www.ccb.jhu.edu/software/hisat2/downloads/evaluation"
    
    # Download the reference human genome, SNPs, and gene annotations
    if not os.path.exists("data"):
        os.mkdir("data")
    os.chdir("data")
    genome_files = ["genome.fa", "genome.fa.fai", "genes.gtf", "snpCommon.txt", "genome.snp", "genome.ss"]
    small_genome_files = ["22.fa", "22.fa.fai", "genes_22.gtf", "22.snp", "22.ss", \
                              "22_20-21M.fa", "22_20-21M.fa.fai", "genes_22_20-21M.gtf", "22_20-21M.snp", "22_20-21M.ss"]
    files = []
    if not small:
        files += genome_files
    files += small_genome_files
    for file in files:
        if os.path.exists(file):
            continue
        wget_cmd = "wget %s/data/%s" % (data_root, file)
        print >> sys.stderr, wget_cmd
        os.system(wget_cmd)
    os.chdir("..")

    # Download indexes
    if not os.path.exists("indexes"):
        os.mkdir("indexes")
    os.chdir("indexes")
    aligners = ["HISAT2", "HISAT", "Bowtie", "STAR", "GSNAP"]
    for genome in ["genome", "22", "22_20-21M"]:
        if small and genome == "genome":
            continue
        for aligner in aligners:
            if genome == "genome":
                aligner_dir = aligner
            else:
                aligner_dir = aligner + "_" + genome
            if os.path.exists(aligner_dir):
                continue
            cmd = "wget %s/indexes/%s.tar.gz; tar xvzf %s.tar.gz; rm %s.tar.gz" % \
                (data_root, aligner_dir, aligner_dir, aligner_dir)
            print >> sys.stderr, cmd
            os.system(cmd)
    os.chdir("..")
            
    
if __name__ == "__main__":
    parser = ArgumentParser(
        description='Get reference genome, annotations, and indexes')
    parser.add_argument('-s', '--small',
                        dest='small',
                        action='store_true',
                        default=False,
                        help='small testset')
    args = parser.parse_args()
    get_data(args.small)
