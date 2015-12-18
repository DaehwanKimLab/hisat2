#!/usr/bin/env python

import sys, os
from argparse import ArgumentParser, FileType

def get_data(small = False):
    data_root = "http://www.ccb.jhu.edu/software/hisat2/downloads/evaluation"
    
    # Download the reference human genome, SNPs, and gene annotations
    if not os.path.exists("data"):
        os.mkdir("data")
    os.chdir("data")
    genome_files = ["genome.fa", "genome.fa.fai", "genome.gtf", "snpCommon.txt", "genome.snp", "genome.ss", "genome.exon"]
    small_genome_files = ["22.fa", "22.fa.fai", "22.gtf", "22.snp", "22.ss", "22.exon", \
                              "22_20-21M.fa", "22_20-21M.fa.fai", "22_20-21M.gtf", "22_20-21M.snp", "22_20-21M.ss", "22_20-21M.exon"]
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

    # Download simulated and real reads
    if not os.path.exists("reads"):
        os.mkdir("reads")
    os.chdir("reads")
    for type in ["simulation", "real"]:
        if small and type == "real":
            continue
        if not os.path.exists(type):
            os.mkdir(type)
        os.chdir(type)
        if type == "simulation":
            files = ["1M_DNA_reads_22",
                     "1M_DNA_mismatch_reads_22",
                     "1M_DNA_snp_reads_22",
                     "1M_DNA_mismatch_snp_reads_22",
                     "1M_RNA_reads_22",
                     "1M_RNA_constant_reads_22",
                     "1M_RNA_mismatch_reads_22",
                     "1M_RNA_snp_reads_22",
                     "1M_RNA_mismatch_snp_reads_22",
                     "1M_RNA_reads_22_20-21M",
                     "20M_DNA_reads_genome",
                     "20M_DNA_snp_reads_genome",
                     "20M_RNA_reads_genome",
                     "20M_RNA_snp_reads_genome"]
        else:
            files = ["108M_RNA_wgEncodeCshlLongRnaSeq",
                     "62M_RNA_SRR353653",
                     "80M_DNA_SRR345300",
                     "5M_DNA_NA12878D"]
        for file in files:
            if small and file.find("20M") != -1:
                continue
            if os.path.exists(file):
                continue
            cmd = "wget %s/reads/%s/%s.tar.gz; tar xvzf %s.tar.gz; rm %s.tar.gz" % \
                (data_root, type, file, file, file)
            print >> sys.stderr, cmd
            os.system(cmd)
        os.chdir("..")
    
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
