#!/usr/bin/env python

import sys, os
use_message = '''
'''

def get_data():
    data_root = "http://www.ccb.jhu.edu/software/hisat2/downloads/evaluation"
    
    # Download the reference human genome, SNPs, and gene annotations
    if not os.path.exists("data"):
        os.mkdir("data")
    os.chdir("data")
    files = ["genome.fa", "genome.fa.fai", "genes.gtf", "snpCommon.txt", "genome.snp", "genome.ss", \
                 "22.fa", "22.fa.fai", "genes_22.gtf", "22.snp", "22.ss", \
                 "22_20-21M.fa", "22_20-21M.fa.fai", "genes_22_20-21M.gtf", "22_20-21M.snp", "22_20-21M.ss"]
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
        for aligner in aligners:
            if genome == "genome":
                aligner_dir = aligner
            else:
                aligner_dir = aligner + "_" + genome
            if os.path.exists(aligner_dir):
                continue
            cmd = "wget %s/indexes/%s/%s.tar.gz; tar xvzf %s.tar.gz" % (data_root, aligner_dir, aligner_dir)
            print >> sys.stderr, cmd
            os.system(cmd)
    os.chdir("..")
            
    
if __name__ == "__main__":
    get_data()
