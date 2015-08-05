#!/usr/bin/env python

import sys, os
use_message = '''
'''

# Assume Mac Platform
def get_data():
    # Build indexes
    if not os.path.exists("indexes"):
        os.mkdir("indexes")
    os.chdir("indexes")
    aligners = ["HISAT", "Bowtie", "STAR", "GSNAP"]
    for aligner in aligners:
        if os.path.exists(aligner):
            continue
        os.mkdir(aligner)
        os.chdir(aligner)
        if aligner == "HISAT2":
            cmd = "../../aligners/bin/hisat2-build ../../data/genome.fa genome"
            cmd = cmd + "; ../../aligners/bin/hisat2-build -p 4 ../../data/genome.fa --snp ../../data/genome.snp genome_snp"
            cmd = cmd + "; ../../aligners/bin/hisat2-build -p 4 ../../data/genome.fa --ss ../../data/genome.ss genome_ss"
            cmd = cmd + "; ../../aligners/bin/hisat2-build -p 4 ../../data/genome.fa --snp ../../data/genome.snp --ss ../../data/genome.ss genome_snp_ss"
        elif aligner == "HISAT":
            cmd = "../../aligners/bin/hisat-build ../../data/genome.fa genome"
        elif aligner == "Bowtie":
            cmd = "../../aligners/bin/bowtie-build ../../data/genome.fa genome"
        elif aligner == "STAR":
            cmd = "../../aligners/bin/STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles ../../data/genome.fa"
        elif aligner == "GSNAP":
            cmd = "../../aligners/bin/gmap_build.pl -D . -d genome ../../data/genome.fa"
        else:
            assert False
        print >> sys.stderr, cmd
        os.system(cmd)
        os.chdir("..")

    os.chdir("..")
            
    
if __name__ == "__main__":
    get_data()
