#!/usr/bin/env python

import sys, os
use_message = '''
'''

def build_indexes():
    # Build indexes
    if not os.path.exists("indexes"):
        os.mkdir("indexes")
    os.chdir("indexes")
    aligners = ["HISAT2", "HISAT", "Bowtie", "STAR", "GSNAP", "BWA"]
    for genome in ["genome", "22", "22_20-21M"]:
        for aligner in aligners:
            if genome == "genome":
                dir = aligner
            else:
                dir = aligner + "_" + genome
            if os.path.exists(dir):
                continue
            os.mkdir(dir)
            os.chdir(dir)
            if aligner == "HISAT2":
                cmd = "../../aligners/bin/hisat2-build ../../data/%s.fa %s" % (genome, genome)
                cmd = cmd + "; ../../aligners/bin/hisat2-build -p 4 ../../data/%s.fa --snp ../../data/%s.snp %s_snp" % (genome, genome, genome)
                cmd = cmd + "; ../../aligners/bin/hisat2-build -p 4 ../../data/%s.fa --ss ../../data/%s.ss --exon ../../data/%s.exon %s_tran" % (genome, genome, genome, genome)
                cmd = cmd + "; ../../aligners/bin/hisat2-build -p 4 ../../data/%s.fa --snp ../../data/%s.snp --ss ../../data/%s.ss --exon ../../data/%s.exon %s_snp_tran" % (genome, genome, genome, genome, genome)
            elif aligner == "HISAT":
                cmd = "../../aligners/bin/hisat-build ../../data/%s.fa %s" % (genome, genome)
                cmd = cmd + "; ../../aligners/bin/tophat -G ../../data/%s.gtf --transcriptome-index=gtf %s; rm -rf tophat_out" % (genome, genome)
            elif aligner == "Bowtie":
                cmd = "../../aligners/bin/bowtie-build ../../data/%s.fa %s" % (genome, genome)
            elif aligner == "STAR":
                cmd = "../../aligners/bin/STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles ../../data/%s.fa" % (genome)
                cmd = cmd + "; mkdir gtf; ../../aligners/bin/STAR --runMode genomeGenerate --genomeDir gtf --genomeFastaFiles ../../data/%s.fa --sjdbGTFfile ../../data/%s.gtf --sjdbOverhang 99 --runThreadN 4" % (genome, genome)
            elif aligner == "GSNAP":
                cmd = "../../aligners/bin/gmap_build -B ../../aligners/bin -D . -d %s ../../data/%s.fa" % (genome, genome)
            elif aligner == "BWA":
                cmd = "../../aligners/bin/bwa -p %s.fa ../../data/%s.fa" % (genome, genome)
            else:
                assert False
            print >> sys.stderr, cmd
            os.system(cmd)
            os.chdir("..")

    os.chdir("..")
            
    
if __name__ == "__main__":
    build_indexes()
