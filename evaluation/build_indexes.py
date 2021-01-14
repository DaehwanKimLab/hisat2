#!/usr/bin/env python3

import sys, os
from argparse import ArgumentParser, FileType

# GRCh38 release 84
def build_indexes(input_aligners,
                  input_genomes,
                  fresh):
    # Build indexes
    if not os.path.exists("indexes"):
        os.mkdir("indexes")
    os.chdir("indexes")
    
    if input_aligners:
        aligners = input_aligners
    else:
        aligners = ["hisat2", "hisat-gt", "hisat", "bowtie", "star", "gsnap", "bwa", "minimap2", "salmon", "kallisto"]
        
    if input_genomes:
        genomes = input_genomes
    else:
        genomes = ["22_20-21M", "22", "genome"]

    for genome in genomes:
        for aligner in aligners:
            if genome == "genome":
                dir = aligner
            else:
                dir = aligner + "_" + genome
            if os.path.exists(dir):
                if fresh:
                    os.system("rm -rf %s" % (dir))
                else:
                    continue
            os.mkdir(dir)
            os.chdir(dir)
            if aligner == "hisat2":
                cmd = "../../aligners/bin/hisat2-build ../../data/%s.fa %s" % (genome, genome)
                cmd = cmd + "; ../../aligners/bin/hisat2-build -p 4 ../../data/%s.fa --snp ../../data/%s.snp --haplotype ../../data/%s.haplotype %s_snp" % (genome, genome, genome, genome)
                cmd = cmd + "; ../../aligners/bin/hisat2-build -p 4 ../../data/%s.fa --ss ../../data/%s.ss --exon ../../data/%s.exon %s_tran" % (genome, genome, genome, genome)
                cmd = cmd + "; ../../aligners/bin/hisat2-build -p 4 ../../data/%s.fa --snp ../../data/%s.snp --haplotype ../../data/%s.haplotype --ss ../../data/%s.ss --exon ../../data/%s.exon %s_snp_tran" % (genome, genome, genome, genome, genome, genome)
            elif aligner == "hisat-gt":
                cmd = "../../aligners/bin/hisat2_extract_transcript_graph.py -g ../../data/{genome}.gtf -r ../../data/{genome}.fa -s ../../data/{genome}.snp -p ../../data/{genome}.haplotype --chrtome -o ../../data/{genome}".format(genome=genome)
                cmd = cmd + ";mv ../../data/{genome}.gt.snp ../../data/{genome}.gt.snp.org; ../../aligners/bin/filter_snp.py ../../data/{genome}.gt.snp.org ../../data/{genome}.gt.fa ../../data/{genome}.gt.snp".format(genome=genome)
                cmd = cmd + ";../../aligners/bin/hisat2-build -p 4 ../../data/%s.gt.fa --ss ../../data/%s.gt.ss %s.gt" % (genome, genome, genome)
                cmd = cmd + ";../../aligners/bin/hisat2-build -p 4 ../../data/{genome}.gt.fa --ss ../../data/{genome}.gt.ss --snp ../../data/{genome}.gt.snp --haplotype ../../data/{genome}.gt.haplotype {genome}_snp.gt".format(genome=genome)
            elif aligner == "hisat":
                cmd = "../../aligners/bin/hisat-build ../../data/%s.fa %s" % (genome, genome)
                cmd = cmd + "; ../../aligners/bin/tophat -G ../../data/%s.gtf --transcriptome-index=gtf %s; rm -rf tophat_out" % (genome, genome)
            elif aligner == "bowtie":
                cmd = "../../aligners/bin/bowtie-build ../../data/%s.fa %s" % (genome, genome)
            elif aligner == "bowtie2":
                cmd = "../../aligners/bin/bowtie2-build --threads 6 ../../data/%s.fa %s" % (genome, genome)
            elif aligner == "star":
                cmd = "../../aligners/bin/STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles ../../data/%s.fa" % (genome)
                cmd = cmd + "; mkdir gtf; ../../aligners/bin/STAR --runMode genomeGenerate --genomeDir gtf --genomeFastaFiles ../../data/%s.fa --sjdbGTFfile ../../data/%s.gtf --sjdbOverhang 99 --runThreadN 4" % (genome, genome)
            elif aligner == "gsnap":
                cmd = "../../aligners/bin/gmap_build -B ../../aligners/bin -D . -d %s ../../data/%s.fa" % (genome, genome)
            elif aligner == "bwa":
                cmd = "../../aligners/bin/bwa index -p %s.fa ../../data/%s.fa" % (genome, genome)
            elif aligner == "minimap2":
                cmd = "../../aligners/bin/minimap2 -x sr -d %s.mmi ../../data/%s.fa" % (genome, genome)
            elif aligner == "vg":
                assert False
            elif aligner == 'salmon':
                cmd = ""
                if genome == 'genome':
                    # download cdna from ensembl-84
                    # wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
                    cmd = ' [[ -e ../../data/genome.cdna.fa ]] || (wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz && gzip -cd Homo_sapiens.GRCh38.cdna.all.fa.gz > ../../data/genome.cdna.fa)'
                    cmd += ' && ../../aligners/bin/salmon index -t ../../data/genome.cdna.fa -i genome -p 5'
            elif aligner == 'kallisto':
                cmd = ""
                if genome == 'genome':
                    # download cdna from ensembl-84
                    # wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
                    cmd = ' [[ -e ../../data/genome.cdna.fa ]] || (wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz && gzip -cd Homo_sapiens.GRCh38.cdna.all.fa.gz > ../../data/genome.cdna.fa)'
                    cmd += ' && ../../aligners/bin/kallisto index -i genome ../../data/genome.cdna.fa'
            else:
                assert False
            print(cmd, file=sys.stderr)
            os.system(cmd)
            os.chdir("..")

    os.chdir("..")
            
    
if __name__ == "__main__":
    parser = ArgumentParser(
        description='Build indexes for HISAT2, HISAT-GT, STAR, Bowtie1/2, GSNAP, BWA-mem, etc.')
    parser.add_argument('--aligner-list',
                        dest='aligner_list',
                        type=str,
                        default="",
                        help='comma-separated list of aligners (e.g. hisat2,hisat-gt,bowtie2,bwa')
    parser.add_argument('--genome-list',
                        dest='genome_list',
                        type=str,
                        default="",
                        help='comma-separated list of genomes (e.g. 22,genome')
    parser.add_argument('--fresh',
                        dest='fresh',
                        action='store_true',
                        help='delete existing alignment related directories (e.g. hisat2_single)')
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()

    aligners = []
    for aligner in args.aligner_list.split(','):
        if aligner == "":
            continue
        aligners.append(aligner)

    genomes = []
    for genome in args.genome_list.split(','):
        if genome == "":
            continue
        genomes.append(genome)

    build_indexes(aligners,
                  genomes,
                  args.fresh)
    
