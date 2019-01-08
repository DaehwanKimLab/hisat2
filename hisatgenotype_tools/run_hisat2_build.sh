#!/bin/bash -l
#SBATCH --job-name=infphio.genotype.hisat2-build
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400G
#SBATCH --partition=lrgmem
#SBATCH --time=168:0:0 
#SBATCH --workdir=/home-1/dkim136@jhu.edu/infphilo/hisat2/hisat2/evaluation/tests/genotype

/home-1/dkim136@jhu.edu/infphilo/hisat2/hisat2/hisat2-build -p 4 --snp genotype_genome.snp --haplotype genotype_genome.haplotype genotype_genome.fa genotype_genome
