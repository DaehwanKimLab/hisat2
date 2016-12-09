#!/bin/bash -l
#SBATCH --job-name=infphio.genotype
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400G
#SBATCH --partition=lrgmem
#SBATCH --time=168:0:0 
#SBATCH --workdir=/home-1/dkim136@jhu.edu/infphilo/hisat2/hisat2/evaluation/tests/HLA_novel

/home-1/dkim136@jhu.edu/infphilo/hisat2/hisat2/evaluation/tests/HLA_novel/hisatgenotype_build_genome.py -p 4 --verbose --commonvar genome.fa genotype_genome
