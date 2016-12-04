#!/bin/bash -l
#SBATCH --job-name=infphio.genotype
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400G
#SBATCH --partition=lrgmem
#SBATCH --time=168:0:0 
#SBATCH --workdir=/home-1/dkim136@jhu.edu/infphilo/hisat2/hisat2/evaluation/tests/genotype

/home-1/dkim136@jhu.edu/infphilo/hisat2/hisat2/evaluation/tests/genotype/hisat2_build_genotype_genome.py --verbose -p 4 genome.fa genotype_genome
