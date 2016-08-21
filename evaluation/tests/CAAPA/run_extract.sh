#!/bin/bash -l
#SBATCH --job-name=infphio.HLA.CP.extract.genome
#SBATCH --nodes=1
#SBATCH --cpus-per-task=44
#SBATCH --mem=440G
#SBATCH --partition=lrgmem
#SBATCH --time=166:0:0 
#SBATCH --workdir=/home-1/dkim136@jhu.edu/infphilo/hisat2/evaluation/tests/HLA_CP_extract_genome_partial
/home-1/dkim136@jhu.edu/infphilo/hisat2/evaluation/tests/HLA_CP_extract_genome_partial/hisat2_cp_extract_reads.py --reference-type genome --partial -p 44 --job-range 8,12
