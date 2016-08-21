#!/bin/bash -l
#SBATCH --job-name=infphio.HLA.CP
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --partition=shared
#SBATCH --time=12:0:0 
#SBATCH --workdir=/home-1/dkim136@jhu.edu/infphilo/hisat2/evaluation/tests/HLA_CP_extract_genome_partial
/home-1/dkim136@jhu.edu/infphilo/hisat2/evaluation/tests/HLA_CP_extract_genome_partial/hisat2_test_HLA_genotyping_CAAPA.py -p 24 > caapa_hla.txt
