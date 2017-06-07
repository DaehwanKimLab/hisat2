#!/bin/bash -l
#SBATCH --job-name=infphio.HLA.ILMN.extract.genome
#SBATCH --nodes=1
#SBATCH --cpus-per-task=17
#SBATCH --mem=120G
#SBATCH --partition=shared
#SBATCH --time=166:0:0 
#SBATCH --workdir=/home-1/dkim136@jhu.edu/infphilo/hisat2/evaluation/tests/HLA_novel

/home-1/dkim136@jhu.edu/infphilo/hisat2/evaluation/tests/HLA_novel/scripts/extract_reads.py --base-fname genotype_genome --reference-type genome --read-dir /home-1/dkim136@jhu.edu/ssalzbe1/users/infphilo/platinum_genomes --out-dir ILMN -p 17

