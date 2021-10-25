# HISAT2 Read Simulator Expanded

## Author Information

The original simulation code was completed in 2015 by Dr. Daehwan Kim. 
Enhancements by Micah Thornton (Kim Lab) beginning 2020. 

## Overview 

The original HISAT2 Simulator for synthesizing reads from a given genetic sequence file will generate reads 
that come from one of several sequencing kinds.  The primary goal of this branch is to update the 
original simulator which does not currently allow for the simulation of BiSulfite reads for analysis 
for use with the original BSmooth statistical pipeline (Kasper Hansen et al 2012-2020) as well as the simulation
of reads of other kinds. 

## Original Simulator Information 
Simulate reads from GENOME (fasta) and GTF files

positional arguments:
  genome_file           input GENOME file
  gtf_file              input GTF file
  snp_file              input SNP file
  base_fname            output base filename

optional arguments:
  -h, --help            show this help message and exit
  -d, --dna             DNA-seq reads (default: RNA-seq reads)
  --single-end          single-end reads (default: paired-end reads)
  -r READ_LEN, --read-length READ_LEN
                        read length (default: 100)
  -f FRAG_LEN, --fragment-length FRAG_LEN
                        fragment length (default: 250)
  -n NUM_FRAG, --num-fragment NUM_FRAG
                        number of fragments (default: 1000000)
  -e EXPR_PROFILE, --expr-profile EXPR_PROFILE
                        expression profile: flux or constant (default: flux)
  --repeat-info REPEAT_FNAME
                        repeat information filename
  --error-rate ERROR_RATE
                        per-base sequencing error rate (%) (default: 0.0)
  --max-mismatch MAX_MISMATCH
                        max mismatches due to sequencing errors (default: 3)
  --random-seed RANDOM_SEED
                        random seeding value (default: 0)
  --snp-prob SNP_PROB   probability of a read including a snp when the read spans the snp ranging from 0.0 to 1.0 (default: 1.0)
  --sanity-check        sanity check
  -v, --verbose         also print some statistics to stderr
  --version             show program's version number and exit

## Updates and Edits 

09-17-2020:  Added a class for methylation status of CpGs to be generated using a Bernoulli trials approach. 

09-18-2020:  In order to verify the proper functionality of the simulation script, one may use the included Makefile by running 'make Tests'; 
the resulting test information will be outputted to a log file called sim-test-res.log. This is a work in progress, for the intention of testing
various simulated sequencing reads. 

