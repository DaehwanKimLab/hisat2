#!/usr/bin/env python
#
# Copyright 2017, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT-genotype.
#
# HISAT-genotype is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT-genotype is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT-genotype.  If not, see <http://www.gnu.org/licenses/>.
#


import sys, os, subprocess, re
import random
import argparse
from datetime import datetime

##################################################
#   Common Arguments
##################################################
def args_common(parser, threads = True, debug = False):
    if threads:
        parser.add_argument("-p", "--threads",
                            dest="threads",
                            type=int,
                            default=1,
                            help="Number of threads")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='Print statistics to stderr')
    if debug:
        parser.add_argument("--debug",
                            dest="debug",
                            type=str,
                            default="",
                            help="Test database or code (options: basic, pair, full, single-end, test_list, test_id)(e.g., test_id:10,basic)")

def args_databases(parser, genome = False):
    if genome:
        parser.add_argument("-x", "--ref-genome",
                            dest="genotype_genome",
                            type=str,
                            default="",
                            help="Base name for genome index if not genotype_genome (default: empty)")
    parser.add_argument("--base", "--base-fname",
                        dest="base_fname",
                        type=str,
                        default="",
                        help="Base file name for index, variants, haplotypes, etc. (e.g. hla, rbg, codis) (default: empty)")
    parser.add_argument("--locus-list",
                        dest="locus_list",
                        type=str,
                        default="",
                        help="A comma-separated list of gene names (default: empty, all genes)")

def args_set_aligner(parser, missmatch = True):
    parser.add_argument("--aligner",
                        dest="aligner",
                        type=str,
                        default="hisat2",
                        help="Set aligner to use (ex. hisat2, bowtie2) (default: hisat2)")
    parser.add_argument("--linear-index",
                        dest="graph_index",
                        action="store_false",
                        help="Use linear index") 
    if missmatch:
        parser.add_argument("--num-mismatch",
                            dest="num_mismatch",
                            type=int,
                            default=0,
                            help="Maximum number of mismatches per read alignment to be considered (default: 0)")

def args_aligner_inputs(parser, keep=False):
    parser.add_argument('-f', '--fasta',
                        dest='fastq',
                        action='store_false',
                        help='FASTA format')
    parser.add_argument("-U",
                        dest="read_fname_U",
                        type=str,
                        default="",
                        help="filename for single-end reads")
    parser.add_argument("-1",
                        dest="read_fname_1",
                        type=str,
                        default="",
                        help="filename for paired-end reads")
    parser.add_argument("-2",
                        dest="read_fname_2",
                        type=str,
                        default="",
                        help="filename for paired-end reads")
    if keep:
        parser.add_argument("--keep-alignment",
                            dest="keep_alignment",
                            action="store_true",
                            help="Keep alignment file")

def args_assembly(parser):
    parser.add_argument("--assembly",
                        dest="assembly",
                        action="store_true",
                        help="Perform assembly")
    parser.add_argument("--assembly-name",
                        dest="output_base",
                        type=str,
                        default="assembly_graph",
                        help="Assembly base file name (default: assembly_graph)")
    parser.add_argument("--assembly-verbose",
                        dest="assembly_verbose",
                        action="store_true",
                        help="Output intermediate assembly information")

def args_input_output(parser, indir = True, outdir = True):
    if indir:
        parser.add_argument('--in-dir',
                            dest="read_dir",
                            type=str,
                            default="",
                            help='Input directory (e.g. read_input) (default: (empty))')
    if outdir:
        parser.add_argument("--out-dir",
                            dest="out_dir",
                            type=str,
                            default="",
                            help='Output directory (default: (empty))')

def args_bamfile(parser):
    parser.add_argument("--bamfile",
                        dest="alignment_fname",
                        type=str,
                        default="",
                        help="BAM file name")

def args_reference_type(parser):
    parser.add_argument("--reference-type",
                        dest="reference_type",
                        type=str,
                        default="gene",
                        help="Reference type: gene, chromosome, and genome (default: gene)")

def args_no_partial(parser):
    parser.add_argument('--no-partial',
                        dest='partial',
                        action='store_false',
                        help='Include partial alleles (e.g. A_nuc.fasta)')

def args_single_end(parser):
    parser.add_argument('--single-end',
                        dest='paired',
                        action='store_false',
                        help='Choose to set files to single-ended reads (unnecessary when using -U -1 -2 options)')
 
def args_var_gaps(parser):
    parser.add_argument("--inter-gap",
                        dest="inter_gap",
                        type=int,
                        default=30,
                        help="Maximum distance for variants to be in the same haplotype")
    parser.add_argument("--intra-gap",
                        dest="intra_gap",
                        type=int,
                        default=50,
                        help="Break a haplotype into several haplotypes")


##################################################
#   Script Specific Arguments
##################################################
def args_extract_reads(parser):
    parser.add_argument("--suffix",
                        dest="suffix",
                        type=str,
                        default="fq.gz",
                        help="Read file suffix (Default: fq.gz)")  
    parser.add_argument('--simulation',
                        dest='simulation',
                        action='store_true',
                        help='Simulated reads (Default: False)')    
    parser.add_argument("--pp", "--threads-aprocess",
                        dest="threads_aprocess",
                        type=int,
                        default=1,
                        help="Number of threads a process")
    parser.add_argument("--max-sample",
                        dest="max_sample",
                        type=int,
                        default=sys.maxint,
                        help="Number of samples to be extracted (default: sys.maxint)")
    parser.add_argument("--job-range",
                        dest="job_range",
                        type=str,
                        default="0,1",
                        help="two numbers (e.g. 1,3)")
    parser.add_argument("--extract-whole",
                        dest="extract_whole",
                        action='store_true',
                        help="Extract all reads")


def args_extract_vars(parser):
    parser.add_argument("--whole-haplotype",
                        dest="whole_haplotype",
                        action="store_true",
                        help="Include partial alleles (e.g. A_nuc.fasta)")
    parser.add_argument("--min-var-freq",
                        dest="min_var_freq",
                        type=float,
                        default=0.0,
                        help="Exclude variants whose freq is below than this value in percentage (Default: 0.0)")    
    parser.add_argument("--ext-seq",
                        dest="ext_seq_len",
                        type=int,
                        default=0,
                        help="Length of extra sequences flanking backbone sequences (Default: 0)")
    parser.add_argument("--leftshift",
                        dest="leftshift",
                        action="store_true",
                        help="Shift deletions to the leftmost")

def args_locus(parser):
    parser.add_argument("--simulate-interval",
                        dest="simulate_interval",
                        type=int,
                        default=10,
                        help="Reads simulated at every these base pairs (default: 10)")
    parser.add_argument("--read-len",
                        dest="read_len",
                        type=int,
                        default=100,
                        help="Length of simulated reads (default: 100)")
    parser.add_argument("--fragment-len",
                        dest="fragment_len",
                        type=int,
                        default=350,
                        help="Length of fragments (default: 350)")
    parser.add_argument("--best-alleles",
                        dest="best_alleles",
                        action='store_true',
                        help="")
    parser.add_argument("--random-seed",
                        dest="random_seed",
                        type=int,
                        default=1,
                        help="A seeding number for randomness (default: 1)")
    parser.add_argument("--num-editdist",
                        dest="num_editdist",
                        type=int,
                        default=2,
                        help="Maximum number of mismatches per read alignment to be considered (default: 2)")
    parser.add_argument("--perbase-errorrate",
                        dest="perbase_errorrate",
                        type=float,
                        default=0.0,
                        help="Per basepair error rate in percentage when simulating reads (default: 0.0)")
    parser.add_argument("--perbase-snprate",
                        dest="perbase_snprate",
                        type=float,
                        default=0.0,
                        help="Per basepair SNP rate in percentage when simulating reads (default: 0.0)")
    parser.add_argument("--skip-fragment-regions",
                        dest="skip_fragment_regions",
                        type=str,
                        default="",
                        help="A comma-separated list of regions from which no reads originate, e.g., 500-600,1200-1400 (default: None).")
    parser.add_argument('--verbose-level',
                        dest='verbose_level',
                        type=int,
                        default=0,
                        help='also print some statistics to stderr (default: 0)')
    parser.add_argument("--no-error-correction",
                        dest="error_correction",
                        action="store_false",
                        help="Correct sequencing errors")
    parser.add_argument("--only-locus-list",
                        dest="only_locus_list",
                        type=str,
                        default="",
                        help="A comma-separated list of genes (default: empty, all genes)")
    parser.add_argument("--discordant",
                        dest="discordant",
                        action="store_true",
                        help="Allow discordantly mapped pairs or singletons")
    parser.add_argument("--type-primary-exons",
                        dest="type_primary_exons",
                        action="store_true",
                        help="Look at primary exons first")
    parser.add_argument("--keep-low-abundance-alleles",
                        dest="remove_low_abundance_alleles",
                        action="store_false",
                        help="Do not remove alleles with low abundance while performing typing")
    parser.add_argument("--display-alleles",
                        dest="display_alleles",
                        type=str,
                        default="",
                        help="A comma-separated list of alleles to display in HTML (default: empty)")


def args_build_genome(parser):
    parser.add_argument("--commonvar",
                        dest="use_commonvar",
                        action="store_true",
                        help="Include common variants from dbSNP")
    parser.add_argument("--clinvar",
                        dest="use_clinvar",
                        action="store_true",
                        help="Include variants from ClinVar database")

def args_locus_samples(parser):
    parser.add_argument("--region-list",
                        dest="region_list",
                        type=str,
                        default="",
                        help="A comma-separated list of regions (overwrites --base and --locus-list) (default: empty)")
    parser.add_argument("--num-editdist",
                        dest="num_editdist",
                        type=int,
                        default=2,
                        help="Maximum number of mismatches per read alignment to be considered (default: 2)")
    parser.add_argument("--max-sample",
                        dest="max_sample",
                        type=int,
                        default=sys.maxint,
                        help="Number of samples to be analyzed (default: sys.maxint)")
    parser.add_argument('--platinum-check',
                        dest='platinum_check',
                        action='store_true',
                        help='Check for concordance of platinum genomes')

def args_HLA_genotyping_PGs(parser, gold_allele_info):
    parser.add_argument("--hla-list",
                        dest="hla_list",
                        type=str,
                        default="A,B,C,DQA1,DQB1,DRB1",
                        help="A comma-separated list of HLA genes (default: A,B,C,DQA1,DQB1,DRB1)")
    genomes_default = ','.join(gold_allele_info.keys())
    parser.add_argument("--genome-list",
                        dest="genome_list",
                        type=str,
                        default=genomes_default,
                        help="A comma-separated list of genomes (default: %s)" % genomes_default)
    parser.add_argument("--exclude-allele-list",
                        dest="exclude_allele_list",
                        type=str,
                        default="",
                        help="A comma-separated list of alleles to be excluded")

def args_hla_cyp(parser):
    parser.add_argument("--reads",
                        dest="read_fname",
                        type=str,
                        default="",
                        help="Fastq read file name")    
    parser.add_argument("--allele-list",
                        dest = "default_allele_list",
                        type=str,
                        default="",
                        help="A comma-separated list of HLA alleles to be tested. Alleles are retrieved from default backbone data (all alleles included in backbone).")
    parser.add_argument('--partial',
                        dest='partial',
                        action='store_true',
                        help='Include partial alleles (e.g. A_nuc.fasta)')
    parser.add_argument("--aligner-list",
                        dest="aligners",
                        type=str,
                        default="",
                        help="A comma-separated list of aligners (Overwrites --aligner option) (e.g. hisat2.graph,hisat2.linear,bowtie2.linear)")
    parser.add_argument("--simulate-interval",
                        dest="simulate_interval",
                        type=int,
                        default=1,
                        help="Reads simulated at every these base pairs (default: 1)")
    parser.add_argument("--coverage",
                        dest="coverage",
                        action='store_true',
                        help="Experimental purpose (assign reads based on coverage)")
    parser.add_argument("--best-alleles",
                        dest="best_alleles",
                        action='store_true',
                        help="")
    parser.add_argument("--exclude-allele-list",
                        dest="exclude_allele_list",
                        type=str,
                        default="",
                        help="A comma-separated list of alleles to be excluded. Enter a number N to randomly select N alleles for exclusion and N non-excluded alleles for testing (2N tested in total).")
    parser.add_argument("--novel_allele_detection",
                        dest="novel_allele_detection",
                        action='store_true',
                        help="Change test to detection of new alleles. Report sensitivity and specificity rate at the end.")

def args_convert_codis(parser):
    parser.add_argument("--min-freq",
                        dest="min_freq",
                        type=float,
                        default=0.0,
                        help="minimum allele frequency (default: 0.0)")    


