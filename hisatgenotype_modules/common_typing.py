#!/usr/bin/env python

import sys, os
import math
import random
from copy import deepcopy


"""
"""
def read_genome(genome_file):
    chr_dic, chr_names, chr_full_names = {}, [], []
    chr_name, chr_full_name, sequence = "", "", ""
    for line in genome_file:
        if line.startswith(">"):
            if chr_name and sequence:
                chr_dic[chr_name] = sequence
                chr_names.append(chr_name)
            chr_full_name = line.strip()[1:]
            chr_name = line.strip().split()[0][1:]
            chr_full_names.append(chr_full_name)
            sequence = ""
        else:
            sequence += line.strip()
    if chr_name and sequence:
        chr_dic[chr_name] = sequence
        chr_names.append(chr_name)
        chr_full_names.append(chr_full_name)
    return chr_dic, chr_names, chr_full_names


"""
"""
def read_allele_sequences(fname):
    allele_seqs = {}
    allele_name, sequence = "", ""
    for line in open(fname):
        if line.startswith(">"):
            if allele_name != "" and allele_name not in allele_seqs:
                allele_seqs[allele_name] = sequence
            allele_name = line.strip()[1:]
            sequence = ""
        else:
            sequence += line.strip()
    if allele_name != "" and allele_name not in allele_seqs:
        allele_seqs[allele_name] = sequence
    return allele_seqs


"""
"""
def read_variants(fname):
    allele_vars = {}
    for line in open(fname):
        var_id, type, allele_name, left, data = line.strip().split()
        left = int(left)
        if type == "deletion":
            data = int(data)
        if allele_name not in allele_vars:
            allele_vars[allele_name] = []
        allele_vars[allele_name].append([left, type, data, var_id])
    return allele_vars


"""
"""
def read_haplotypes(fname):
    allele_haplotypes = {}
    for line in open(fname):
        haplotype_id, allele_name, left, right, vars = line.strip().split()
        vars = vars.split(',')
        left, right = int(left), int(right)
        if allele_name not in allele_haplotypes:
            allele_haplotypes[allele_name] = []
        allele_haplotypes[allele_name].append([left, right, vars])
    return allele_haplotypes


"""
"""
def read_links(fname):
    links = []
    for line in open(fname):
        var_id, allele_names = line.strip().split('\t')
        links.append([var_id, allele_names])
    return links


"""
Compare two variants
"""
def compare_vars(a, b):
    a_pos, a_type, a_data = a[:3]
    b_pos, b_type, b_data = b[:3]

    if a_pos != b_pos:
        return a_pos - b_pos
    if a_type != b_type:
         if a_type == 'I':
             return -1
         elif b_type == 'I':
             return 1
         if a_type == 'S':
             return -1
         else:
             return 1
    if a_data < b_data:
        return -1
    elif a_data > b_data:
        return 1
    else:
        return 0


"""
"""
def check_files(fnames):
    for fname in fnames:
        if not os.path.exists(fname):
            return False
    return True
   
    
"""
Download GRCh38 human reference and HISAT2 indexes
"""
def download_genome_and_index(ex_path):
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]
    if not check_files(HISAT2_fnames):
        os.system("wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz; tar xvzf grch38.tar.gz; rm grch38.tar.gz")
        hisat2_inspect = os.path.join(ex_path, "hisat2-inspect")
        os.system("%s grch38/genome > genome.fa" % hisat2_inspect)
        os.system("samtools faidx genome.fa")


