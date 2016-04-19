#!/usr/bin/env python

#
# Copyright 2016, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
#


import os, sys, subprocess, re
import inspect
from argparse import ArgumentParser, FileType


"""
"""
def read_genome(genome_file):
    chr_dic, chr_names = {}, []
    chr_name, sequence = "", ""
    for line in genome_file:
        if line.startswith(">"):
            if chr_name and sequence:
                chr_dic[chr_name] = sequence
                chr_names.append(chr_name)
            chr_name = line.strip().split()[0][1:]
            sequence = ""
        else:
            sequence += line.strip()
    if chr_name and sequence:
        chr_dic[chr_name] = sequence
        chr_names.append(chr_name)
    return chr_dic, chr_names


"""
"""
def read_sequences(fname):
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
def build_genotype_genome(reference,
                          base_fname,                          
                          gap,
                          split,
                          verbose):    
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(build_genotype_genome))
    ex_path = os.path.dirname(curr_script)

    def check_files(fnames):
        for fname in fnames:
            if not os.path.exists(fname):
                return False
        return True

    # Download HISAT2 index
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]
    if not check_files(HISAT2_fnames):
        os.system("wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz; tar xvzf grch38.tar.gz; rm grch38.tar.gz")
        hisat2_inspect = os.path.join(ex_path, "hisat2-inspect")
        os.system("%s grch38/genome > genome.fa" % hisat2_inspect)
        os.system("samtools faidx genome.fa")

    # Load genomic sequences
    chr_dic, chr_names = read_genome(open(reference))

    # Genes to be genotyped
    genotype_genes = {}

    # Clone a git repository, IMGTHLA
    if not os.path.exists("IMGTHLA"):
        os.system("git clone https://github.com/jrob119/IMGTHLA.git")

    # Extract HLA variants, backbone sequence, and other sequeces
    HLA_fnames = ["hla_backbone.fa",
                  "hla.ref",
                  "hla.snp",
                  "hla.haplotype",
                  "hla.link"]

    if not check_files(HLA_fnames):
        extract_hla_script = os.path.join(ex_path, "hisat2_extract_HLA_vars.py")
        extract_cmd = [extract_hla_script]
        # if partial:
        #    extract_cmd += ["--partial"]
        extract_cmd += ["--gap", "30",
                        "--split", "50"]
        if verbose:
            print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
        proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()
        if not check_files(HLA_fnames):
            print >> sys.stderr, "Error: extract_HLA_vars failed!"
            sys.exit(1)

    # Read HLA genes
    if os.path.exists("hla.ref"):
        for line in open("hla.ref"):
            HLA_name, chr, left, right, length, exon_str = line.strip().split()
            left, right = int(left), int(right)
            length = int(length)
            if chr not in chr_names:
                continue
            if chr not in genotype_genes:
                genotype_genes[chr] = []
            genotype_genes[chr].append([left, right, length, HLA_name, "hla"])

    # Write genotype genome
    var_num, haplotype_num = 0, 0
    genome_out_file = open("%s.fa" % base_fname, 'w')
    var_out_file = open("%s.snp" % base_fname, 'w')
    haplotype_out_file = open("%s.haplotype" % base_fname, 'w')
    link_out_file = open("%s.link" % base_fname, 'w')
    for chr in chr_names:
        assert chr in chr_dic
        chr_seq = chr_dic[chr]
        chr_len = len(chr_seq)
        if chr not in genotype_genes:
            continue

        chr_genes = genotype_genes[chr]
        def gene_cmp(a, b):
            a_left, a_right, a_length = a[:3]
            b_left, b_right, b_length = b[:3]
            if a_left != b_left:
                return a_left - b_left
            if a_right != b_right:
                return a_right - b_right
            return a_lenght - b_length
        chr_genes = sorted(chr_genes, cmp=gene_cmp)

        out_chr_seq = ""
        
        off = 0
        prev_right = 0
        for gene in chr_genes:
            left, right, length, name, family = gene
            # Read HLA backbone sequences
            allele_seqs = read_sequences("%s_backbone.fa" % family)

            # Read HLA variants
            allele_vars = read_variants("hla.snp")

            # Read HLA haplotypes
            allele_haplotypes = read_haplotypes("hla.haplotype")

            # Read HLA link information between haplotypes and variants
            links = read_links("hla.link")

            if name not in allele_seqs or \
                    name not in allele_vars or \
                    name not in allele_haplotypes:
                continue
            allele_seq = allele_seqs[name]
            vars = allele_vars[name]
            haplotypes = allele_haplotypes[name]
            assert length == len(allele_seq)
            assert left < chr_len and right < chr_len
            # Skipping overlapping genes
            if left < prev_right:
                print >> sys.stderr, "Warning: skipping %s ..." % (name)
                continue

            varID2htID = {}

            assert left < right
            prev_length = right - left + 1
            assert prev_length <= length

            if prev_right < left:
                out_chr_seq += chr_seq[prev_right:left]
            out_chr_seq += allele_seq

            # Output variants (genotype_genome.snp)
            for var in vars:
                var_left, var_type, var_data, var_id = var
                new_var_id = "hv%d" % var_num
                varID2htID[var_id] = new_var_id
                print >> var_out_file, "%s\t%s\t%s\t%d\t%s" % \
                    (new_var_id, var_type, chr, var_left + left + off, var_data)
                var_num += 1

            # Output haplotypes (genotype_genome.haplotype)
            for haplotype in haplotypes:
                ht_left, ht_right, ht_vars = haplotype
                new_ht_vars = []
                for var_id in ht_vars:
                    assert var_id in varID2htID
                    new_ht_vars.append(varID2htID[var_id])
                print >> haplotype_out_file, "ht%d\t%s\t%d\t%d\t%s" % \
                    (haplotype_num, chr, ht_left + left + off, ht_right + left + off, ','.join(new_ht_vars))
                haplotype_num += 1

            # Output link information between alleles and variants (genotype_genome.link)
            for link in links:
                var_id, allele_names = link
                if var_id not in varID2htID:
                    continue
                new_var_id = varID2htID[var_id]
                print >> link_out_file, "%s\t%s" % (new_var_id, allele_names)
                
            off += (length - prev_length)

            prev_right = right + 1
            
                
        out_chr_seq += chr_seq[prev_right:]

        assert len(out_chr_seq) == len(chr_seq) + off

        # Output chromosome sequence
        print >> genome_out_file, ">%s" % (chr)
        line_width = 60
        for s in range(0, len(out_chr_seq), line_width):
            print >> genome_out_file, out_chr_seq[s:s+line_width]

    genome_out_file.close()
    var_out_file.close()
    haplotype_out_file.close()
    link_out_file.close()

        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Extract HLA variants from HLA multiple sequence alignments")
    parser.add_argument("reference",
                        nargs='?',
                        type=str,
                        help="Reference genome")
    parser.add_argument("base_fname",
                        nargs='?',
                        type=str,
                        help="base filename for genotype genome")
    parser.add_argument("-g", "--gap",
                        dest="gap",
                        type=int,
                        default=30,
                        help="Maximum distance for variants to be in the same haplotype")
    parser.add_argument("-s", "--split",
                        dest="split",
                        type=int,
                        default=50,
                        help="Break a haplotype into several haplotypes")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="also print some statistics to stderr")

    args = parser.parse_args()
    if not args.reference or not args.base_fname:
        parser.print_help()
        sys.exit(1)
    if args.gap > args.split:
        print >> sys.stderr, "Error: -g/--gap (%d) must be smaller than -s/--split (%d)" % (args.gap, args.split)
        sys.exit(1)
    build_genotype_genome(args.reference,
                          args.base_fname,
                          args.gap,
                          args.split,
                          args.verbose)
