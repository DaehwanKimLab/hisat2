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


import sys, os, subprocess, re
import inspect, random
import math
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
def genotype_test(reference_type,
                  base_fname,
                  threads,
                  simulate_interval,
                  num_mismatch,
                  verbose,
                  daehwan_debug):
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(genotype_test))
    ex_path = os.path.dirname(curr_script)

    def check_files(fnames):
        for fname in fnames:
            if not os.path.exists(fname):
                return False
        return True

    # variants, backbone sequence, and other sequeces
    genotype_fnames = ["%s.fa" % base_fname,
                       "%s.gene" % base_fname,
                       "%s.snp" % base_fname,
                       "%s.haplotype" % base_fname,
                       "%s.link" % base_fname,
                       "%s.coord" % base_fname,
                       "%s.clnsig" % base_fname]
    # hisat2 graph index files
    genotype_fnames += ["%s.%d.ht2" % (base_fname, i+1) for i in range(8)]
    # hla files
    genotype_fnames += ["hla_sequences.fa"]                       

    if not check_files(genotype_fnames):
        print >> sys.stderr, "Error: some files are missing!"
        sys.exit(1)

    # Load genomic sequences
    chr_dic, chr_names = read_genome(open("%s.fa" % base_fname))

    # Read variants with clinical significance
    clnsigs = {}
    for line in open("%s.clnsig" % base_fname):
        var_id, var_gene, var_clnsig = line.strip().split('\t')
        clnsigs[var_id] = [var_gene, var_clnsig]

    vars = {}
    for line in open("%s.snp" % base_fname):
        var_id, type, chr, left, data = line.strip().split()
        left = int(left)
        if type == "deletion":
            data = int(data)
        vars[var_id] = [chr, left, type, data]

    # Read HLA alleles (names and sequences)
    HLAs = {}
    def read_HLA_alleles(fname, HLAs):
        for line in open(fname):
            if line.startswith(">"):
                HLA_name = line.strip().split()[0][1:]
                HLA_gene = HLA_name.split('*')[0]
                if not HLA_gene in HLAs:
                    HLAs[HLA_gene] = {}
                if not HLA_name in HLAs[HLA_gene]:
                    HLAs[HLA_gene][HLA_name] = ""
            else:
                HLAs[HLA_gene][HLA_name] += line.strip()
        return HLAs
    read_HLA_alleles("hla_sequences.fa", HLAs)

    # Test HLA genotyping
    test_list = []
    basic_test, pair_test = True, False
    if daehwan_debug:
        if "basic_test" in daehwan_debug:
            basic_test, pair_test = True, False
        else:
            basic_test, pair_test = False, True

    # daehwan - for debugging purposes
    basic_test, pair_test = False, True

    test_passed = {}
    test_list = []
    genes = HLAs.keys()
    if basic_test:
        for gene in genes:
            HLA_gene_alleles = HLAs[gene]
            for HLA_name in HLA_gene_alleles:
                if HLA_name.find("BACKBONE") != -1:
                    continue
                test_list.append([[HLA_name]])
    if pair_test:
        test_size = 500
        allele_count = 2
        for test_i in range(test_size):
            test_pairs = []
            for gene in genes:
                HLA_gene_alleles = []
                for allele in HLAs[gene]:
                    if allele.find("BACKBONE") != -1:
                        continue
                    HLA_gene_alleles.append(allele)
                nums = [i for i in range(len(HLA_gene_alleles))]
                random.shuffle(nums)
                test_pairs.append(sorted([HLA_gene_alleles[nums[i]] for i in range(allele_count)]))
            test_list.append(test_pairs)

    # DK - QC measures are necessary - see the Omixon paper (2013) at http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0078410
    # (1) read length of at least 75 bp
    # (2) exons 2 and 3 covered by 70% or higher as these exons of HLA-A,B,C are highly polymorphic
    for test_i in range(len(test_list)):
        if "test_id" in daehwan_debug:
            daehwan_test_ids = daehwan_debug["test_id"].split('-')
            if str(test_i + 1) not in daehwan_test_ids:
                continue

        print >> sys.stderr, "Test %d" % (test_i + 1)
        test_HLA_list = test_list[test_i]

        # daehwan - for debugging purposes
        # test_HLA_list = [["A*11:50Q", "A*11:01:01:01", "A*01:01:01:01"]]
        partial_alleles = []
        for test_HLA_names in test_HLA_list:
            for test_HLA_name in test_HLA_names:
                gene = test_HLA_name.split('*')[0]
                test_HLA_seq = HLAs[gene][test_HLA_name]
                seq_type = "partial" if test_HLA_name in partial_alleles else "full"
                print >> sys.stderr, "\t%s - %d bp (%s sequence)" % (test_HLA_name, len(test_HLA_seq), seq_type)

        HLA_reads_1, HLA_reads_2 = [], []

        # Simulate reads from two HLA alleles
        def simulate_reads(seq, simulate_interval = 1, frag_len = 250, read_len = 100):
            comp_table = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
            reads_1, reads_2 = [], []
            for i in range(0, len(seq) - frag_len + 1, simulate_interval):
                reads_1.append(seq[i:i+read_len])
                tmp_read_2 = reversed(seq[i+frag_len-read_len:i+frag_len])
                read_2 = ""
                for s in tmp_read_2:
                    if s in comp_table:
                        read_2 += comp_table[s]
                    else:
                        read_2 += s
                reads_2.append(read_2)
            return reads_1, reads_2

        var_ids = clnsigs.keys()
        random.shuffle(var_ids)
        add_count = 0
        for var_id in var_ids:
            assert var_id in vars
            var_chr, var_left, var_type, var_data = vars[var_id]
            assert var_chr in chr_dic
            assert var_id in clnsigs
            var_gene, var_clnsig = clnsigs[var_id]
            print >> sys.stderr, "\t%s %s: %s:%d %s %s (%s)" % \
                (var_gene, var_id, var_chr, var_left, var_type, var_data, var_clnsig)

            chr_seq = chr_dic[var_chr]
            flank_len = 195
            var_seq = chr_seq[var_left-flank_len:var_left]
            if var_type in ["single", "insertion"]:
                var_seq += var_data
                if var_type == "single":
                    var_seq += chr_seq[var_left+1:var_left+1+flank_len]
                else:
                    var_seq += chr_seq[var_left:var_left+flank_len]
            else:
                var_seq += chr_seq[var_left+var_data:var_left+var_data+1+flank_len]
            tmp_reads_1, tmp_reads_2 = simulate_reads(var_seq, simulate_interval)
            HLA_reads_1 += tmp_reads_1
            HLA_reads_2 += tmp_reads_2

            tmp_reads_1, tmp_reads_2 = simulate_reads(chr_seq[var_left-flank_len:var_left+flank_len], simulate_interval)
            HLA_reads_1 += tmp_reads_1
            HLA_reads_2 += tmp_reads_2

            add_count += 1
            if add_count >= 20:
                break

        for test_HLA_names in test_HLA_list:
            gene = test_HLA_names[0].split('*')[0]

            for test_HLA_name in test_HLA_names:
                HLA_seq = HLAs[gene][test_HLA_name]
                tmp_reads_1, tmp_reads_2 = simulate_reads(HLA_seq, simulate_interval)
                HLA_reads_1 += tmp_reads_1
                HLA_reads_2 += tmp_reads_2

        # Write reads into a fasta read file
        def write_reads(reads, idx):
            read_file = open('hla_input_%d.fa' % idx, 'w')
            for read_i in range(len(reads)):
                print >> read_file, ">%d" % (read_i + 1)
                print >> read_file, reads[read_i]
            read_file.close()
        write_reads(HLA_reads_1, 1)
        write_reads(HLA_reads_2, 2)

        print >> sys.stderr, "\n\t\thisat2 graph on %s" % (reference_type)

        hisat_genotype = os.path.join(ex_path, "hisat2_genotype.py")
        genotype_cmd = [hisat_genotype,
                        "-f",
                        "--base-name", "genotype_genome",
                        "-1", "hla_input_1.fa",
                        "-2", "hla_input_2.fa"]
        if verbose:
            print >> sys.stderr, ' '.join(genotype_cmd)
        genotype_proc = subprocess.Popen(genotype_cmd)
        genotype_proc.communicate()



        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='test HLA genotyping')
    parser.add_argument("--reference-type",
                        dest="reference_type",
                        type=str,
                        default="gene",
                        help="Reference type: gene, chromosome, and genome (default: gene)")
    parser.add_argument("--base-name",
                        dest="base_fname",
                        type=str,
                        default="genotype_genome",
                        help="base filename for genotype genome")
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    parser.add_argument("--simulate-interval",
                        dest="simulate_interval",
                        type=int,
                        default=1,
                        help="Reads simulated at every these base pairs (default: 1)")
    parser.add_argument("--num-mismatch",
                        dest="num_mismatch",
                        type=int,
                        default=0,
                        help="Maximum number of mismatches per read alignment to be considered (default: 0)")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')
    parser.add_argument("--daehwan-debug",
                        dest="daehwan_debug",
                        type=str,
                        default="",
                        help="e.g., test_id:10,read_id:10000,basic_test")

    args = parser.parse_args()
    if not args.reference_type in ["gene", "chromosome", "genome"]:
        print >> sys.stderr, "Error: --reference-type (%s) must be one of gene, chromosome, and genome." % (args.reference_type)
        sys.exit(1)
    daehwan_debug = {}
    if args.daehwan_debug != "":
        for item in args.daehwan_debug.split(','):
            if ':' in item:
                key, value = item.split(':')
                daehwan_debug[key] = value
            else:
                daehwan_debug[item] = 1

    random.seed(1)
    genotype_test(args.reference_type,
                  args.base_fname,
                  args.threads,
                  args.simulate_interval,
                  args.num_mismatch,
                  args.verbose,
                  daehwan_debug)
