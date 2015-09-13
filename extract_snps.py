#!/usr/bin/env python

#
# Copyright 2015, Daehwan Kim <infphilo@gmail.com>
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


import sys, re
from collections import defaultdict as dd, Counter
from argparse import ArgumentParser, FileType

"""
"""
def reverse_complement(seq):
    result = ""
    for nt in seq:
        base = nt
        if nt == 'A':
            base = 'T'
        elif nt == 'a':
            base = 't'
        elif nt == 'C':
            base = 'G'
        elif nt == 'c':
            base = 'g'
        elif nt == 'G':
            base = 'C'
        elif nt == 'g':
            base = 'c'
        elif nt == 'T':
            base = 'A'
        elif nt == 't':
            base = 'a'

        result = base + result

    return result


"""
"""
def read_genome(genome_file):
    chr_dic = {}

    chr_name, sequence = "", ""
    for line in genome_file:
        if line[0] == ">":
            if chr_name and sequence:
                chr_dic[chr_name] = sequence

            chr_name = line[1:-1]
            sequence = ""
        else:
            sequence += line[:-1]

    if chr_name and sequence:
        chr_dic[chr_name] = sequence

    return chr_dic


"""
"""
def extract_snps(genome_file, snp_file, verbose = False, testset = False):
    # load genomic sequences
    chr_dic = read_genome(genome_file)

    if testset:
        testset_bname = genome_file.name.split(".")[0]
        ref_testset_file = open(testset_bname + ".ref.testset.fa", "w")
        alt_testset_file = open(testset_bname + ".alt.testset.fa", "w")

    # load SNPs
    snp_list = []
    for line in snp_file:
        if not line or line.startswith('#'):
            continue

        line = line.strip()
        try:
            """
            id, chr, start, end, rs_id, score, strand, refNCBI, refUCSC, observed, molType, classType, valid, \
                avHet, avHetSE, func, locType, weight, exceptions, submitterCount, submitters, \
                alleleFreqCount, alleles, alleleNs, alleleFreqs, bitfields = line.split("\t")
            """
            id, chr, start, end, rs_id, score, strand, refNCBI, refUCSC, observed, molType, classType = line.split('\t')[:12]
        except ValueError:
            continue

        start, end = int(start), int(end)
        score = int(score)

        if molType != "genomic":
            continue

        if classType not in ["single", "deletion", "insertion"]:
            continue

        if classType == "single":
            if start + 1 != end:
                continue
        elif classType == "deletion":
            assert start < end
        else:
            assert classType == "insertion"
            if start != end:
                continue
            
        if chr not in chr_dic:
            continue
        chr_seq = chr_dic[chr]
        chr_len = len(chr_seq)

        if start >= len(chr_seq):
            continue

        # daehwan - for debugging purposes
        # """
        if len(snp_list) > 0:
            _, _, last_chr, last_start, _, _  = snp_list[-1]
            if chr == last_chr and abs(start - last_start) <= 20:
                continue
        # """

        observed = observed.upper()
        allele_list = observed.split("/")
        # Reverse complement alleles if strand is negative
        if strand == "-":
            tmp_allele_list = []
            for allele in allele_list:
                tmp_allele_list.append(reverse_complement(allele))
            allele_list = tmp_allele_list
            
        if classType == "single":
            ref_base = chr_seq[start].upper()
            if ref_base not in allele_list:
                continue
            for allele in allele_list:
                if allele not in "ACGT" or len(allele) != 1:
                    continue
                if allele == ref_base:
                    continue
                snp_list.append([rs_id, classType, chr, start, end, allele])

                if testset:
                    ref_seq = chr_seq[start-50:start+50]
                    alt_seq = chr_seq[start-50:start] + allele + chr_seq[start+1:start+50]
                    print >> ref_testset_file, ">%s_single_%d" % (rs_id, start - 50)
                    print >> ref_testset_file, ref_seq
                    print >> alt_testset_file, ">%s_single_%d_%s" % (rs_id, start - 50, ref_seq)
                    print >> alt_testset_file, alt_seq
                
        elif classType == "deletion":
            snp_list.append([rs_id, classType, chr, start, end, ""])
            delLen = end -start
            if testset and delLen > 0 and delLen <= 10:
                ref_seq = chr_seq[start-50:start+50]
                alt_seq = chr_seq[start-50:start] + chr_seq[start+delLen:start+50+delLen]
                print >> ref_testset_file, ">%s_deletion_%d" % (rs_id, start - 50)
                print >> ref_testset_file, ref_seq
                print >> alt_testset_file, ">%s_deletion_%d_%s" % (rs_id, start - 50, ref_seq)
                print >> alt_testset_file, alt_seq
        else:
            assert classType == "insertion"
            for allele in allele_list:
                if allele == "-" or len(allele) <= 0:
                    continue
                if re.match('^[ACGT]+$', allele):
                    snp_list.append([rs_id, classType, chr, start, end, allele])
                    insLen = len(allele)
                    if testset and insLen > 0 and insLen <= 10:
                        ref_seq = chr_seq[start-50:start+50]
                        alt_seq = chr_seq[start-50:start] + allele + chr_seq[start:start+50-insLen]
                        print >> ref_testset_file, ">%s_insertion_%d" % (rs_id, start - 50)
                        print >> ref_testset_file, ref_seq
                        print >> alt_testset_file, ">%s_insertion_%d_%s" % (rs_id, start - 50, ref_seq)
                        print >> alt_testset_file, alt_seq

    if testset:
        ref_testset_file.close()
        alt_testset_file.close()

                    
    # Sort SNPs (snp_list) according to chromosomes, genomic coordinates, types, and alleles
    def snp_cmp(a, b):
        _, classType1, chr1, start1, end1, allele1 = a
        _, classType2, chr2, start2, end2, allele2 = b
        if chr1 < chr2:
            return -1
        elif chr2 < chr1:
            return 1
        
        if start1 < start2:
            return -1
        elif start2 < start1:
            return 1
        
        if end1 < end2:
            return -1
        elif end2 < end1:
            return 1
        
        if classType1 != classType2:
            if classType1 == "single":
                return -1
            elif classType2 == "single":
                return 1
            elif classType1 == "deletion":
                return -1
            else:
                assert classType1 == "insertion" and classType2 == "deletion"
                return 1

        if len(allele1) != len(allele2):
            if len(allele1) < len(allele2):
                return -1
            else:
                return 1
            
        if allele1 < allele2:
            return -1
        elif allele2 < allele1:
            return 1
        else:
            return 0

    snp_list = sorted(snp_list, cmp=snp_cmp)
    assert len(snp_list) > 0
    tmp_snp_list = snp_list[:1]
    for snp in snp_list[1:]:
        if cmp(tmp_snp_list[-1], snp) != 0:
            tmp_snp_list.append(snp)
    snp_list = tmp_snp_list


    print >> sys.stderr, "Number of SNPs: %d" % (len(snp_list))
    for snp in snp_list:
        id, classType, chr, start, end, allele = snp
        if classType in ["single", "insertion"]:
            out = "%s\t%s\t%s\t%d\t%s" % (id, classType, chr, start, allele)
        elif classType in ["deletion"]:
            out = "%s\t%s\t%s\t%d\t%d" % (id, classType, chr, start, end - start)
        print out
                     

if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract SNPs from a SNP file')
    parser.add_argument('genome_file',
        nargs='?',
        type=FileType('r'),
        help='input genome file')
    parser.add_argument('snp_file',
        nargs='?',
        type=FileType('r'),
        help='input snp file')
    parser.add_argument('-v', '--verbose',
        dest='verbose',
        action='store_true',
        help='also print some statistics to stderr')
    parser.add_argument('--testset',
        dest='testset',
        action='store_true',
        help='print test reads')

    args = parser.parse_args()
    if not args.snp_file:
        parser.print_help()
        exit(1)
    extract_snps(args.genome_file, args.snp_file, args.verbose, args.testset)
