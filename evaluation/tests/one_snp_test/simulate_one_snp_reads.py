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

import sys, math, random, re
from collections import defaultdict, Counter
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
            chr_name = line.strip().split()[0][1:]
            sequence = ""
        else:
            sequence += line[:-1]

    if chr_name and sequence:
        chr_dic[chr_name] = sequence

    return chr_dic


"""
"""
def read_snp(snp_file):
    snps = defaultdict(list)
    for line in snp_file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        try:
            snpID, type, chr, pos, data = line.split('\t')
            chr = chr.split()[0]
        except ValueError:
            continue

        assert type in ["single", "deletion", "insertion"]
        if type == "deletion":
            data = int(data)
        snps[chr].append([snpID, type, int(pos), data])

    return snps


"""
"""
def getSamAlignment(chr_seq, read_len, snp):
    snp_id, snp_type, snp_pos, snp_data = snp

    # Define MD, XM, NM, Zs, read_seq
    read_seq, MD, Zs  = "", "", ""
    left_read_len = read_len / 2
    pos = snp_pos - left_read_len
    assert pos >= 0
    if snp_type == "single":
        cigar_str = "100M"
        read_seq = chr_seq[pos:pos+left_read_len] + snp_data + chr_seq[pos+left_read_len+1:pos+read_len]
        MD = "%d%s%d" % (left_read_len, chr_seq[pos+left_read_len], read_len - left_read_len - 1)
        Zs = "%d|S|%s" % (left_read_len, snp_id)
    elif snp_type == "deletion":
        del_len = int(snp_data)
        cigar_str = "%dM%dD%dM" % (left_read_len, del_len, read_len - left_read_len)
        read_seq = chr_seq[pos:pos+left_read_len] + chr_seq[pos+left_read_len+del_len:pos+read_len+del_len]
        MD = "%d^%s%d" % (left_read_len, chr_seq[pos+left_read_len:pos+left_read_len+del_len], read_len - left_read_len)
        Zs = "%d|D|%s" % (left_read_len, snp_id)
    else:
        assert snp_type == "insertion"
        ins_len = len(snp_data)
        assert ins_len < read_len
        cigar_str = "%dM%dI%dM" % (left_read_len, ins_len, read_len - left_read_len - ins_len)
        read_seq = chr_seq[pos:pos+left_read_len] + snp_data + chr_seq[pos+left_read_len:pos+read_len-ins_len]
        MD = "%d" % (read_len - ins_len)
        Zs = "%d|I|%s" % (left_read_len, snp_id)

    if len(read_seq) != read_len:
        print >> sys.stderr, "read length differs:", len(read_seq), "vs.", read_len
        print >> sys.stderr, pos, cigar_str, MD, Zs
        assert False

    ref_read_seq = chr_seq[pos:pos+read_len]
    return pos, cigar_str, MD, Zs, read_seq, ref_read_seq


"""
"""
cigar_re = re.compile('\d+\w')
def samRepOk(genome_seq, read_seq, chr, pos, cigar, MD, Zs):
    assert chr in genome_seq
    chr_seq = genome_seq[chr]
    assert pos < len(chr_seq)

    # Calculate XM and NM based on Cigar and Zs
    cigars = cigar_re.findall(cigar)
    cigars = [[int(cigars[i][:-1]), cigars[i][-1]] for i in range(len(cigars))]
    ref_pos, read_pos = pos, 0
    ann_ref_seq, ann_ref_rel, ann_read_seq, ann_read_rel = [], [], [], []
    for i in range(len(cigars)):
        cigar_len, cigar_op = cigars[i]
        if cigar_op == "M":
            partial_ref_seq = chr_seq[ref_pos:ref_pos+cigar_len]
            partial_read_seq = read_seq[read_pos:read_pos+cigar_len]
            assert len(partial_ref_seq) == len(partial_read_seq)
            ann_ref_seq += list(partial_ref_seq)
            ann_read_seq += list(partial_read_seq)
            for j in range(len(partial_ref_seq)):
                if partial_ref_seq[j] == partial_read_seq[j]:
                    ann_ref_rel.append("=")
                    ann_read_rel.append("=")
                else:
                    ann_ref_rel.append("X")
                    ann_read_rel.append("X")
            ref_pos += cigar_len
            read_pos += cigar_len
        elif cigar_op == "D":
            partial_ref_seq = chr_seq[ref_pos:ref_pos+cigar_len]
            ann_ref_rel += list(partial_ref_seq)
            ann_ref_seq += list(partial_ref_seq)
            ann_read_rel += (["-"] * cigar_len)
            ann_read_seq += (["-"] * cigar_len)
            ref_pos += cigar_len
        elif cigar_op == "I":
            partial_read_seq = read_seq[read_pos:read_pos+cigar_len]
            ann_ref_rel += (["-"] * cigar_len)
            ann_ref_seq += (["-"] * cigar_len)
            ann_read_rel += list(partial_read_seq)
            ann_read_seq += list(partial_read_seq) 
            read_pos += cigar_len
        elif cigar_op == "N":
            ref_pos += cigar_len
        else:
            assert False
    
    assert len(ann_ref_seq) == len(ann_read_seq)
    assert len(ann_ref_seq) == len(ann_ref_rel)
    assert len(ann_ref_seq) == len(ann_read_rel)
    ann_Zs_seq = ["0" for i in range(len(ann_ref_seq))]

    Zss, Zs_i, snp_pos_add = [], 0, 0
    if Zs != "":
        Zss = Zs.split(',')
        Zss = [zs.split('|') for zs in Zss]

    ann_read_pos = 0
    for zs in Zss:
        zs_pos, zs_type, zs_id = zs
        zs_pos = int(zs_pos)
        for i in range(zs_pos):
            while ann_read_rel[ann_read_pos] == '-':
                ann_read_pos += 1
            ann_read_pos += 1
        if zs_type == "S":
            ann_Zs_seq[ann_read_pos] = "1"
            ann_read_pos += 1
        elif zs_type == "D":
            while ann_read_rel[ann_read_pos] == '-':
                ann_Zs_seq[ann_read_pos] = "1"
                ann_read_pos += 1
        elif zs_type == "I":
            while ann_ref_rel[ann_read_pos] == '-':
                ann_Zs_seq[ann_read_pos] = "1"
                ann_read_pos += 1
        else:
            assert False

    tMD, tXM, tNM = "", 0, 0
    match_len = 0
    i = 0
    while i < len(ann_ref_seq):
        if ann_ref_rel[i] == "=":
            assert ann_read_rel[i] == "="
            match_len += 1
            i += 1
            continue
        assert ann_read_rel[i] != "="
        if ann_ref_rel[i] == "X" and ann_read_rel[i] == "X":
            if match_len > 0:
                tMD += ("{}".format(match_len))
                match_len = 0
            tMD += ann_ref_seq[i]
            if ann_Zs_seq[i] == "0":
                tXM += 1
                tNM += 1
            i += 1
        else:
            assert ann_ref_rel[i] == "-" or ann_read_rel[i] == "-"
            if ann_ref_rel[i] == '-':
                while ann_ref_rel[i] == '-':
                    if ann_Zs_seq[i] == "0":
                        tNM += 1
                    i += 1
            else:
                assert ann_read_rel[i] == '-'
                del_seq = ""
                while  ann_read_rel[i] == '-':
                    del_seq += ann_ref_seq[i]
                    if ann_Zs_seq[i] == "0":
                        tNM += 1
                    i += 1
                if match_len > 0:
                    tMD += ("{}".format(match_len))
                    match_len = 0
                tMD += ("^{}".format(del_seq))

    if match_len > 0:
        tMD += ("{}".format(match_len))

    if tMD != MD:
        print >> sys.stderr, chr, pos, cigar, MD, Zs
        print >> sys.stderr, tMD
        assert False
        
        
"""
"""
def simulate_reads(genome_file, snp_file, base_fname, \
                       paired_end, read_len, frag_len, \
                       num_frag, sanity_check, verbose):
    if read_len > frag_len:
        frag_len = read_len

    genome_seq = read_genome(genome_file)
    snps = read_snp(snp_file)
    chr_ids = genome_seq.keys()

    sam_file = open(base_fname + ".sam", "w")

    # Write SAM header
    print >> sam_file, "@HD\tVN:1.0\tSO:unsorted"
    for chr in genome_seq.keys():
        print >> sam_file, "@SQ\tSN:%s\tLN:%d" % (chr, len(genome_seq[chr]))
    
    read_file = open(base_fname + "_snp_1.fa", "w")
    ref_read_file = open(base_fname + "_ref_1.fa", "w")
    if paired_end:
        read2_file = open(base_fname + "_snp_2.fa", "w")
        ref_read2_file = open(base_fname + "_ref_2.fa", "w")

    cur_read_id = 1
    for chr in snps:
        assert chr in genome_seq
        chr_seq = genome_seq[chr]
        chr_len = len(chr_seq)
        if chr in snps:
            chr_snps = snps[chr]
        else:
            chr_snps = []
        for snp in chr_snps:
            # SAM specification (v1.4)
            # http://samtools.sourceforge.net/
            flag, flag2 = 99, 163  # 83, 147
            pos, cigar_str, MD, Zs, read_seq, ref_read_seq = getSamAlignment(chr_seq, read_len, snp)
            # pos2, cigar2_str, MD2, Zs2, read2_seq = getSamAlignment(chr_seq, frag_pos+frag_len-read_len, read_len, snp)
            if sanity_check:
                samRepOk(genome_seq, read_seq, chr, pos, cigar_str, MD, Zs)
                # samRepOk(genome_seq, read2_seq, chr, pos2, cigar2_str, MD2, Zs2)

            if Zs != "":
                Zs = ("\tZs:Z:{}".format(Zs))
            # if Zs2 != "":
            #    Zs2 = ("\tZs:Z:{}".format(Zs2))

            read_id_str = "{}_{}_{}_{}".format(cur_read_id, chr, pos, cigar_str)
            print >> read_file, ">{}".format(read_id_str)
            print >> read_file, read_seq
            print >> sam_file, "{}\t{}\t{}\t{}\t255\t{}\t{}\t{}\t0\t{}\t*\tXM:i:0\tNM:i:0\tMD:Z:{}{}".format(read_id_str, flag, chr, pos + 1, cigar_str, chr, pos + 1, read_seq, MD, Zs)

            print >> ref_read_file, ">{}_{}_{}_100M".format(cur_read_id, chr, pos)
            print >> ref_read_file, ref_read_seq
            """
            if paired_end:
                print >> read2_file, ">{}".format(cur_read_id)
                print >> read2_file, reverse_complement(read2_seq)
                print >> sam_file, "{}\t{}\t{}\t{}\t255\t{}\t{}\t{}\t0\t{}\t*\tXM:i:0\tNM:i:0\tMD:Z:{}{}".format(cur_read_id, flag2, chr, pos2 + 1, cigar2_str, chr, pos + 1, read2_seq, MD2, Zs2)
            """

            cur_read_id += 1
            
    sam_file.close()
    read_file.close()
    ref_read_file.close()
    if paired_end:
        read2_file.close()
        ref_read2_file.close()


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Simulate reads from GENOME (fasta) and GTF files')
    parser.add_argument('genome_file',
                        nargs='?',
                        type=FileType('r'),
                        help='input GENOME file')
    parser.add_argument('snp_file',
                        nargs='?',
                        type=FileType('r'),
                        help='input SNP file')
    parser.add_argument('base_fname',
                        nargs='?',
                        type=str,
                        help='output base filename')
    parser.add_argument('--paired-end',
                        dest='paired_end',
                        action='store_true',
                        help='single-end reads (default: paired-end reads)')
    parser.add_argument('-r', '--read-length',
                        dest='read_len',
                        action='store',
                        type=int,
                        default=100,
                        help='read length (default: 100)')
    parser.add_argument('-f', '--fragment-length',
                        dest='frag_len',
                        action='store',
                        type=int,
                        default=250,
                        help='fragment length (default: 250)')
    parser.add_argument('-n', '--num-fragment',
                        dest='num_frag',
                        action='store',
                        type=int,
                        default=1000000,
                        help='number of fragments (default: 1000000)')
    parser.add_argument('--sanity-check',
                        dest='sanity_check',
                        action='store_true',
                        help='sanity check')
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')
    parser.add_argument('--version', 
                        action='version',
                        version='%(prog)s 2.1.0')
    args = parser.parse_args()
    if not args.genome_file or not args.snp_file:
        parser.print_help()
        exit(1)
    simulate_reads(args.genome_file, args.snp_file, args.base_fname, \
                       args.paired_end, args.read_len, args.frag_len, \
                       args.num_frag, args.sanity_check, args.verbose)
