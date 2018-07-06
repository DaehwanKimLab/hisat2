#!/usr/bin/python
import sys, subprocess
import re
from argparse import ArgumentParser, FileType
from collections import defaultdict, Counter

flag_include_N = True

"""
"""
def read_genome(genome_file):
    chr_dic = {}
    chr_name, sequence = "", ""
    for line in genome_file:
        if line.startswith(">"):
            if chr_name and sequence:
                chr_dic[chr_name] = sequence
            chr_name = line.strip().split()[0][1:]
            sequence = ""
        else:
            line = line.strip()
            if not flag_include_N:
                # remove N-bases
                line = line.replace('N', '')

            sequence += line;    

    if chr_name and sequence:
        chr_dic[chr_name] = sequence
    return chr_dic

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
def read_snp(snp_file):
    snps = defaultdict(dict)
    for line in snp_file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        try:
            snpID, type, chr, pos, data = line.split('\t')
        except ValueError:
            continue

        assert type in ["single", "deletion", "insertion"]
        if type == "deletion":
            data = int(data)

        snps[chr][snpID] = [snpID, type, int(pos), data];

    return snps


def indelCount(snp_list, snp_id_list):
    indel = 0

    for snp_id in snp_id_list:
        snp = snp_list[snp_id]

        if snp[1] == 'deletion':
            indel -= int(snp[3])
        elif snp[1] == 'insertion':
            indel += len(snp[3])

    return indel

def applySNPs(snp_list, ref_sqn, snp_id_list, base_pos):

    ref_pos = 0
    read_pos = 0
    read = "" 

    for snp_id in snp_id_list:
        snp = snp_list[snp_id]

        pos = snp[2] - base_pos;

        while ref_pos < pos:
            read += ref_sqn[ref_pos]
            ref_pos += 1


        if snp[1] == 'single':
            read += snp[3]
            ref_pos += 1
        elif snp[1] == 'deletion':
            ref_pos += int(snp[3])
        elif snp[1] == 'insertion':
            read += snp[3]

        #print snp_id, snp_list[snp_id]

    while ref_pos < len(ref_sqn):
        read += ref_sqn[ref_pos]
        ref_pos += 1

    return read


def main(genome_file, rpt_name):
    # load genome sequeuce
    chr_dic = read_genome(genome_file)

    rpt_fa_name = rpt_name + ".rep.fa"
    rpt_info_name = rpt_name + ".rep.info"
    rpt_snp_name = rpt_name + ".rep.snp"

    # load repeat sequence
    fp = open(rpt_fa_name, 'r')
    rpt_dic = read_genome(fp)
    fp.close()

    # load repeat snp
    fp = open(rpt_snp_name, 'r')
    rpt_snps = read_snp(fp)
    fp.close()

    # Validates
    # load repeat info
    fp = open(rpt_info_name, 'r')
    repeat_sequence = ""
    repeat_length = 0
    snp_cnt = 0
    indel = 0
    snp_id_list = []

    for line in fp:
        line = line.strip()

        if line.startswith('>'):
            line = line[1:]
            fields = line.split()

            #print fields

            name, rpt_seq_name, rpt_pos, rpt_len, pos_cnt, snp_cnt = fields[0:6]

            snp_cnt = int(snp_cnt)
            rpt_pos = int(rpt_pos)
            rpt_len = int(rpt_len)

            if snp_cnt > 0:
                snp_id_list = fields[6].split(',')
            else:
                snp_id_list = []

            #print name, snp_cnt, snp_list

            # make repeat_sequence (with snp)
            

            repeat_sequence = rpt_dic[rpt_seq_name][rpt_pos:rpt_pos + rpt_len]
            indel = 0

            if snp_cnt > 0:
                # apply snps
                repeat_sequence = applySNPs(rpt_snps[rpt_seq_name], repeat_sequence, snp_id_list, rpt_pos)
                # in/del count
                indel = indelCount(rpt_snps[rpt_seq_name], snp_id_list)

            #repeat_length = rpt_len + indel
            repeat_length = len(repeat_sequence)

            #print repeat_sequence

        else:
            coords = line.split()
            for coord in coords:
                chr, pos, strand = coord.split(':')
                pos = int(pos)

                # get string
                seq = chr_dic[chr][pos:pos + repeat_length]
                if strand == '-':
                    seq = reverse_complement(seq)

                if seq != repeat_sequence:
                    print 'Mismatch', seq, repeat_sequence, snp_cnt, coord, snp_id_list, repeat_length
                    
    fp.close()


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Validate repeat files')

    parser.add_argument('genome_file',
                        nargs='?',
                        type=FileType('r'),
                        help='input genome file (e.g. genome.fa)')

    parser.add_argument('-r', 
                        dest='rpt_name',
                        type=str,
                        help='Repeat Name')

    args = parser.parse_args()
    if not args.genome_file or not args.rpt_name:
        parser.print_help()
        exit(1)

    main(args.genome_file, args.rpt_name)
