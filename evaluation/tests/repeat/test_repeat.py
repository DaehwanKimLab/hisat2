#!/usr/bin/env python
import sys, os, subprocess, random
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
def read_genome(genome_filename):
    chr_dic = {}
    genome_file = open(genome_filename, "r")

    chr_name, sequence = "", ""
    for line in genome_file:
        if line[0] == ">":
            if chr_name and sequence:
                chr_dic[chr_name] = sequence

            chr_name = line[1:-1].split()[0]
            sequence = ""
        else:
            sequence += line[:-1]

    if chr_name and sequence:
        chr_dic[chr_name] = sequence

    genome_file.close()

    print >> sys.stderr, "genome is loaded"
    
    return chr_dic


"""
"""
def generate_random_seq(seq_len):
    assert seq_len > 0
    random_seq = ""
    for i in xrange(seq_len):
        random_seq += "ACGT"[random.randint(0, 3)]
    return random_seq


"""
"""
def test_repeat(verbose):
    random.seed(1)
    
    backbone_seq = generate_random_seq(500)
    mm_seq = backbone_seq[:]
    mm_seq = mm_seq[:50] + ("A" if mm_seq[50] != "A" else "C") + mm_seq[51:]
    mm_seq2 = backbone_seq[:]
    mm_seq2 = mm_seq2[:450] + ("A" if mm_seq2[450] != "A" else "C") + mm_seq2[451:]
    del_seq = backbone_seq[:]
    del_seq = del_seq[:50] + del_seq[52:150] + del_seq[152:]
    del_seq2 = backbone_seq[:]
    del_seq2 = del_seq2[:350] + del_seq2[352:450] + del_seq2[452:]
    indel_seq = backbone_seq[:]
    indel_seq = indel_seq[:30] + indel_seq[32:130] + "AAA" + indel_seq[130:]
    indel_seq2 = backbone_seq[:]
    indel_seq2 = indel_seq2[:30] + "AAA" + indel_seq2[30:130] + indel_seq2[132:]
        
    seqs = [
        # dummy_seq,
        ["bb01", backbone_seq],
        ["bb02", backbone_seq],
        ["bb03", backbone_seq],
        ["bb04", backbone_seq],
        ["bb05", backbone_seq],
        ["mm01", mm_seq],
        ["mm02", mm_seq],
        ["dd01", del_seq],
        ["dd02", del_seq],
        ["dd03", del_seq2],
        ["dd04", del_seq2],
        ["id01", indel_seq],
        ["id02", indel_seq],
        ["id03", indel_seq],
        ["id04", indel_seq],
        ["id05", indel_seq],
        ["id06", indel_seq],
        ["id07", indel_seq2],
    ]
    
    for id, seq in seqs:
        print ">%s" % id
        print generate_random_seq(20)
        print seq
        print generate_random_seq(20)


"""
"""
if __name__ == "__main__":
    parser = ArgumentParser(
        description='')
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()
    test_repeat(args.verbose)

    
