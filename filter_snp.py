#!/usr/bin/env python3

from sys import stderr, exit
from collections import defaultdict, Counter
from argparse import ArgumentParser, FileType
import itertools
import pprint
import bisect

bVerbose = False
bDebug = False


def read_genome(genome_file, chr_filter=None):
    chr_dic = {}
    # chr_filter = [str(x) for x in list(range(1, 23)) + ['X', 'Y']]

    def check_chr_filter(chr_name, chr_filter):
        if chr_filter is None:
            return True
        if chr_name in chr_filter:
            return True
        return False

    chr_name, sequence = "", ""
    for line in genome_file:
        if line[0] == ">":
            if chr_name and sequence and check_chr_filter(chr_name, chr_filter):
                chr_dic[chr_name] = sequence
            chr_name = line.strip().split()[0][1:]
            sequence = ""
        else:
            sequence += line[:-1]

    if chr_name and sequence and check_chr_filter(chr_name, chr_filter):
        chr_dic[chr_name] = sequence

    if bVerbose:
        lengths = [len(value) for value in chr_dic.values()]
        print('Number of Chromosomes: {}'.format(len(chr_dic)))
        print('Total length: {}'.format(sum(lengths)))

    return chr_dic


def read_snps(snp_file):
    snps = defaultdict(list)
    snps_by_chr = defaultdict(list)

    if not snp_file:
        return snps, snps_by_chr

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

        snps[snpID] = [snpID, type, int(pos), data, chr]

    for snp_id, value in snps.items():
        chr_name = value[4]
        snps_by_chr[chr_name].append(value)

    # sort by pos
    for chrname in snps_by_chr:
        snps_by_chr[chrname].sort(key=lambda x: x[2])

    return snps, snps_by_chr


def snp_check2(chrname, ref, snps, out_fp, stats, flank=10):
    if bVerbose:
        print('Checking....', chrname, flank,  file=stderr)

    pairlist = list()
    pairset = set()

    for i, snp in enumerate(snps):
        if snp[1] != 'single':
            continue

        pos = snp[2]
        snpbase = snp[3]
        refbase = ref[pos]

        if pos < flank:
            if bVerbose:
                print('Skip...', snp, file=stderr)
            continue

        if pos > (len(ref) - flank):
            if bVerbose:
                print('Skip...', snp, file=stderr)
            continue

        flankseq = ref[pos - flank:pos] + refbase + ref[pos + 1: pos + 1 + flank]

        pairlist.append([(refbase, snpbase), flankseq, pos])
        pairset.add((refbase, snpbase, flankseq))

    flankdict = defaultdict(list)

    for x in pairlist:
        pair = x[0]
        flankseq = x[1]

        flankdict[(pair, flankseq)].append(x)

    count = 0
    filtered_pos = dict()
    for x, v in flankdict.items():
        if len(v) > 1:
            count += len(v)

            for i in v:
                filtered_pos[i[2]] = 1

    # do_filter
    for snp in snps:
        pos = snp[2]

        if snp[1] == 'single' and pos in filtered_pos:
            continue

        print('\t'.join([snp[0], snp[1], snp[4], str(snp[2]), str(snp[3])]), file=out_fp)

    count = len(filtered_pos)
    stats['Removed SNPs'] += count
    stats['{} Removed SNPs'.format(chrname)] = count


def filter_SNPs(snp_file, genome_file, out_fname, flank_len):
    stats = defaultdict(lambda: 0)

    genome = read_genome(genome_file)
    snps, snps_by_chr = read_snps(snp_file)

    stats['Number of Chromosomes'] = len(genome)
    stats['Number of SNPs'] = len(snps)

    out_file = open(out_fname, 'w')

    for chrname in genome:
        snp_check2(chrname, genome[chrname], snps_by_chr[chrname], out_file, stats, flank_len)

    out_file.close()

    print(' Chromosomes:', stats['Number of Chromosomes'])
    print('        SNPs:', stats['Number of SNPs'])
    print('Removed SNPs:', stats['Removed SNPs'])

    if bVerbose:
        for chrname in genome:
            print('{} Removed SNPs:'.format(chrname), stats['{} Removed SNPs'.format(chrname)])

    if bDebug:
        pprint.pprint(stats)


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Filter SNPs')

    parser.add_argument('snp_file', metavar='SNP_FILE', type=FileType('r'),
                        help='input SNP file')

    parser.add_argument('genome_file', metavar='GENOME_FILE', type=FileType('r'),
                        help='reference genome file')

    parser.add_argument('out_fname', metavar='OUTPUT_FILE', type=str,
                        help='output SNP filename')

    parser.add_argument('--flank-len', dest='flank_len',
                        type=int,
                        action='store',
                        default=10,
                        help='Flank sequence length')

    parser.add_argument('--debug',
                        dest='bDebug',
                        action='store_true',
                        help='Run in debug mode')

    parser.add_argument('--verbose',
                        dest='bVerbose',
                        action='store_true',
                        help='Show more messages')

    args = parser.parse_args()

    if args.bDebug is not None:
        bDebug = args.bDebug

    if args.bVerbose is not None:
        bVerbose = args.bVerbose

    filter_SNPs(args.snp_file, args.genome_file, args.out_fname, args.flank_len)
