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
            
            chr_name = line[1:-1]
            sequence = ""
        else:
            sequence += line[:-1]

    if chr_name and sequence:
        chr_dic[chr_name] = sequence
    
    return chr_dic


"""
"""
def read_transcript(gtf_file, frag_len):
    genes = defaultdict(list)
    transcripts = {}

    # Parse valid exon lines from the GTF file into a dict by transcript_id
    for line in gtf_file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '#' in line:
            line = line.split('#')[0].strip()
        try:
            chrom, source, feature, left, right, score, \
                strand, frame, values = line.split('\t')
        except ValueError:
            continue
        # Zero-based offset
        left, right = int(left) - 1, int(right) - 1
        if feature != 'exon' or left >= right:
            continue

        values_dict = {}
        for attr in values.split(';')[:-1]:
            attr, _, val = attr.strip().partition(' ')
            values_dict[attr] = val.strip('"')

        if 'gene_id' not in values_dict or \
                'transcript_id' not in values_dict:
            continue

        transcript_id = values_dict['transcript_id']
        if transcript_id not in transcripts:
            transcripts[transcript_id] = [chrom, strand, [[left, right]]]
            genes[values_dict['gene_id']].append(transcript_id)
        else:
            transcripts[transcript_id][2].append([left, right])

    # Sort exons and merge where separating introns are <=5 bps
    for tran, [chr, strand, exons] in transcripts.items():
            exons.sort()
            tmp_exons = [exons[0]]
            for i in range(1, len(exons)):
                if exons[i][0] - tmp_exons[-1][1] <= 5:
                    tmp_exons[-1][1] = exons[i][1]
                else:
                    tmp_exons.append(exons[i])
            transcripts[tran] = [chr, strand, tmp_exons]

    tmp_transcripts = {}
    for tran, [chr, strand, exons] in transcripts.items():
        exon_lens = [e[1] - e[0] + 1 for e in exons]
        transcript_len = sum(exon_lens)
        if transcript_len >= frag_len:
            tmp_transcripts[tran] = [chr, strand, transcript_len, exons]

    transcripts = tmp_transcripts

    return genes, transcripts
    

"""
"""
def read_snp(snp_file):
    snps = []
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
        snps.append([snpID, type, chr, int(pos), data])

    return snps


"""
"""
def sanity_check_input(genome_seq, genes, transcripts, snps, frag_len):
    num_canon_ss, num_ss = 0, 0
    for transcript, [chr, strand, transcript_len, exons] in transcripts.items():
        assert transcript_len >= frag_len
        if len(exons) <= 1:
            continue
        if chr not in genome_seq:
            continue
        chr_seq = genome_seq[chr]
        for i in range(len(exons) - 1):
            left1, right1 = exons[i]
            assert left1 < right1
            left2, right2 = exons[i+1]
            assert left2 < right2
            assert left1 < left2 and right1 < right2
            donor = chr_seq[right1+1:right1+3]
            acceptor = chr_seq[left2-2:left2]
            if strand == "-":
                donor, acceptor = reverse_complement(acceptor), reverse_complement(donor)
            if donor == "GT" and acceptor == "AG":
                num_canon_ss += 1
            num_ss += 1

    print >> sys.stderr, "GT/AG splice sites: {}/{} ({:.2%})".format(num_canon_ss, num_ss, (float(num_canon_ss) / num_ss))

    num_alt_single, num_single = 0, 0
    for snp in snps:
        snpID, type, chr, pos, data = snp
        if type != "single":
            continue
        if chr not in genome_seq:
            continue
        chr_seq = genome_seq[chr]
        assert pos < len(chr_seq)
        if chr_seq[pos] != data:
            num_alt_single += 1
        num_single += 1

    print >> sys.stderr, "Alternative bases: {}/{} ({:.2%})".format(num_alt_single, num_single, (float(num_alt_single) / num_single))


"""
"""
def generate_expr_profile(expr_profile_type, num_transcripts = 10000):
    # Modelling and simulating generic RNA-Seq experiments with the flux simulator
    # http://nar.oxfordjournals.org/content/suppl/2012/06/29/gks666.DC1/nar-02667-n-2011-File002.pdf
    def calc_expr(x, a):
        x, a, b = float(x), 9500.0, 9500.0
        k = -0.6
        return (x**k) * math.exp(x/a * (x/b)**2)
    
    expr_profile = [0.0] * num_transcripts
    for i in range(len(expr_profile)):
        if expr_profile_type == "flux":
            expr_profile[i] = calc_expr(i + 1, num_transcripts)
        elif expr_profile_type == "constant":
            expr_profile[i] = 1.0
        else:
            assert False

    expr_sum = sum(expr_profile)
    expr_profile = [expr_profile[i] / expr_sum for i in range(len(expr_profile))]
    assert abs(sum(expr_profile) - 1.0) < 0.001
    return expr_profile


"""
"""
cigar_re = re.compile('\d+\w')
def mapRepOk(genome_seq, read_seq, chr, pos, cigar, XM, NM, MD):
    assert chr in genome_seq
    chr_seq = genome_seq[chr]
    assert pos < len(chr_seq)

    tXM, tNM, tMD = 0, 0, ""
    cigars = cigar_re.findall(cigar)
    cigars = [[int(cigars[i][:-1]), cigars[i][-1]] for i in range(len(cigars))]
    read_pos, ref_pos = 0, pos
    for i in range(len(cigars)):
        cigar_len, cigar_op = cigars[i]
        if cigar_op == "M":
            read_partial_seq = read_seq[read_pos:read_pos+cigar_len]
            ref_partial_seq = chr_seq[ref_pos:ref_pos+cigar_len]
            assert len(read_partial_seq) == len(ref_partial_seq)
            for j in range(len(read_partial_seq)):
                if read_partial_seq[j] != ref_partial_seq[j]:
                    tXM += 1
                    tNM += 1

        # if cigar_op == "N":
        #    left, right = ref_pos - 1, ref_pos + cigar_len

        if cigar_op in "MND":
            ref_pos += cigar_len
        if cigar_op in "MI":
            read_pos += cigar_len

    assert tXM == XM and tNM == NM

    
"""
"""
def simulate_reads(genome_file, gtf_file, snp_file, base_fname, \
                       rna, paired_end, read_len, frag_len, \
                       num_frag, expr_profile_type, random_seed, \
                       sanity_check, verbose):
    random.seed(random_seed)
    if read_len > frag_len:
        frag_len = read_len

    genome_seq = read_genome(genome_file)
    genes, transcripts = read_transcript(gtf_file, frag_len)
    snps = read_snp(snp_file)

    if sanity_check:
        sanity_check_input(genome_seq, genes, transcripts, snps, frag_len)

    num_transcripts = min(len(transcripts), 10000)
    expr_profile = generate_expr_profile(expr_profile_type, num_transcripts)
    expr_profile = [int(expr_profile[i] * num_frag) for i in range(len(expr_profile))]

    assert num_frag >= sum(expr_profile)
    expr_profile[0] += (num_frag - sum(expr_profile))
    assert num_frag == sum(expr_profile)

    transcript_ids = transcripts.keys()
    random.shuffle(transcript_ids)

    sam_file = open(base_fname + ".sam", "w")
    read_file = open(base_fname + "_1.fa", "w")
    if paired_end:
        read2_file = open(base_fname + "_2.fa", "w")

    assert len(transcript_ids) >= len(expr_profile)
    cur_read_id = 1
    for t in range(len(expr_profile)):
        transcript_id = transcript_ids[t]
        chr, strand, transcript_len, exons = transcripts[transcript_id]
        t_num_frags = expr_profile[t]
        t_seq = ""
        assert chr in genome_seq
        chr_seq = genome_seq[chr]
        for e in exons:
            assert e[0] < e[1]
            t_seq += chr_seq[e[0]:e[1]+1]

        assert len(t_seq) == transcript_len

        for f in range(t_num_frags):
            frag_pos = random.randint(0, transcript_len - frag_len)
            assert frag_pos + frag_len <= transcript_len

            def getMapInfo(exons, frag_pos, read_len):
                pos, cigar = exons[0][0], ""
                e_pos = 0
                prev_e = None
                for e_i in range(len(exons)):
                    e = exons[e_i]
                    if prev_e:
                        i_len = e[0] - prev_e[1] - 1
                        pos += i_len
                    e_len = e[1] - e[0] + 1
                    if e_len <= frag_pos:
                        frag_pos -= e_len
                        pos += e_len
                    else:
                        pos += frag_pos
                        e_pos = frag_pos
                        break                        
                    prev_e = e
                    
                assert e_i < len(exons)
                assert e_pos < exons[e_i][1] - exons[e_i][0] + 1
                prev_e = None
                for e in exons[e_i:]:
                    if prev_e:
                        i_len = e[0] - prev_e[1] - 1
                        cigar += ("{}N".format(i_len))
                    e_len = e[1] - e[0] + 1
                    if e_len < e_pos + read_len:
                        m_len = e_len - e_pos
                        e_pos = 0
                        read_len -= m_len
                        cigar += ("{}M".format(m_len))
                    else:
                        cigar += ("{}M".format(read_len))
                        read_len = 0
                        break
                    prev_e = e

                return pos, cigar
            
            flag, flag2 = 99, 163  # 83, 147
            pos, cigar = getMapInfo(exons, frag_pos, read_len)
            pos2, cigar2 = getMapInfo(exons, frag_pos+frag_len-read_len, read_len)
            read_seq, read2_seq  = t_seq[frag_pos:frag_pos+read_len], t_seq[frag_pos+frag_len-read_len:frag_pos+frag_len]
            XM, XM2 = 0, 0
            NM, NM2 = 0, 0
            MD, MD2 = "100", "100"

            swapped = False
            if paired_end:
                if random.randint(0, 1) == 1:
                    swapped = True
                if swapped:
                    flag, flag2 = flag2 - 16, flag - 16
                    pos, pos2 = pos2, pos
                    cigar, cigar2 = cigar2, cigar
                    read_seq, read2_seq = read2_seq, read_seq
                    XM, XM2 = XM2, XM
                    NM, NM2 = NM2, NM
                    MD, MD2 = MD2, MD
                
            if sanity_check:
                mapRepOk(genome_seq, read_seq, chr, pos, cigar, XM, NM, MD)
                mapRepOk(genome_seq, read2_seq, chr, pos2, cigar2, XM2, NM2, MD2)

            print >> read_file, ">{}".format(cur_read_id)
            if swapped:
                print >> read_file, reverse_complement(read_seq)
            else:
                print >> read_file, read_seq
            print >> sam_file, "{}\t{}\t{}\t{}\t255\t{}\t{}\t{}\t0\t{}\t*\tXM:i:{}\tNM:i:{}\tMD:Z:{}\tTI:Z:{}".format(cur_read_id, flag, chr, pos, cigar, chr, pos2, read_seq, XM, NM, MD, transcript_id)
            if paired_end:
                print >> read2_file, ">{}".format(cur_read_id)
                if swapped:
                    print >> read2_file, read2_seq
                else:
                    print >> read2_file, reverse_complement(read2_seq)
                print >> sam_file, "{}\t{}\t{}\t{}\t255\t{}\t{}\t{}\t0\t{}\t*\tXM:i:{}\tNM:i:{}\tMD:Z:{}\tTI:Z:{}".format(cur_read_id, flag2, chr, pos2, cigar2, chr, pos, read2_seq, XM2, NM2, MD2, transcript_id)

            cur_read_id += 1

            
    sam_file.close()
    read_file.close()
    if paired_end:
        read2_file.close()


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Simulate reads from GENOME (fasta) and GTF files')
    parser.add_argument('genome_file',
                        nargs='?',
                        type=FileType('r'),
                        help='input GENOME file')
    parser.add_argument('gtf_file',
                        nargs='?',
                        type=FileType('r'),
                        help='input GTF file')
    parser.add_argument('snp_file',
                        nargs='?',
                        type=FileType('r'),
                        help='input SNP file')
    parser.add_argument('base_fname',
                        nargs='?',
                        type=str,
                        help='output base filename')
    parser.add_argument('-d', '--dna',
                        dest='rna',
                        action='store_false',
                        default=True,
                        help='DNA-seq reads (default: RNA-seq reads)')
    parser.add_argument('--single-end',
                        dest='paired_end',
                        action='store_false',
                        default=True,
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
    parser.add_argument('-e', '--expr-profile',
                        dest='expr_profile',
                        action='store',
                        type=str,
                        default='flux',
                        help='expression profile: flux or constant (default: flux)')
    parser.add_argument('--random-seed',
                        dest='random_seed',
                        action='store',
                        type=int,
                        default=0,
                        help='random seeding value (default: 0)')
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
                        version='%(prog)s 2.0.0-alpha')
    args = parser.parse_args()
    if not args.gtf_file:
        parser.print_help()
        exit(1)
    simulate_reads(args.genome_file, args.gtf_file, args.snp_file, args.base_fname, \
                       args.rna, args.paired_end, args.read_len, args.frag_len, \
                       args.num_frag, args.expr_profile, args.random_seed, \
                       args.sanity_check, args.verbose)
