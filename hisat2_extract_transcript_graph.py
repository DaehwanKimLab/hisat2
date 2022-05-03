#!/usr/bin/env python3
#
# Copyright 2020, Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>
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

from sys import stderr, exit
from collections import defaultdict, Counter
from argparse import ArgumentParser, FileType
import itertools
import pprint
import bisect

bDebug = False
bVerbose = False
bUniqueSNP = False
bChrTome = False


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


def write_fasta(fp, ref_id, seq, width=60):
    print('>{}'.format(ref_id), file=fp)

    for i in range(0, len(seq), width):
        print(seq[i:i + width], file=fp)

    return


def load_transcript(genome_seq, gtf_file):
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
            if bVerbose:
                print("Warning: Can't parse line:", line)
            continue

        if chrom not in genome_seq:
            if bVerbose:
                print('chr {} is not in genome: {}'.format(chrom, line))
            continue

        # Zero-based offset
        left, right = int(left) - 1, int(right) - 1
        if left >= right:
            if bDebug:
                print('Wrong position', line)
            continue

        values_dict = {}
        for attr in values.split(';')[:-1]:
            attr, _, val = attr.strip().partition(' ')
            values_dict[attr] = val.strip('"')

        if feature == 'gene':
            gene_id = values_dict['gene_id']
            if gene_id in genes:
                print('Duplicated gene:', gene_id)
                continue

            gene_biotype = values_dict['gene_biotype']
            gene_name = values_dict['gene_name']
            genes[gene_id] = [list(), gene_biotype, gene_name, chrom, list()]

        elif feature == 'exon':
            transcript_id = values_dict['transcript_id']
            gene_id = values_dict['gene_id']

            if gene_id not in genes:
                gene_name = values_dict['gene_name']
                gene_biotype = ''
                genes[gene_id] = [list(), gene_biotype, gene_name, chrom, list()]

            if transcript_id not in transcripts:
                transcripts[transcript_id] = [chrom, strand, [[left, right]], gene_id]
                genes[gene_id][0].append(transcript_id)
            else:
                transcripts[transcript_id][2].append([left, right])
        else:
            continue

    return genes, transcripts


def read_transcript(genome_seq, gtf_file, min_transcript_len=0):
    genes, transcripts = load_transcript(genome_seq, gtf_file)

    merged_exon = 0
    min_transcript_len_filtered = 0

    # Sort exons and merge where separating introns are <=5 bps
    for tran, [chrom, strand, exons, gene_id] in transcripts.items():
        exons.sort()
        tmp_exons = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - tmp_exons[-1][1] <= 5:
                merged_exon += 1
                tmp_exons[-1][1] = exons[i][1]
            else:
                tmp_exons.append(exons[i])
        transcripts[tran] = [chrom, strand, tmp_exons, gene_id]

    tmp_transcripts = {}
    for tran, [chrom, strand, exons, gene_id] in transcripts.items():
        exon_lens = [e[1] - e[0] + 1 for e in exons]
        transcript_len = sum(exon_lens)
        if transcript_len >= min_transcript_len:
            tmp_transcripts[tran] = [chrom, strand, transcript_len, exons, gene_id, list()]
        else:
            min_transcript_len_filtered += 1

    transcripts = tmp_transcripts

    if bVerbose:
        num_exons = [len(values[3]) for values in transcripts.values()]
        total_exons = sum(num_exons)
        avg_exons = total_exons / len(num_exons)

        num_protein_coding = 0
        for gene in genes.values():
            if gene[1] == 'protein_coding':
                num_protein_coding += 1

        assert len(num_exons) == len(transcripts)
        print("Number of genes: {}, protein_coding_genes: {}".format(len(genes), num_protein_coding))
        print("Number of transcripts: {}".format(len(transcripts)))
        print("merged_exon: {}, min_transcript_len_filtered: {}".format(merged_exon, min_transcript_len_filtered))
        print("total exons: {}, average exons: {:.2f}".format(total_exons, avg_exons))

    return genes, transcripts


def read_snps(snp_file):
    print('Read SNPs')
    snps = defaultdict(list)

    if not snp_file:
        return snps

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

        snps[chr].append([snpID, type, int(pos), data, chr])

    return snps


def read_haplotypes(hap_file):
    print('Read Haplotypes')

    haplotypes = defaultdict(list)


    if not hap_file:
        return haplotypes

    for line in hap_file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        # ht-id chrname left right snplist(,)

        try:
            htid, chrname, left, right, snplist = line.split('\t')

            haplotypes[chrname].append([htid, chrname, int(left), int(right), snplist.split(',')])
        except ValueError:
            continue

    return haplotypes

def find_snp_gene(chr_snp_list, chr_exon_list, snps_list):

    # chr_snp_list and chr_exon_list are sorted by position
    # for each snp, find gene_id and added to snps_list[gene_id]

    snp_i = 0
    exon_j = 0

    while snp_i < len(chr_snp_list) and exon_j < len(chr_exon_list):

        # snp position(inclusive)
        snp_start = chr_snp_list[snp_i][2]
        snp_end = snp_start
        snp_type = chr_snp_list[snp_i][1]

        if snp_type == 'deletion':
            snp_end = snp_start + chr_snp_list[snp_i][3] - 1

        # find gene that includes [snp_start, snp_end] region

        j = exon_j
        while j < len(chr_exon_list):
            ex = chr_exon_list[j]

            if snp_end < ex[0]:
                break

            if ex[0] <= snp_start and snp_end <= ex[1]:
                gene_id = ex[2]
                snps_list[gene_id].append(chr_snp_list[snp_i])
                break

            j += 1

        exon_j = j
        snp_i += 1

    return


def make_gene_snps(snps, genes, exons_list=None):
    # snps_list[gene] = [list of snps within the gene]
    snps_list = defaultdict(list)

    # sort snps by position
    for chrname in snps:
        snps[chrname].sort(key=lambda x: x[2])
        if bVerbose:
            print(chrname, len(snps[chrname]))

    # create temporary list of exons per each chromosome
    chr_exons_list = defaultdict(list)
    chr_gene_list = groupby_gene(genes)
    for chrname, gene_list in chr_gene_list.items():
        for gene_id in gene_list:
            ex = genes[gene_id][4]
            chr_exons_list[chrname] += [[e[0], e[1], gene_id] for e in ex]

    for chrname in chr_exons_list:
        chr_exons_list[chrname].sort()

    for chrname, chr_snp_list in snps.items():
        find_snp_gene(chr_snp_list, chr_exons_list[chrname], snps_list)

    return snps_list


def merge_range(A, B):
    assert A[0] <= B[0]

    if (A[1] + 1) < B[0]:
        return None

    return [A[0], max(A[1], B[1])]


def make_consensus_exons(exons_list):
    assert len(exons_list) >= 1

    exons_list.sort()

    assert exons_list[0][0] <= exons_list[0][1]
    cons_exon_list = [exons_list[0]]

    for i in range(1, len(exons_list)):
        last_exon = cons_exon_list[-1]
        curr_exon = exons_list[i]

        assert curr_exon[0] <= curr_exon[1]

        # check two ranges
        m = merge_range(last_exon, curr_exon)
        if m is None:
            cons_exon_list.append([*curr_exon])
        else:
            cons_exon_list[-1] = m

    accum_len = 0
    for e in cons_exon_list:
        e.append(accum_len)
        accum_len += e[1] - e[0] + 1

    return cons_exon_list


def get_seq(ref, exons):
    seq = ''

    for e in exons:
        seq += ref[e[0]:e[1] + 1]

    return seq


def write_transcripts_seq(fp, gene_id, chr_name, genome_seq, exon_list):
    chr_seq = genome_seq[chr_name]
    seq = get_seq(chr_seq, exon_list)

    write_fasta(fp, gene_id, seq)

    return


def map_to_exons(old_exon, consensus_exon_list):
    for idx, cur_exon in enumerate(consensus_exon_list):
        if old_exon[0] >= cur_exon[0] and old_exon[1] <= cur_exon[1]:
            offset = old_exon[0] - cur_exon[0]

            return [cur_exon[2] + offset, cur_exon[2] + offset + (old_exon[1] - old_exon[0])], idx

    print("Warning: Can't find old_exon", old_exon, file=stderr)

    return [-1, -1], -1


def write_transcripts_snp(fp, gene_id, trans_ids, transcripts, exon_list, real_gene_id='', offset=0):
    snps = set()

    for trans_id in trans_ids:
        trans = transcripts[trans_id]

        if bDebug:
            pprint.pprint(trans)

        old_exons = trans[3]

        assert len(old_exons) >= 1

        new_exons = list()
        for old_exon in old_exons:
            ne, idx = map_to_exons(old_exon, exon_list)

            if bDebug:
                print(ne, idx)

            if idx == -1:
                continue

            new_exons.append(ne)

        for i in range(1, len(new_exons)):
            gap = new_exons[i][0] - new_exons[i - 1][1] - 1
            if gap > 0:
                if bUniqueSNP:
                    snps.add(('', new_exons[i - 1][1] + 1 + offset, gap, ''))
                else:
                    snps.add((trans[0], new_exons[i - 1][1] + 1 + offset, gap, trans_id))

    snps_list = sorted(list(snps))

    snp_count = 0

    snp_name = gene_id
    if real_gene_id:
        snp_name = real_gene_id
    for s in snps_list:
        print('{}.{}\t{}\t{}\t{}\t{}'.format(snp_name, snp_count, 'deletion', gene_id, s[1], s[2]), file=fp)
        snp_count += 1

    return


def write_transcripts_ss(fp, gene_id, gene, trans_ids, transcripts, exon_list, offset=0):
    ss = set()

    for trans_id in trans_ids:
        trans = transcripts[trans_id]

        if bDebug:
            pprint.pprint(trans)

        old_exons = trans[3]

        assert len(old_exons) >= 1

        new_exons = list()
        for idx in trans[5]:
            new_start = exon_list[idx][2]
            new_len = exon_list[idx][1] - exon_list[idx][0] + 1

            new_exons.append([new_start, new_start + new_len - 1])

        for i in range(1, len(new_exons)):
            gap = new_exons[i][0] - new_exons[i - 1][1] - 1
            if gap > 0:
                ss.add((new_exons[i - 1][1], new_exons[i][0], trans[1]))

    ss_list = sorted(list(ss))

    for s in ss_list:
        print('{}\t{}\t{}\t{}'.format(gene_id, s[0] + offset, s[1] + offset, s[2]), file=fp)

    return


def write_transcripts_gene_snp(fp, new_chrname, trans_ids, transcripts, snp_list, exon_list, real_gene_id='', offset=0):
    #print(new_chrname, real_gene_id)
    for snp in snp_list:

        new_rsid = "{}.{}".format(snp[0], real_gene_id)

        old_pos = [snp[2], snp[2]]
        new_pos, idx = map_to_exons(old_pos, exon_list)

        if bDebug:
            print(old_pos, new_pos, idx)

        if idx == -1:
            print('Cannot map', snp)
            continue

        print('\t'.join([new_rsid, snp[1], new_chrname, str(new_pos[0] + offset), str(snp[3])]), file=fp)
        # update new position
        snp[2] = new_pos[0] + offset

    return


def write_transcripts_map(fp, gene_id, chr_name, exon_list, real_gene_id='', offset=0, tome_len=0):
    header = '>{}\t{}'.format(gene_id, chr_name)

    if real_gene_id:
        header += '\t{}\t{}\t{}'.format(real_gene_id, offset, tome_len)

    print(header, file=fp)

    for i, e in enumerate(exon_list):
        fp.write('{}:{}:{}'.format(chr_name, e[0], e[1] - e[0] + 1))
        if i % 4 == 3:
            fp.write('\n')
        else:
            fp.write('\t')

    if len(exon_list) % 4 != 0:
        fp.write('\n')

    return


def write_transcripts_tidmap(fp, chrtome_name, chr_name, gene_id, tid_map):
    num_bit_len = tid_map[0]
    num_trans = len(tid_map[1])
    header = '>{}\t{}\t{}\t{}\t{}'.format(chrtome_name, chr_name, gene_id, num_bit_len, num_trans)

    print(header, file=fp)

    for t in tid_map[1]:
        fp.write('{:X}\t{}\n'.format(t[0], t[1]))

    return


def write_transcripts_haplotype(trhap_file, chrtome_name, old_chrname, haplotypes, tmp_snp_list):
    # haplotype file
    # ht-id chrname left right snplist(,)

    snp_dict = {t[0]: t for t in tmp_snp_list}

    # htid, chrname, int(left), int(right), snplist.split(',')])
    for ht in haplotypes:
        """
        if ht[0] == 'ht7472786':
            print(ht)
        """
        snps = []
        new_snps = []
        for snpid in ht[4]:
            if snpid in snp_dict:
                snps.append(snpid)
                new_snps.append('{}.{}'.format(snpid, snp_dict[snpid][5]))

        if not snps:
            # print('Skip:', ht)
            continue

        def get_left_right(s):
            l = s[2]
            r = s[2] if s[1] != 'deletion' else s[2] + s[3] - 1
            return l, r

        new_left, new_right = get_left_right(snp_dict[snps[0]])

        for s in snps[1:]:
            l, r = get_left_right(snp_dict[s])

            if l < new_left:
                new_left = l
            if r > new_right:
                new_right = r

        print('\t'.join([ht[0], chrtome_name, str(new_left), str(new_right), ','.join(new_snps)]), file=trhap_file)


    # print('tmp')
    return


def groupby_gene(genes):
    chr_gene_list = dict()

    for k, g in itertools.groupby(genes.items(), key=lambda x : x[1][3]):
        chr_gene_list[k] = [ x[0] for x in g ]

    return chr_gene_list


def make_chrtome_seq(exon_list, out_seq, ref_seq):
    out_seq = get_seq(ref_seq, exon_list)
    return


def generate_chr_tome(genes, transcripts, gene_ids, snps_list, haplotypes, genome_seq, tid_map,
                      trseq_file, trsnp_file, trmap_file, trss_file, trhap_file, trtid_file):

    chr_gene_list = groupby_gene(genes)

    for chrname, chrgene in chr_gene_list.items():
        print(chrname, len(chrgene))
        offset = 0
        seq = ''

        chrtome_name = '{}_tome'.format(chrname)

        tmp_snp_list = list()

        for g in chrgene:
            gene = genes[g]
            exons = gene[4]
            tmp_seq = get_seq(genome_seq[chrname], exons)
            seq += tmp_seq + 'N'

            tran_ids = gene[0]

            # write ss, snp, map
            write_transcripts_ss(trss_file, chrtome_name, gene, tran_ids, transcripts, exons, offset)
            write_transcripts_gene_snp(trsnp_file, chrtome_name, tran_ids, transcripts, snps_list[g], exons, g, offset)
            write_transcripts_map(trmap_file, chrtome_name, chrname, exons, g, offset, len(tmp_seq))

            write_transcripts_tidmap(trtid_file, chrtome_name, chrname, g, tid_map[g])

            offset += len(tmp_seq) + 1

            tmp_snp_list += [[*snp, g] for snp in snps_list[g]]

        # write new haplotypes
        write_transcripts_haplotype(trhap_file, chrtome_name, chrname, haplotypes[chrname], tmp_snp_list)

        write_fasta(trseq_file, chrtome_name, seq)

    return


def make_unique_range(ranges):

    tmp_range = list()
    for r in ranges:
        tmp_range.append([r[0], -1])
        tmp_range.append([r[1], 1])

    tmp_range.sort()

    range_list = list()
    opened = 0
    last_pos = -1

    for r in tmp_range:

        if opened == 0:
            assert r[1] < 0

            last_pos = r[0]
            opened += 1
            continue

        if r[1] < 0:
            opened += 1
            if last_pos == r[0]:
                continue

            assert last_pos <= r[0] - 1
            range_list.append([last_pos, r[0] - 1])
            last_pos = r[0]
        else:
            opened -= 1
            if last_pos == r[0] + 1:
                continue

            assert last_pos <= r[0]
            range_list.append([last_pos, r[0]])
            last_pos = r[0] + 1

    assert opened == 0

    if len(range_list) > 0:
        last_pos = range_list[0][1]
        for r in range_list[1:]:
            assert last_pos < r[0]
            last_pos = r[1]

    return range_list


def generate_common_exons(genes, transcripts):
    if bVerbose:
        print('Genes:', len(genes))
        print('Transcripts:', len(transcripts))

    for gene_id, gene_value in genes.items():
        if bDebug:
            print(gene_id)

        all_exons = []

        for tid in gene_value[0]:
            trans = transcripts[tid]
            assert gene_id == trans[4]
            all_exons += trans[3]

        unit_exons = make_unique_range(all_exons)

        tmp_common_exons = list()
        offset = 0
        for r in unit_exons:
            tmp_common_exons.append([r[0], r[1], offset])
            offset += (r[1] - r[0] + 1)

        if bDebug:
            print(len(tmp_common_exons), tmp_common_exons)


        gene_value[4] = tmp_common_exons

    return


def map_new_exons(new_exons, old_exons):
    id_list = list()

    new_i = 0
    new_j = 0

    for old_e in old_exons:

        while new_i < len(new_exons):
            if new_exons[new_i][0] == old_e[0]:
                break
            new_i += 1

        assert new_i < len(new_exons)

        new_j = new_i

        while new_j < len(new_exons):
            if new_exons[new_j][1] == old_e[1]:
                break
            new_j += 1

        assert new_j < len(new_exons)
        assert new_i <= new_j

        tmp_len = 0
        for i in range(new_i, new_j + 1):
            id_list.append(i)
            tmp_len += (new_exons[i][1] - new_exons[i][0] + 1)

        assert tmp_len == (old_e[1] - old_e[0] + 1)

        new_i = new_j

    return id_list


def translate_to_common_exons(genes, transcripts):
    # print(transcripts)

    for tran_id, tran_item in transcripts.items():
        old_exons = tran_item[3]
        gene_id = tran_item[4]
        new_exons = genes[gene_id][4]

        id_list = map_new_exons(new_exons, old_exons)
        # print(tran_item)
        tran_item[5] = id_list
        # print(gene_id, old_exons, new_exons, id_list)

    return


def generate_exon_tid_map(genes, transcripts):
    # for each gene, make a list of (exon bit map, tid)

    tid_map = dict()

    # gene [list(), gene_biotype, gene_name, chrom, list()]
    for gene_id, gene in genes.items():

        tid_map[gene_id] = [len(gene[4]), list()]

        for tid in gene[0]:
            tran = transcripts[tid]

            # tran[5] is a list of exon index
            bitmap = 0
            for i in tran[5]:
                bitmap += 2 ** i

            # print(len(tran[5]), max(tran[5]), len(gene[4]), bitmap, gene_id, tid)

            tid_map[gene_id][1].append([bitmap, tid])

    return tid_map


def extract_transcript_graph(genome_file, gtf_file, snp_file, hap_file, base_fname):
    genome_seq = read_genome(genome_file)
    genes, transcripts = read_transcript(genome_seq, gtf_file)
    snps = read_snps(snp_file)
    haplotypes = read_haplotypes(hap_file)

    """
    if bDebug:
        pprint.pprint(genes)
        #pprint.pprint(transcripts)

        return
    """

    # get sorted id
    gene_ids = sorted(list(genes.keys()))
    # exons_list = {}

    generate_common_exons(genes, transcripts)

    translate_to_common_exons(genes, transcripts)

    tid_map = generate_exon_tid_map(genes, transcripts)

    # get snps only in the genes
    snps_list = make_gene_snps(snps, genes)

    with open(base_fname + ".gt.fa", "w") as trseq_file, \
            open(base_fname + ".gt.snp", "w") as trsnp_file, \
            open(base_fname + ".gt.map", "w") as trmap_file, \
            open(base_fname + ".gt.ss", "w") as trss_file, \
            open(base_fname + ".gt.haplotype", "w") as trhap_file, \
            open(base_fname + ".gt.tmap", "w") as trtid_file:

        generate_chr_tome(genes, transcripts, gene_ids, snps_list, haplotypes, genome_seq, tid_map,
                          trseq_file, trsnp_file, trmap_file, trss_file, trhap_file, trtid_file)


    return


def write_exon_info(fp, tid, chr_name, strand, transcript_len, exons, gene_id):
    # FIXME: sync with graph (chrtome-mode)
    print('>{}\t{}\t{}\t{}\t{}'.format(tid, chr_name, strand, transcript_len, gene_id), file=fp)

    for i, e in enumerate(exons):
        fp.write('{}:{}:{}'.format(chr_name, e[0], e[1] - e[0] + 1))
        if i % 4 == 3:
            fp.write('\n')
        else:
            fp.write('\t')

    if len(exons) % 4 != 0:
        fp.write('\n')


def extract_transcript_linear(genome_file, gtf_file, snp_file, hap_file, base_fname):
    genome_seq = read_genome(genome_file)
    genes, transcripts = read_transcript(genome_seq, gtf_file)

    transcript_ids = sorted(list(transcripts.keys()))

    with open(base_fname + ".gt.fa", "w") as trseq_file, \
            open(base_fname + ".gt.map", "w") as trmap_file:
        for tid in transcript_ids:
            chrom, strand, transcript_len, exons, gene_id = transcripts[tid]

            seq = get_seq(genome_seq[chrom], exons)

            write_fasta(trseq_file, tid, seq)

            write_exon_info(trmap_file, tid, *transcripts[tid])

    return


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract transcripts')

    parser.add_argument('-g', '--gtf-file',
                        dest='gtf_file',
                        type=FileType('r'),
                        help='input GTF file')

    parser.add_argument('-s', '--snp-file',
                        dest='snp_file',
                        type=FileType('r'),
                        help='input SNP file')

    parser.add_argument('-p', '--haplotype',
                        dest='hap_file',
                        type=FileType('r'),
                        help='input Haplotype file')

    parser.add_argument('-r', '--genome',
                        dest='genome_file',
                        type=FileType('r'),
                        help='reference genome file')

    parser.add_argument('-o', '--output',
                        dest='out_fname',
                        type=str,
                        help='output filename prefix')

    parser.add_argument('--linear',
                        dest='bLinear',
                        action='store_true',
                        default=False,
                        help='Create linear-transcripts files')

    parser.add_argument('--chrtome',
                        dest='bChrTome',
                        action='store_true',
                        default=False,
                        help='Create transcriptome for chromosome')

    parser.add_argument('--unique-snp',
                        dest='bUniqueSNP',
                        action='store_true',
                        default=False,
                        help='Create SNP file with unique SNP')

    parser.add_argument('--debug',
                        dest='bDebug',
                        action='store_true',
                        help='Run in debug mode')

    parser.add_argument('--verbose',
                        dest='bVerbose',
                        action='store_true',
                        help='Show more messages')

    args = parser.parse_args()
    if not args.gtf_file or not args.genome_file: # or not args.out_fname or not args.hap_file:
        parser.print_help()
        exit(1)

    if args.bDebug is not None:
        bDebug = args.bDebug

    if args.bVerbose is not None:
        bVerbose = args.bVerbose

    if args.bUniqueSNP is not None:
        bUniqueSNP = args.bUniqueSNP

    if args.bChrTome is not None:
        bChrTome = args.bChrTome

    if args.bLinear:
        extract_transcript_linear(
            args.genome_file,
            args.gtf_file,
            args.snp_file,
            args.hap_file,
            args.out_fname)
    else:
        extract_transcript_graph(
            args.genome_file,
            args.gtf_file,
            args.snp_file,
            args.hap_file,
            args.out_fname)
