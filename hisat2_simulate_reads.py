#!/usr/bin/env python3
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
#######################################################################
#                        Update Log                                   #
# 09-17-2020:  Micah Thornton  --  Added Class for Methylation        #
#               Pattern Based on Bernoulli Trials                     #  


import os, sys, math, random, re, copy
from collections import defaultdict, Counter
import itertools
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
python2 style randint
"""

def myrandint(m, x):
    s = x - m + 1
    return m + int(random.random() * s)

"""
Random source for sequencing errors
"""
class ErrRandomSource:
    def __init__(self, prob = 0.0, size = 1 << 20):
        self.size = size
        self.rands = []
        for i in range(self.size):
            if random.random() < prob:
                self.rands.append(1)
            else:
                self.rands.append(0)
        self.cur = 0
        
    def getRand(self):
        assert self.cur < len(self.rands)
        rand = self.rands[self.cur]
        self.cur = (self.cur + 1) % len(self.rands)
        return rand

def strsub(char, i, stri):
        if i == 0:
            return(char + stri[1:]);
        elif i == len(stri)-1:
            return(stri[:i]+char);
        else:
            return(stri[:i] + char + stri[(i+1):]);
'''
A Bernoulli Random Source for Methylation status at CpG locations
'''
#class MethylationPatternBernoulliCpG: 
#    def __init__(self, prob = 0.0, nts, seed=0xFEED): 
#        self.size = len(nts); 
#        self.methylated = []; 
#        self.cpgs = []; 
#        random.seed(seed);
#        for i in range(self.size)-1: 
#            if nts[i] == 'C' and nts[i+1] == 'G': 
#                self.cpgs.append(i); 
#        for i in self.cpgs: 
#            if random.random() <= prob: 
#                self.methylated.append(i);  
#        self.cur = 0; 
#
#    def getMethStat(self): 
#        assert self.cur < self.size; 
#        if self.cur in self.methylated:
#            ret = 1; 
#        else: 
#            ret = 0;
#        self.cur = (self.cur + 1) % self.size
#        return ret;
#
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


    chr_filter = [str(x) for x in list(range(1, 23)) + ['X', 'Y']]
    #chr_filter = None

    if chr_filter:
        for chr_id, chr_seq in chr_dic.items():
            if not chr_id in chr_filter: 
                chr_dic.pop(chr_id, None)
    
    return chr_dic


"""
"""
def read_transcript(genome_seq,
                    gtf_file,
                    frag_len,
                    verbose = False):
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
        if chrom not in genome_seq:
            continue
        
        # Zero-based offset
        left, right = int(left) - 1, int(right) - 1
        if left >= right:
            continue

        values_dict = {}
        for attr in values.split(';')[:-1]:
            attr, _, val = attr.strip().partition(' ')
            values_dict[attr] = val.strip('"')

        if feature == 'gene':
            gene_id = values_dict['gene_id']
            if gene_id in genes:
                print(gene_id, 'already exists', file=sys.stderr)
                continue

            gene_biotype = values_dict['gene_biotype']
            gene_name = values_dict['gene_name']
            genes[gene_id] = [list(), gene_biotype, gene_name, chrom]

        elif feature == 'exon':
            transcript_id = values_dict['transcript_id']
            gene_id = values_dict['gene_id']

            if transcript_id not in transcripts:
                transcripts[transcript_id] = [chrom, strand, [[left, right]], gene_id]
                genes[gene_id][0].append(transcript_id)
            else:
                transcripts[transcript_id][2].append([left, right])
        else:
            continue


    merged_exon = 0
    min_transcript_len_filtered = 0

    # Sort exons and merge where separating introns are <=5 bps
    for tran, [chr, strand, exons, gene_id] in transcripts.items():
        exons.sort()
        tmp_exons = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - tmp_exons[-1][1] <= 5:
                merged_exon += 1
                tmp_exons[-1][1] = exons[i][1]
            else:
                tmp_exons.append(exons[i])
        transcripts[tran] = [chr, strand, tmp_exons, gene_id]

    tmp_transcripts = {}
    for tran, [chr, strand, exons, gene_id] in transcripts.items():
        exon_lens = [e[1] - e[0] + 1 for e in exons]
        transcript_len = sum(exon_lens)
        if transcript_len >= frag_len:
            tmp_transcripts[tran] = [chr, strand, transcript_len, exons, gene_id]
        else:
            min_transcript_len_filtered += 1
            genes[gene_id][0].remove(tran)
            # Remove a gene with no transcripts
            if len(genes[gene_id][0]) == 0:
                del genes[gene_id]

    transcripts = tmp_transcripts

    if verbose:
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
        print("merged exons: {}, filtered transcripts (minimum length: {}): {}".format(merged_exon, frag_len, min_transcript_len_filtered))
        print("total exons: {}, average exons: {:.2f}".format(total_exons, avg_exons))

    return genes, transcripts
    

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
        except ValueError:
            continue

        assert type in ["single", "deletion", "insertion"]
        if type == "deletion":
            data = int(data)
        snps[chr].append([snpID, type, int(pos), data])

    return snps


"""
"""
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

            j += 1

        exon_j = j
        snp_i += 1

    return


"""
"""
def groupby_gene(genes):
    chr_genes = dict()

    for k, g in itertools.groupby(genes.items(), key=lambda x : x[1][3]):
        chr_genes[k] = [ x[0] for x in g ]

    return chr_genes


"""
"""
def make_gene_snps(genome_snps, genes, super_transcripts):
    # snps[gene] = [list of snps within the gene]
    snps = defaultdict(list)

    # Sort SNPs by position
    for chr in genome_snps:
        genome_snps[chr].sort(key=lambda x: x[2])

    # Create a temporary list of exons per chromosome
    chr_exons = defaultdict(list)
    chr_genes = groupby_gene(genes)
    for chr, genes in chr_genes.items():
        for gene_id in genes:
            if gene_id not in super_transcripts:
                continue
            
            super_transcript = super_transcripts[gene_id][3]
            chr_exons[chr] += [[e[0], e[1], gene_id] for e in super_transcript]

    for chr in chr_exons:
        chr_exons[chr].sort()

    for chr, chr_snps in genome_snps.items():
        find_snp_gene(chr_snps, chr_exons[chr], snps)

    return snps


"""
"""
def make_super_transcripts(exons_list):
    assert len(exons_list) >= 1

    exons_list.sort()

    assert exons_list[0][0] <= exons_list[0][1]
    cons_exon_list = [exons_list[0]]

    for i in range(1, len(exons_list)):
        last_exon = cons_exon_list[-1]
        curr_exon = exons_list[i]

        assert curr_exon[0] <= curr_exon[1]

        def merge_range(A, B):
            assert A[0] <= B[0]
            if (A[1] + 1) < B[0]:
                return None

            return [A[0], max(A[1], B[1])]

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


"""
"""
def sanity_check_input(genome_seq, genes, transcripts, snps, frag_len):
    num_canon_ss, num_ss = 0, 0
    for transcript, [chr, strand, transcript_len, exons, gene_id] in transcripts.items():
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

    if num_ss > 0:
        print("GT/AG splice sites: {}/{} ({:.2%})".format(num_canon_ss, num_ss, (float(num_canon_ss) / num_ss)), file=sys.stderr)

    num_alt_single, num_single = 0, 0
    for chr, chr_snps in snps.items():
        if chr not in genome_seq:
            continue
        chr_seq = genome_seq[chr]
        prev_snp = None
        for snp in chr_snps:
            snpID, type, pos, data = snp
            if prev_snp:
                assert prev_snp[2] <= pos
            prev_snp = snp
            if type != "single":
                continue
            assert pos < len(chr_seq)
            if chr_seq[pos] != data:
                num_alt_single += 1
            num_single += 1

    if num_single > 0:
        print("Alternative bases: {}/{} ({:.2%})".format(num_alt_single, num_single, (float(num_alt_single) / num_single)), file=sys.stderr)


"""
"""
def generate_rna_expr_profile(expr_profile_type, num_transcripts = 10000):
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
def generate_dna_expr_profile(genome_seq):
    expr_profile = []
    for chr_id, chr_seq in genome_seq.items():
        expr_profile.append(len(chr_seq))
    expr_sum = float(sum(expr_profile))
    expr_profile = [expr_profile[i] / expr_sum for i in range(len(expr_profile))]
    assert abs(sum(expr_profile) - 1.0) < 0.001
    return expr_profile


"""
"""
def getSNPs(ref_snps, left, right):
    low, high = 0, len(ref_snps)
    while low < high:
        mid = (low + high) // 2
        snpID, type, pos, data = ref_snps[mid]
        if pos < left:
            low = mid + 1
        else:
            high = mid - 1

    snps = []
    for snp in ref_snps[low:]:
        snpID, type, pos, data = snp
        pos2 = pos
        if type == "deletion":
            pos2 += data
        if pos2 >= right:
            break
        if pos >= left:
            if len(snps) > 0:
                _, prev_type, prev_pos, prev_data = snps[-1]
                assert prev_pos <= pos
                prev_pos2 = prev_pos
                if prev_type == "deletion":
                    prev_pos2 += prev_data
                if pos <= prev_pos2:
                    continue
            snps.append(snp)

    return snps


"""
"""
def getSAMAlignment(ref_seq, frag_pos, read_len, ref_snps, err_rand_src, max_mismatch):
    # Find the genomic position for frag_pos
    pos, cigars, cigar_descs = frag_pos, [], []
    tmp_read_len = read_len

    ref_len = len(ref_seq)
    match_len = 0
    mismatch, remain_chr_len = 0, ref_len - (frag_pos + read_len)
    assert remain_chr_len >= 0

    # Retreive SNPs
    snps = getSNPs(ref_snps, frag_pos, frag_pos + read_len)

    # Simulate mismatches due to sequencing errors
    mms = []
    for i in range(pos, min(ref_len, pos + tmp_read_len - 1)):
        if err_rand_src.getRand() == 1:
            assert i < len(chr_seq)
            err_base = "A"
            rand = myrandint(0, 2)
            if chr_seq[i] == "A":
                err_base = "GCT"[rand]
            elif chr_seq[i] == "C":
                err_base = "AGT"[rand]
            elif chr_seq[i] == "G":
                err_base = "ACT"[rand]
            else:
                err_base = "ACG"[rand]                    
            mms.append(["", "single", i, err_base])

    tmp_diffs = snps + mms
    tmp_diffs = sorted(tmp_diffs, key=lambda t: t[2])
    diffs = []
    if len(tmp_diffs) > 0:
        diffs = tmp_diffs[:1]
        for diff in tmp_diffs[1:]:
            _, tmp_type, tmp_pos, tmp_data = diff
            _, prev_type, prev_pos, prev_data = diffs[-1]
            if prev_type == "deletion":
                prev_pos += prev_data
            if tmp_pos <= prev_pos:
                continue
            diffs.append(diff)

    cigar_left, desc_left = pos, pos
    cigar_descs.append([])
    prev_diff = None
    for diff in diffs:
        diff_id, diff_type, diff_pos, diff_data = diff
        if prev_diff:
            prev_diff_id, prev_diff_type, prev_diff_pos, prev_diff_data = prev_diff
            if prev_diff_type == "deletion":
                prev_diff_pos += prev_diff_data
            assert prev_diff_pos < diff_pos
        diff_pos2 = diff_pos
        if diff_type == "deletion":
            diff_pos2 += diff_data
        if cigar_left + tmp_read_len - 1 < diff_pos2:
            break

        if diff_type == "single":
            if mismatch + 1 > max_mismatch:
                continue
            cigar_descs[-1].append([diff_pos - desc_left, diff_data, diff_id])
            desc_left = diff_pos + 1
            mismatch += 1
            
        elif diff_type == "deletion":
            del_len = diff_data
            if mismatch + del_len > max_mismatch:
                continue
            if len(cigars) <= 0 and diff_pos - cigar_left <= 0:
                continue                
            if remain_chr_len < del_len:
                continue
            remain_chr_len -= del_len
            if diff_pos - cigar_left > 0:
                cigars.append("{}M".format(diff_pos - cigar_left))
                cigar_descs[-1].append([diff_pos - desc_left, "", ""])
                cigar_descs.append([])
            cigars.append("{}D".format(del_len))
            cigar_descs[-1].append([0, del_len, diff_id])
            cigar_descs.append([])
            tmp_read_len -= (diff_pos - cigar_left)
            cigar_left = desc_left = diff_pos + del_len

        elif diff_type == "insertion":
            ins_len = len(diff_data)
            if mismatch + ins_len > max_mismatch:
                continue
            if len(cigars) <= 0 and diff_pos - cigar_left <= 0:
                continue
            if cigar_left + tmp_read_len - 1 < diff_pos + ins_len:
                break
            if diff_pos - cigar_left > 0:
                cigars.append("{}M".format(diff_pos - cigar_left))
                cigar_descs[-1].append([diff_pos - desc_left, "", ""])
                cigar_descs.append([])
            cigars.append("{}I".format(ins_len))
            cigar_descs[-1].append([0, diff_data, diff_id])
            cigar_descs.append([])
            tmp_read_len -= (diff_pos - cigar_left)
            tmp_read_len -= ins_len
            cigar_left = desc_left = diff_pos

        else:
            assert False
            
        prev_diff = diff

    cigar_right = min(ref_len, cigar_left + tmp_read_len - 1)
    cigar_len = cigar_right - cigar_left + 1
    remain_len = cigar_right - desc_left + 1
    if remain_len > 0:
        cigar_descs[-1].append([remain_len, "", ""])
    tmp_read_len -= cigar_len
    cigars.append(("{}M".format(cigar_len)))

    # Define MD, XM, NM, Zs, read_seq
    MD, XM, NM, Zs, read_seq = "", 0, 0, "", ""
    assert len(cigars) == len(cigar_descs)
    MD_match_len, Zs_match_len = 0, 0
    cur_pos = frag_pos
    for c in range(len(cigars)):
        cigar = cigars[c]
        cigar_len, cigar_op = int(cigar[:-1]), cigar[-1]
        cigar_desc = cigar_descs[c]
        if cigar_op == 'N':
            continue
        if cigar_op == 'M':
            for add_match_len, alt_base, snp_id in cigar_desc:
                MD_match_len += add_match_len
                Zs_match_len += add_match_len
                assert cur_pos + add_match_len <= ref_len
                read_seq += ref_seq[cur_pos:cur_pos+add_match_len]
                cur_pos += add_match_len
                if alt_base != "":
                    if MD_match_len > 0:
                        MD += ("{}".format(MD_match_len))
                        MD_match_len = 0
                    MD += ref_seq[cur_pos]
                    if snp_id != "":
                        if Zs != "":
                            Zs += ","
                        Zs += ("{}|S|{}".format(Zs_match_len, snp_id))
                        Zs_match_len = 0
                    else:
                        Zs_match_len += 1
                    if snp_id == "":
                        XM += 1
                        NM += 1
                    read_seq += alt_base
                    cur_pos += 1
        elif cigar_op == 'D':
            assert len(cigar_desc) == 1
            add_match_len, del_len, snp_id = cigar_desc[0]
            MD_match_len += add_match_len
            Zs_match_len += add_match_len
            if MD_match_len > 0:
                MD += ("{}".format(MD_match_len))
                MD_match_len = 0
            MD += ("^{}".format(ref_seq[cur_pos:cur_pos+cigar_len]))
            read_seq += ref_seq[cur_pos:cur_pos+add_match_len]
            if Zs != "":
                Zs += ","
            Zs += ("{}|D|{}".format(Zs_match_len, cigar_desc[0][-1]))
            Zs_match_len = 0
            cur_pos += cigar_len
        elif cigar_op == 'I':
            assert len(cigar_desc) == 1
            add_match_len, ins_seq, snp_id = cigar_desc[0]
            ins_len = len(ins_seq)
            MD_match_len += add_match_len
            Zs_match_len += add_match_len
            read_seq += ref_seq[cur_pos:cur_pos+add_match_len]
            read_seq += ins_seq
            if Zs != "":
                Zs += ","
            Zs += ("{}|I|{}".format(Zs_match_len, cigar_desc[0][-1]))
            Zs_match_len = 0
        else:
            assert False

    if MD_match_len > 0:
        MD += ("{}".format(MD_match_len))

    if len(read_seq) != read_len:
        print("read length differs:", len(read_seq), "vs.", read_len, file=sys.stderr)
        print(pos, "".join(cigars), cigar_descs, MD, XM, NM, Zs, file=sys.stderr)
        assert False

    return pos, cur_pos, cigars, cigar_descs, MD, XM, NM, Zs, read_seq


"""
"""
def convertCoordFromTranscriptToGenome(transcript, pos):
    # Convert a transcript-based coordinate to a chromosome-based coordinate
    gpos = 0
    for ei, e in enumerate(transcript):
        e_left, e_right = e[:2]
        e_len = e_right - e_left + 1
        if pos < e_len:
            gpos = e_left + pos
            return gpos, ei
        else:
            pos -= e_len

    return gpos, -1


"""
"""
def convertCoordFromTranscriptToGenome_two(transcript, tleft, tright):
    gleft, gleft_i = convertCoordFromTranscriptToGenome(transcript, tleft)
    gright, gright_i = convertCoordFromTranscriptToGenome(transcript, tleft)

    return gleft, gright, gleft_i, gright_i


"""
"""
def convertCoordFromGenomeToTranscript(transcript, gpos):
    # Convert a chromosome-based coordinate to a transcript-based coordinate
    tpos = 0
    for ei, e in enumerate(transcript):
        e_left, e_right = e[:2]
        e_len = e_right - e_left + 1
        if gpos < e_left:
            return -1, -1
        if gpos <= e_right:
            tpos += (gpos - e_left)
            return tpos, ei
        else:
            tpos += e_len

    return -1, -1


"""
"""
def convertCoordFromGenomeToTranscript_two(transcript, gleft, gright):
    tleft, tleft_i = convertCoordFromGenomeToTranscript(transcript, gleft)
    tright, tright_i = convertCoordFromGenomeToTranscript(transcript, gleft)

    return tleft, tright, tleft_i, tright_i


"""
"""
def convertSAMAlignmentFromTranscriptToGenome(transcript, tpos, tcigars):
    gleft, gpos_i = convertCoordFromTranscriptToGenome(transcript, tpos)
    gright = gleft
    gcigars = []
    for cigar in tcigars:
        cigar_len, cigar_op = int(cigar[:-1]), cigar[-1]

        if cigar_op == 'M' or cigar_op == 'D':
            while cigar_len > 0:
                assert gpos_i < len(transcript)
                e_left, e_right = transcript[gpos_i][:2]
                assert e_left <= gright and gright <= e_right
                e_remain_len = e_right - gright + 1
                if cigar_len < e_remain_len:
                    gcigars.append("%d%s" % (cigar_len, cigar_op))
                    gright += cigar_len
                    break
                else:
                    gcigars.append("%d%s" % (e_remain_len, cigar_op))
                    cigar_len -= e_remain_len
                    gpos_i += 1
                    if gpos_i == len(transcript):
                        assert cigar_len == 0
                        break

                    assert gpos_i < len(transcript)
                    gright = transcript[gpos_i][0]

                    if cigar_len > 0:
                        gcigars.append("%dN" % (gright - e_right - 1))        
        else:
            assert cigar_op == 'I'
            gcigars.append(cigar)

    return gleft, gcigars

def convertSAMAlignmentFromTranscriptToReal(transcript, tpos, tcigars):
    return convertSAMAlignmentFromTranscriptToGenome(transcript, tpos, tcigars)


"""
"""
def convertSAMAlignmentFromTranscriptToGene(transcript, tpos, tcigars, super_transcript):
    def convertTranscript(super_transcript, transcript, tpos):
        assert len(transcript) > 0

        tmp_exons = []
        s_i, t_i = 0, 0
        while t_i < len(transcript):
            se_left, se_right, se_pos = super_transcript[s_i]
            te_left, te_right = transcript[t_i]
            assert se_left <= te_left
            if te_left > se_right:
                s_i += 1
                continue

            tmp_exons.append([te_left - se_left + se_pos, te_right - se_left + se_pos])
            t_i += 1

        exons = [tmp_exons[0]]
        for te in tmp_exons[1:]:
            if exons[-1][1] + 1 == te[0]:
                exons[-1][1] = te[1]
            else:
                exons.append(te)

        return exons
        
    exons = convertTranscript(super_transcript, transcript, tpos)
    s_pos, s_cigars = convertSAMAlignmentFromTranscriptToReal(exons, tpos, tcigars)

    return s_pos, s_cigars
    

"""
"""
cigar_re = re.compile('\d+\w')
def SAMRepOk(ref_seq, read_seq, chr, pos, cigar, XM, NM, MD, Zs, max_mismatch):
    assert pos < len(ref_seq)

    # Calculate XM and NM based on Cigar and Zs
    cigars = cigar_re.findall(cigar)
    cigars = [[int(cigars[i][:-1]), cigars[i][-1]] for i in range(len(cigars))]
    ref_pos, read_pos = pos, 0
    ann_ref_seq, ann_ref_rel, ann_read_seq, ann_read_rel = [], [], [], []
    for i in range(len(cigars)):
        cigar_len, cigar_op = cigars[i]
        if cigar_op == "M":
            partial_ref_seq = ref_seq[ref_pos:ref_pos+cigar_len]
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
            partial_ref_seq = ref_seq[ref_pos:ref_pos+cigar_len]
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

    if tMD != MD or tXM != XM or tNM != NM or XM > max_mismatch or XM != NM:
        print(chr, pos, cigar, MD, XM, NM, Zs, file=sys.stderr)
        print(tMD, tXM, tNM, file=sys.stderr)
        assert False


"""
"""
def simulate_DNA_reads(base_fname, repeat_fname, genome_seq, snps,
                       paired_end, read_len, frag_len, num_frag,
                       error_rate, err_rand_src, max_mismatch, snp_prob,
                       bisulfite, meth_prob,
                       allele_specific,
                       sanity_check, verbose):
    expr_profile = generate_dna_expr_profile(genome_seq)
    expr_profile = [int(expr_profile[i] * num_frag) for i in range(len(expr_profile))]
    assert num_frag >= sum(expr_profile)
    while sum(expr_profile) < num_frag:
        for i in range(min(num_frag - sum(expr_profile), len(expr_profile))):
            expr_profile[i] += 1
    assert num_frag == sum(expr_profile)
    
    repeat_loci = {}
    if repeat_fname != "" and os.path.exists(repeat_fname):
        for line in open(repeat_fname):
            if line.startswith('>'):
                continue
            coords = line.strip().split()
            for coord in coords:
                chr, pos, strand = coord.split(':')
                if chr not in repeat_loci:
                    repeat_loci[chr] = []
                repeat_loci[chr].append([int(pos), strand])

    chr_ids = list(genome_seq.keys())

    sam_file = open(base_fname + ".sam", "w")
    stat_file = open(base_fname + ".stat", "w")

    # Write SAM header
    print("@HD\tVN:1.0\tSO:unsorted", file=sam_file)
    for chr in genome_seq.keys():
        print("@SQ\tSN:%s\tLN:%d" % (chr, len(genome_seq[chr])), file=sam_file)
    
    read_file = open(base_fname + "_1.fa", "w")
    if paired_end:
        read2_file = open(base_fname + "_2.fa", "w")

    cur_read_id = 1
    for t in range(len(expr_profile)):
        t_num_frags = expr_profile[t]
        chr = chr_ids[t]
        print(chr, t_num_frags, file=sys.stderr)

        if allele_specific:
            ref_num_frags = int(t_num_frags / 2)
            alt_num_frags = t_num_frags - ref_num_frags
            print("%s_ref\t%d" % (chr, ref_num_frags), file=stat_file)
            print("%s_alt\t%d" % (chr, ref_num_frags), file=stat_file)

        assert chr in genome_seq
        chr_seq = genome_seq[chr]
        chr_len = len(chr_seq)
        if chr in repeat_loci:
            chr_repeat_loci = repeat_loci[chr]
        else:
            chr_repeat_loci = []
            
        t_seq = chr_seq
        exons = [[0, chr_len - 1]]

        if chr in snps:
            chr_snps = snps[chr]
        else:
            chr_snps = []

        # Retreive SNPs
        snps = chr_snps
        if snp_prob < 1.0 and len(snps) > 0:
            snps_ = []
            for snp in snps:
                if random.random() <= snp_prob:
                    snps_.append(snp)
            snps = snps_

        for f in range(t_num_frags):
            while True:
                if len(chr_repeat_loci):
                    locus_id = myrandint(0, len(chr_repeat_loci) - 1)
                    frag_pos = chr_repeat_loci[locus_id][0]
                else:
                    frag_pos = myrandint(0, chr_len - frag_len)
                if 'N' not in chr_seq[frag_pos:frag_pos + frag_len]:
                    break

            # DK - 50%/50% alleles
            if allele_specific:
                if f < t_num_frags / 2:
                    t_snp_prob = 0.0
                else:
                    t_snp_prob = 1.0
            else:
                t_snp_prob = snp_prob

            # SAM specification (v1.4)
            # http://samtools.sourceforge.net/
            flag, flag2 = 99, 163  # 83, 147
            pos, rpos, cigars, cigar_descs, MD, XM, NM, Zs, read_seq = getSAMAlignment(chr_seq, frag_pos, read_len, snps, err_rand_src, max_mismatch)
            pos2, rpos2, cigars2, cigar2_descs, MD2, XM2, NM2, Zs2, read2_seq = getSAMAlignment(chr_seq, frag_pos+frag_len-read_len, read_len, snps, err_rand_src, max_mismatch)
            swapped = False
            if paired_end:
                #if random.randint(0, 1) == 1:
                if myrandint(0, 1) == 1:
                    swapped = True
                if swapped:
                    flag, flag2 = flag - 16, flag2 - 16
                    pos, pos2 = pos2, pos
                    rpos, rpos2 = rpos2, rpos
                    cigars, cigars2 = cigars2, cigars
                    cigar_descs, cigar2_descs = cigar2_descs, cigar_descs
                    read_seq, read2_seq = read2_seq, read_seq
                    XM, XM2 = XM2, XM
                    NM, NM2 = NM2, NM
                    MD, MD2 = MD2, MD
                    Zs, Zs2 = Zs2, Zs

            cigar_str, cigar2_str = "".join(cigars), "".join(cigars2)
            if sanity_check:
                SAMRepOk(chr_seq, read_seq, chr, pos, cigar_str, XM, NM, MD, Zs, max_mismatch)
                SAMRepOk(chr_seq, read2_seq, chr, pos2, cigar2_str, XM2, NM2, MD2, Zs2, max_mismatch)

            if Zs != "":
                Zs = ("\tZs:Z:{}".format(Zs))
            if Zs2 != "":
                Zs2 = ("\tZs:Z:{}".format(Zs2))

            if bisulfite:
                assert meth_prob >= 0 and meth_prob <= 1; 
                for i in range(len(read_seq)): 
                    if i == len(read_seq)-1: 
                        if read_seq[i] == 'C' and genome_seq[chr][pos+len(read_seq)] == "G": 
                            if random.random() <= meth_prob: 
                                read_seq = strsub("T", i, read_seq); 
                    else: 
                        if read_seq[i] == 'C' and genome_seq[chr][pos+len(read_seq)] == "G": 
                            if random.random() <= meth_prob: 
                                read_seq = strsub("T", i, read_seq); 

            print(">{}".format(cur_read_id), file=read_file)
            if swapped:
                print(reverse_complement(read_seq), file=read_file)
            else:
                print(read_seq, file=read_file)
                
            print("{}\t{}\t{}\t{}\t255\t{}\t{}\t{}\t0\t{}\t*\tXM:i:{}\tNM:i:{}\tMD:Z:{}{}".format(cur_read_id, flag, chr, pos + 1, cigar_str, chr, pos2 + 1, read_seq, XM, NM, MD, Zs), file=sam_file)
            if paired_end:
                print(">{}".format(cur_read_id), file=read2_file)
                if swapped:
                    print(read2_seq, file=read2_file)
                else:
                    print(reverse_complement(read2_seq), file=read2_file)
                print("{}\t{}\t{}\t{}\t255\t{}\t{}\t{}\t0\t{}\t*\tXM:i:{}\tNM:i:{}\tMD:Z:{}{}".format(cur_read_id, flag2, chr, pos2 + 1, cigar2_str, chr, pos + 1, read2_seq, XM2, NM2, MD2, Zs2), file=sam_file)            

            cur_read_id += 1
            
    sam_file.close()
    stat_file.close()
    read_file.close()
    if paired_end:
        read2_file.close()


"""
"""
def simulate_RNA_reads(base_fname, genome_seq, genes, transcripts, genome_snps,
                       paired_end, read_len, frag_len, num_frag,
                       error_rate, err_rand_src, max_mismatch, snp_prob, expr_profile_type,
                       bisulfite, meth_prob,
                       allele_specific, num_genes, num_transcripts,
                       sanity_check, verbose):


    num_genes = min(len(genes), num_genes)
    gene_ids = list(genes.keys())
    random.shuffle(gene_ids)
    gene_ids = gene_ids[:num_genes]

    # DK - debugging purposes
    # gene_ids = set(["ENSG00000100095"])

    transcript_ids = []
    for gene_id in gene_ids:
        transcript_ids += genes[gene_id][0]

    num_transcripts = min(len(transcript_ids), num_transcripts)

    expr_profile = generate_rna_expr_profile(expr_profile_type, num_transcripts)
    expr_profile = [int(expr_profile[i] * num_frag) for i in range(len(expr_profile))]
    
    assert num_frag >= sum(expr_profile)
    while sum(expr_profile) < num_frag:
        for i in range(min(num_frag - sum(expr_profile), len(expr_profile))):
            expr_profile[i] += 1
    assert num_frag == sum(expr_profile)
    
    random.shuffle(transcript_ids)
    assert len(transcript_ids) >= len(expr_profile)

    transcript_exprs = {}
    for e, num_frags in enumerate(expr_profile):
        transcript_id = transcript_ids[e]
        transcript_exprs[transcript_id] = num_frags

    # Construct a super transcript for each gene in which transcripts are included
    super_transcripts = {}
    for gene_id in gene_ids:
        tran_ids = genes[gene_id][0]
        chr, strand = transcripts[tran_ids[0]][:2]
        tmp_exons = []

        # Collect all exons of the gene
        for tid in tran_ids:
            if tid not in transcript_ids:
                continue
            tmp_exons += copy.deepcopy(transcripts[tid][3])

        exons = make_super_transcripts(tmp_exons)

        def length_of_transcript(exons):
            tran_len = 0
            for e in exons:
                e_l, e_r = e[:2]
                tran_len += (e_r - e_l + 1)
            return tran_len

        transcript_len = length_of_transcript(exons)
        super_transcripts[gene_id] = [chr, strand, transcript_len, exons, gene_id]

    # Get snps in the super transcripts
    snps = make_gene_snps(genome_snps, genes, super_transcripts)

    def find_compatible_transcripts(gene_id, transcript_id, left, right):
        assert left <= right
        
        compatible_transcripts = []
        def find_position(exons, left, right, texons):
            gleft, gright, gleft_i, gright_i = convertCoordFromTranscriptToGenome_two(exons, left, right)
            assert gleft >= 0 and gleft <= gright
            assert gleft_i >= 0 and gleft_i <= gright_i

            left2, right2, left2_i, right2_i = convertCoordFromGenomeToTranscript_two(texons, gleft, gright)
            if left2_i < 0 or right2_i < 0:
                return -1
            assert left2_i <= right2_i

            # Check whether the two transcript have the same number of exons
            if gright_i - gleft_i != right2_i - left2_i:
                return -1

            numexons = gright_i - gleft_i + 1
            for i in range(numexons):
                e, te = exons[gleft_i + i], texons[left2_i + i]
                if i > 0 and e[0] != te[0]:
                    return -1
                if i < numexons - 1 and e[1] != te[1]:
                    return -1
                            
            return left2

        exons = transcripts[transcript_id][3]
        for tid in genes[gene_id][0]:
            if tid == transcript_id:
                continue

            if tid not in transcript_exprs:
                continue

            texons = transcripts[tid][3]
            tleft = find_position(exons, left, right, texons)
            if tleft >= 0:
                compatible_transcripts.append([tid, tleft])
                
        return compatible_transcripts

    sam_file = open(base_fname + ".sam", "w")
    sam_transcript_file = open(base_fname + ".tran.sam", "w")
    sam_gene_file = open(base_fname + ".gene.sam", "w")
    stat_file = open(base_fname + ".stat", "w")

    # Write SAM header
    print("@HD\tVN:1.0\tSO:unsorted", file=sam_file)
    for chr in genome_seq.keys():
        print("@SQ\tSN:%s\tLN:%d" % (chr, len(genome_seq[chr])), file=sam_file)

    print("@HD\tVN:1.0\tSO:unsorted", file=sam_transcript_file)
    for transcript_id in transcript_ids:
        transcript_len = transcripts[transcript_id][2]
        print("@SQ\tSN:%s\tLN:%d" % (transcript_id, transcript_len), file=sam_transcript_file)

    print("@HD\tVN:1.0\tSO:unsorted", file=sam_gene_file)
    for gene_id in gene_ids:
        transcript_len = super_transcripts[gene_id][2]
        print("@SQ\tSN:%s\tLN:%d" % (gene_id, transcript_len), file=sam_gene_file)
    
    read_file = open(base_fname + "_1.fa", "w")
    if paired_end:
        read2_file = open(base_fname + "_2.fa", "w")

    cur_read_id = 1
    for gene_id in gene_ids:
        chr, strand, super_transcript_len, super_exons = super_transcripts[gene_id][:4]

        assert chr in genome_seq
        chr_seq = genome_seq[chr]
        chr_len = len(chr_seq)
            
        super_transcript_seq = ""
        for e in super_exons:
            assert e[0] < e[1]
            super_transcript_seq += chr_seq[e[0]:e[1]+1]
        assert len(super_transcript_seq) == super_transcript_len

        for transcript_id in genes[gene_id][0]:
            if transcript_id not in transcript_exprs:
                continue

            t_num_frags = transcript_exprs[transcript_id]
            print(gene_id, transcript_id, t_num_frags, file=sys.stderr)
            transcript_len, exons = transcripts[transcript_id][2:4]

            if allele_specific:
                ref_num_frags = int(t_num_frags / 2)
                alt_num_frags = t_num_frags - ref_num_frags
                print("%s\t%s_ref\t%d" % (gene_id, transcript_id, ref_num_frags), file=stat_file)
                print("%s\t%s_alt\t%d\tsnps" % (gene_id, transcript_id, ref_num_frags), file=stat_file)

            transcript_seq = ""
            for e in exons:
                assert e[0] < e[1]
                transcript_seq += chr_seq[e[0]:e[1]+1]
            assert len(transcript_seq) == transcript_len

            if gene_id in snps:
                gene_snps = copy.deepcopy(snps[gene_id])
            else:
                gene_snps = []

            # Retreive SNPs
            if snp_prob < 1.0 and len(gene_snps) > 0:
                gene_snps_ = []
                for snp in snps:
                    if random.random() <= snp_prob:
                        snps_.append(snp)
                gene_snps = snps_

            # Convert SNP coordinate from Genome to Transcript
            transcript_snps = []
            for g, snp in enumerate(gene_snps):
                snp_id, snp_type, snp_pos, snp_data = snp
                tpos, tpos_i = convertCoordFromGenomeToTranscript(exons, snp_pos)
                if tpos < 0:
                    continue

                transcript_snps.append([snp_id, snp_type, tpos, snp_data])

            for f in range(t_num_frags):
                frag_pos = myrandint(0, transcript_len - frag_len)

                # DK - 50%/50% alleles
                if allele_specific:
                    if f < t_num_frags / 2:
                        t_snp_prob = 0.0
                    else:
                        t_snp_prob = 1.0
                else:
                    t_snp_prob = snp_prob

                # SAM specification (v1.4)
                # http://samtools.sourceforge.net/
                flag, flag2 = 99, 163  # 83, 147
                pos, rpos, cigars, cigar_descs, MD, XM, NM, Zs, read_seq = getSAMAlignment(transcript_seq, frag_pos, read_len, transcript_snps, err_rand_src, max_mismatch)
                pos2, rpos2, cigars2, cigar2_descs, MD2, XM2, NM2, Zs2, read2_seq = getSAMAlignment(transcript_seq, frag_pos+frag_len-read_len, read_len, transcript_snps, err_rand_src, max_mismatch)
                swapped = False
                if paired_end:
                    if myrandint(0, 1) == 1:
                        swapped = True
                    if swapped:
                        flag, flag2 = flag - 16, flag2 - 16
                        pos, pos2 = pos2, pos
                        rpos, rpos2 = rpos2, rpos
                        cigars, cigars2 = cigars2, cigars
                        cigar_descs, cigar2_descs = cigar2_descs, cigar_descs
                        read_seq, read2_seq = read2_seq, read_seq
                        XM, XM2 = XM2, XM
                        NM, NM2 = NM2, NM
                        MD, MD2 = MD2, MD
                        Zs, Zs2 = Zs2, Zs

                cigar_str, cigar2_str = "".join(cigars), "".join(cigars2)
                if sanity_check:
                    SAMRepOk(transcript_seq, read_seq, transcript_id, pos, cigar_str, XM, NM, MD, Zs, max_mismatch)
                    SAMRepOk(transcript_seq, read2_seq, transcript_id, pos2, cigar2_str, XM2, NM2, MD2, Zs2, max_mismatch)

                Zs_str, Zs2_str = "", ""
                if Zs != "":
                    Zs_str = ("\tZs:Z:{}".format(Zs))
                if Zs2 != "":
                    Zs2_str = ("\tZs:Z:{}".format(Zs2))

                XS = "\tXS:A:{}".format(strand)
                TI = "\tTI:Z:{}-{}".format(transcript_id, pos + 1)
                TI2 = "\tTI:Z:{}-{}".format(transcript_id, pos2 + 1)

                def get_TO(compatible_transcripts):
                    TO = ""
                    if compatible_transcripts:
                        TO = "\tTO:Z:"
                        for i, [tid, tleft] in enumerate(compatible_transcripts):
                            if i > 0:
                                TO += "|"
                            TO += "{}-{}".format(tid, tleft + 1)
                    return TO

                compatible_transcripts = find_compatible_transcripts(gene_id, transcript_id, pos, rpos)
                TO = get_TO(compatible_transcripts)
                compatible_transcripts2 = find_compatible_transcripts(gene_id, transcript_id, pos2, rpos2)
                TO2 = get_TO(compatible_transcripts2)

                # C to T conversion for BS-seq simulation
                if bisulfite:
                    assert meth_prob >= 0 and meth_prob <= 1; 
                    for i in range(len(read_seq)): 
                        if i == len(read_seq)-1: 
                            if read_seq[i] == 'C' and genome_seq[chr][pos+len(read_seq)] == "G": 
                                if random.random() <= meth_prob: 
                                    read_seq = strsub("T", i, read_seq); 
                        else: 
                            if read_seq[i] == 'C' and genome_seq[chr][pos+len(read_seq)] == "G": 
                                if random.random() <= meth_prob: 
                                    read_seq = strsub("T", i, read_seq); 

                # Write left and right reads
                print(">{}".format(cur_read_id), file=read_file)
                if swapped:
                    print(reverse_complement(read_seq), file=read_file)
                else:
                    print(read_seq, file=read_file)
                if paired_end:
                    print(">{}".format(cur_read_id), file=read2_file)
                    if swapped:
                        print(read2_seq, file=read2_file)
                    else:
                        print(reverse_complement(read2_seq), file=read2_file)
                
                # Transcript-based SAM
                print("{}\t{}\t{}\t{}\t255\t{}\t{}\t{}\t0\t{}\t*\tXM:i:{}\tNM:i:{}\tMD:Z:{}{}{}{}".format(cur_read_id, flag, transcript_id, pos + 1, cigar_str, transcript_id, pos2 + 1, read_seq, XM, NM, MD, Zs_str, XS, TO), file=sam_transcript_file)
                if paired_end:
                    print("{}\t{}\t{}\t{}\t255\t{}\t{}\t{}\t0\t{}\t*\tXM:i:{}\tNM:i:{}\tMD:Z:{}{}{}{}".format(cur_read_id, flag2, transcript_id, pos2 + 1, cigar2_str, transcript_id, pos + 1, read2_seq, XM2, NM2, MD2, Zs2_str, XS, TO2), file=sam_transcript_file)

                # Gene/Supertranscript-based SAM
                spos, scigars = convertSAMAlignmentFromTranscriptToGene(exons, pos, cigars, super_exons)
                spos2, scigars2 = convertSAMAlignmentFromTranscriptToGene(exons, pos2, cigars2, super_exons)
                scigar_str, scigar2_str = "".join(scigars), "".join(scigars2)
                if sanity_check:
                    SAMRepOk(super_transcript_seq, read_seq, gene_id, spos, scigar_str, XM, NM, MD, Zs, max_mismatch)
                    SAMRepOk(super_transcript_seq, read2_seq, gene_id, spos2, scigar2_str, XM2, NM2, MD2, Zs2, max_mismatch)

                print("{}\t{}\t{}\t{}\t255\t{}\t{}\t{}\t0\t{}\t*\tXM:i:{}\tNM:i:{}\tMD:Z:{}{}{}{}{}".format(cur_read_id, flag, gene_id, spos + 1, scigar_str, gene_id, spos2 + 1, read_seq, XM, NM, MD, Zs_str, XS, TI, TO), file=sam_gene_file)
                if paired_end:
                    print("{}\t{}\t{}\t{}\t255\t{}\t{}\t{}\t0\t{}\t*\tXM:i:{}\tNM:i:{}\tMD:Z:{}{}{}{}{}".format(cur_read_id, flag2, gene_id, spos2 + 1, scigar2_str, gene_id, spos + 1, read2_seq, XM2, NM2, MD2, Zs2_str, XS, TI2, TO2), file=sam_gene_file)

                # Genome-based SAM
                gpos, gcigars = convertSAMAlignmentFromTranscriptToGenome(exons, pos, cigars)
                gpos2, gcigars2 = convertSAMAlignmentFromTranscriptToGenome(exons, pos2, cigars2)
                gcigar_str, gcigar2_str = "".join(gcigars), "".join(gcigars2)
                if sanity_check:
                    SAMRepOk(chr_seq, read_seq, chr, gpos, gcigar_str, XM, NM, MD, Zs, max_mismatch)
                    SAMRepOk(chr_seq, read2_seq, chr, gpos2, gcigar2_str, XM2, NM2, MD2, Zs2, max_mismatch)

                print("{}\t{}\t{}\t{}\t255\t{}\t{}\t{}\t0\t{}\t*\tXM:i:{}\tNM:i:{}\tMD:Z:{}{}{}{}{}".format(cur_read_id, flag, chr, gpos + 1, gcigar_str, chr, gpos2 + 1, read_seq, XM, NM, MD, Zs_str, XS, TI, TO), file=sam_file)
                if paired_end:
                    print("{}\t{}\t{}\t{}\t255\t{}\t{}\t{}\t0\t{}\t*\tXM:i:{}\tNM:i:{}\tMD:Z:{}{}{}{}{}".format(cur_read_id, flag2, chr, gpos2 + 1, gcigar2_str, chr, gpos + 1, read2_seq, XM2, NM2, MD2, Zs2_str, XS, TI2, TO2), file=sam_file)

                cur_read_id += 1
            
    sam_file.close()
    sam_transcript_file.close()
    sam_gene_file.close()
    stat_file.close()
    read_file.close()
    if paired_end:
        read2_file.close()
        
        
"""
"""
def simulate_reads(genome_file, gtf_file, snp_file, base_fname,
                   rna, paired_end, read_len, frag_len,
                   num_frag, expr_profile_type, repeat_fname,
                   error_rate, max_mismatch,
                   random_seed, snp_prob,
                   bisulfite, meth_prob,
                   allele_specific, num_gene, num_tran,
                   sanity_check, verbose):
    random.seed(random_seed, version=1)
    err_rand_src = ErrRandomSource(error_rate / 100.0)
    
    if read_len > frag_len:
        frag_len = read_len

    genome_seq = read_genome(genome_file)
    snps = read_snp(snp_file)

    genes, transcripts = {}, {}
    if rna:
        genes, transcripts = read_transcript(genome_seq, gtf_file, frag_len, verbose)

    if sanity_check:
        sanity_check_input(genome_seq, genes, transcripts, snps, frag_len)

    if rna:
        simulate_RNA_reads(base_fname, genome_seq, genes, transcripts, snps,
                           paired_end, read_len, frag_len, num_frag,
                           error_rate, err_rand_src, max_mismatch, snp_prob, expr_profile_type,
                           bisulfite, meth_prob,
                           allele_specific, num_gene, num_tran,
                           sanity_check, verbose)
    else:
        simulate_DNA_reads(base_fname, repeat_fname, genome_seq, snps,
                           paired_end, read_len, frag_len, num_frag,
                           error_rate, err_rand_src, max_mismatch, snp_prob,
                           bisulfite, meth_prob,
                           allele_specific,
                           sanity_check, verbose)
        

"""
"""
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
    parser.add_argument('--repeat-info',
                        dest='repeat_fname',
                        action='store',
                        type=str,
                        default='',
                        help='repeat information filename')
    parser.add_argument('--error-rate',
                        dest='error_rate',
                        action='store',
                        type=float,
                        default=0.0,
                        help='per-base sequencing error rate (%%) (default: 0.0)')
    parser.add_argument('--max-mismatch',
                        dest='max_mismatch',
                        action='store',
                        type=int,
                        default=3,
                        help='max mismatches due to sequencing errors (default: 3)')
    parser.add_argument('--random-seed',
                        dest='random_seed',
                        action='store',
                        type=int,
                        default=0,
                        help='random seeding value (default: 0)')
    parser.add_argument('--snp-prob',
                        dest='snp_prob',
                        action='store',
                        type=float,
                        default=1.0,
                        help='probability of a read including a snp when the read spans the snp ranging from 0.0 to 1.0 (default: 1.0)')
    parser.add_argument('--bisulfite', 
                        dest='bisulfite',
                        action='store_true',
                        default=False,
                        help='Allows for the simulation of Bisulfite Sequencing Reads (WIP - For now, just randomly flips Cyotosines to Thymines If they are followed by a Guanine)')
    parser.add_argument('--meth-prob', 
                        dest='meth_prob',
                        action='store',
                        type=float, 
                        default=0.0,
                        help='The CpG location specification for methylation probability across all reads (default: 0.0)')
    parser.add_argument('-a', '--allele-specific', 
                        dest='allele_specific',
                        action='store_true',
                        default=False,
                        help='Allows for allele-specific gene expressions')
    parser.add_argument('--num-gene',
                        dest='num_gene',
                        action='store',
                        type=int,
                        default=5000,
                        help='number of genes (default: 5000)')
    parser.add_argument('--num-transcript',
                        dest='num_tran',
                        action='store',
                        type=int,
                        default=10000,
                        help='number of transcripts (default: 10000)')
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
    if not args.genome_file or not args.gtf_file or not args.snp_file:
        parser.print_help()
        exit(1)
    if not args.rna:
        args.expr_profile = "constant"
        
    simulate_reads(args.genome_file, args.gtf_file, args.snp_file, args.base_fname,
                   args.rna, args.paired_end, args.read_len, args.frag_len,
                   args.num_frag, args.expr_profile, args.repeat_fname,
                   args.error_rate, args.max_mismatch,
                   args.random_seed, args.snp_prob, 
                   args.bisulfite, args.meth_prob,
                   args.allele_specific, args.num_gene, args.num_tran,
                   args.sanity_check, args.verbose)

