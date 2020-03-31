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

import os, sys, math, random, re
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
def read_transcript(genome_seq, gtf_file, frag_len):
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
def getSNPs(chr_snps, left, right):
    low, high = 0, len(chr_snps)
    while low < high:
        mid = (low + high) // 2
        snpID, type, pos, data = chr_snps[mid]
        if pos < left:
            low = mid + 1
        else:
            high = mid - 1

    snps = []
    for snp in chr_snps[low:]:
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
def getSamAlignment(rna, exons, chr_seq, trans_seq, frag_pos, read_len, chr_snps, snp_prob, err_rand_src, max_mismatch):
    # Find the genomic position for frag_pos and exon number
    tmp_frag_pos, tmp_read_len = frag_pos, read_len
    pos, cigars, cigar_descs = exons[0][0], [], []
    e_pos = 0
    prev_e = None
    for e_i in range(len(exons)):
        e = exons[e_i]
        if prev_e:
            i_len = e[0] - prev_e[1] - 1
            pos += i_len
        e_len = e[1] - e[0] + 1
        if e_len <= tmp_frag_pos:
            tmp_frag_pos -= e_len
            pos += e_len
        else:
            pos += tmp_frag_pos
            e_pos = tmp_frag_pos
            break                        
        prev_e = e

    # Define Cigar and its descriptions
    assert e_i < len(exons)
    e_len = exons[e_i][1] - exons[e_i][0] + 1
    assert e_pos < e_len
    cur_pos = pos
    match_len = 0
    prev_e = None
    mismatch, remain_trans_len = 0, len(trans_seq) - (frag_pos + read_len)
    assert remain_trans_len >= 0
    for e_i in range(e_i, len(exons)):
        e = exons[e_i]
        if prev_e:
            i_len = e[0] - prev_e[1] - 1
            cur_pos += i_len
            cigars.append(("{}N".format(i_len)))
            cigar_descs.append([])
        tmp_e_left = e_left = e[0] + e_pos
        e_pos = 0

        # Retreive SNPs
        if rna:
            snps = getSNPs(chr_snps, e_left, e[1])
        else:
            snps = getSNPs(chr_snps, frag_pos, frag_pos + read_len)

        if snp_prob < 1.0 and len(snps) > 0:
            snps_ = []
            for snp in snps:
                if random.random() <= snp_prob:
                    snps_.append(snp)
            snps = snps_
            
        # Simulate mismatches due to sequencing errors
        mms = []
        for i in range(e_left, min(e[1], e_left + tmp_read_len - 1)):
            if err_rand_src.getRand() == 1:
                assert i < len(chr_seq)
                err_base = "A"
                #rand = random.randint(0, 2)
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
#        def diff_sort(a , b):
#            return a[2] - b[2]

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
            if e_left + tmp_read_len - 1 < diff_pos2 or e[1] < diff_pos2:
                break
            
            if diff_type == "single":
                if mismatch + 1 > max_mismatch:
                    continue
                cigar_descs[-1].append([diff_pos - tmp_e_left, diff_data, diff_id])
                tmp_e_left = diff_pos + 1
                mismatch += 1
            elif diff_type == "deletion":
                del_len = diff_data
                if mismatch + del_len > max_mismatch:
                    continue
                if len(cigars) <= 0 and diff_pos - e_left <= 0:
                    continue                
                if remain_trans_len < del_len:
                    continue
                remain_trans_len -= del_len
                if diff_pos - e_left > 0:
                    cigars.append("{}M".format(diff_pos - e_left))
                    cigar_descs[-1].append([diff_pos - tmp_e_left, "", ""])
                    cigar_descs.append([])
                cigars.append("{}D".format(del_len))
                cigar_descs[-1].append([0, del_len, diff_id])
                cigar_descs.append([])
                tmp_read_len -= (diff_pos - e_left)
                e_left = tmp_e_left = diff_pos + del_len

            elif diff_type == "insertion":
                ins_len = len(diff_data)
                if mismatch + ins_len > max_mismatch:
                    continue
                if len(cigars) <= 0 and diff_pos - e_left <= 0:
                    continue
                if e_left + tmp_read_len - 1 < diff_pos + ins_len:
                    break
                if diff_pos - e_left > 0:
                    cigars.append("{}M".format(diff_pos - e_left))
                    cigar_descs[-1].append([diff_pos - tmp_e_left, "", ""])
                    cigar_descs.append([])
                cigars.append("{}I".format(ins_len))
                cigar_descs[-1].append([0, diff_data, diff_id])
                cigar_descs.append([])
                tmp_read_len -= (diff_pos - e_left)
                tmp_read_len -= ins_len
                e_left = tmp_e_left = diff_pos

            else:
                assert False
            prev_diff = diff

        e_right = min(e[1], e_left + tmp_read_len - 1)
        e_len = e_right - e_left + 1
        remain_e_len = e_right - tmp_e_left + 1
        if remain_e_len > 0:
            cigar_descs[-1].append([remain_e_len, "", ""])
        if e_len < tmp_read_len:
            tmp_read_len -= e_len
            cigars.append(("{}M".format(e_len)))
        else:
            assert e_len == tmp_read_len
            cigars.append(("{}M".format(tmp_read_len)))
            tmp_read_len = 0
            break
        prev_e = e

    # Define MD, XM, NM, Zs, read_seq
    MD, XM, NM, Zs, read_seq = "", 0, 0, "", ""
    assert len(cigars) == len(cigar_descs)
    MD_match_len, Zs_match_len = 0, 0
    cur_trans_pos = frag_pos
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
                assert cur_trans_pos + add_match_len <= len(trans_seq)
                read_seq += trans_seq[cur_trans_pos:cur_trans_pos+add_match_len]
                cur_trans_pos += add_match_len
                if alt_base != "":
                    if MD_match_len > 0:
                        MD += ("{}".format(MD_match_len))
                        MD_match_len = 0
                    MD += trans_seq[cur_trans_pos]
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
                    cur_trans_pos += 1
        elif cigar_op == 'D':
            assert len(cigar_desc) == 1
            add_match_len, del_len, snp_id = cigar_desc[0]
            MD_match_len += add_match_len
            Zs_match_len += add_match_len
            if MD_match_len > 0:
                MD += ("{}".format(MD_match_len))
                MD_match_len = 0
            MD += ("^{}".format(trans_seq[cur_trans_pos:cur_trans_pos+cigar_len]))
            read_seq += trans_seq[cur_trans_pos:cur_trans_pos+add_match_len]
            if Zs != "":
                Zs += ","
            Zs += ("{}|D|{}".format(Zs_match_len, cigar_desc[0][-1]))
            Zs_match_len = 0
            cur_trans_pos += cigar_len
        elif cigar_op == 'I':
            assert len(cigar_desc) == 1
            add_match_len, ins_seq, snp_id = cigar_desc[0]
            ins_len = len(ins_seq)
            MD_match_len += add_match_len
            Zs_match_len += add_match_len
            read_seq += trans_seq[cur_trans_pos:cur_trans_pos+add_match_len]
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

    return pos, cigars, cigar_descs, MD, XM, NM, Zs, read_seq


"""
"""
cigar_re = re.compile('\d+\w')
def samRepOk(genome_seq, read_seq, chr, pos, cigar, XM, NM, MD, Zs, max_mismatch):
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

    if tMD != MD or tXM != XM or tNM != NM or XM > max_mismatch or XM != NM:
        print(chr, pos, cigar, MD, XM, NM, Zs, file=sys.stderr)
        print(tMD, tXM, tNM, file=sys.stderr)
        assert False
        
        
"""
"""
def simulate_reads(genome_file, gtf_file, snp_file, base_fname,
                   rna, paired_end, read_len, frag_len,
                   num_frag, expr_profile_type, repeat_fname,
                   error_rate, max_mismatch,
                   random_seed, snp_prob, sanity_check, verbose):
    random.seed(random_seed, version=1)
    err_rand_src = ErrRandomSource(error_rate / 100.0)
    
    if read_len > frag_len:
        frag_len = read_len

    genome_seq = read_genome(genome_file)
    if rna:
        genes, transcripts = read_transcript(genome_seq, gtf_file, frag_len)
    else:
        genes, transcripts = {}, {}
    snps = read_snp(snp_file)

    if sanity_check:
        sanity_check_input(genome_seq, genes, transcripts, snps, frag_len)

    if rna:
        num_transcripts = min(len(transcripts), 10000)
        expr_profile = generate_rna_expr_profile(expr_profile_type, num_transcripts)
    else:
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

    if rna:
        transcript_ids = sorted(list(transcripts.keys()))
        random.shuffle(transcript_ids, random=random.random)
        assert len(transcript_ids) >= len(expr_profile)
    else:
        chr_ids = list(genome_seq.keys())

    sam_file = open(base_fname + ".sam", "w")

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
        if rna:
            transcript_id = transcript_ids[t]
            chr, strand, transcript_len, exons = transcripts[transcript_id]
            print(transcript_id, t_num_frags, file=sys.stderr)
        else:
            chr = chr_ids[t]
            print(chr, t_num_frags, file=sys.stderr)

        assert chr in genome_seq
        chr_seq = genome_seq[chr]
        chr_len = len(chr_seq)
        if chr in repeat_loci:
            chr_repeat_loci = repeat_loci[chr]
        else:
            chr_repeat_loci = []
            
        if rna:
            t_seq = ""
            for e in exons:
                assert e[0] < e[1]
                t_seq += chr_seq[e[0]:e[1]+1]
            assert len(t_seq) == transcript_len
        else:
            t_seq = chr_seq
            exons = [[0, chr_len - 1]]

        if chr in snps:
            chr_snps = snps[chr]
        else:
            chr_snps = []

        for f in range(t_num_frags):
            if rna:
                #frag_pos = random.randint(0, transcript_len - frag_len)
                frag_pos = myrandint(0, transcript_len - frag_len)
            else:
                while True:
                    if len(chr_repeat_loci):
                        #locus_id = random.randint(0, len(chr_repeat_loci) - 1)
                        locus_id = myrandint(0, len(chr_repeat_loci) - 1)
                        frag_pos = chr_repeat_loci[locus_id][0]
                    else:
                        #frag_pos = random.randint(0, chr_len - frag_len)
                        frag_pos = myrandint(0, chr_len - frag_len)
                    if 'N' not in chr_seq[frag_pos:frag_pos + frag_len]:
                        break

            # SAM specification (v1.4)
            # http://samtools.sourceforge.net/
            flag, flag2 = 99, 163  # 83, 147
            pos, cigars, cigar_descs, MD, XM, NM, Zs, read_seq = getSamAlignment(rna, exons, chr_seq, t_seq, frag_pos, read_len, chr_snps, snp_prob, err_rand_src, max_mismatch)
            pos2, cigars2, cigar2_descs, MD2, XM2, NM2, Zs2, read2_seq = getSamAlignment(rna, exons, chr_seq, t_seq, frag_pos+frag_len-read_len, read_len, chr_snps, snp_prob, err_rand_src, max_mismatch)
            swapped = False
            if paired_end:
                #if random.randint(0, 1) == 1:
                if myrandint(0, 1) == 1:
                    swapped = True
                if swapped:
                    flag, flag2 = flag - 16, flag2 - 16
                    pos, pos2 = pos2, pos
                    cigars, cigars2 = cigars2, cigars
                    cigar_descs, cigar2_descs = cigar2_descs, cigar_descs
                    read_seq, read2_seq = read2_seq, read_seq
                    XM, XM2 = XM2, XM
                    NM, NM2 = NM2, NM
                    MD, MD2 = MD2, MD
                    Zs, Zs2 = Zs2, Zs

            cigar_str, cigar2_str = "".join(cigars), "".join(cigars2)
            if sanity_check:
                samRepOk(genome_seq, read_seq, chr, pos, cigar_str, XM, NM, MD, Zs, max_mismatch)
                samRepOk(genome_seq, read2_seq, chr, pos2, cigar2_str, XM2, NM2, MD2, Zs2, max_mismatch)

            if Zs != "":
                Zs = ("\tZs:Z:{}".format(Zs))
            if Zs2 != "":
                Zs2 = ("\tZs:Z:{}".format(Zs2))

            if rna:
                XS = "\tXS:A:{}".format(strand)
                TI = "\tTI:Z:{}".format(transcript_id)
            else:
                XS, TI = "", ""                

            print(">{}".format(cur_read_id), file=read_file)
            if swapped:
                print(reverse_complement(read_seq), file=read_file)
            else:
                print(read_seq, file=read_file)
            print("{}\t{}\t{}\t{}\t255\t{}\t{}\t{}\t0\t{}\t*\tXM:i:{}\tNM:i:{}\tMD:Z:{}{}{}{}".format(cur_read_id, flag, chr, pos + 1, cigar_str, chr, pos2 + 1, read_seq, XM, NM, MD, Zs, XS, TI), file=sam_file)
            if paired_end:
                print(">{}".format(cur_read_id), file=read2_file)
                if swapped:
                    print(read2_seq, file=read2_file)
                else:
                    print(reverse_complement(read2_seq), file=read2_file)
                print("{}\t{}\t{}\t{}\t255\t{}\t{}\t{}\t0\t{}\t*\tXM:i:{}\tNM:i:{}\tMD:Z:{}{}{}{}".format(cur_read_id, flag2, chr, pos2 + 1, cigar2_str, chr, pos + 1, read2_seq, XM2, NM2, MD2, Zs2, XS, TI), file=sam_file)

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
                   args.random_seed, args.snp_prob, args.sanity_check, args.verbose)
