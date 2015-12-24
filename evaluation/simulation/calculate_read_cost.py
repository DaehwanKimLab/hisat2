#!/usr/bin/env python

import sys, os, subprocess
import multiprocessing
import string, re
import platform
from datetime import datetime, date, time
import copy

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

            chr_name = line[1:-1]
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
def extract_splice_sites(gtf_fname):
    trans = {}

    gtf_file = open(gtf_fname)
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
        left, right = int(left), int(right)

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
        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, [[left, right]]]
        else:
            trans[transcript_id][2].append([left, right])

    gtf_file.close()
    
    # Sort exons and merge where separating introns are <=5 bps
    for tran, [chrom, strand, exons] in trans.items():
            exons.sort()
            tmp_exons = [exons[0]]
            for i in range(1, len(exons)):
                if exons[i][0] - tmp_exons[-1][1] <= 5:
                    tmp_exons[-1][1] = exons[i][1]
                else:
                    tmp_exons.append(exons[i])
            trans[tran] = [chrom, strand, tmp_exons]

    # Calculate and print the unique junctions
    junctions = set()
    for chrom, strand, exons in trans.values():
        for i in range(1, len(exons)):
            junctions.add(to_junction_str([chrom, exons[i-1][1], exons[i][0]]))

    return junctions


cigar_re = re.compile('\d+\w')


def to_junction_str(junction):
    return "%s-%d-%d" % (junction[0], junction[1], junction[2])


def to_junction(junction_str):
    fields = junction_str.split("-")
    if len(fields) > 3:
        chr, left, right = "-".join(fields[:-2]), fields[-2], fields[-1]
    else:
        assert len(fields) == 3
        chr, left, right = fields

    return [chr, int(left), int(right)]


def junction_cmp(a, b):
    if a[0] != b[0]:
        if a[0] < b[0]:
            return -1
        else:
            return 1

    if a[1] != b[1]:
        if a[1] < b[1]:
            return -1
        else:
            return 1

    if a[2] != b[2]:
        if a[2] < b[2]:
            return -1
        else:
            return 1

    return 0


# chr and pos are assumed to be integers
def get_junctions(chr, pos, cigar_str, min_anchor_len = 0, read_len = 100):
    junctions = []    
    right_pos = pos
    cigars = cigar_re.findall(cigar_str)
    cigars = [[int(cigars[i][:-1]), cigars[i][-1]] for i in range(len(cigars))]

    left_anchor_lens = []
    cur_left_anchor_len = 0
    for i in range(len(cigars)):
        length, cigar_op = cigars[i]
        if cigar_op in "MI":
            cur_left_anchor_len += length
        elif cigar_op == "N":
            assert cur_left_anchor_len > 0
            left_anchor_lens.append(cur_left_anchor_len)
            cur_left_anchor_len = 0
        
    for i in range(len(cigars)):
        length, cigar_op = cigars[i]
        if cigar_op == "N":
            left, right = right_pos - 1, right_pos + length

            if i > 0 and cigars[i-1][1] in "ID":
                if cigars[i-1][1] == "I":
                    left += cigars[i-1][0]
                else:
                    left -= cigars[i-1][0]
            if i + 1 < len(cigars) and cigars[i+1][1] in "ID":
                if cigars[i+1][1] == "I":
                    right -= cigars[i+1][0]
                else:
                    right += cigars[i+1][0]

            junction_idx = len(junctions)
            assert junction_idx < len(left_anchor_lens)
            left_anchor_len = left_anchor_lens[junction_idx]
            assert left_anchor_len > 0 and left_anchor_len < read_len
            right_anchor_len = read_len - left_anchor_len
            if left_anchor_len >= min_anchor_len and right_anchor_len >= min_anchor_len:
                junctions.append([chr, left, right])
                
        if cigar_op in "MND":
            right_pos += length

    return junctions

def get_right(pos, cigars):
    right_pos = pos
    cigars = cigar_re.findall(cigars)
    for cigar in cigars:
        length = int(cigar[:-1])
        cigar_op = cigar[-1]
        if cigar_op in "MDN":
            right_pos += length

    return right_pos

def get_cigar_chars(cigars):
    cigars = cigar_re.findall(cigars)
    cigar_chars = ""
    for cigar in cigars:
        cigar_op = cigar[-1]
        cigar_chars += cigar_op

    return cigar_chars

def get_cigar_chars_MN(cigars):
    cigars = cigar_re.findall(cigars)
    cigar_chars = ""
    for cigar in cigars:
        cigar_op = cigar[-1]
        if cigar_op in "MN":
            if cigar_chars == "" or cigar_chars[-1] != cigar_op:
                cigar_chars += cigar_op

    return cigar_chars

def is_small_anchor_junction_read(cigars):
    cigar_list = []
    for cigar in cigar_re.findall(cigars):
        cigar_op = cigar[-1]
        cigar_len = int(cigar[:-1])
        cigar_list.append([cigar_op, cigar_len])

    if len(cigar_list) < 3:
        return False

    if cigar_list[0][0] != 'M' or cigar_list[-1][0] != 'M':
        return False

    if cigar_list[0][1] > 10 and cigar_list[-1][1] > 10:
        return False

    if cigar_list[1][0] != 'N' or cigar_list[-2][0] != 'N':
        return False

    return True


def is_canonical_junction(chr_dic, junction):
    chr, left, right = junction
    donor = chr_dic[chr][left:left+2]
    acceptor = chr_dic[chr][right-3:right-1]
    rev_donor = reverse_complement(acceptor)
    rev_acceptor = reverse_complement(donor)

    if (donor == "GT" and acceptor == "AG") or \
            (rev_donor == "GT" and rev_acceptor == "AG"):
        return True

    return False


def is_small_exon_junction_read(cigars, min_exon_len = 23):
    cigars = cigar_re.findall(cigars)
    for i in range(1, len(cigars) - 1):
        cigar = cigars[i]
        cigar_op = cigar[-1]
        cigar_len = int(cigar[:-1])

        prev_op = cigars[i-1][-1]
        next_op = cigars[i+1][-1]
        
        if prev_op == 'N' and cigar_op == 'M' and next_op == 'N':
            if cigar_len <= min_exon_len:
                return True

    return False


def extract_single(infilename, outfilename, chr_dic, aligner):
    infile = open(infilename)
    outfile = open(outfilename, "w")
    for line in infile:
        if line[0] == '@':
            continue

        cols = line[:-1].split()
        read_id, flag, chr, pos, mapQ, cigar_str = cols[:6]
        read_seq = cols[9]
        if len(read_id) >= 3 and read_id[-2] == "/":
            read_id = read_id[:-2]

        if read_id.find("seq.") == 0:
            read_id = read_id[4:]

        if aligner == "gsnap":
            chr = chr.replace("_", ":")

        flag = int(flag)
        pos = int(pos)

        if flag & 0x4 != 0:
            continue

        NM = ""
        for i in range(11, len(cols)):
            col = cols[i]
            # "nM" from STAR
            if col[:2] == "NM" or col[:2] == "nM":
                NM = col
        assert NM != ""
        NM = int(NM[5:])

        read_pos, right_pos = 0, pos - 1
        cigars = cigar_re.findall(cigar_str)
        cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
        for i in range(len(cigars)):
            cigar_op, length = cigars[i]
            
            if cigar_op == "S":
                assert i == 0 or i == len(cigars) - 1
                if i == 0:
                    assert cigars[i+1][0] == "M"
                    ref_seq = chr_dic[chr][right_pos-length:right_pos]
                else:
                    assert cigars[i-1][0] == "M"
                    ref_seq = chr_dic[chr][right_pos:right_pos+length]

                ref_seq = ref_seq.upper()
                if length == len(ref_seq):
                    for j in range(length):
                        if ref_seq[j] != "N" and read_seq[read_pos+j] != ref_seq[j]:
                            NM += 1
                else:
                    NM += length

            if cigar_op in "MND":
                right_pos += length

            if cigar_op in "MIS":
                read_pos += length

        if cigars[0][0] == "S":
            assert cigars[1][0] == "M"
            pos -= cigars[0][1]
            cigars[1][1] += cigars[0][1]
            cigars = cigars[1:]
        if cigars[-1][0] == "S":
            assert cigars[-2][0] == "M"
            cigars[-2][1] += cigars[-1][1]
            cigars = cigars[:-1]

        cigar_str = ""
        for cigar in cigars:
            cigar_op, length = cigar
            cigar_str += ("%d%s" % (length, cigar_op))

        p_str = "%s\t%s\t%d\t%s\tNM:i:%d" % \
            (read_id, chr, pos, cigar_str, NM)
        
        print >> outfile, p_str

    outfile.close()
    infile.close()
    

def extract_pair(infilename, outfilename, chr_dic, aligner):
    read_dic = {}
    pair_reported = set()

    infile = open(infilename)
    outfile = open(outfilename, "w")
    for line in infile:
        if line[0] == '@':
            continue

        cols = line[:-1].split()
        read_id, flag, chr1, pos1, mapQ, cigar1_str, chr2, pos2 = cols[:8]
        read_seq = cols[9]
        if len(read_id) >= 3 and read_id[-2] == "/":
            read_id = read_id[:-2]

        if read_id.find("seq.") == 0:
            read_id = read_id[4:]

        if aligner == "gsnap":
            chr1 = chr1.replace("_", ":")
            chr2 = chr2.replace("_", ":")

        flag = int(flag)
        canonical_pos1, canonical_pos2 = int(pos1), int(pos2)
        pos1 = canonical_pos1

        if flag & 0x4 != 0:
            continue

        NM1 = ""
        for i in range(11, len(cols)):
            col = cols[i]
            # "nM" from STAR
            if col[:2] == "NM" or col[:2] == "nM":
                NM1 = col
        assert NM1 != ""
        NM1 = int(NM1[5:])

        read_pos, right_pos = 0, pos1 - 1
        cigars = cigar_re.findall(cigar1_str)
        cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
        for i in range(len(cigars)):
            cigar_op, length = cigars[i]
            
            if cigar_op == "S":
                assert i == 0 or i == len(cigars) - 1
                if i == 0:
                    assert cigars[i+1][0] == "M"
                    ref_seq = chr_dic[chr1][right_pos-length:right_pos]
                else:
                    assert cigars[i-1][0] == "M"
                    ref_seq = chr_dic[chr1][right_pos:right_pos+length]

                ref_seq = ref_seq.upper()
                if length == len(ref_seq):
                    for j in range(length):
                        if ref_seq[j] != "N" and read_seq[read_pos+j] != ref_seq[j]:
                            NM1 += 1
                else:
                    NM1 += length

            if cigar_op in "MND":
                right_pos += length

            if cigar_op in "MIS":
                read_pos += length

        if cigars[0][0] == "S":
            assert cigars[1][0] == "M"
            pos1 -= cigars[0][1]
            cigars[1][1] += cigars[0][1]
            cigars = cigars[1:]
        if cigars[-1][0] == "S":
            assert cigars[-2][0] == "M"
            cigars[-2][1] += cigars[-1][1]
            cigars = cigars[:-1]

        cigar1_str = ""
        for cigar in cigars:
            cigar_op, length = cigar
            cigar1_str += ("%d%s" % (length, cigar_op))

        if chr2 == '*':
            continue

        if chr2 == '=':
            chr2 = chr1

        me = "%s\t%s\t%d" % (read_id, chr1, canonical_pos1)
        partner = "%s\t%s\t%d" % (read_id, chr2, canonical_pos2)
        if partner in read_dic:
            maps = read_dic[partner]
            for map in maps:
                if map[0] == me:
                    cigar2_str, NM2, pos2 = map[1:4]
                    if int(pos2) > int(pos1):
                        p_str = "%s\t%s\t%d\t%s\t%s\t%d\t%s\tNM:i:%d\tNM:i:%d" % \
                                (read_id, chr1, pos1, cigar1_str, chr2, pos2, cigar2_str, NM1, NM2)
                    else:
                        p_str = "%s\t%s\t%d\t%s\t%s\t%d\t%s\tNM:i:%d\tNM:i:%d" % \
                                (read_id, chr2, pos2, cigar2_str, chr1, pos1, cigar1_str, NM2, NM1)

                    if p_str not in pair_reported:
                        pair_reported.add(p_str)
                        print >> outfile, p_str

        if not me in read_dic:
            read_dic[me] = []

        read_dic[me].append([partner, cigar1_str, NM1, pos1])

    outfile.close()
    infile.close()


def is_junction_read(junctions_dic, chr, pos, cigar_str):
    result_junctions = []
    junctions = get_junctions(chr, pos, cigar_str)
    for junction in junctions:
        junction_str = to_junction_str(junction)
        result_junctions.append([junction_str, junction_str in junctions_dic])

    return result_junctions


def is_junction_pair(junctions_dic, chr, pos, cigar_str, mate_chr, mate_pos, mate_cigar_str):
    junctions = is_junction_read(junctions_dic, chr, pos, cigar_str)
    mate_junctions = is_junction_read(junctions_dic, mate_chr, mate_pos, mate_cigar_str)
    junctions += mate_junctions
    return junctions


def find_in_gtf_junctions(chr_dic, gtf_junctions, junction, relax_dist = 5):
    def find_in_gtf_junctions(gtf_junctions, junction):
        l, u = 0, len(gtf_junctions)
        while l < u:
            m = (l + u) / 2
            assert m >= 0 and m < len(gtf_junctions)
            cmp_result = junction_cmp(junction, gtf_junctions[m])
            if cmp_result == 0:
                return m
            elif cmp_result < 0:
                u = m
            else:
                l = m + 1

        return l

    chr, left, right = junction
    gtf_index = find_in_gtf_junctions(gtf_junctions, [chr, left - relax_dist, right - relax_dist])

    if gtf_index >= 0:
        i = gtf_index
        while i < len(gtf_junctions):
            chr2, left2, right2 = gtf_junctions[i]
            if chr2 > chr or \
                    left2 - left > relax_dist or \
                    right2 - right > relax_dist:
                break

            if abs(left - left2) <= relax_dist and left - left2 == right - right2:
                test_small = ":" in chr
                if is_canonical_junction(chr_dic, gtf_junctions[i]):
                    if left == left2:
                        return i
                    else:
                        return -1
                else:
                    return i
            i += 1

    return -1

def compare_single_sam(RNA, reference_sam, query_sam, mapped_fname, chr_dic, gtf_junctions, gtf_junctions_set, ex_gtf_junctions):
    aligned, multi_aligned = 0, 0
    db_dic, db_junction_dic = {}, {}
    mapped_file = open(mapped_fname, "w")
    file = open(reference_sam, "r")
    junction_read_dic = {}
    for line in file:
        if line[0] == '@':
            continue

        read_name, chr, pos, cigar, NM = line[:-1].split()
        pos = int(pos)
                
        if read_name.find("seq.") == 0:
            read_name = read_name[4:]

        if len(read_name) > 2 and read_name[-2] == '/':
            read_name = read_name[:-2]

        multi_aligned += 1
        if read_name not in db_dic:
            db_dic[read_name] = []
            aligned += 1

        pos2 = get_right(pos, cigar)
        db_dic[read_name].append([chr, pos, pos2, cigar])

        read_junctions = is_junction_read(gtf_junctions_set, chr, pos, cigar)
        if len(read_junctions) > 0:
            if read_name not in db_junction_dic:
                db_junction_dic[read_name] = []

            for junction_str, is_gtf_junction in read_junctions:
                db_junction_dic[read_name].append([junction_str, is_gtf_junction])

                if junction_str not in junction_read_dic:
                    junction_read_dic[junction_str] = []
                junction_read_dic[junction_str].append(line[:-1])

    file.close()

    temp_junctions, temp_gtf_junctions = set(), set()
    for read_name, can_junctions in db_junction_dic.items():
        if len(can_junctions) <= 0:
            continue

        # daehwan - for debugging purposes
        # 1. select the best candidate among spliced alignments if multiple

        def pickup_junction(can_junctions):
            junctions = [can_junctions[0]]
            for i in range(1, len(can_junctions)):
                def get_intron_len(can_junction):
                    chr, left, right = to_junction(can_junction)
                    return right - left - 1

                intron, intron_cmp = get_intron_len(junctions[0][0]), get_intron_len(can_junctions[i][0])

                if intron > intron_cmp:
                    junctions = [can_junctions[i]]
                elif intron == intron_cmp:
                    junctions.append(can_junctions[i])

            return junctions

        # can_junctions = pickup_junction(can_junctions)

        for can_junction in can_junctions:
            found_junction_str = None
            junction_str, is_gtf_junction = can_junction
            if is_gtf_junction:
                found_junction_str = junction_str

            if not found_junction_str:
                junction = to_junction(junction_str)
                gtf_index = find_in_gtf_junctions(chr_dic, gtf_junctions, junction)

                if gtf_index >= 0:
                    is_gtf_junction = True
                    found_junction_str = to_junction_str(gtf_junctions[gtf_index])

            if found_junction_str:
                temp_gtf_junctions.add(found_junction_str)
                temp_junctions.add(found_junction_str)
            else:
                if junction_str not in temp_junctions:
                    None
                    # assert junction_str in junction_read_dic
                    # daehwan - for debugging purposes
                    """
                    if len(junction_read_dic[junction_str]) <= 2:
                        canonical = is_canonical_junction(chr_dic, to_junction(junction_str))
                        if not canonical:
                            print >> sys.stdout, read_name, junction_str, len(junction_read_dic[junction_str]), can_junctions
                            for line in junction_read_dic[junction_str]:
                                print >> sys.stdout, "\t", line
                    """
                temp_junctions.add(junction_str)


    temp2_junctions = []
    for junction in temp_junctions:
        temp2_junctions.append(to_junction(junction))
    temp_junctions = sorted(list(temp2_junctions), cmp=junction_cmp)
    temp2_junctions = []
    for can_junction in temp_junctions:
        if len(temp2_junctions) <= 0:
            temp2_junctions.append(can_junction)
        else:
            chr, left, right = temp2_junctions[-1]
            chr2, left2, right2 = can_junction
            if chr == chr2 and \
                    abs(left - left2) == abs(right - right2) and \
                    abs(left - left2) <= 10 and \
                    not to_junction_str(can_junction) in temp_gtf_junctions:
                continue
            temp2_junctions.append(can_junction)

    temp_junctions = set()
    for junction in temp2_junctions:
        temp_junctions.add(to_junction_str(junction))

    file = open(query_sam)
    mapped, unmapped, unique_mapped, mapping_point = 0, 0, 0, 0.0
    for line in file:
        if line[0] == '@':
            continue
        if RNA:
            read_name, chr, pos, cigar, trans_id, NM = line[:-1].split()
        else:
            read_name, chr, pos, cigar, NM = line[:-1].split()
        pos = int(pos)
        pos2 = get_right(pos, cigar)
        if read_name not in db_dic:
            unmapped += 1
            continue

        maps = db_dic[read_name]
        found = False
        if [chr, pos, pos2, cigar] in maps:
            found = True

        if not found:
            for map in maps:
                if chr == map[0] and \
                       pos == map[1] and \
                       pos2 == map[2] and \
                       get_cigar_chars(cigar) == get_cigar_chars(map[3]):

                    read_junctions = is_junction_read(gtf_junctions_set, map[0], map[1], map[3])
                    if True:
                        found_list = [False for i in range(len(read_junctions))]
                        for j in  range(len(read_junctions)):
                            junction_str, is_gtf_junction = read_junctions[j]
                            junction = to_junction(junction_str)
                            gtf_index = find_in_gtf_junctions(chr_dic, gtf_junctions, junction)
                            if gtf_index >= 0:
                                found_list[j] = True
                        found = not (False in found_list)
                    else:
                        found = False
                    break

        if found:
            print >> mapped_file, read_name
            mapped += 1
            if len(maps) == 1:
                unique_mapped += 1
            mapping_point += (1.0 / len(maps))
        else:
            unmapped += 1
            
    file.close()
    mapped_file.close()

    # daehwan - for debugging purposes
    false_can_junctions, false_noncan_junctions = 0, 0
    for junction_str in temp_junctions:
        if junction_str in temp_gtf_junctions:
            continue
        if junction_str in ex_gtf_junctions:
            continue
        if is_canonical_junction(chr_dic, to_junction(junction_str)):
            false_can_junctions += 1
        else:
            false_noncan_junctions += 1
    print >> sys.stderr, "\t\t\tfalse junctions: %d (canonical), %d (non-canonical)" % (false_can_junctions, false_noncan_junctions)
    
    return mapped, unique_mapped, unmapped, aligned, multi_aligned, len(temp_junctions), len(temp_gtf_junctions), mapping_point


def compare_paired_sam(RNA, reference_sam, query_sam, mapped_fname, chr_dic, gtf_junctions, gtf_junctions_set, ex_gtf_junctions):
    aligned, multi_aligned = 0, 0
    db_dic, db_junction_dic, junction_pair_dic = {}, {}, {}
    mapped_file = open(mapped_fname, "w")
    file = open(reference_sam, "r")
    for line in file:
        if line[0] == '@':
            continue
        read_name, chr, pos, cigar, chr2, pos2, cigar2, NM, NM2 = line[:-1].split()
        pos, pos2 = int(pos), int(pos2)

        if read_name.find("seq.") == 0:
            read_name = read_name[4:]

        if len(read_name) > 2 and read_name[-2] == '/':
            read_name = read_name[:-2]

        multi_aligned += 1        
        if read_name not in db_dic:
            db_dic[read_name] = []
            aligned += 1

        pos_right, pos2_right = get_right(pos, cigar), get_right(pos2, cigar2)
        db_dic[read_name].append([chr, pos, pos_right, cigar, pos2, pos2_right, cigar2])

        pair_junctions = is_junction_pair(gtf_junctions_set, chr, pos, cigar, chr2, pos2, cigar2)
        if len(pair_junctions) > 0:
            if read_name not in db_junction_dic:
                db_junction_dic[read_name] = []

            for junction_str, is_gtf_junction in pair_junctions:
                db_junction_dic[read_name].append([junction_str, is_gtf_junction])

                # daehwan - for debugging purposes
                if junction_str not in junction_pair_dic:
                    junction_pair_dic[junction_str] = []
                junction_pair_dic[junction_str].append(line[:-1])

    file.close()

    temp_junctions, temp_gtf_junctions = set(), set()
    for read_name, can_junctions in db_junction_dic.items():
        if len(can_junctions) <= 0:
            continue

        # daehwan - for debugging purposes
        # 1. select the best candidate among spliced alignments if multiple

        def pickup_junction(can_junctions):
            junctions = [can_junctions[0]]
            for i in range(1, len(can_junctions)):
                def get_intron_len(can_junction):
                    chr, left, right = to_junction(can_junction)
                    return right - left - 1

                intron, intron_cmp = get_intron_len(junctions[0][0]), get_intron_len(can_junctions[i][0])

                if intron > intron_cmp:
                    junctions = [can_junctions[i]]
                elif intron == intron_cmp:
                    junctions.append(can_junctions[i])

            return junctions

        # can_junctions = pickup_junction(can_junctions)

        for can_junction in can_junctions:
            found_junction_str = None
            junction_str, is_gtf_junction = can_junction

            # daehwan - for debugging purposes
            assert junction_str in junction_pair_dic
            if len(junction_pair_dic[junction_str]) <= 5:
                continue

            if is_gtf_junction:
                found_junction_str = junction_str

            if not found_junction_str:
                junction = to_junction(junction_str)
                gtf_index = find_in_gtf_junctions(chr_dic, gtf_junctions, junction)

                if gtf_index >= 0:
                    is_gtf_junction = True
                    found_junction_str = to_junction_str(gtf_junctions[gtf_index])

            if found_junction_str:
                temp_gtf_junctions.add(found_junction_str)
                temp_junctions.add(found_junction_str)
            else:
                if junction_str not in temp_junctions:
                    None
                    # assert junction_str in junction_read_dic
                    # print >> sys.stdout, read_name, junction_str, len(junction_read_dic[junction_str])
                    # for line in junction_read_dic[junction_str]:
                    #     print >> sys.stdout, "\t", line

                temp_junctions.add(junction_str)

    # daehwan - for debugging purposes
    filter_junction_db = {}

    temp2_junctions = []
    for junction in temp_junctions:
        temp2_junctions.append(to_junction(junction))
    temp_junctions = sorted(list(temp2_junctions), cmp=junction_cmp)
    temp2_junctions = []
    for can_junction in temp_junctions:
        if len(temp2_junctions) <= 0:
            temp2_junctions.append(can_junction)

            # daehwan - for debugging purposes
            # assert to_junction_str(can_junction) in junction_pair_dic
            # filter_junction_db[to_junction_str(can_junction)] = len(junction_pair_dic[to_junction_str(can_junction)])
        else:
            chr, left, right = temp2_junctions[-1]
            chr2, left2, right2 = can_junction
            if chr == chr2 and \
                    abs(left - left2) == abs(right - right2) and \
                    abs(left - left2) <= 10 and \
                    not to_junction_str(can_junction) in temp_gtf_junctions:
                
                # daehwan - for debugging purposes
                # assert to_junction_str(temp2_junctions[-1]) in junction_pair_dic
                # assert to_junction_str(temp2_junctions[-1]) in filter_junction_dic
                # filter_junction_db[to_junction_str(temp2_junctions[-1])] += len(junction_pair_dic[to_junction_str(temp2_junctions[-1])])
                
                continue

            temp2_junctions.append(can_junction)
            
            # daehwan - for debugging purposes
            # assert to_junction_str(can_junction) in junction_pair_dic
            # filter_junction_db[to_junction_str(can_junction)] = len(junction_pair_dic[to_junction_str(can_junction)])
            
    temp_junctions = set()
    for junction in temp2_junctions:
        # daehwan - for debugging purposes
        # assert to_junction_str(junction) in filter_junction_dic
        # if filter_junction_dic[to_junction_str(junction)] <= 5:
        #    continue
        
        temp_junctions.add(to_junction_str(junction))

    file = open(query_sam)
    mapped, unique_mapped, unmapped, mapping_point = 0, 0, 0, 0.0
    for line in file:
        if line[0] == '@':
            continue

        if RNA:
            read_name, chr, pos, cigar, chr2, pos2, cigar2, trans_id, NM, NM2 = line[:-1].split()
        else:
            read_name, chr, pos, cigar, chr2, pos2, cigar2, NM, NM2 = line[:-1].split()
        pos, pos2 = int(pos), int(pos2)
        pos_right, pos2_right = get_right(pos, cigar), get_right(pos2, cigar2)

        if read_name not in db_dic:
            unmapped += 1
            continue

        maps = db_dic[read_name]

        found = False
        if [chr, pos, pos_right, cigar, pos2, pos2_right, cigar2] in maps:
            found = True

        if not found:
            for map in maps:
                if chr == map[0] and \
                       pos == map[1] and \
                       pos_right == map[2] and \
                       get_cigar_chars(cigar) == get_cigar_chars(map[3]) and \
                       pos2 == map[4] and \
                       pos2_right == map[5] and \
                       get_cigar_chars(cigar2) == get_cigar_chars(map[6]):

                    pair_junctions = is_junction_pair(gtf_junctions_set, map[0], map[1], map[3], map[0], map[4], map[6])
                    if True:
                        found_list = [False for i in range(len(pair_junctions))]
                        for j in  range(len(pair_junctions)):
                            junction_str, is_gtf_junction = pair_junctions[j]
                            junction = to_junction(junction_str)
                            gtf_index = find_in_gtf_junctions(chr_dic, gtf_junctions, junction)
                            if gtf_index >= 0:
                                found_list[j] = True
                        found = not (False in found_list)
                    else:
                        found = False
                    break

        if found:
            print >> mapped_file, read_name
            mapped += 1
            if len(maps) == 1:
                unique_mapped += 1
            mapping_point += (1.0 / len(maps))
        else:
            unmapped += 1
            
    file.close()
    mapped_file.close()

    # daehwan - for debugging purposes
    false_can_junctions, false_noncan_junctions = 0, 0
    for junction_str in temp_junctions:
        if junction_str in temp_gtf_junctions:
            continue
        if is_canonical_junction(chr_dic, to_junction(junction_str)):
            false_can_junctions += 1
        else:
            false_noncan_junctions += 1
    print >> sys.stderr, "\t\t\tfalse junctions: %d (canonical), %d (non-canonical)" % (false_can_junctions, false_noncan_junctions)
        
    
    return mapped, unique_mapped, unmapped, aligned, multi_aligned, len(temp_junctions), len(temp_gtf_junctions), mapping_point


def extract_mapped_unmapped(read_fname, mapped_id_fname, mapped_fname, unmapped_fname, read2_fname = "", mapped2_fname = "", unmapped2_fname = ""):
    mapped_ids = set()
    mapped_id_file = open(mapped_id_fname)
    for line in mapped_id_file:
        read_id = int(line[:-1])
        mapped_ids.add(read_id)
                              
    mapped_id_file.close()

    def write_reads(read_fname, mapped_fname, unmapped_fname):
        mapped_file = open(mapped_fname, "w")
        unmapped_file = open(unmapped_fname, "w")
        read_file = open(read_fname)
        write = False
        for line in read_file:
            if line[0] == "@":
                read_id = int(line[1:-1])
                write = read_id in mapped_ids

            if write:
                print >> mapped_file, line[:-1]
            else:
                print >> unmapped_file, line[:-1]

        read_file.close()
        mapped_file.close()
        unmapped_file.close()

    write_reads(read_fname, mapped_fname, unmapped_fname)
    if read2_fname != "":
        assert mapped2_fname != ""
        assert unmapped2_fname != ""
        write_reads(read2_fname, mapped2_fname, unmapped2_fname)


def sql_execute(sql_db, sql_query):
    sql_cmd = [
        "sqlite3", sql_db,
        "-separator", "\t",
        "%s;" % sql_query
        ]
    # print >> sys.stderr, sql_cmd
    sql_process = subprocess.Popen(sql_cmd, stdout=subprocess.PIPE)
    output = sql_process.communicate()[0][:-1]
    return output


def create_sql_db(sql_db):
    if os.path.exists(sql_db):
        print >> sys.stderr, sql_db, "already exists!"
        return
    
    columns = [
        ["id", "integer primary key autoincrement"],
        ["genome", "text"],
        ["head", "text"],
        ["end_type", "text"],
        ["type", "text"],
        ["aligner", "text"],
        ["version", "text"],
        ["num_reads", "integer"],
        ["mapped_reads", "integer"],
        ["unique_mapped_reads", "integer"],
        ["unmapped_reads", "integer"],
        ["mapping_point", "real"],
        ["time", "real"],
        ["true_gtf_junctions", "integer"],
        ["temp_junctions", "integer"],
        ["temp_gtf_junctions", "integer"],
        ["host", "text"],
        ["created", "text"],
        ["cmd", "text"]
        ]
    
    sql_create_table = "CREATE TABLE ReadCosts ("
    for i in range(len(columns)):
        name, type = columns[i]
        if i != 0:
            sql_create_table += ", "
        sql_create_table += ("%s %s" % (name, type))
    sql_create_table += ");"
    sql_execute(sql_db, sql_create_table)


def write_analysis_data(sql_db, genome_name, database_name):
    if not os.path.exists(sql_db):
        return
    
    aligners = []
    sql_aligners = "SELECT aligner FROM ReadCosts GROUP BY aligner"
    output = sql_execute(sql_db, sql_aligners)
    aligners = output.split()

    can_read_types = ["all", "M", "2M_gt_15", "2M_8_15", "2M_1_7", "gt_2M"]    
    tmp_read_types = []
    sql_types = "SELECT type FROM ReadCosts GROUP BY type"
    output = sql_execute(sql_db, sql_types)
    tmp_read_types = output.split()

    read_types = []
    for read_type in can_read_types:
        if read_type in tmp_read_types:
            read_types.append(read_type)

    for paired in [False, True]:
        database_fname = genome_name + "_" + database_name
        if paired:
            end_type = "paired"
            database_fname += "_paired"
        else:
            end_type = "single"
            database_fname += "_single"
        database_fname += ".analysis"
        database_file = open(database_fname, "w")
        print >> database_file, "end_type\ttype\taligner\tnum_reads\ttime\tmapped_reads\tunique_mapped_reads\tunmapped_reads\tmapping_point\ttrue_gtf_junctions\ttemp_junctions\ttemp_gtf_junctions"
        for aligner in aligners:
            for read_type in read_types:
                sql_row = "SELECT end_type, type, aligner, num_reads, time, mapped_reads, unique_mapped_reads, unmapped_reads, mapping_point, true_gtf_junctions, temp_junctions, temp_gtf_junctions FROM ReadCosts"
                sql_row += " WHERE genome = '%s' and head = '%s' and aligner = '%s' and type = '%s' and end_type = '%s' ORDER BY created DESC LIMIT 1" % (genome_name, database_name, aligner, read_type, end_type)
                output = sql_execute(sql_db, sql_row)
                if output:
                    print >> database_file, output

        database_file.close()


def calculate_read_cost():
    sql_db_name = "analysis.db"
    if not os.path.exists(sql_db_name):
        create_sql_db(sql_db_name)

    num_cpus = multiprocessing.cpu_count()
    if num_cpus > 8:
        num_threads = min(8, num_cpus)
        desktop = False
    else:
        num_threads = min(3, num_cpus)
        desktop = True

    data_base = "sim"
    test_large_index = False
    verbose = False
    just_runtime = False
    sql_write = True
    aligners = [
        ["hisat2", "", "", ""],
        ["hisat2", "", "snp", ""],
        ["hisat2", "", "snp_tran", ""],
        ["hisat", "x1", "", ""],
        ["hisat", "", "", ""],
        ["hisat", "x2", "", ""],
        ["tophat2", "", "", ""],
        ["star", "", "", ""],
        ["star", "x2", "", ""],
        ["gsnap", "", "", ""],
        ["bowtie", "", "", ""],
        ["bowtie2", "", "", ""],
        ["bwa", "mem", "", ""]
        ]
    readtypes = ["all", "M", "2M_gt_15", "2M_8_15", "2M_1_7", "gt_2M"]
     
    aligners = [
        # ["hisat", "", "", ""],        
        # ["hisat2", "", "", ""],
        # ["hisat2", "x2", "", ""],
        # ["hisat2", "x1", "tran", ""],
        # ["hisat2", "", "tran", ""],
        # ["hisat2", "", "", "201b"],
        # ["hisat2", "", "", ""],
        # ["hisat2", "x1", "tran", "201b"],
        # ["hisat2", "x1", "tran", ""],
        # ["hisat2", "", "snp", "201b"],
        # ["hisat2", "", "snp", ""],
        ["hisat2", "x1", "snp_tran", "201b"],
        ["hisat2", "x1", "snp_tran", ""],
        # ["hisat2", "x1", "snp_tran_ercc", ""],
        # ["tophat2", "gtfonly", "", ""],
        # ["tophat2", "gtf", "", ""],
        # ["star", "", "", ""],
        # ["star", "gtf", "", ""],
        # ["bowtie", "", "", ""],
        # ["bowtie2", "", "", ""],
        # ["gsnap", "", "", ""],
        # ["hisat2", "", "snp", ""],
        # ["hisat2", "", "tran", ""],
        # ["hisat2", "", "snp_tran", ""],
        ]
    readtypes = ["all"]
    verbose = True
    debug = False
    # sql_write = False
    # just_runtime = True

    cwd = os.getcwd()
    if len(cwd.split("reads_")) > 1:
        genome = cwd.split("reads_")[1]
    else:
        genome = "genome"
    RNA = (cwd.find("RNA") != -1)

    test_small = (genome != "genome")

    if just_runtime:
        verbose = True

    chr_dic = read_genome("../../data/" + genome + ".fa")
    gtf_junctions = extract_splice_sites("../../data/genome.gtf")
    align_stat = []
    # for paired in [False, True]:
    for paired in [False]:
        for readtype in readtypes:
            if paired:
                base_fname = data_base + "_paired"
                type_sam_fname = base_fname + "_" + readtype + ".sam"
                type_read1_fname = base_fname +  "_1_" + readtype + ".fa"
                type_read2_fname = base_fname +  "_2_" + readtype + ".fa"
                type_junction_fname = base_fname + "_" + readtype + ".junc"
            else:
                base_fname = data_base + "_single"
                type_sam_fname = base_fname + "_" + readtype + ".sam"
                type_read1_fname = base_fname + "_" + readtype + ".fa"
                type_read2_fname = ""
                type_junction_fname = base_fname + "_" + readtype + ".junc"

            assert os.path.exists(type_sam_fname) and os.path.exists(type_junction_fname)
            numreads = 0
            type_sam_file = open(type_sam_fname)
            for line in type_sam_file:
                numreads += 1
            type_sam_file.close()
            if numreads <= 0:
                continue
            print >> sys.stderr, "%s\t%d" % (readtype, numreads)

            junctions, junctions_set = [], set()
            type_junction_file = open(type_junction_fname)
            for line in type_junction_file:
                chr, left, right = line[:-1].split()
                left, right = int(left), int(right)
                junctions.append([chr, left, right])
                junctions_set.add(to_junction_str([chr, left, right]))
                
            type_junction_file.close()

            aligner_bin_base = "../../../aligners/bin"
            def get_aligner_version(aligner, version):
                version = ""
                if aligner == "hisat2" or \
                        aligner == "hisat" or \
                        aligner == "bowtie" or \
                        aligner == "bowtie2":
                    if version:
                        cmd = ["%s/%s_%s/%s" % (aligner_bin_base, aligner, version, aligner)]
                    else:
                        cmd = ["%s/%s" % (aligner_bin_base, aligner)]
                    cmd += ["--version"]                    
                    cmd_process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                    version = cmd_process.communicate()[0][:-1].split("\n")[0]
                    version = version.split()[-1]
                elif aligner == "tophat2":
                    cmd = ["%s/tophat" % (aligner_bin_base)]
                    cmd += ["--version"]
                    cmd_process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                    version = cmd_process.communicate()[0][:-1].split()[-1]
                elif aligner in ["star", "starx2"]:
                    version = "2.4.2a"
                elif aligner == "gsnap":
                    cmd = ["%s/gsnap" % (aligner_bin_base)]
                    cmd_process = subprocess.Popen(cmd, stderr=subprocess.PIPE)
                    version = cmd_process.communicate()[1][:-1].split("\n")[0]
                    version = version.split()[2]
                elif aligner == "bwa":
                    cmd = ["%s/bwa" % (aligner_bin_base)]
                    cmd_process = subprocess.Popen(cmd, stderr=subprocess.PIPE)
                    version = cmd_process.communicate()[1][:-1].split("\n")[2]
                    version = version.split()[1]
                    
                return version

            index_base = "../../../indexes"
            index_add = ""
            if genome != "genome":
                index_add = "_" + genome
            def get_aligner_cmd(RNA, aligner, type, index_type, version, read1_fname, read2_fname, out_fname, cmd_idx = 0):
                cmd = []
                if aligner == "hisat2":
                    if version:
                        cmd = ["%s/hisat2_%s/hisat2" % (aligner_bin_base, version)]
                    else:
                        cmd = ["%s/hisat2" % (aligner_bin_base)]
                    if num_threads > 1:
                        cmd += ["-p", str(num_threads)]
                    cmd += ["-f"]
                    # cmd += ["-k", "5"]
                    # cmd += ["--score-min", "C,-18"]
                    # cmd += ["--pen-cansplice", "0"]
                    # cmd += ["--pen-noncansplice", "12"]
                    # cmd += ["--pen-intronlen", "G,-8,1"]
                    # cmd += ["--metrics", "1",
                    #        "--metrics-file", "metrics.out"]

                    if not RNA:
                        cmd += ["--no-spliced-alignment"]

                    if type in ["x1", "x2"]:
                        cmd += ["--no-temp-splicesite"]

                    # cmd += ["--no-anchorstop"]
                    # cmd += ["--sp", "1000,1000"]

                    # daehwan - for debugging purposes
                    # cmd += ["--dta"]
                    # cmd += ["--dta-cufflinks"]

                    if type == "x2":
                        if cmd_idx == 0:
                            cmd += ["--novel-splicesite-outfile"]
                        else:
                            cmd += ["--novel-splicesite-infile"]
                        cmd += ["splicesites.txt"]
                    
                    # "--novel-splicesite-infile",
                    # "../splicesites.txt",
                    # "--rna-strandness",
                    # "FR",
                    index_cmd = "%s/HISAT2%s/" % (index_base, index_add) + genome
                    if index_type:
                        index_cmd += ("_" + index_type)
                    cmd += [index_cmd]
                    if paired:
                        cmd += ["-1", read1_fname,
                                "-2", read2_fname]
                    else:
                        cmd += [read1_fname]                        
                elif aligner == "hisat":
                    cmd = ["%s/hisat" % (aligner_bin_base)]
                    if num_threads > 1:
                        cmd += ["-p", str(num_threads)]
                    cmd += ["-f"]
                    # cmd += ["-k", "5"]
                    # cmd += ["--score-min", "C,-18"]
                    if version != "":
                        version = int(version)
                    else:
                        version = sys.maxint

                    if not RNA:
                        cmd += ["--no-spliced-alignment"]

                    if type in ["x1", "x2"]:
                        cmd += ["--no-temp-splicesite"]

                    """
                    cmd += ["--rdg", "100,100",
                            "--rfg", "100,100"]
                    """

                    if type == "x2":
                        if cmd_idx == 0:
                            cmd += ["--novel-splicesite-outfile"]
                        else:
                            cmd += ["--novel-splicesite-infile"]
                        cmd += ["splicesites.txt"]
                    
                    # "--novel-splicesite-infile",
                    # "../splicesites.txt",
                    # "--rna-strandness",
                    # "FR",
                    cmd += ["%s/HISAT%s/" % (index_base, index_add) + genome]
                    if paired:
                        cmd += ["-1", read1_fname,
                                "-2", read2_fname]
                    else:
                        cmd += [read1_fname]                        
                elif aligner == "tophat2":
                    cmd = ["%s/tophat" % (aligner_bin_base)]
                    if num_threads > 1:
                        cmd += ["-p", str(num_threads)]
                    if type.find("gtf") != -1:
                        cmd += ["--transcriptome-index=%s/HISAT%s/gtf/%s" % (index_base, index_add, genome)]
                    if type == "gtfonly":
                        cmd += ["--transcriptome-only"]
                    cmd += ["--read-edit-dist", "3"]
                    cmd += ["--no-sort-bam"]
                    cmd += ["--read-realign-edit-dist", "0"]
                    cmd += ["--keep-tmp",
                            "%s/HISAT%s/" % (index_base, index_add) + genome,
                            read1_fname]
                    if paired:
                        cmd += [read2_fname]
                elif aligner == "star":
                    cmd = ["%s/STAR" % (aligner_bin_base)]
                    if num_threads > 1:
                        cmd += ["--runThreadN", str(num_threads)]
                    cmd += ["--genomeDir"]
                    if type == "gtf":
                        cmd += ["%s/STAR%s/gtf" % (index_base, index_add)]
                    else:
                        cmd += ["%s/STAR%s" % (index_base, index_add)]
                    if desktop:
                        cmd += ["--genomeLoad", "NoSharedMemory"]
                    else:
                        cmd += ["--genomeLoad", "LoadAndKeep"]
                    if type == "x2":
                        if cmd_idx == 1:
                            cmd += ["--alignSJDBoverhangMin", "1"]
                    cmd += ["--readFilesIn",
                            read1_fname]
                    if paired:
                        cmd += [read2_fname]
                    if paired:
                        cmd += ["--outFilterMismatchNmax", "6"]
                    else:
                        cmd += ["--outFilterMismatchNmax", "3"]
                elif aligner == "bowtie":
                    cmd = ["%s/bowtie" % (aligner_bin_base)]
                    if num_threads > 1:
                        cmd += ["-p", str(num_threads)]
                    cmd += ["-f", 
                            "--sam",
                            "-k", "10"]
                    cmd += ["-n", "3"]
                    if paired:
                        cmd += ["-X", "500"]
                    cmd += ["%s/Bowtie%s/" % (index_base, index_add) + genome]
                    if paired:
                        cmd += ["-1", read1_fname,
                                "-2", read2_fname]
                    else:
                        cmd += [read1_fname]
                elif aligner == "bowtie2":
                    cmd = ["%s/bowtie2" % (aligner_bin_base)]
                    if num_threads > 1:
                        cmd += ["-p", str(num_threads)]
                    cmd += ["-f", "-k", "10"]
                    cmd += ["--score-min", "C,-18"]
                    cmd += ["-x %s/HISAT%s/" % (index_base, index_add) + genome]
                    if paired:
                        cmd += ["-1", read1_fname,
                                "-2", read2_fname]
                    else:
                        cmd += [read1_fname]
                elif aligner == "gsnap":
                    cmd = ["%s/gsnap" % (aligner_bin_base),
                           "-A",
                           "sam"]
                    if num_threads > 1:
                        cmd += ["-t", str(num_threads)]
                    cmd += ["--max-mismatches=3",
                            "-D", "%s/GSNAP%s" % (index_base, index_add),
                            "-N", "1",
                            "-d", genome,
                            read1_fname]
                    if paired:
                        cmd += [read2_fname]
                elif aligner == "bwa":
                    cmd = ["%s/bwa" % (aligner_bin_base)]
                    if type in ["mem", "aln"]:
                        cmd += [type]
                    elif type == "sw":
                        cmd += ["bwa" + type]
                    if num_threads > 1:
                        cmd += ["-t", str(num_threads)]
                    cmd += ["%s/BWA%s/%s.fa" % (index_base, index_add, genome)]
                    cmd += [read1_fname]
                    if paired:
                        cmd += [read2_fname]
                    else:
                        assert False

                return cmd

            if test_small:
                init_time = {"hisat2" : 0.2, "hisat" : 0.1, "bowtie" : 0.1, "bowtie2" : 0.2, "gsnap" : 0.0, "bwa" : 0.0}
                if desktop:
                    init_time["star"] = 1.7
                else:
                    init_time["star"] = 0.0
            else:
                if desktop:
                    init_time = {"hisat2" : 3.0, "hisat" : 3.0, "bowtie" : 1.3, "bowtie2" : 1.9, "star" : 27.0, "gsnap" : 12, "bwa" : 1.3}
                else:
                    init_time = {"hisat2" : 9.5, "hisat" : 9.5, "bowtie" : 3.3, "bowtie2" : 4.1, "star" : 1.7, "gsnap" : 0.1, "bwa" : 3.3}
                    
            init_time["tophat2"] = 0.0
            for aligner, type, index_type, version in aligners:
                aligner_name = aligner + type
                if version:
                    aligner_name += ("_%s" % version)
                if aligner == "hisat2" and index_type != "":
                    aligner_name += ("_" + index_type)
                two_step = (aligner == "tophat2" or type == "x2" or (aligner in ["hisat2", "hisat"] and type == ""))
                if RNA and readtype != "M":
                    if aligner in ["bowtie", "bowtie2", "bwa"]:
                        continue
                if readtype != "all":
                    if two_step:
                        continue
                if not RNA and readtype != "all":
                    continue

                print >> sys.stderr, "\t%s\t%s" % (aligner_name, str(datetime.now()))
                if paired:
                    aligner_dir = aligner_name + "_paired"
                else:
                    aligner_dir = aligner_name + "_single"
                if not os.path.exists(aligner_dir):
                    os.mkdir(aligner_dir)
                os.chdir(aligner_dir)

                out_fname = base_fname + "_" + readtype + ".sam"
                out_fname2 = out_fname + "2"
                duration = -1.0
                if not os.path.exists(out_fname):
                    if not os.path.exists("../one.fa") or not os.path.exists("../two.fa"):
                        os.system("head -400 ../%s_1.fa > ../one.fa" % (data_base))
                        os.system("head -400 ../%s_2.fa > ../two.fa" % (data_base))

                    if just_runtime:
                        out_fname = "/dev/null"
                        out_fname2 = "/dev/null"

                    if not two_step:
                        align_stat.append([readtype, aligner_name])

                    # dummy commands for caching
                    if aligner != "tophat2":
                        dummy_cmd = get_aligner_cmd(RNA, aligner, type, index_type, version, "../one.fa", "../two.fa", "/dev/null")
                        if verbose:
                            print >> sys.stderr, "\t", datetime.now(), " ".join(dummy_cmd)
                        if aligner in ["hisat2", "hisat", "bowtie", "bowtie2", "gsnap", "bwa"]:
                            proc = subprocess.Popen(dummy_cmd, stdout=open("/dev/null", "w"), stderr=subprocess.PIPE)
                        else:
                            proc = subprocess.Popen(dummy_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        proc.communicate()
                        if verbose:
                            print >> sys.stderr, "\t", datetime.now(), "finished"

                    # align all reads
                    aligner_cmd = get_aligner_cmd(RNA, aligner, type, index_type, version, "../" + type_read1_fname, "../" + type_read2_fname, out_fname)
                    start_time = datetime.now()
                    if verbose:
                        print >> sys.stderr, "\t", start_time, " ".join(aligner_cmd)
                    if aligner in ["hisat2", "hisat", "bowtie", "bowtie2", "gsnap", "bwa"]:
                        proc = subprocess.Popen(aligner_cmd, stdout=open(out_fname, "w"), stderr=subprocess.PIPE)
                    else:
                        proc = subprocess.Popen(aligner_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    proc.communicate()
                    finish_time = datetime.now()
                    duration = finish_time - start_time
                    assert aligner in init_time
                    duration = duration.total_seconds() - init_time[aligner]
                    if duration < 0.1:
                        duration = 0.1
                    if verbose:
                        print >> sys.stderr, "\t", finish_time, "finished:", duration

                    if debug and aligner == "hisat2":
                        os.system("cat metrics.out")
                        print >> sys.stderr, "\ttime: %.4f" % (duration)

                    if aligner == "star" and type in ["", "gtf"]:
                        os.system("mv Aligned.out.sam %s" % out_fname)
                    elif aligner in ["hisat2", "hisat"] and type == "x2":
                        aligner_cmd = get_aligner_cmd(RNA, aligner, type, index_type, version, "../" + type_read1_fname, "../" + type_read2_fname, out_fname, 1)
                        start_time = datetime.now()
                        if verbose:
                            print >> sys.stderr, "\t", start_time, " ".join(aligner_cmd)
                        proc = subprocess.Popen(aligner_cmd, stdout=open(out_fname, "w"), stderr=subprocess.PIPE)
                        proc.communicate()
                        finish_time = datetime.now()
                        duration += (finish_time - start_time).total_seconds()
                        assert aligner in init_time
                        duration -= init_time[aligner]
                        if duration < 0.1:
                            duration = 0.1
                        if verbose:
                            print >> sys.stderr, "\t", finish_time, "finished:", duration
                    elif aligner == "star" and type == "x2":
                        assert os.path.exists("SJ.out.tab")
                        os.system("awk 'BEGIN {OFS=\"\t\"; strChar[0]=\".\"; strChar[1]=\"+\"; strChar[2]=\"-\";} {if($5>0){print $1,$2,$3,strChar[$4]}}' SJ.out.tab > SJ.out.tab.Pass1.sjdb")
                        for file in os.listdir("."):
                            if file in ["SJ.out.tab.Pass1.sjdb", "genome.fa"]:
                                continue
                            os.remove(file)
                        star_index_cmd = "STAR --genomeDir ./ --runMode genomeGenerate --genomeFastaFiles genome.fa --sjdbFileChrStartEnd SJ.out.tab.Pass1.sjdb --sjdbOverhang 99 --runThreadN %d" % (num_threads)
                        if verbose:
                            print >> sys.stderr, "\t", datetime.now(), star_index_cmd
                        os.system(star_index_cmd)
                        if verbose:
                            print >> sys.stderr, "\t", datetime.now(), " ".join(dummy_cmd)
                        proc = subprocess.Popen(dummy_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        proc.communicate()
                        if verbose:
                            print >> sys.stderr, "\t", datetime.now(), "finished"
                        aligner_cmd = get_aligner_cmd(RNA, aligner, type, index_type, version, "../" + type_read1_fname, "../" + type_read2_fname, out_fname, 1)
                        start_time = datetime.now()
                        if verbose:
                            print >> sys.stderr, "\t", start_time, " ".join(aligner_cmd)
                        proc = subprocess.Popen(aligner_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        proc.communicate()
                        finish_time = datetime.now()
                        duration += (finish_time - start_time).total_seconds()
                        assert aligner in init_time
                        duration -= init_time[aligner]
                        if duration < 0.1:
                            duration = 0.1
                        if verbose:
                            print >> sys.stderr, "\t", finish_time, "finished:", duration
                        os.system("mv Aligned.out.sam %s" % out_fname)
                    elif aligner == "tophat2":
                        os.system("samtools sort -n tophat_out/accepted_hits.bam accepted_hits; samtools view -h accepted_hits.bam > %s" % out_fname)

                    if aligner in ["gsnap", "tophat2"]:
                        os.system("tar cvzf %s.tar.gz %s &> /dev/null" % (out_fname, out_fname))

                if just_runtime:
                    os.chdir("..")
                    continue

                if not os.path.exists(out_fname2):
                    if paired:
                        extract_pair(out_fname, out_fname2, chr_dic, aligner)
                    else:
                        extract_single(out_fname, out_fname2, chr_dic, aligner)

                for readtype2 in readtypes:
                    if not two_step and readtype != readtype2:
                        continue

                    type_sam_fname2 = base_fname + "_" + readtype2 + ".sam"
                    if os.path.exists(type_sam_fname2 + ".done"):
                        continue
                    
                    if paired:
                        type_read_fname2 = base_fname + "_1_" + readtype2 + ".fa"
                    else:
                        type_read_fname2 = base_fname + "_" + readtype2 + ".fa"
                    mapped_id_fname = base_fname + "_" + readtype2 + ".read_id"
                    if paired:
                        mapped, unique_mapped, unmapped, aligned, multi_aligned, temp_junctions, temp_gtf_junctions, mapping_point = compare_paired_sam(RNA, out_fname2, "../" + type_sam_fname2, mapped_id_fname, chr_dic, junctions, junctions_set, gtf_junctions)
                    else:
                        mapped, unique_mapped, unmapped, aligned, multi_aligned, temp_junctions, temp_gtf_junctions, mapping_point = compare_single_sam(RNA, out_fname2, "../" + type_sam_fname2, mapped_id_fname, chr_dic, junctions, junctions_set, gtf_junctions)
                    proc = subprocess.Popen(["wc", "-l", "../" + type_read_fname2], stdout=subprocess.PIPE)
                    out = proc.communicate()[0]
                    numreads = int(out.split()[0]) / 2
                    assert mapped + unmapped == numreads
                    
                    if two_step:
                        print >> sys.stderr, "\t\t%s" % readtype2
                    print >> sys.stderr, "\t\taligned: %d, multi aligned: %d" % (aligned, multi_aligned)
                    print >> sys.stderr, "\t\tcorrectly mapped: %d (%.2f%%) mapping_point: %.2f" % (mapped, float(mapped) * 100.0 / numreads, mapping_point * 100.0 / numreads)
                    print >> sys.stderr, "\t\tuniquely and correctly mapped: %d (%.2f%%)" % (unique_mapped, float(unique_mapped) * 100.0 / numreads)
                    if readtype == readtype2:
                        print >> sys.stderr, "\t\t\t%d reads per sec (all)" % (numreads / max(1.0, duration))
                    print >> sys.stderr, "\t\tjunc. sensitivity %d / %d (%.2f%%), junc. accuracy: %d / %d (%.2f%%)" % \
                        (temp_gtf_junctions, len(junctions), float(temp_gtf_junctions) * 100.0 / max(1, len(junctions)), \
                             temp_gtf_junctions, temp_junctions, float(temp_gtf_junctions) * 100.0 / max(1, temp_junctions))

                    if duration > 0.0:
                        if sql_write and os.path.exists("../" + sql_db_name):
                            if paired:
                                end_type = "paired"
                            else:
                                end_type = "single"
                            sql_insert = "INSERT INTO \"ReadCosts\" VALUES(NULL, '%s', '%s', '%s', '%s', '%s', '%s', %d, %d, %d, %d, %f, %f, %d, %d, %d, '%s', datetime('now', 'localtime'), '%s');" % \
                                (genome, data_base, end_type, readtype2, aligner_name, get_aligner_version(aligner, version), numreads, mapped, unique_mapped, unmapped, mapping_point, duration, len(junctions), temp_junctions, temp_gtf_junctions, platform.node(), " ".join(aligner_cmd))
                            sql_execute("../" + sql_db_name, sql_insert)     

                        if two_step:
                            align_stat.append([readtype2, aligner_name, numreads, duration, mapped, unique_mapped, unmapped, mapping_point, len(junctions), temp_junctions, temp_gtf_junctions])
                        else:
                            align_stat[-1].extend([numreads, duration, mapped, unique_mapped, unmapped, mapping_point, len(junctions), temp_junctions, temp_gtf_junctions])

                    os.system("touch %s.done" % type_sam_fname2)                    

                os.chdir("..")

    print >> sys.stdout, "\t".join(["type", "aligner", "all", "all_time", "mapped", "unique_mapped", "unmapped", "mapping point", "true_gtf_junctions", "temp_junctions", "temp_gtf_junctions"])
    for line in align_stat:
        outstr = ""
        for item in line:
            if outstr != "":
                outstr += "\t"
            outstr += str(item)
        print >> sys.stdout, outstr

    if os.path.exists(sql_db_name):
        write_analysis_data(sql_db_name, genome, data_base)
        


if __name__ == "__main__":
    calculate_read_cost()
