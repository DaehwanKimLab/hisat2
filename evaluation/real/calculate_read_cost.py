#!/usr/bin/env python

import sys, os, subprocess
import multiprocessing
import platform
import string
import re
from datetime import datetime, date, time

MAX_EDIT = 51

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
    chr, left, right = junction_str.split("-")
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


cigar_re = re.compile('\d+\w')

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

def is_non_canonical_junction_read(chr_dic, chr, left, cigars, canonical_junctions = [["GT", "AG"], ["GC", "AG"], ["AT", "AC"]]):
    pos = left
    for cigar in cigar_re.findall(cigars):
        cigar_op = cigar[-1]
        cigar_len = int(cigar[:-1])

        if cigar_op in 'MD':
            pos += cigar_len
        elif cigar_op == 'N':
            right = pos + cigar_len

            donor = chr_dic[chr][pos-1:pos+1]
            acceptor = chr_dic[chr][right-3:right-1]

            rev_donor = reverse_complement(acceptor)
            rev_acceptor = reverse_complement(donor)

            # print donor, acceptor
            # print rev_donor, rev_acceptor

            if [donor, acceptor] not in canonical_junctions and [rev_donor, rev_acceptor] not in canonical_junctions:
                return True

            pos = right
            
    return False

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

def is_junction_read(chr_dic, gtf_junctions, chr, pos, cigar_str):
    result_junctions = []
    junctions = get_junctions(chr, pos, cigar_str, 0, 101)
    for junction in junctions:
        junction_str = to_junction_str(junction)
        is_gtf_junction = False
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

        # allow small (<= 5bp) discrepancy for non-canonical splice sites.
        relaxed_junction_dist = 5
        chr, left, right = junction
        gtf_index = find_in_gtf_junctions(gtf_junctions, [chr, left - relaxed_junction_dist, right - relaxed_junction_dist])
        if gtf_index >= 0:
            i = gtf_index
            while i < len(gtf_junctions):
                chr2, left2, right2 = gtf_junctions[i]
                if chr2 > chr or \
                        left2 - left > relaxed_junction_dist or \
                        right2 - right > relaxed_junction_dist:
                    break

                if abs(left - left2) <= relaxed_junction_dist and left - left2 == right - right2:
                    canonical = is_canonical_junction(chr_dic, gtf_junctions[i])
                    if left == left2 or not canonical:
                        is_gtf_junction = True
                        break

                i += 1

        result_junctions.append([junction_str, is_gtf_junction])

    is_gtf_junction_read = False
    if len(result_junctions) > 0:
        is_gtf_junction_read = True
        for junction_str, is_gtf_junction in result_junctions:
            if not is_gtf_junction:
                is_gtf_junction_read = False
                break
        
    return result_junctions, len(result_junctions) > 0, is_gtf_junction_read


def is_junction_pair(chr_dic, gtf_junctions, chr, pos, cigar_str, mate_chr, mate_pos, mate_cigar_str):
    junctions, junction_read, gtf_junction_read = is_junction_read(chr_dic, gtf_junctions, chr, pos, cigar_str)
    mate_junctions, mate_junction_read, mate_gtf_junction_read = is_junction_read(chr_dic, gtf_junctions, mate_chr, mate_pos, mate_cigar_str)
    junctions += mate_junctions
    junction_pair = len(junctions) > 0
    if junction_pair:
        gtf_junction_pair = True
        if junction_read and not gtf_junction_read:
            gtf_junction_pair = False
        if mate_junction_read and not mate_gtf_junction_read:
            gtf_junction_pair = False
    else:
        gtf_junction_pair = False
        
    return junctions, junction_pair, gtf_junction_pair


def extract_reads_and_pairs(chr_dic, sam_filename, read_filename, pair_filename, unmapped_read_1_fq_name, unmapped_read_2_fq_name):
    temp_read_filename, temp_pair_filename = read_filename + ".temp", pair_filename + ".temp"
    temp_read_file, temp_pair_file = open(temp_read_filename, "w"), open(temp_pair_filename, "w")

    unmapped_read_1_fq, unmapped_read_2_fq = open(unmapped_read_1_fq_name, "w"), open(unmapped_read_2_fq_name, "w")
    hisat2 = read_filename.find("hisat2") != -1 or pair_filename.find("hisat2") != -1    

    read_dic = {}
    prev_read_id = ""
    sam_file = open(sam_filename, "r")
    for line in sam_file:
        if line[0] == "@":
            continue

        fields = line[:-1].split()
        read_id, flag, chr, pos, mapQ, cigar_str, mate_chr, mate_pos, template_len, read_seq, read_qual = fields[:11]

        flag, pos, mate_pos = int(flag), int(pos), int(mate_pos)
        read_seq = read_seq.upper()

        if flag & 0x04 != 0 or \
               chr == "*" or \
               cigar_str == "*":
            """
            if flag & 0x80 != 0:
                print >> unmapped_read_2_fq, "@%s\n%s\n+%s\n%s" % (read_id, read_seq, read_id, read_qual)
            else:
                print >> unmapped_read_1_fq, "@%s\n%s\n+%s\n%s" % (read_id, read_seq, read_id, read_qual)
            """
            continue

        if mate_chr == '=':
            mate_chr = chr

        if len(read_id) >= 3 and read_id[-2] == "/":
            read_id = read_id[:-2]

        if read_id.find("seq.") == 0:
            read_id = read_id[4:]

        if read_id != prev_read_id:
            read_dic = {}

        HISAT2_XM, HISAT2_NM = 0, 0
        if hisat2:
            for field in fields[11:]:
                if field[:5] == "XM:i:":
                    HISAT2_XM = int(field[5:])
                elif field[:5] == "NM:i:":
                    HISAT2_NM = int(field[5:])

        prev_read_id = read_id
        
        XM, gap = 0, 0
        read_pos, right_pos = 0, pos - 1,
        junction_read = False
        cigars = cigar_re.findall(cigar_str)
        for i in range(len(cigars)):
            cigar = cigars[i]
            length = int(cigar[:-1])
            cigar_op = cigar[-1]

            if cigar_op == "S":
                if i != 0 and i != len(cigars) - 1:
                    print >> sys.stderr, "S is located at %dth out of %d %s" % (i+1, len(cigars), cigar_str)

            if cigar_op in "MS":
                if cigar_op == "S" and i == 0:
                    ref_seq = chr_dic[chr][right_pos-length:right_pos]
                else:
                    ref_seq = chr_dic[chr][right_pos:right_pos+length]

                ref_seq = ref_seq.upper()
                if length == len(ref_seq):
                    for j in range(length):
                        if ref_seq[j] != "N" and read_seq[read_pos+j] != ref_seq[j]:
                            XM += 1
                            if hisat2 and cigar_op == "S":
                                HISAT2_XM += 1
                                HISAT2_NM += 1
                else:
                    XM += length

            if cigar_op in "MND":
                right_pos += length

            if cigar_op in "MIS":
                read_pos += length
                
            if cigar_op in "ID":
                gap += length

            if cigar_op == "N":
                junction_read = True
                
        NM = XM + gap
        if hisat2:
            XM, NM = HISAT2_XM, HISAT2_NM
        if NM < MAX_EDIT:
            print >> temp_read_file, "%s\t%d\t%s\t%s\t%s\tXM:i:%d\tNM:i:%d" % \
                  (read_id, flag, chr, pos, cigar_str, XM, NM)

            found = False
            me = "%s\t%s\t%d" % (read_id, chr, pos)
            partner = "%s\t%s\t%d" % (read_id, mate_chr, mate_pos)
            if partner in read_dic:
                maps = read_dic[partner]
                for map in maps:
                    if map[0] == me:
                        mate_flag, mate_cigar_str, mate_XM, mate_NM = map[1:]
                        if mate_pos > pos:
                            flag, chr, pos, cigar_str, XM, NM, mate_flag, mate_chr_str, mate_pos, mate_cigar_str, mate_XM, mate_NM = \
                                  mate_flag, mate_chr, mate_pos, mate_cigar_str, mate_XM, mate_NM, flag, chr, pos, cigar_str, XM, NM

                        print >> temp_pair_file, "%s\t%d\t%s\t%d\t%s\tXM:i:%d\tNM:i:%d\t%d\t%s\t%d\t%s\tXM:i:%d\tNM:i:%d" % \
                              (read_id, mate_flag, mate_chr, mate_pos, mate_cigar_str, mate_XM, mate_NM, flag, chr, pos, cigar_str, XM, NM)
                        found = True
                        break

            if not found:
                if not me in read_dic:
                    read_dic[me] = []

                read_dic[me].append([partner, flag, cigar_str, XM, NM])

    sam_file.close()

    temp_read_file.close()
    temp_pair_file.close()

    unmapped_read_1_fq.close()
    unmapped_read_2_fq.close()


    sort = False
    if sort:
        command = "sort %s | uniq > %s; rm %s" % (temp_read_filename, read_filename, temp_read_filename)
        os.system(command)
        
        command = "sort %s | uniq > %s; rm %s" % (temp_pair_filename, pair_filename, temp_pair_filename)
        os.system(command)
    else:
        command = "mv %s %s; mv %s %s" % (temp_read_filename, read_filename, temp_pair_filename, pair_filename)
        os.system(command)


def remove_redundant_junctions(junctions):
    temp_junctions = []
    for junction in junctions:
        temp_junctions.append(to_junction(junction))
    junctions = sorted(list(temp_junctions), cmp=junction_cmp)
    temp_junctions = []
    for can_junction in junctions:
        if len(temp_junctions) <= 0:
            temp_junctions.append(can_junction)
        else:
            chr, left, right = temp_junctions[-1]
            chr2, left2, right2 = can_junction
            if chr == chr2 and \
                    abs(left - left2) == abs(right - right2) and \
                    abs(left - left2) <= 10:
                continue
            temp_junctions.append(can_junction)
    junctions = set()
    for junction in temp_junctions:
        junctions.add(to_junction_str(junction))

    return junctions



def read_stat(read_filename, gtf_junctions, chr_dic = None, debug = False):
    read_stat = [[0, 0, 0] for i in range(MAX_EDIT)]
    temp_junctions = [set() for i in range(MAX_EDIT)]
    temp_gtf_junctions = [set() for i in range(MAX_EDIT)]

    alignment = []
    prev_read_id = ""
    read_file = open(read_filename, "r")
    for line in read_file:
        read_id, flag, chr, pos, cigar_str, XM, NM = line[:-1].split()
        flag, pos = int(flag), int(pos)
        XM, NM = int(XM[5:]), int(NM[5:])

        read_junctions, junction_read, gtf_junction_read = \
            is_junction_read(chr_dic, gtf_junctions, chr, pos, cigar_str)

        if junction_read:
            for junction_str, is_gtf_junction in read_junctions:
                if NM < len(temp_junctions):
                    temp_junctions[NM].add(junction_str)

                    if is_gtf_junction:
                        temp_gtf_junctions[NM].add(junction_str)

        if read_id != prev_read_id:
            if prev_read_id != "":
                NM2, junction_read2, gtf_junction_read2 = alignment
                if NM2 < len(read_stat):
                    read_stat[NM2][0] += 1
                    
                    if junction_read2:
                        read_stat[NM2][1] += 1

                        if gtf_junction_read2:
                            read_stat[NM2][2] += 1

            alignment = []

        prev_read_id = read_id

        if not alignment:
            alignment = [NM, junction_read, gtf_junction_read]
        elif alignment[0] > NM or \
                (alignment[0] == NM and not alignment[2] and junction_read):
            alignment = [NM, junction_read, gtf_junction_read]

    read_file.close()

    for i in range(len(read_stat)):
        temp_junctions[i] = remove_redundant_junctions(temp_junctions[i])
        temp_gtf_junctions[i] = remove_redundant_junctions(temp_gtf_junctions[i])

    for i in range(len(read_stat)):
        read_stat[i].append(len(temp_junctions[i]))
        read_stat[i].append(len(temp_gtf_junctions[i]))

    if alignment:
        NM2, junction_read2, gtf_junction_read2 = alignment
        if NM2 < len(read_stat):
            read_stat[NM2][0] += 1

            if junction_read2:
                read_stat[NM2][1] += 1

                if gtf_junction_read2:
                    read_stat[NM2][2] += 1

    return read_stat


def pair_stat(pair_filename, gtf_junctions, chr_dic):
    pair_stat = [[0, 0, 0] for i in range(MAX_EDIT)]
    dis_pair_stat = [0 for i in range(MAX_EDIT)]
    temp_junctions = [set() for i in range(MAX_EDIT)]
    temp_gtf_junctions = [set() for i in range(MAX_EDIT)]

    alignment, dis_alignments = [], []
    prev_read_id = ""
    pair_filename = open(pair_filename, "r")
    for line in pair_filename:
        read_id, flag, chr, pos, cigar_str, XM, NM, mate_flag, mate_chr, mate_pos, mate_cigar_str, mate_XM, mate_NM = line[:-1].split()
        flag, pos, XM, NM, mate_flag, mate_pos, mate_XM, mate_NM = \
             int(flag), int(pos), int(XM[5:]), int(NM[5:]), int(mate_flag), int(mate_pos), int(mate_XM[5:]), int(mate_NM[5:])

        pair_XM = XM + mate_XM
        pair_NM = NM + mate_NM

        pair_junctions, junction_pair, gtf_junction_pair = \
            is_junction_pair(chr_dic, gtf_junctions, chr, pos, cigar_str, mate_chr, mate_pos, mate_cigar_str)

        if junction_pair:
            for junction_str, is_gtf_junction in pair_junctions:
                if pair_NM < len(temp_junctions):
                    temp_junctions[pair_NM].add(junction_str)

                    if is_gtf_junction:
                        temp_gtf_junctions[pair_NM].add(junction_str)

        if read_id != prev_read_id:
            if prev_read_id != "":
                NM2, junction_read2, gtf_junction_read2 = alignment
                if NM2 < len(pair_stat):
                    pair_stat[NM2][0] += 1

                    if junction_read2:
                        pair_stat[NM2][1] += 1
                        if gtf_junction_read2:
                            pair_stat[NM2][2] += 1

            for NM2 in dis_alignments:
                if NM2 < len(dis_pair_stat):
                    dis_pair_stat[NM2] += 1
    
            alignment = []
            dis_alignment = []

        prev_read_id = read_id

        if not alignment:
            alignment = [pair_NM, junction_pair, gtf_junction_pair]
        elif alignment[0] > pair_NM or \
                (alignment[0] == pair_NM and not alignment[2] and junction_pair):
            alignment = [pair_NM, junction_pair, gtf_junction_pair]

        if mate_chr != chr or ((flag & 0x10) != 0 or (mate_flag & 0x10) == 0):
            if len(dis_alignments) == 0:
                dis_alignments = [pair_NM]
            elif dis_alignments[0] > pair_NM:
                dis_alignments = [pair_NM]

    pair_filename.close()

    if alignment:
        NM2, junction_read2, gtf_junction_read2 = alignment
        if NM2 < len(pair_stat):
            pair_stat[NM2][0] += 1

            if junction_read2:
                pair_stat[NM2][1] += 1
                if gtf_junction_read2:
                    pair_stat[NM2][2] += 1

    assert len(dis_alignments) <= 1
    for NM2 in dis_alignments:
        if NM2 < len(dis_pair_stat):
            dis_pair_stat[NM2] += 1

    for i in range(len(pair_stat)):
        temp_junctions[i] = remove_redundant_junctions(temp_junctions[i])
        temp_gtf_junctions[i] = remove_redundant_junctions(temp_gtf_junctions[i])

    for i in range(len(pair_stat)):
        pair_stat[i].append(len(temp_junctions[i]))
        pair_stat[i].append(len(temp_gtf_junctions[i]))                            

    return pair_stat, dis_pair_stat


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
        ["reads", "text"],
        ["genome", "text"],
        ["end_type", "text"],
        ["aligner", "text"],
        ["version", "test"],
        ["use_annotation", "text"],
        ["edit_distance", "integer"],
        ["mapped_reads", "integer"],
        ["junction_reads", "integer"],
        ["gtf_junction_reads", "integer"],
        ["junctions", "integer"],
        ["gtf_junctions", "integer"],
        ["runtime", "real"],
        ["host", "text"],
        ["created", "text"],
        ["cmd", "text"]
        ]
    
    sql_create_table = "CREATE TABLE Mappings ("
    for i in range(len(columns)):
        name, type = columns[i]
        if i != 0:
            sql_create_table += ", "
        sql_create_table += ("%s %s" % (name, type))
    sql_create_table += ");"
    sql_execute(sql_db, sql_create_table)


def write_analysis_data(sql_db, database_name, paired):
    if not os.path.exists(sql_db):
        return

    if paired:
        paired = "paired"
    else:
        paired = "single"
    
    aligners = []
    sql_aligners = "SELECT aligner FROM Mappings WHERE end_type = '%s' GROUP BY aligner" % (paired)
    output = sql_execute(sql_db, sql_aligners)
    aligners = output.split()

    database_fname = database_name + "_" + paired + ".analysis"
    database_file = open(database_fname, "w")

    print >> database_file, "aligner\tuse_annotation\tend_type\tedit_distance\tmapped_reads\tjunction_reads\tgtf_junction_reads\tjunctions\tgtf_junctions\truntime"
    for aligner in aligners:
        for edit_distance in range(MAX_EDIT):
            sql_row = "SELECT aligner, use_annotation, end_type, edit_distance, mapped_reads, junction_reads, gtf_junction_reads, junctions, gtf_junctions, runtime FROM Mappings"
            sql_row += " WHERE reads = '%s' and aligner = '%s' and edit_distance = %d and end_type = '%s' ORDER BY created DESC LIMIT 1" % (database_name, aligner, edit_distance, paired)
            output = sql_execute(sql_db, sql_row)
            if output:
                print >> database_file, output
            
    database_file.close()


def calculate_read_cost():
    sql_db_name = "analysis.db"
    if not os.path.exists(sql_db_name):
        create_sql_db(sql_db_name)

    full_workdir = os.getcwd()
    workdir = full_workdir.split("/")[-1]
    
    num_cpus = multiprocessing.cpu_count()
    if num_cpus > 8:
        num_threads = min(8, num_cpus)
        desktop = False
    else:
        num_threads = min(3, num_cpus)
        desktop = True

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
        ["bwa", "", "", ""]
        ]
    is_large_file = False
    assert os.path.exists("1.fq")
    if os.path.getsize("1.fq") > (2 * 1024 * 1024 * 1024):
        is_large_file = True

    aligners = [
        # ["hisat", "", "", ""],
        # ["hisat2", "", "", ""],
        # ["hisat2", "", "snp", ""],
        # ["hisat2", "", "tran", ""],
        ["hisat2", "x1", "snp", ""],
        ["bwa", "mem", "", ""],
        # ["bwa", "sw", "", ""],
        # ["hisat", "", "", ""],
        # ["star", "", "", ""],
        # ["star", "x2", "", ""],
        # ["tophat2", "", "", ""],
        # ["bowtie", "", "", ""],
        # ["bowtie2", "", "", ""],
        # ["hisat2", "x1", "", ""],
        # ["hisat2", "x2", "", ""],
        # ["hisat2", "", "tran", ""],
        # ["hisat2", "", "snp_tran", ""],
        ]

    # sql_write = False
    verbose = True
    debug = False

    genome = "genome"
    cwd = os.getcwd()    
    RNA = (cwd.find("RNA") != -1)

    chr_dic = read_genome("../../../data/" + genome + ".fa")
    gtf_junction_strs = extract_splice_sites("../../../data/genome.gtf")
    gene = "no"
    gtf_junctions = []
    for junction_str in gtf_junction_strs:
        junction = to_junction(junction_str)
        gtf_junctions.append(junction)
    gtf_junctions = sorted(gtf_junctions, cmp=junction_cmp)            

    print >> sys.stderr, "aligner\tuse_annotation\tend_type\tedit_distance\tmapped_reads\tjunction_reads\tgtf_junction_reads\tjunctions\tgtf_junctions\truntime"
    
    # for paired in [False, True]:
    for paired in [True]:
        type_read1_fname = "1.fq"
        if paired:
            type_read2_fname = "2.fq"
        else:
            type_read2_fname = ""

        aligner_bin_base = "../../../../aligners/bin"
        def get_aligner_version(aligner):
            version = ""
            if aligner == "hisat2" or \
                    aligner == "hisat" or \
                    aligner == "bowtie" or \
                    aligner == "bowtie2":
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

        index_base = "../../../../indexes"
        index_add = ""
        if genome != "genome":
            index_add = "_" + genome
        def get_aligner_cmd(RNA, aligner, type, index_type, version, read1_fname, read2_fname, out_fname, cmd_idx = 0):
            cmd = []
            if aligner == "hisat2":
                cmd = ["%s/hisat2" % (aligner_bin_base)]
                if num_threads > 1:
                    cmd += ["-p", str(num_threads)]
                    
                # cmd += ["-k", "5"]
                # cmd += ["--score-min", "C,-18"]

                # daehwan - for debugging purposes
                # cmd += ["--score-min", "C,-50"]
                    
                if version != "":
                    version = int(version)
                else:
                    version = sys.maxint
                # cmd += ["--pen-cansplice", "0"]
                # cmd += ["--pen-noncansplice", "12"]
                # cmd += ["--pen-intronlen", "G,-8,1"]
                # cmd += ["--metrics", "1",
                #         "--metrics-file", "metrics.out"]

                if not RNA:
                    cmd += ["--no-spliced-alignment"]

                if type in ["x1", "x2"] or not RNA:
                    cmd += ["--no-temp-splicesite"]

                # daehwan - for debugging purposes
                # cmd += ["--dta"]
                """
                if index_type == "tran":
                    cmd += ["--no-anchorstop"]
                    cmd += ["-k", "100"]
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
                # cmd += ["-k", "5"]
                # cmd += ["--score-min", "C,-18"]
                if version != "":
                    version = int(version)
                else:
                    version = sys.maxint

                if not RNA:
                    cmd += ["--no-spliced-alignment"]

                if type in ["x1", "x2"] or not RNA:
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
                if type == "x2" and cmd_idx == 1:
                    cmd += ["--genomeDir", "."]
                else:
                    cmd += ["--genomeDir", "%s/STAR%s" % (index_base, index_add)]
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
                cmd += ["--sam",
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
                cmd += ["-k", "10"]
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

        init_time = {"hisat2" : 0.0, "hisat" : 0.0, "bowtie" : 0.0, "bowtie2" : 0.0, "star" : 0.0, "olego" : 0.0, "gsnap" : 0.0, "tophat2" : 0.0, "bwa" : 0.0}
        if not is_large_file:
            if desktop:
                init_time = {"hisat2" : 3.0, "hisat" : 3.0, "bowtie" : 1.3, "bowtie2" : 1.9, "star" : 27.0, "gsnap" : 1, "bwa" : 1.3}
            else:
                init_time = {"hisat2" : 9.5, "hisat" : 9.5, "bowtie" : 3.3, "bowtie2" : 4.1, "star" : 1.7, "gsnap" : 0.1, "bwa" : 3.3}
        init_time["tophat2"] = 0.0
                    
        for aligner, type, index_type, version in aligners:
            aligner_name = aligner + type + version
            if aligner == "hisat2" and index_type != "":
                aligner_name += ("_" + index_type)
            if paired:
                aligner_dir = aligner_name + "_paired"
            else:
                aligner_dir = aligner_name + "_single"
            if not os.path.exists(aligner_dir):
                os.mkdir(aligner_dir)
            os.chdir(aligner_dir)

            out_fname = "accepted.sam"
            aligner_cmd = get_aligner_cmd(RNA, aligner, type, index_type, version, "../" + type_read1_fname, "../" + type_read2_fname, out_fname)
            duration = 0.1
            if not os.path.exists(out_fname):
                if not os.path.exists("../one.fq") or not os.path.exists("../two.fq"):
                    os.system("head -400 ../1.fq > ../one.fq")
                    os.system("head -400 ../2.fq > ../two.fq")

                # dummy commands for caching index
                if aligner not in ["tophat2"]:
                    for i in range(2):
                        dummy_cmd = get_aligner_cmd(RNA, aligner, type, index_type, version, "../one.fq", "../two.fq", "/dev/null")
                        start_time = datetime.now()
                        if verbose:
                            print >> sys.stderr, start_time, "\t", " ".join(dummy_cmd)
                        if aligner in ["hisat2", "hisat", "bowtie", "bowtie2", "gsnap", "bwa"]:
                            proc = subprocess.Popen(dummy_cmd, stdout=open("/dev/null", "w"), stderr=subprocess.PIPE)
                        else:
                            proc = subprocess.Popen(dummy_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        proc.communicate()
                        finish_time = datetime.now()
                        duration = finish_time - start_time
                        duration = duration.total_seconds()
                        if verbose:
                            print >> sys.stderr, finish_time, "duration:", duration

                # align all reads
                if paired:
                    sweep_read_cmd = "cat ../%s ../%s > /dev/null" % (type_read1_fname, type_read2_fname)
                else:
                    sweep_read_cmd = "cat ../%s > /dev/null" % (type_read1_fname)
                print >> sys.stderr, datetime.now(), "\t", sweep_read_cmd
                os.system(sweep_read_cmd)

                skip_alignment = False
                if paired and aligner == "olego" and os.path.exists(out_fname + "1"):
                    skip_alignment = True

                if not skip_alignment:
                    aligner_cmd = get_aligner_cmd(RNA, aligner, type, index_type, version, "../" + type_read1_fname, "../" + type_read2_fname, out_fname)
                    start_time = datetime.now()
                    if verbose:
                        print >> sys.stderr, start_time, "\t", " ".join(aligner_cmd)
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
                        print >> sys.stderr, finish_time, "duration:", duration

                    if debug and aligner == "hisat" and type == "x1":
                        os.system("cat metrics.out")
                        print >> sys.stderr, "\ttime: %.4f" % (duration)
                        # break

                if aligner == "star" and type in ["", "gtf"]:
                    os.system("mv Aligned.out.sam %s" % out_fname)
                elif aligner in ["hisat2", "hisat"] and type == "x2":
                    aligner_cmd = get_aligner_cmd(RNA, aligner, type, index_type, version, "../" + type_read1_fname, "../" + type_read2_fname, out_fname, 1)
                    if verbose:
                        print >> sys.stderr, start_time, "\t", " ".join(aligner_cmd)
                    start_time = datetime.now()
                    proc = subprocess.Popen(aligner_cmd, stdout=open(out_fname, "w"), stderr=subprocess.PIPE)
                    proc.communicate()
                    finish_time = datetime.now()
                    duration += (finish_time - start_time).total_seconds()
                    assert aligner in init_time
                    duration -= init_time[aligner]
                    if duration < 0.1:
                        duration = 0.1
                    if verbose:
                        print >> sys.stderr, finish_time, "duration:", duration
                elif aligner == "star" and type == "x2":
                    assert os.path.exists("SJ.out.tab")
                    os.system("awk 'BEGIN {OFS=\"\t\"; strChar[0]=\".\"; strChar[1]=\"+\"; strChar[2]=\"-\";} {if($5>0){print $1,$2,$3,strChar[$4]}}' SJ.out.tab > SJ.out.tab.Pass1.sjdb")
                    for file in os.listdir("."):
                        if file in ["SJ.out.tab.Pass1.sjdb", "genome.fa"]:
                            continue
                        os.remove(file)
                    star_index_cmd = "STAR --genomeDir ./ --runMode genomeGenerate --genomeFastaFiles ../../../../data/genome.fa --sjdbFileChrStartEnd SJ.out.tab.Pass1.sjdb --sjdbOverhang 100 --runThreadN %d" % (num_threads)
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

                os.system("echo %s %s %f >> runtime" % (str(datetime.now()), aligner, duration))

                if aligner in ["star", "tophat2", "gsnap"]:
                    os.system("tar cvzf %s.tar.gz %s &> /dev/null" % (out_fname, out_fname))

            if just_runtime:
                os.chdir("..")
                continue

            suffix = aligner
            read_sam, pair_sam = suffix + ".read.sam", suffix + ".pair.sam"
            unmapped_read_1_fq, unmapped_read_2_fq = suffix + ".unmapped.1.fq", suffix + ".unmapped.2.fq"
            if not os.path.exists(read_sam) or not os.path.exists(pair_sam):
                extract_reads_and_pairs(chr_dic, out_fname, read_sam, pair_sam, unmapped_read_1_fq, unmapped_read_2_fq)

            done_filename = suffix + ".done"
            if not os.path.exists(done_filename):
                done_file = open(done_filename, "w")
                if paired:
                    sum, dis_sum = [0, 0, 0, 0, 0], 0
                    stat, dis_stat = pair_stat(pair_sam, gtf_junctions, chr_dic)
                    output = ""
                    for i in range(len(stat)):
                        for j in range(len(sum)):
                            sum[j] += stat[i][j]

                        dis_sum += dis_stat[i]
                        mapped_reads, junction_reads, gtf_junction_reads, num_junctions, num_gtf_junctions = sum
                        output += "%s\t%s\tpaired\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n" % (aligner_name, gene, i, mapped_reads, junction_reads, gtf_junction_reads, num_junctions, num_gtf_junctions, duration)

                        if sql_write and os.path.exists("../" + sql_db_name):
                            sql_insert = "INSERT INTO \"Mappings\" VALUES(NULL, '%s', '%s', '%s', '%s', '%s', '%s', %d, %d, %d, %d, %d, %d, %f, '%s', datetime('now', 'localtime'), '%s');" % \
                                    (workdir, genome, "paired", aligner_name, get_aligner_version(aligner), "no", i, mapped_reads, junction_reads, gtf_junction_reads, num_junctions, num_gtf_junctions, duration, platform.node(), " ".join(aligner_cmd))
                            sql_execute("../" + sql_db_name, sql_insert)     
                    

                    print >> sys.stderr, output,
                    print >> done_file, output
                else:
                    sum = [0, 0, 0, 0, 0]
                    stat = read_stat(read_sam, gtf_junctions, chr_dic)
                    output = ""
                    for i in range(len(stat)):
                        for j in range(len(sum)):
                            sum[j] += stat[i][j]

                        mapped_reads, junction_reads, gtf_junction_reads, num_junctions, num_gtf_junctions = sum
                        output += "%s\t%s\tsingle\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n" % (aligner_name, gene, i, mapped_reads, junction_reads, gtf_junction_reads, num_junctions, num_gtf_junctions, duration)

                        if sql_write and os.path.exists("../" + sql_db_name):
                            sql_insert = "INSERT INTO \"Mappings\" VALUES(NULL, '%s', '%s', '%s', '%s', '%s', '%s', %d, %d, %d, %d, %d, %d, %f, '%s', datetime('now', 'localtime'), '%s');" % \
                                    (workdir, genome, "single", aligner_name, get_aligner_version(aligner), "no", i, mapped_reads, junction_reads, gtf_junction_reads, num_junctions, num_gtf_junctions, duration, platform.node(), " ".join(aligner_cmd))
                            sql_execute("../" + sql_db_name, sql_insert)                    
                        
                    print >> sys.stderr, output,
                    print >> done_file, output
                    
                done_file.close()


            os.chdir("..")

        if os.path.exists(sql_db_name):
            write_analysis_data(sql_db_name, workdir, paired)



if __name__ == "__main__":
    calculate_read_cost()
