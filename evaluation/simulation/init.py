#!/usr/bin/env python

import sys, os
import string, re

use_message = '''
'''

def extract_pair(RNA):
    read_dic = {}
    pair_reported = set()

    out_file = open("sim_paired.sam", "w")
    hits_file = open("sim.sam")
    for line in hits_file:
        if line[0] == '@':
            continue

        cols = line[:-1].split()
        read_id, flag, chr1, pos1, mapQ, cigar1, chr2, pos2 = cols[:8]
        if len(read_id) >= 3 and read_id[-2] == "/":
            read_id = read_id[:-2]

        if read_id.find("seq.") == 0:
            read_id = read_id[4:]

        flag = int(flag)
        pos1, pos2 = int(pos1), int(pos2)

        if flag & 0x4 != 0:
            continue

        TI, NM1 = "", ""
        for i in range(11, len(cols)):
            col = cols[i]
            if col[:2] == "TI":
                TI = "\t" + col[5:]
            # "nM" from STAR
            elif col[:2] == "NM" or col[:2] == "nM":
                NM1 = col
        assert NM1 != ""
        
        if chr2 == '*':
            continue

        if chr2 == '=':
            chr2 = chr1

        me = "%s\t%s\t%d" % (read_id, chr1, pos1)
        partner = "%s\t%s\t%d" % (read_id, chr2, pos2)
        if partner in read_dic:
            maps = read_dic[partner]
            for map in maps:
                if map[0] == me:
                    cigar2, NM2 = map[1:3]
                    if int(pos2) > int(pos1):
                        p_str = "%s\t%s\t%d\t%s\t%s\t%d\t%s%s\t%s\t%s" % \
                                (read_id, chr1, pos1, cigar1, chr2, pos2, cigar2, TI, NM1, NM2)
                    else:
                        p_str = "%s\t%s\t%d\t%s\t%s\t%d\t%s%s\t%s\t%s" % \
                                (read_id, chr2, pos2, cigar2, chr1, pos1, cigar1, TI, NM2, NM1)

                    if p_str not in pair_reported:
                        pair_reported.add(p_str)
                        print >> out_file, p_str

        if not me in read_dic:
            read_dic[me] = []

        read_dic[me].append([partner, cigar1, NM1])

        
    hits_file.close()
    out_file.close()


def init_reads(read_dir, RNA):
    if RNA:
        sam_cmd = "awk '{TI = $NF; NM = $13; if ($2 < 128) print $1\"\\t\"$3\"\\t\"$4\"\\t\"$6\"\\t\"substr(TI,6)\"\\t\"NM}' sim.sam > sim_1.sam"
    else:
        sam_cmd = "awk '{NM = $13; if ($2 < 128) print $1\"\\t\"$3\"\\t\"$4\"\\t\"$6\"\\t\"NM}' sim.sam > sim_1.sam"
    os.system(sam_cmd)
    if RNA:
        sam_cmd = "awk '{TI = $NF; NM = $13; if ($2 >= 128) print $1\"\\t\"$3\"\\t\"$4\"\\t\"$6\"\\t\"substr(TI,6)\"\\t\"NM}' sim.sam > sim_2.sam"
    else:
        sam_cmd = "awk '{NM = $13; if ($2 >= 128) print $1\"\\t\"$3\"\\t\"$4\"\\t\"$6\"\\t\"NM}' sim.sam > sim_2.sam"
    os.system(sam_cmd)
    extract_pair(RNA)


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

def classify_reads(RNA):
    if RNA:
        readtypes = ["all", "M", "2M_gt_15", "2M_8_15", "2M_1_7", "gt_2M"]
    else:
        readtypes = ["all" ,"M"]
        
    readtype_order = {}
    for i in range(len(readtypes)):
        readtype = readtypes[i]
        assert readtype not in readtype_order
        readtype_order[readtype] = i

    for paired in [False, True]:
        for readtype in readtypes:
            if paired:
                base_fname = "sim_paired"
                type_sam_fname = base_fname + "_" + readtype + ".sam"
                type_read1_fname = base_fname +  "_1_" + readtype + ".fa"
                type_read2_fname = base_fname +  "_2_" + readtype + ".fa"
                type_junction_fname = base_fname + "_" + readtype + ".junc"
            else:
                base_fname = "sim_single"
                type_sam_fname = base_fname + "_" + readtype + ".sam"
                type_read1_fname = base_fname + "_" + readtype + ".fa"
                type_read2_fname = ""
                type_junction_fname = base_fname + "_" + readtype + ".junc"
                
            if os.path.exists(type_sam_fname) and \
                    os.path.exists(type_junction_fname):
                continue
            
            read_ids = set()
            junctions = []
            type_sam_file = open(type_sam_fname, "w")
            if paired:
                sam_file = open("sim_paired.sam")
            else:
                sam_file = open("sim_1.sam")

            cigar_re = re.compile('\d+\w')
            def get_read_type(cigar):
                cigars = cigar_re.findall(cigar)

                def get_cigar_chars_MN(cigars):
                    cigars = cigar_re.findall(cigars)
                    cigar_chars = ""
                    for cigar in cigars:
                        cigar_op = cigar[-1]
                        if cigar_op in "MN":
                            if cigar_chars == "" or cigar_chars[-1] != cigar_op:
                                cigar_chars += cigar_op

                    return cigar_chars

                cigar_str = get_cigar_chars_MN(cigar)
                if cigar_str == "M":
                    return cigar_str
                elif cigar_str == "MNM":
                    assert len(cigars) >= 3
                    left_anchor = 0
                    for ci in range(len(cigars)):
                        c = cigars[ci][-1]
                        c_len = int(cigars[ci][:-1])
                        if c in "MI":
                            left_anchor += c_len
                        else:
                            break
                    assert left_anchor > 0
                    right_anchor = 0
                    for ci in reversed(range(len(cigars))):
                        c = cigars[ci][-1]
                        c_len = int(cigars[ci][:-1])
                        if c in "MI":
                            right_anchor += c_len
                        else:
                            break
                    assert right_anchor > 0                            
                    min_anchor = min(left_anchor, right_anchor)
                    if min_anchor > 15:
                        return "2M_gt_15"
                    elif min_anchor >= 8 and min_anchor <= 15:
                        return "2M_8_15"
                    elif min_anchor >= 1 and min_anchor <= 7:
                        return "2M_1_7"
                    else:
                        assert False
                else:
                    assert cigar_str not in ["M", "MNM"]
                    return "gt_2M"

            # chr and pos are assumed to be integers
            def get_junctions(chr, pos, cigar_str):
                junctions = []    
                right_pos = pos
                cigars = cigar_re.findall(cigar_str)
                cigars = [[int(cigars[i][:-1]), cigars[i][-1]] for i in range(len(cigars))]
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

                        junctions.append([chr, left, right])

                    if cigar_op in "MND":
                        right_pos += length

                return junctions

            for line in sam_file:
                if paired:
                    if RNA:
                        read_id, chr, pos, cigar, chr2, pos2, cigar2, trans_id, NM, NM2 = line[:-1].split()
                    else:
                        read_id, chr, pos, cigar, chr2, pos2, cigar2, NM, NM2 = line[:-1].split()
                else:
                    if RNA:
                        read_id, chr, pos, cigar, trans_id, NM = line[:-1].split()
                    else:
                        read_id, chr, pos, cigar, NM = line[:-1].split()

                read_id = int(read_id)
                readtype2 = get_read_type(cigar)
                if paired:
                    readtype3 = get_read_type(cigar2)
                    assert readtype2 in readtype_order
                    assert readtype3 in readtype_order
                    if readtype_order[readtype2] < readtype_order[readtype3]:
                        readtype2 = readtype3

                if readtype == "all" or readtype == readtype2:
                    read_ids.add(read_id)
                    print >> type_sam_file, line[:-1]
                    junctions += get_junctions(chr, int(pos), cigar)
                    if paired:
                        junctions += get_junctions(chr2, int(pos2), cigar2)

            sam_file.close()
            type_sam_file.close()

            # make this set non-redundant
            junction_set = set()
            for junction in junctions:
                junction_set.add(to_junction_str(junction))

            junctions = []
            for junction_str in junction_set:
                junctions.append(to_junction(junction_str))

            # sort the list of junctions
            junctions = sorted(junctions, cmp=junction_cmp)

            # write the junctions into a file
            type_junction_file = open(type_junction_fname, "w")
            for junction in junctions:
                print >> type_junction_file, "%s\t%d\t%d" % (junction[0], junction[1], junction[2])
            type_junction_file.close()

            def write_reads(read_fname, type_read_fname):
                type_read_file = open(type_read_fname, "w")
                read_file = open(read_fname)

                write = False
                for line in read_file:
                    if line[0] == ">":
                        read_id = int(line[1:-1])
                        write = read_id in read_ids

                    if write:
                        print >> type_read_file, line[:-1]

                read_file.close()
                type_read_file.close()

            if paired:
                write_reads("sim_1.fa", type_read1_fname)
                write_reads("sim_2.fa", type_read2_fname)
            else:
                write_reads("sim_1.fa", type_read1_fname)


def init():
    read_dir_base = "../reads/simulation/"
    read_dirs = os.listdir(read_dir_base)
    for read_dir in read_dirs:
        if os.path.exists(read_dir):
            continue
        
        if not os.path.exists(read_dir_base + read_dir + "/sim.sam") or \
                not os.path.exists(read_dir_base + read_dir + "/sim_1.fa") or \
                not os.path.exists(read_dir_base + read_dir + "/sim_2.fa"):
            continue

        print >> sys.stderr, "Processing", read_dir, "..."

        os.mkdir(read_dir)
        os.chdir(read_dir)
        os.system("ln -s ../%s%s/* ." % (read_dir_base, read_dir))
        os.system("ln -s sim_1.fa 1.fa")
        os.system("ln -s sim_2.fa 2.fa")

        RNA = (read_dir.find("RNA") != -1)
        init_reads(read_dir, RNA)
        classify_reads(RNA)            
        os.system("ln -s ../calculate_read_cost.py .")

        os.chdir("..")

    
if __name__ == "__main__":
    init()
