#!/usr/bin/env python

import sys, os
import string, re
from hisat_env import *
from analyze_common import *
use_message = '''
'''


def init_reads(read_prefix, num_reads):
    full_workdir = os.getcwd()
    workdir = full_workdir.split("/")[-1]

    files = [
        "sim.indels",
        "sim.mismatches",
        "sim.sam",
        "sim.simexpr",
        "sim.transfrags.gtf",
        "sim_1.fq",
        "sim_2.fq"
        ]

    for file in files:
        if os.path.exists("reads/" + file):
            continue        
        scp_cmd = "scp %s/simulation/%s/reads/%s reads" % (hisat_data_location, workdir, file)
        print >> sys.stderr, scp_cmd
        os.system(scp_cmd)

    os.chdir("reads")
    if not os.path.exists(read_prefix + "_1.fq") or not os.path.exists(read_prefix + "_2.fq"):
        head_cmd = "head -%d sim_1.fq > %s_1.fq" % (num_reads * 4, read_prefix)
        print >> sys.stderr, head_cmd
        os.system(head_cmd)
        head_cmd = "head -%d sim_2.fq > %s_2.fq" % (num_reads * 4, read_prefix)
        print >> sys.stderr, head_cmd
        os.system(head_cmd)

        shuffle_reads_cmd = "../../shuffle_reads.py %s_1.fq %s_2.fq" % (read_prefix, read_prefix)
        shuffle_reads_cmd += "; mv %s_1.fq.shuffle %s_1.fq" % (read_prefix, read_prefix)
        shuffle_reads_cmd += "; mv %s_2.fq.shuffle %s_2.fq" % (read_prefix, read_prefix)
        print >> sys.stderr, shuffle_reads_cmd
        os.system(shuffle_reads_cmd)

    if not os.path.exists(read_prefix + ".sam"):
        awk_cmd = "awk '{if($1 <= %d) print}' sim.sam > %s.sam" % (num_reads, read_prefix)
        print >> sys.stderr, awk_cmd
        os.system(awk_cmd)

    if not os.path.exists(read_prefix + "_paired.sam"):
        extract_pair_cmd = "../../extract_paired.py %s.sam > %s_paired.sam" % (read_prefix, read_prefix)
        print >> sys.stderr, extract_pair_cmd
        os.system(extract_pair_cmd)

    if not os.path.exists(read_prefix + "_1.sam"):
        sam_cmd = "awk '{if ($12 ~ /^XS/) {TI = $19; NM = $18} else {TI = $18; NM = $17} if ($2 %% 128 >= 64) print $1\"\\t\"$3\"\\t\"$4\"\\t\"$6\"\\t\"substr(TI,6)\"\\t\"NM}' %s.sam > %s_1.sam" % (read_prefix, read_prefix)
        print >> sys.stderr, sam_cmd
        os.system(sam_cmd)
    if not os.path.exists(read_prefix + "_2.sam"):
        sam_cmd = "awk '{if ($12 ~ /^XS/) {TI = $19; NM = $18} else {TI = $18; NM = $17} if ($2 %% 128 >= 128) print $1\"\\t\"$3\"\\t\"$4\"\\t\"$6\"\\t\"substr(TI,6)\"\\t\"NM}' %s.sam > %s_2.sam" % (read_prefix, read_prefix)
        print >> sys.stderr, sam_cmd
        os.system(sam_cmd)
    os.chdir("..")
    
    ln_cmd = "ln -s reads/%s_* ." % (read_prefix)
    print >> sys.stderr, ln_cmd
    os.system(ln_cmd)
    ln_cmd = "ln -s %s_1.fq 1.fq; ln -s %s_2.fq 2.fq" % (read_prefix, read_prefix)
    print >> sys.stderr, ln_cmd
    os.system(ln_cmd)
    

def init_indexes():
    aligners = ["STAR", "GSNAP", "OLego"]
    for aligner in aligners:
        if os.path.exists(aligner):
            continue
        ln_cmd = "ln -s ../../%s ." % (aligner)
        print >> sys.stderr, ln_cmd
        os.system(ln_cmd)

    ln_cmd = "ln -s ../genome* ../genes.gtf ."
    print >> sys.stderr, ln_cmd
    os.system(ln_cmd)
    

def link_python_codes():
    python_codes = ["calculate_read_cost.py"]
    for python_code in python_codes:
        ln_cmd = "ln -s ../%s ." % (python_code)
        print >> sys.stderr, ln_cmd
        os.system(ln_cmd)


def classify_reads(read_prefix):    
    readtypes = ["all", "M", "2M_gt_15", "2M_8_15", "2M_1_7", "gt_2M"]
    readtype_order = {}
    for i in range(len(readtypes)):
        readtype = readtypes[i]
        assert readtype not in readtype_order
        readtype_order[readtype] = i


    for paired in [False, True]:
    # for paired in [False]:
        for readtype in readtypes:
            if paired:
                base_fname = read_prefix + "_paired"
                type_sam_fname = base_fname + "_" + readtype + ".sam"
                type_read1_fname = base_fname +  "_1_" + readtype + ".fq"
                type_read2_fname = base_fname +  "_2_" + readtype + ".fq"
                type_junction_fname = base_fname + "_" + readtype + ".junc"
            else:
                base_fname = read_prefix + "_single"
                type_sam_fname = base_fname + "_" + readtype + ".sam"
                type_read1_fname = base_fname + "_" + readtype + ".fq"
                type_read2_fname = ""
                type_junction_fname = base_fname + "_" + readtype + ".junc"
                
            if os.path.exists(type_sam_fname) and \
                    os.path.exists(type_junction_fname):
                continue
            
            read_ids = set()
            junctions = []
            type_sam_file = open(type_sam_fname, "w")
            if paired:
                sam_file = open(read_prefix + "_paired.sam")
            else:
                sam_file = open(read_prefix + "_1.sam")

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
                    read_id, chr, pos, cigar, chr2, pos2, cigar2, trans_id, NM, NM2 = line[:-1].split()
                else:
                    read_id, chr, pos, cigar, trans_id, NM = line[:-1].split()

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
                    if line[0] == "@":
                        read_id = int(line[1:-1])
                        write = read_id in read_ids

                    if write:
                        print >> type_read_file, line[:-1]

                read_file.close()
                type_read_file.close()

            if paired:
                write_reads(read_prefix + "_1.fq", type_read1_fname)
                write_reads(read_prefix + "_2.fq", type_read2_fname)
            else:
                write_reads(read_prefix + "_1.fq", type_read1_fname)
