#!/usr/bin/env python
#
# Copyright 2017, Daehwan Kim <infphilo@gmail.com> and Bongsoo Park <genomicspark@gmail.com>
#
# This file is part of HISAT-bisulfite.
#
# HISAT-bisulfite is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT-bisulfite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT-bisulfite.  If not, see <http://www.gnu.org/licenses/>.
#

import sys, os, subprocess
import multiprocessing
import string, re
import platform
from datetime import datetime, date, time
import copy
from argparse import ArgumentParser, FileType
import hisatgenotype_typing_common as typing_common


def simulate_reads(genome,
                   sim_name,
                   cpg,
                   mismatch,
                   numreads):
    if not os.path.exists("%s.sam" % sim_name) or \
            not os.path.exists("%s_1.fa" % sim_name) or \
            not os.path.exists("%s_2.fa" % sim_name):
        if cpg:
            cpg_fname = "%s.sense.snp" % genome
        else:
            cpg_fname = "/dev/null"
        cmd_add = ""
        if mismatch:
            cmd_add += "--error-rate 0.5 "
        cmd = "hisat2_simulate_reads.py --sanity-check --dna %s --num-fragment %d %s.fa /dev/null %s %s" % \
            (cmd_add, numreads, genome, cpg_fname, sim_name)
        print >> sys.stderr, "\tSimulating reads:", cmd
        os.system(cmd)
        assert os.path.exists("%s.sam" % sim_name) and \
            os.path.exists("%s_1.fa" % sim_name) and \
            os.path.exists("%s_2.fa" % sim_name)
        
    sam_cmd = "awk '{NM = $13; if ($2 < 128) print $1\"\\t\"$3\"\\t\"$4\"\\t\"$6\"\\t\"NM}' %s.sam > %s_1.sam" % (sim_name, sim_name)
    os.system(sam_cmd)
    sam_cmd = "awk '{NM = $13; if ($2 >= 128) print $1\"\\t\"$3\"\\t\"$4\"\\t\"$6\"\\t\"NM}' %s.sam > %s_2.sam" % (sim_name, sim_name)
    os.system(sam_cmd)

    # Extract pairs
    out_file = open("%s_paired.sam" % sim_name, 'w')
    hits_file = open("%s.sam" % sim_name)
    def write_reads(file, reads):
        assert len(reads) == 2
        chr1, pos1, cigar1, NM1 = reads[0]
        chr2, pos2, cigar2, NM2 = reads[1]
        if pos2 > pos1:
            p_str = "%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s" % \
                (prev_read_id, chr1, pos1, cigar1, chr2, pos2, cigar2, NM1, NM2)
        else:
            p_str = "%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s" % \
                (prev_read_id, chr2, pos2, cigar2, chr1, pos1, cigar1, NM2, NM1)
        print >> file, p_str
        
    prev_read_id, reads = None, []
    for line in hits_file:
        if line[0] == '@':
            continue
        cols = line[:-1].split()
        read_id, flag, chr, pos, mapQ, cigar = cols[:6]
        pos = int(pos)
        NM = ""
        for col in cols:
            if col[:2] == "NM":
                NM = col[5:]
        assert NM != ""
        if prev_read_id != None and prev_read_id != read_id:
            write_reads(out_file, reads)
            reads = []
        prev_read_id = read_id            
        reads.append([chr, pos, cigar, NM])
    if len(reads) > 0:
        write_reads(out_file, reads)
    out_file.close()


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


"""
"""
def get_cigar_chars_MN(cigars):
    cigars = cigar_re.findall(cigars)
    cigar_chars = ""
    for cigar in cigars:
        cigar_op = cigar[-1]
        if cigar_op in "MN":
            if cigar_chars == "" or cigar_chars[-1] != cigar_op:
                cigar_chars += cigar_op

    return cigar_chars


"""
"""
def extract_single(infilename, outfilename, chr_dic, aligner, version, debug_dic):
    infile = open(infilename)
    outfile = open(outfilename, "w")
    prev_read_id = ""
    num_reads, num_aligned_reads, num_ualigned_reads = 0, 0, 0
    if aligner == "hisat2":
        prev_NH, NH_real = 0, 0
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

        if read_id != prev_read_id:
            num_reads += 1

        flag, pos, mapQ = int(flag), int(pos), int(mapQ)
        if flag & 0x4 != 0:
            prev_read_id = read_id
            continue

        NH = ""
        NM = ""
        for i in range(11, len(cols)):
            col = cols[i]
            # "nM" from STAR
            if col.startswith("NM") or col.startswith("nM"):
                NM = col
            elif col.startswith("NH"):
                NH = col
        assert NM != ""
        NM = int(NM[5:])
        if NH != "":
            NH = int(NH[5:])
            if aligner == "hisat2":
                if prev_read_id == read_id:
                    assert prev_NH == NH
                if NH == 1 or mapQ == 60:
                    assert NH == 1 and mapQ == 60
                    
        if read_id != prev_read_id:
            num_aligned_reads += 1
            if aligner == "hisat2" and \
               NH == 1:
                num_ualigned_reads += 1
                
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

        if aligner == "hisat2":
            if prev_read_id != read_id:
                if prev_read_id != "":
                    assert prev_NH == NH_real
                NH_real = 1
            else:
                NH_real += 1
            prev_NH = NH
        prev_read_id = read_id

    if aligner == "hisat2":
        if prev_read_id != "":
            assert prev_NH == NH_real

    outfile.close()
    infile.close()


"""
"""
def extract_pair(infilename, outfilename, chr_dic, aligner, version, debug_dic):
    read_dic = {}
    pair_reported = set()

    infile = open(infilename)
    outfile = open(outfilename, "w")
    num_pairs, num_conc_aligned_pairs, num_conc_ualigned_pairs, num_disc_aligned_pairs = 0, 0, 0, 0
    num_aligned_reads, num_ualigned_reads = 0, 0
    
    prev_read_id, pair_list = "", set()
    if aligner == "hisat2":
        prev_NH1, prev_NH2 = 0, 0
        NH1_real, NH2_real = 0, 0
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

        if read_id != prev_read_id:
            num_pairs += 1
            pair_list = set()
        
        flag = int(flag)
        canonical_pos1, canonical_pos2 = int(pos1), int(pos2)
        left_read = (flag & 0x40 != 0)
        pos1 = canonical_pos1
        mapQ = int(mapQ)
        if flag & 0x4 != 0:
            prev_read_id, is_prev_read_left = read_id, left_read
            continue

        concordant = (flag & 0x2 != 0)        
        NH, NM1, YT = "", "", ""
        for i in range(11, len(cols)):
            col = cols[i]
            # "nM" from STAR
            if col.startswith("NM") or col.startswith("nM"):
                NM1 = col
            elif col.startswith("NH"):
                NH = col
            elif col.startswith("YT"):
                YT = col[5:]
        assert NM1 != ""
        NM1 = int(NM1[5:])        

        if aligner == "hisat2":
            assert NH != ""
            NH = int(NH[5:])
            if prev_read_id == read_id:
                if left_read:
                    assert prev_NH1 == 0 or prev_NH1 == NH
                else:
                    assert prev_NH2 == 0 or prev_NH2 == NH

        unpaired = (flag & 0x8 != 0) or (YT in ["UU", "UP"])
        if unpaired:
            if left_read not in pair_list:
                pair_list.add(left_read)
                num_aligned_reads += 1
                if aligner == "hisat2" and NH == 1:
                    num_ualigned_reads += 1
                    assert mapQ == 60
        else:
            if read_id != prev_read_id:
                if concordant:
                    num_conc_aligned_pairs += 1
                    if aligner == "hisat2" and NH == 1:
                        num_conc_ualigned_pairs += 1                        
                else:
                    if aligner == "hisat2":
                        assert YT == "DP"
                    num_disc_aligned_pairs += 1
                        
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
        if aligner == "hisat2":
            if prev_read_id != read_id:
                if prev_read_id != "":
                    assert prev_NH1 == NH1_real
                    assert prev_NH2 == NH2_real
                    prev_NH1, prev_NH2 = 0, 0
                if left_read:
                    NH1_real, NH2_real = 1, 0
                else:
                    NH1_real, NH2_real = 0, 1
            else:
                if left_read:
                    NH1_real += 1
                else:
                    NH2_real += 1
            if left_read:
                prev_NH1 = NH
            else:
                prev_NH2 = NH
        prev_read_id = read_id

    if aligner == "hisat2":
        if prev_read_id != "":
            assert prev_NH1 == NH1_real
            assert prev_NH2 == NH2_real

    outfile.close()
    infile.close()


"""
"""
def compare_single_sam(reference_sam, query_sam, mapped_fname, chr_dic):
    aligned, multi_aligned = 0, 0
    db_dic = {}
    mapped_file = open(mapped_fname, "w")
    file = open(reference_sam, "r")
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

    file = open(query_sam)
    mapped, unmapped, unique_mapped, mapping_point = 0, 0, 0, 0.0
    for line in file:
        if line[0] == '@':
            continue
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
                    break

        if found:
            print >> mapped_file, read_name
            mapped += 1
            if len(maps) == 1:
                unique_mapped += 1
            mapping_point += (1.0 / len(maps))
        else:
            unmapped += 1
            
    mapped_file.close()
    return mapped, unique_mapped, unmapped, aligned, multi_aligned, mapping_point


"""
"""
def compare_paired_sam(reference_sam, query_sam, mapped_fname, chr_dic):
    aligned, multi_aligned = 0, 0
    db_dic = {}
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

    file = open(query_sam)
    mapped, unique_mapped, unmapped, mapping_point = 0, 0, 0, 0.0
    for line in file:
        if line[0] == '@':
            continue
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
                    break

        if found:
            print >> mapped_file, read_name
            mapped += 1
            if len(maps) == 1:
                unique_mapped += 1
            mapping_point += (1.0 / len(maps))
        else:
            unmapped += 1
            
    mapped_file.close()
    return mapped, unique_mapped, unmapped, aligned, multi_aligned, mapping_point


"""
"""
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


"""
"""
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


"""
"""
def create_sql_db(sql_db):
    if os.path.exists(sql_db):
        print >> sys.stderr, sql_db, "already exists!"
        return
    
    columns = [
        ["id", "integer primary key autoincrement"],
        ["genome", "text"],
        ["head", "text"],
        ["end_type", "text"],
        ["aligner", "text"],
        ["version", "text"],
        ["num_reads", "integer"],
        ["mapped_reads", "integer"],
        ["unique_mapped_reads", "integer"],
        ["unmapped_reads", "integer"],
        ["mapping_point", "real"],
        ["time", "real"],
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


"""
"""
def write_analysis_data(sql_db, genome_name):
    if not os.path.exists(sql_db):
        return
    
    aligners = []
    sql_aligners = "SELECT aligner FROM ReadCosts GROUP BY aligner"
    output = sql_execute(sql_db, sql_aligners)
    aligners = output.split()

    for paired in [False, True]:
        database_fname = genome_name
        if paired:
            end_type = "paired"
            database_fname += "_paired"
        else:
            end_type = "single"
            database_fname += "_single"
        database_fname += ".analysis"
        database_file = open(database_fname, "w")
        print >> database_file, "end_type\taligner\tnum_reads\ttime\tmapped_reads\tunique_mapped_reads\tunmapped_reads\tmapping_point"
        for aligner in aligners:
            sql_row = "SELECT end_type, aligner, num_reads, time, mapped_reads, unique_mapped_reads, unmapped_reads, mapping_point FROM ReadCosts"
            sql_row += " WHERE genome = '%s' and head = '%s' and aligner = '%s' and end_type = '%s' ORDER BY created DESC LIMIT 1" % (genome_name, database_name, aligner, end_type)
            output = sql_execute(sql_db, sql_row)
            if output:
                print >> database_file, output

        database_file.close()


"""
"""
def eval(aligners,
         ref_genome,
         region,
         verbose):

    # Download HISAT2 index
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]
    if not typing_common.check_files(HISAT2_fnames):
        typing_common.download_genome_and_index()
    
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
    verbose = True
    just_runtime = False
    sql_write = False
    """
    aligners = [
        ["hisat2", "", "", ""],
        ["bowtie2", "", "", ""],
        ["bwa", "mem", "", ""]
        ]
    """

    if len(region) > 0:
        if len(region) == 1:
            genome = "genome_" + region[0]
        else:
            genome = "genome_%s:%d-%d" % (region[0], region[1], region[2])
    else:
        genome = "genome"

    # genomic sequence, variants, and haplotypes                                                                                                                                                                                                         
    epigenome_fnames = ["%s.fa" % genome,
                        "%s.snp" % genome,
                        "%s.haplotype" % genome]

    for aligner, _, index_type, _ in aligners:
        if aligner == "hisat2":
            index_fnames = ["%s.%s.%d.ht2" % (genome, index_type, i+1) for i in range(8)]
        elif aligner == "bowtie2":
            index_fnames = ["%s.%d.bt2" % (genome, i+1) for i in range(4)]
            index_fnames = ["%s.rev.%d.bt2" % (genome, i+1) for i in range(2)]
        else:
            assert aligner == "bwamem"
            assert False
            
        if not typing_common.check_files(epigenome_fnames + index_fnames):
            build_cmd = ["hisatbisulfite_build_epigenome.py",
                         "--base-fname", "genome",
                         "--abundance-cutoff", "0.05",
                         "--sample-cutoff", "5",
                         "--aligner", aligner]
            if len(region) > 0:
                region_str = region[0]
                if len(region) > 1:
                    region_str = "%s:%d-%d" % (region[0], region[1], region[2])
                build_cmd += ["--region", region_str]

            if index_type == "graph":
                build_cmd += ["--index-suffix", "graph"]
            else:
                build_cmd += ["--linear-index"]
                if aligner == "hisat2":
                    build_cmd += ["--index-suffix", "linear"]
            print >> sys.stderr, "\tBuilding a %s's %s index:" % (aligner, index_type), ' '.join(build_cmd)
            proc = subprocess.Popen(build_cmd,
                                    stderr=open("/dev/null", 'w'))
            proc.communicate()
            assert typing_common.check_files(epigenome_fnames + index_fnames)

    # Load genomic sequences
    chr_dic, chr_names, chr_full_names = typing_common.read_genome(open("%s.fa" % genome))

    cpg, mismatch, numreads = True, False, 1000000

    align_stat = []
    for paired in [False, True]:
        sim_name = genome
        if cpg:
            sim_name += "_cpg"
        if mismatch:
            sim_name += "_mismatch"
        sim_name = "%s_%dM" % (sim_name, numreads / 1000000)    
        if paired:
            type_sam_fname = sim_name + "_paired.sam"
            type_read1_fname = sim_name +  "_1.fa"
            type_read2_fname = sim_name +  "_2.fa"
        else:
            type_sam_fname = sim_name + "_1.sam"
            type_read1_fname = sim_name + "_1.fa"
            type_read2_fname = ""

        if not os.path.exists(type_sam_fname):
            simulate_reads(genome,
                           sim_name,
                           cpg,
                           mismatch,
                           numreads)
        numreads = 0
        type_sam_file = open(type_sam_fname)
        for line in type_sam_file:
            numreads += 1
        type_sam_file.close()
        if numreads <= 0:
            continue
        print >> sys.stderr, "%d" % numreads

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
                    # cmd = ["%s/%s" % (aligner_bin_base, aligner)]
                    cmd = [aligner]
                cmd += ["--version"]                    
                cmd_process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                version = cmd_process.communicate()[0][:-1].split("\n")[0]
                version = version.split()[-1]
            elif aligner in ["star", "starx2"]:
                version = "2.4.2a"
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
        def get_aligner_cmd(aligner, type, index_type, version, read1_fname, read2_fname, out_fname, cmd_idx = 0):
            cmd = []
            if aligner == "hisat2":
                if version:
                    cmd = ["%s/hisat2_%s/hisat2" % (aligner_bin_base, version)]
                else:
                    # cmd = ["%s/hisat2" % (aligner_bin_base)]
                    cmd = [aligner]
                if num_threads > 1:
                    cmd += ["-p", str(num_threads)]
                cmd += ["-f"]
                # cmd += ["-k", "10"]
                # cmd += ["--max-seeds", "100"]
                # cmd += ["--score-min", "C,-18"]
                cmd += ["--no-spliced-alignment"]

                # cmd += ["--no-anchorstop"]
                if version == "" or \
                   (version != "" and int(version) >= 210):
                    cmd += ["--new-summary",
                            "--summary-file", out_fname + ".summary"]

                if index_type == "linear":
                    cmd += ["--score-min", "C,-50"]

                # index_cmd = "%s/HISAT2%s/" % (index_base, index_add) + genome
                # if index_type:
                #     index_cmd += ("_" + index_type)
                index_cmd = "../%s.%s" % (genome, index_type)
                cmd += ["-x", index_cmd]
                if paired:
                    cmd += ["-1", read1_fname,
                            "-2", read2_fname]
                else:
                    cmd += [read1_fname]                        
            elif aligner == "star":
                cmd = ["%s/STAR" % (aligner_bin_base)]
                if num_threads > 1:
                    cmd += ["--runThreadN", str(num_threads)]
                cmd += ["--genomeDir"]
                if cmd_idx == 0:
                    if type == "gtf":
                        cmd += ["%s/STAR%s/gtf" % (index_base, index_add)]
                    else:
                        cmd += ["%s/STAR%s" % (index_base, index_add)]
                else:
                    assert cmd_idx == 1
                    cmd += ["."]

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
                # cmd = ["%s/bowtie2" % (aligner_bin_base)]
                cmd = [aligner]
                if num_threads > 1:
                    cmd += ["-p", str(num_threads)]
                cmd += ["-f"]
                # cmd += ["-k", "10"]
                # cmd += ["--score-min", "C,-18"]
                # cmd += ["-x", "%s/HISAT%s/" % (index_base, index_add) + genome]
                cmd += ["-x", "../%s" % genome]
                if paired:
                    cmd += ["-1", read1_fname,
                            "-2", read2_fname]
                else:
                    cmd += [read1_fname]
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

        if len(region) > 0:
            init_time = {"hisat2" : 0.2, "bowtie" : 0.1, "bowtie2" : 0.2, "bwa" : 0.0}
            if desktop:
                init_time["star"] = 1.7
            else:
                init_time["star"] = 0.0
        else:
            if desktop:
                init_time = {"hisat2" : 3.0, "bowtie" : 1.3, "bowtie2" : 1.9, "star" : 27.0, "bwa" : 1.3}
            else:
                init_time = {"hisat2" : 9.5, "bowtie" : 3.3, "bowtie2" : 4.1, "star" : 1.7, "bwa" : 3.3}

        for aligner, type, index_type, version in aligners:
            aligner_name = aligner + type
            if version != "":
                aligner_name += ("_%s" % version)
            if aligner == "hisat2" and index_type != "":
                aligner_name += ("_" + index_type)
            print >> sys.stderr, "\t%s\t%s" % (aligner_name, str(datetime.now()))
            if paired:
                aligner_dir = aligner_name + "_paired"
            else:
                aligner_dir = aligner_name + "_single"
            if not os.path.exists(aligner_dir):
                os.mkdir(aligner_dir)
            else:
                continue
            os.chdir(aligner_dir)

            out_fname = sim_name + ".sam"
            out_fname2 = out_fname + "2"
            duration = -1.0
            if not os.path.exists(out_fname):
                if not os.path.exists("../one.fa") or not os.path.exists("../two.fa"):
                    os.system("head -400 ../%s_1.fa > ../one.fa" % sim_name)
                    os.system("head -400 ../%s_2.fa > ../two.fa" % sim_name)

                if just_runtime:
                    out_fname = "/dev/null"
                    out_fname2 = "/dev/null"

                align_stat.append(["all", aligner_name])

                # dummy commands for caching index and simulated reads
                dummy_cmd = get_aligner_cmd(aligner, type, index_type, version, "../one.fa", "../two.fa", "/dev/null")
                if verbose:
                    print >> sys.stderr, "\t", datetime.now(), " ".join(dummy_cmd)
                if aligner in ["hisat2", "bowtie", "bowtie2", "bwa"]:
                    proc = subprocess.Popen(dummy_cmd, stdout=open("/dev/null", "w"), stderr=subprocess.PIPE)
                else:
                    proc = subprocess.Popen(dummy_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                proc.communicate()
                if verbose:
                    print >> sys.stderr, "\t", datetime.now(), "finished"

                # Align all reads
                aligner_cmd = get_aligner_cmd(aligner, type, index_type, version, "../" + type_read1_fname, "../" + type_read2_fname, out_fname)
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
                if aligner == "star":
                    os.system("mv Aligned.out.sam %s" % out_fname)
                elif aligner in ["hisat2", "hisat"] and type == "x2":
                    aligner_cmd = get_aligner_cmd(aligner, type, index_type, version, "../" + type_read1_fname, "../" + type_read2_fname, out_fname, 1)
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

            if just_runtime:
                os.chdir("..")
                continue

            if not os.path.exists(out_fname2):
                debug_dic = {}
                if paired:
                    extract_pair(out_fname, out_fname2, chr_dic, aligner, version, debug_dic)
                else:
                    extract_single(out_fname, out_fname2, chr_dic, aligner, version, debug_dic)

            mapped_id_fname = sim_name+ ".read_id"
            if paired:
                mapped, unique_mapped, unmapped, aligned, multi_aligned, mapping_point = \
                    compare_paired_sam(out_fname2, "../" + type_sam_fname, mapped_id_fname, chr_dic)
            else:
                mapped, unique_mapped, unmapped, aligned, multi_aligned, mapping_point = \
                    compare_single_sam(out_fname2, "../" + type_sam_fname, mapped_id_fname, chr_dic)
            proc = subprocess.Popen(["wc", "-l", "../" + type_read1_fname], stdout=subprocess.PIPE)
            out = proc.communicate()[0]
            numreads = int(out.split()[0]) / 2
            assert mapped + unmapped == numreads

            print >> sys.stderr, "\t\taligned: %d, multi aligned: %d" % (aligned, multi_aligned)
            print >> sys.stderr, "\t\tcorrectly mapped: %d (%.2f%%) mapping_point: %.2f" % (mapped, float(mapped) * 100.0 / numreads, mapping_point * 100.0 / numreads)
            print >> sys.stderr, "\t\tuniquely and correctly mapped: %d (%.2f%%)" % (unique_mapped, float(unique_mapped) * 100.0 / numreads)
            print >> sys.stderr, "\t\t\t%d reads per sec (all)" % (numreads / max(1.0, duration))

            if duration > 0.0:
                if sql_write and os.path.exists("../" + sql_db_name):
                    if paired:
                        end_type = "paired"
                    else:
                        end_type = "single"
                    sql_insert = "INSERT INTO \"ReadCosts\" VALUES(NULL, '%s', '%s', '%s', '%s', '%s', %d, %d, %d, %d, %f, %f, '%s', datetime('now', 'localtime'), '%s');" % \
                        (genome, data_base, end_type, aligner_name, get_aligner_version(aligner, version), numreads, mapped, unique_mapped, unmapped, mapping_point, duration, platform.node(), " ".join(aligner_cmd))
                    sql_execute("../" + sql_db_name, sql_insert)     
                align_stat[-1].extend([numreads, duration, mapped, unique_mapped, unmapped, mapping_point])

            os.system("touch %s.done" % type_sam_fname)
            os.chdir("..")

    print >> sys.stdout, "\t".join(["type", "aligner", "all", "all_time", "mapped", "unique_mapped", "unmapped", "mapping point"])
    for line in align_stat:
        outstr = ""
        for item in line:
            if outstr != "":
                outstr += "\t"
            outstr += str(item)
        print >> sys.stdout, outstr

    if os.path.exists(sql_db_name):
        write_analysis_data(sql_db_name, genome)
        


"""
"""
if __name__ == "__main__":
    parser = ArgumentParser(
        description='test HISAT-bisulfite, and compare HISAT-bisulfite with other popular aligners such as Bowtie1/2, BWA-mem, Bismark, and BSmooth etc.')
    parser.add_argument("--aligner-list",
                        dest="aligner_list",
                        type=str,
                        default="hisat2.graph,hisat2.linear,bowtie2",
                        help="A comma separated list of aligners (default: hisat2.graph,hisat2.linear,bowtie2)")
    parser.add_argument("--genome",
                        dest="genome",
                        type=str,
                        default="human",
                        help="Reference genome (default: human)")
    parser.add_argument("--region",
                        dest="region",
                        type=str,
                        default="",
                    help="Genomic region (1-offset) instead of the whole genome (default: (empty), example: 7 and 10:1000000-2000000)")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()
    aligners = []
    for aligner_index in args.aligner_list.split(','):
        version, index_type = "", "linear"
        if aligner_index.find('.') != -1:
            aligner, index_type = aligner_index.split('.')
        else:
            aligner = aligner_index
        if aligner.find('_') != -1:
            aligner, version = aligner.split('_')
        if aligner not in ["hisat2", "bowtie2", "bwamem"]:
            print >> sys.stderr, "Error: --aligner-list should be among hisat2, bowtie2, and bwamem"
            sys.exit(1)
        aligners.append([aligner, "", index_type, version])

    region = []
    if args.region != "":
        try:
            region_chr, region_range = args.region.split(':')
            if region_range != "":
                region_left, region_right = region_range.split('-')
                region_left, region_right = int(region_left), int(region_right)
                region = [region_chr, region_left, region_right]
        except ValueError:
            print >> sys.stderr, "Error: --region %s is ill-formatted." % args.region
            sys.exit(1)

    eval(aligners,
         args.genome,
         region,
         args.verbose)
