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


import sys, os, subprocess, re
import inspect, random
import math
from datetime import datetime, date, time
from argparse import ArgumentParser, FileType
from copy import deepcopy
from hisatgenotype_modules import typing_common, HLA_typing, assembly_graph



"""
Align reads, and sort the alignments into a BAM file
"""
def align_reads(ex_path,
                aligner,
                simulation,
                index_type,
                read_fname,
                fastq,
                threads,
                out_fname,
                verbose):
    if aligner == "hisat2":
        hisat2 = os.path.join(ex_path, "hisat2")
        aligner_cmd = [hisat2, "--mm"]
        if not simulation:
            aligner_cmd += ["--no-unal"]            
        DNA = True
        if DNA:
            aligner_cmd += ["--no-spliced-alignment"] # no spliced alignment
            aligner_cmd += ["-X", "1000"] # max fragment length
        if index_type == "linear":
            aligner_cmd += ["-k", "10"]
        else:
            aligner_cmd += ["--max-altstried", "64"]
            # DK - debugging purposes
            # aligner_cmd += ["--max-seeds", "100"]
        aligner_cmd += ["-x", "hla.%s" % index_type]
    elif aligner == "bowtie2":
        aligner_cmd = [aligner,
                       "--no-unal",
                       "-k", "10",
                       "-x", "hla"]
    else:
        assert False
    assert len(read_fname) in [1,2]
    aligner_cmd += ["-p", str(threads)]
    if not fastq:
        aligner_cmd += ["-f"]
    if len(read_fname) == 1:
        aligner_cmd += ["-U", read_fname[0]]
    else:
        aligner_cmd += ["-1", "%s" % read_fname[0],
                        "-2", "%s" % read_fname[1]]
    if verbose >= 1:
        print >> sys.stderr, ' '.join(aligner_cmd)
    align_proc = subprocess.Popen(aligner_cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=open("/dev/null", 'w'))

    sambam_cmd = ["samtools",
                  "view",
                  "-bS",
                  "-"]
    sambam_proc = subprocess.Popen(sambam_cmd,
                                   stdin=align_proc.stdout,
                                   stdout=open(out_fname + ".unsorted", 'w'),
                                   stderr=open("/dev/null", 'w'))
    sambam_proc.communicate()
    if index_type == "graph":
        bamsort_cmd = ["samtools",
                       "sort",
                       out_fname + ".unsorted",
                       "-o", out_fname]
        bamsort_proc = subprocess.Popen(bamsort_cmd,
                                        stderr=open("/dev/null", 'w'))
        bamsort_proc.communicate()

        bamindex_cmd = ["samtools",
                        "index",
                        out_fname]
        bamindex_proc = subprocess.Popen(bamindex_cmd,
                                         stderr=open("/dev/null", 'w'))
        bamindex_proc.communicate()

    os.system("rm %s" % (out_fname + ".unsorted"))            


"""
""" 
def normalize(prob):
    total = sum(prob.values())
    for allele, mass in prob.items():
        prob[allele] = mass / total

        
"""
"""
def prob_diff(prob1, prob2):
    diff = 0.0
    for allele in prob1.keys():
        if allele in prob2:
            diff += abs(prob1[allele] - prob2[allele])
        else:
            diff += prob1[allele]
    return diff


"""
"""
def HLA_prob_cmp(a, b):
    if a[1] != b[1]:
        if a[1] < b[1]:
            return 1
        else:
            return -1
    assert a[0] != b[0]
    if a[0] < b[0]:
        return -1
    else:
        return 1


"""
"""
def single_abundance(HLA_cmpt,
                     HLA_length):
    def normalize2(prob, length):
        total = 0
        for allele, mass in prob.items():
            assert allele in length
            total += (mass / length[allele])
        for allele, mass in prob.items():
            assert allele in length
            prob[allele] = mass / length[allele] / total

    HLA_prob, HLA_prob_next = {}, {}
    for cmpt, count in HLA_cmpt.items():
        alleles = cmpt.split('-')
        for allele in alleles:
            if allele not in HLA_prob:
                HLA_prob[allele] = 0.0
            HLA_prob[allele] += (float(count) / len(alleles))

    normalize(HLA_prob)
    def next_prob(HLA_cmpt, HLA_prob, HLA_length):
        HLA_prob_next = {}
        for cmpt, count in HLA_cmpt.items():
            alleles = cmpt.split('-')
            alleles_prob = 0.0
            for allele in alleles:
                if allele not in HLA_prob:
                    continue
                alleles_prob += HLA_prob[allele]
            if alleles_prob <= 0.0:
                continue
            for allele in alleles:
                if allele not in HLA_prob:
                    continue
                if allele not in HLA_prob_next:
                    HLA_prob_next[allele] = 0.0
                HLA_prob_next[allele] += (float(count) * HLA_prob[allele] / alleles_prob)
        normalize(HLA_prob_next)
        return HLA_prob_next

    diff, iter = 1.0, 0
    while diff > 0.0001 and iter < 1000:
        HLA_prob_next = next_prob(HLA_cmpt, HLA_prob, HLA_length)
        diff = prob_diff(HLA_prob, HLA_prob_next)
        HLA_prob = HLA_prob_next

        # Accelerate convergence
        """
        if iter >= 10:
            HLA_prob2 = {}
            for allele, prob in HLA_prob.items():
                if prob >= 0.005:
                    HLA_prob2[allele] = prob
            HLA_prob = HLA_prob2
        """

        # DK - debugging purposes
        if iter % 10 == 0 and False:
            print "iter", iter
            for allele, prob in HLA_prob.items():
                if prob >= 0.01:
                    print >> sys.stderr, "\t", iter, allele, prob, str(datetime.now())
        
        iter += 1
        
    """
    for allele, prob in HLA_prob.items():
        allele_len = HLA_length[allele]
        HLA_prob[allele] /= float(allele_len)
    """
    
    # normalize(HLA_prob)
    normalize2(HLA_prob, HLA_length)
    HLA_prob = [[allele, prob] for allele, prob in HLA_prob.items()]
    HLA_prob = sorted(HLA_prob, cmp=HLA_prob_cmp)
    return HLA_prob

    
"""
"""
def joint_abundance(HLA_cmpt,
                    HLA_length):
    allele_names = set()
    for cmpt in HLA_cmpt.keys():
        allele_names |= set(cmpt.split('-'))
    
    HLA_prob, HLA_prob_next = {}, {}
    for cmpt, count in HLA_cmpt.items():
        alleles = cmpt.split('-')
        for allele1 in alleles:
            for allele2 in allele_names:
                if allele1 < allele2:
                    allele_pair = "%s-%s" % (allele1, allele2)
                else:
                    allele_pair = "%s-%s" % (allele2, allele1)
                if not allele_pair in HLA_prob:
                    HLA_prob[allele_pair] = 0.0
                HLA_prob[allele_pair] += (float(count) / len(alleles))

    if len(HLA_prob) <= 0:
        return HLA_prob

    # Choose top allele pairs
    def choose_top_alleles(HLA_prob):
        HLA_prob_list = [[allele_pair, prob] for allele_pair, prob in HLA_prob.items()]
        HLA_prob_list = sorted(HLA_prob_list, cmp=HLA_prob_cmp)
        HLA_prob = {}
        best_prob = HLA_prob_list[0][1]
        for i in range(len(HLA_prob_list)):
            allele_pair, prob = HLA_prob_list[i]
            if prob * 2 <= best_prob:
                break                        
            HLA_prob[allele_pair] = prob
        normalize(HLA_prob)
        return HLA_prob
    HLA_prob = choose_top_alleles(HLA_prob)

    def next_prob(HLA_cmpt, HLA_prob):
        HLA_prob_next = {}
        for cmpt, count in HLA_cmpt.items():
            alleles = cmpt.split('-')
            prob = 0.0
            for allele in alleles:
                for allele_pair in HLA_prob.keys():
                    if allele in allele_pair:
                        prob += HLA_prob[allele_pair]
            for allele in alleles:
                for allele_pair in HLA_prob.keys():
                    if not allele in allele_pair:
                        continue
                    if allele_pair not in HLA_prob_next:
                        HLA_prob_next[allele_pair] = 0.0
                    HLA_prob_next[allele_pair] += (float(count) * HLA_prob[allele_pair] / prob)
        normalize(HLA_prob_next)
        return HLA_prob_next

    diff, iter = 1.0, 0
    while diff > 0.0001 and iter < 1000:
        HLA_prob_next = next_prob(HLA_cmpt, HLA_prob)
        diff = prob_diff(HLA_prob, HLA_prob_next)
        HLA_prob = HLA_prob_next
        HLA_prob = choose_top_alleles(HLA_prob)
        iter += 1

    HLA_prob = [[allele_pair, prob] for allele_pair, prob in HLA_prob.items()]
    HLA_prob = sorted(HLA_prob, cmp=HLA_prob_cmp)
    return HLA_prob


"""
"""
def lower_bound(Var_list, pos):
    low, high = 0, len(Var_list)
    while low < high:
        m = (low + high) / 2
        m_pos = Var_list[m][0]
        if m_pos < pos:
            low = m + 1
        elif m_pos > pos:
            high = m
        else:
            assert m_pos == pos
            while m > 0:
                if Var_list[m-1][0] < pos:
                    break
                m -= 1
            return m
    return low


"""
   var: ['single', 3300, 'G']
   exons: [[301, 373], [504, 822], [1084, 1417], [2019, 2301], [2404, 2520], [2965, 2997], [3140, 3187], [3357, 3361]]
"""
def var_in_exon(var, exons):
    exonic = False
    var_type, var_left, var_data = var
    var_right = var_left
    if var_type == "deletion":
        var_right = var_left + int(var_data) - 1
    for exon_left, exon_right in exons:
        if var_left >= exon_left and var_right <= exon_right:
            return True
    return False


"""
Report variant IDs whose var is within exonic regions
"""
def get_exonic_vars(Vars, exons):
    vars = set()
    for var_id, var in Vars.items():
        var_type, var_left, var_data = var
        var_right = var_left
        if var_type == "deletion":
            var_right = var_left + int(var_data) - 1
        for exon_left, exon_right in exons:
            if var_left >= exon_left and var_right <= exon_right:
                vars.add(var_id)
    return vars


"""
Get representative alleles among those that share the same exonic sequences
"""
def get_rep_alleles(Links, exon_vars):
    allele_vars = {}
    for var, alleles in Links.items():
        if var not in exon_vars:
            continue
        for allele in alleles:
            if allele not in allele_vars:
                allele_vars[allele] = set()
            allele_vars[allele].add(var)

    allele_groups = {}
    for allele, vars in allele_vars.items():
        vars = '-'.join(vars)
        if vars not in allele_groups:
            allele_groups[vars] = []
        allele_groups[vars].append(allele)

    allele_reps = {} # allele representatives
    allele_rep_groups = {} # allele groups by allele representatives
    for allele_members in allele_groups.values():
        assert len(allele_members) > 0
        allele_rep = allele_members[0]
        allele_rep_groups[allele_rep] = allele_members
        for allele_member in allele_members:
            assert allele_member not in allele_reps
            allele_reps[allele_member] = allele_rep

    return allele_reps, allele_rep_groups
    

"""
Identify alternative alignments
"""
def get_alternatives(ref_seq, Vars, Var_list, verbose):
    # Check deletions' alternatives
    def get_alternatives_recur(ref_seq,
                               Vars,
                               Var_list,
                               Alts,
                               var_id,
                               left,
                               alt_list,
                               var_j,
                               latest_pos,
                               debug = False):
        def add_alt(Alts, alt_list, var_id, j_id):
            if j_id.isdigit():
                if var_id not in Alts:
                    Alts[var_id] = [["1"]]
                else:
                    if Alts[var_id][-1][-1].isdigit():
                        Alts[var_id][-1][-1] = str(int(Alts[var_id][-1][-1]) + 1)
                    else:
                        Alts[var_id][-1].append("1")
            else:
                if var_id not in Alts:
                    Alts[var_id] = [[j_id]]
                else:
                    if Alts[var_id][-1][-1].isdigit():
                        Alts[var_id][-1][-1] = j_id
                    else:
                        Alts[var_id][-1].append(j_id)
                Alts[var_id][-1].append("0")
                        
            if not j_id.isdigit():
                alt_list.append(j_id)
                alts = '-'.join(alt_list)
                if alts not in Alts:
                    Alts[alts] = [[var_id]]
                else:
                    Alts[alts].append([var_id])
                
        var_type, var_pos, var_data = Vars[var_id]
        if left: # Look in left direction
            if var_j < 0:
                return
            j_pos, j_id = Var_list[var_j]
            alt_del = []
            if var_id != j_id and j_pos < var_pos + del_len:
                # Check bases between SNPs
                while latest_pos > j_pos:
                    if debug: print latest_pos - 1, ref_seq[latest_pos - 1], latest_pos - 1 - del_len, ref_seq[latest_pos - 1 - del_len]
                    if ref_seq[latest_pos - 1] != ref_seq[latest_pos - 1 - del_len]:
                        break
                    latest_pos -= 1
                    add_alt(Alts, alt_list, var_id, str(latest_pos))
                if latest_pos - 1 > j_pos:
                    return
                if j_pos == latest_pos - 1:
                    j_type, _, j_data = Vars[j_id]
                    if j_type == "single":
                        if debug: print Vars[j_id]
                        off = var_pos + del_len - j_pos
                        if debug: print var_pos - off, ref_seq[var_pos - off]
                        if debug: print j_pos, ref_seq[j_pos]
                        if j_data == ref_seq[var_pos - off]:
                            add_alt(Alts, alt_list, var_id, j_id)
                            latest_pos = j_pos
                    elif j_type == "deletion":
                        j_del_len = int(j_data)
                        if var_pos < j_pos and var_pos + del_len >= j_pos + j_del_len:
                            alt_list2 = alt_list[:] + [j_id]
                            latest_pos2 = j_pos
                            alt_del = [alt_list2, latest_pos2]
                
            get_alternatives_recur(ref_seq,
                                   Vars,
                                   Var_list,
                                   Alts,
                                   var_id,
                                   left,
                                   alt_list,
                                   var_j - 1,
                                   latest_pos,
                                   debug)

            if alt_del:
                alt_list2, latest_pos2 = alt_del
                if var_id not in Alts:
                    Alts[var_id] = [alt_list2[:]]
                else:
                    Alts[var_id].append(alt_list2[:])
                alt_idx = len(Alts[var_id]) - 1
                get_alternatives_recur(ref_seq,
                                       Vars,
                                       Var_list,
                                       Alts,
                                       var_id,
                                       left,
                                       alt_list2,
                                       var_j - 1,
                                       latest_pos2,
                                       debug)
                # Remove this Deletion if not supported by additional bases?
                assert alt_idx < len(Alts[var_id])
                # DK - for debugging purposes
                if Alts[var_id][alt_idx][-1] == j_id:
                    Alts[var_id] = Alts[var_id][:alt_idx] + Alts[var_id][alt_idx+1:]
              
        else: # Look in right direction
            if var_j >= len(Var_list):
                return
            j_pos, j_id = Var_list[var_j]
            alt_del = []
            if var_id != j_id and j_pos >= var_pos:
                # Check bases between SNPs
                while latest_pos < j_pos:
                    if ref_seq[latest_pos + 1] != ref_seq[var_pos + del_len - 1 - (latest_pos - var_pos)]:
                        break
                    latest_pos += 1
                    add_alt(Alts, alt_list, var_id, str(latest_pos))
                if latest_pos + 1 < j_pos:
                    return
                if j_pos == latest_pos + 1:
                    j_type, _, j_data = Vars[j_id]
                    if j_type == "single":
                        if debug: print Vars[j_id]
                        off = j_pos - var_pos
                        if debug: print var_pos + off, ref_seq[var_pos + off]
                        if debug: print var_pos + del_len + off, ref_seq[var_pos + del_len + off]

                        # DK - for debugging purposes
                        if var_pos + del_len + off >= len(ref_seq):
                            print >> sys.stderr, var_id, var
                            print >> sys.stderr, "var_pos: %d, del_len: %d, off: %d" % (var_pos, del_len, off)
                            print >> sys.stderr, "ref_seq: %d, %d" % (len(ref_seq), var_pos + del_len + off)
                            sys.exit(1)
                        
                        if j_data == ref_seq[var_pos + del_len + off]:
                            add_alt(Alts, alt_list, var_id, j_id)
                            latest_pos = j_pos
                    elif j_type == "deletion":
                        j_del_len = int(j_data)
                        if j_pos + j_del_len < var_pos + del_len:
                            alt_list2 = alt_list[:] + [j_id]
                            latest_pos2 = j_pos + j_del_len - 1
                            alt_del = [alt_list2, latest_pos2]

            get_alternatives_recur(ref_seq,
                                   Vars,
                                   Var_list,
                                   Alts,
                                   var_id,
                                   left,
                                   alt_list,
                                   var_j + 1,
                                   latest_pos,
                                   debug)

            if alt_del:
                alt_list2, latest_pos2 = alt_del
                if var_id not in Alts:
                    Alts[var_id] = [alt_list2[:]]
                else:
                    Alts[var_id].append(alt_list2[:])
                alt_idx = len(Alts[var_id]) - 1
                get_alternatives_recur(ref_seq,
                                       Vars,
                                       Var_list,
                                       Alts,
                                       var_id,
                                       left,
                                       alt_list2,
                                       var_j + 1,
                                       latest_pos2,
                                       debug)
                # Remove this Deletion if not supported by additional bases?
                assert alt_idx < len(Alts[var_id])
                if Alts[var_id][alt_idx][-1] == j_id:
                    Alts[var_id] = Alts[var_id][:alt_idx] + Alts[var_id][alt_idx+1:]

    # Check deletions' alternatives
    Alts_left, Alts_right = {}, {}
    for var_i, var_id in Var_list:
        var_type, var_pos, var_data = var = Vars[var_id]
        if var_type != "deletion" or var_pos == 0:
            continue
        del_len = int(var_data)
        if var_pos + del_len >= len(ref_seq):
            assert var_pos + del_len == len(ref_seq)
            continue
        debug = (var_id == "hv1096a")
        if debug:
            print Vars[var_id]

        alt_list = []
        var_j = lower_bound(Var_list, var_pos + del_len - 1)
        latest_pos = var_pos + del_len
        if var_j < len(Var_list):
            get_alternatives_recur(ref_seq,
                                   Vars,
                                   Var_list,
                                   Alts_left,
                                   var_id,
                                   True, # left
                                   alt_list,
                                   var_j,
                                   latest_pos,
                                   debug)
        alt_list = []
        var_j = lower_bound(Var_list, var_pos)
        latest_pos = var_pos - 1
        assert var_j >= 0
        get_alternatives_recur(ref_seq,
                               Vars,
                               Var_list,
                               Alts_right,
                               var_id,
                               False, # right
                               alt_list,
                               var_j,
                               latest_pos,
                               debug)

        if debug:
            print "DK :-)"
            sys.exit(1)

    def debug_print_alts(Alts, dir):
        for alt_list1, alt_list2 in Alts.items():
            print "\t", dir, ":", alt_list1, alt_list2
            out_str = "\t\t"
            alt_list1 = alt_list1.split('-')
            for i in range(len(alt_list1)):
                alt = alt_list1[i]
                var_type, var_pos, var_data = Vars[alt]
                out_str += ("%s-%d-%s" % (var_type, var_pos, var_data))
                if i + 1 < len(alt_list1):
                    out_str += " "
            for i in range(len(alt_list2)):
                alt_list3 = alt_list2[i]
                out_str += "\t["
                for j in range(len(alt_list3)):
                    alt = alt_list3[j]
                    if alt.isdigit():
                        out_str += alt
                    else:
                        var_type, var_pos, var_data = Vars[alt]
                        out_str += ("%s-%d-%s" % (var_type, var_pos, var_data))
                    if j + 1 < len(alt_list3):
                        out_str += ", "
                out_str += "]"
            print out_str
    if verbose >= 2: debug_print_alts(Alts_left, "left")
    if verbose >= 2: debug_print_alts(Alts_right, "right")

    return Alts_left, Alts_right


"""
Identify ambigious differences that may account for other alleles,
  given a list of differences (cmp_list) between a read and a potential allele   
"""
def identify_ambigious_diffs(Vars, Alts_left, Alts_right, cmp_list, verbose):
    cmp_left, cmp_right = 0, len(cmp_list) - 1
    i = 0
    while i < len(cmp_list):
        cmp_i = cmp_list[i]
        type, pos, length = cmp_i[:3]
        # Check alternative alignments
        if type in ["mismatch", "deletion"]:
            var_id = cmp_i[3]
            if var_id == "unknown":
                i += 1
                continue
            
            # Left direction
            id_str = var_id
            total_del_len = length if type == "deletion" else 0
            for j in reversed(range(0, i)):
                cmp_j = cmp_list[j]
                j_type, j_pos, j_len = cmp_j[:3]
                if j_type != "match":
                    if len(cmp_j) < 4:
                        continue
                    j_var_id = cmp_j[3]
                    id_str += ("-%s" % j_var_id)
                    if j_type == "deletion":
                        total_del_len += j_len
            last_type, last_pos, last_len = cmp_list[0][:3]
            assert last_type in ["match", "mismatch"]
            left_pos = last_pos + total_del_len
            if id_str in Alts_left:
                orig_alts = id_str.split('-')
                alts_list = Alts_left[id_str]
                for alts in alts_list:
                    if alts[-1].isdigit():
                        assert type == "deletion"
                        assert len(orig_alts) == 1
                        alts_id_str = '-'.join(alts[:-1])
                        alt_left_pos = pos
                        alt_total_del_len = 0
                        for alt in alts[:-1]:
                            assert alt in Vars
                            alt_type, alt_pos, alt_data = Vars[alt]
                            alt_left_pos = alt_pos - 1
                            if alt_type == "deletion":
                                alt_total_del_len += int(alt_data)
                        alt_left_pos = alt_left_pos + alt_total_del_len - int(alts[-1]) + 1
                    else:
                        alts_id_str = '-'.join(alts)
                        assert alts_id_str in Alts_left
                        for back_alts in Alts_left[alts_id_str]:
                            back_id_str = '-'.join(back_alts)
                            if back_id_str.find(id_str) != 0:
                                continue
                            assert len(orig_alts) < len(back_alts)
                            assert back_alts[-1].isdigit()
                            alt_left_pos = pos
                            alt_total_del_len = 0
                            for alt in back_alts[:len(orig_alts) + 1]:
                                if alt.isdigit():
                                    alt_left_pos = alt_left_pos - int(alt) + 1
                                else:
                                    assert alt in Vars
                                    alt_type, alt_pos, alt_data = Vars[alt]
                                    alt_left_pos = alt_pos - 1
                                    if alt_type == "deletion":
                                        alt_total_del_len += int(alt_data)
                            alt_left_pos += alt_total_del_len
                        if left_pos >= alt_left_pos:
                            if verbose >= 2:
                                print "LEFT:", cmp_list
                                print "\t", type, "id_str:", id_str, "=>", alts_id_str, "=>", back_alts, "left_pos:", left_pos, "alt_left_pos:", alt_left_pos
                            cmp_left = i + 1
                            break
    
            # Right direction
            if cmp_right + 1 == len(cmp_list):
                id_str = var_id
                total_del_len = length if type == "deletion" else 0
                for j in range(i + 1, len(cmp_list)):
                    cmp_j = cmp_list[j]
                    j_type, j_pos, j_len = cmp_j[:3]
                    if j_type != "match":
                        if len(cmp_j) < 4:
                            continue
                        j_var_id = cmp_j[3]
                        id_str += ("-%s" % j_var_id)
                        if j_type == "deletion":
                            total_del_len += j_len                        
                last_type, last_pos, last_len = cmp_list[-1][:3]
                assert last_type in ["match", "mismatch"]
                right_pos = last_pos + last_len - 1 - total_del_len
                if id_str in Alts_right:
                    orig_alts = id_str.split('-')
                    alts_list = Alts_right[id_str]
                    for alts in alts_list:
                        if alts[-1].isdigit():
                            assert type == "deletion"
                            assert len(orig_alts) == 1
                            alts_id_str = '-'.join(alts[:-1])
                            alt_right_pos = pos
                            alt_total_del_len = 0
                            for alt in alts[:-1]:
                                assert alt in Vars
                                alt_type, alt_pos, alt_data = Vars[alt]
                                alt_right_pos = alt_pos
                                if alt_type == "single":
                                    alt_right_pos += 1
                                else:
                                    assert alt_type == "deletion"
                                    alt_del_len = int(alt_data)
                                    alt_right_pos += alt_del_len
                                    alt_total_del_len += alt_del_len
                            alt_right_pos = alt_right_pos - alt_total_del_len + int(alts[-1]) - 1
                        else:
                            alts_id_str = '-'.join(alts)
                            assert alts_id_str in Alts_right
                            for back_alts in Alts_right[alts_id_str]:
                                back_id_str = '-'.join(back_alts)
                                if back_id_str.find(id_str) != 0:
                                    continue
                                assert len(orig_alts) < len(back_alts)
                                assert back_alts[-1].isdigit()
                                alt_right_pos = pos
                                alt_total_del_len = 0
                                for alt in back_alts[:len(orig_alts) + 1]:
                                    if alt.isdigit():
                                        alt_right_pos = alt_right_pos + int(alt) - 1
                                    else:
                                        assert alt in Vars
                                        alt_type, alt_pos, alt_data = Vars[alt]
                                        alt_right_pos = alt_pos
                                        if alt_type == "single":
                                            alt_right_pos += 1
                                        else:
                                            assert alt_type == "deletion"
                                            alt_del_len = int(alt_data)
                                            alt_right_pos += alt_del_len
                                            alt_total_del_len += alt_del_len
                                alt_right_pos -= alt_total_del_len
                                    
                        if right_pos <= alt_right_pos:
                            if verbose >= 2:
                                print "RIGHT:", cmp_list
                                print "\t", type, "id_str:", id_str, "=>", alts_id_str, "right_pos:", right_pos, "alt_right_pos:", alt_right_pos
                            cmp_right = i - 1
                            break
        i += 1

    return cmp_left, cmp_right


"""
Example,
   gene_name, allele_name (input): A, A*32:01:01
   allele (output): single-136-G-hv47,deletion-285-1-hv57, ... ,single-3473-T-hv1756,deletion-3495-1-hv1763,single-3613-C-hv1799 
"""
def get_allele(gene_name, allele_name, Vars, Var_list, Links):    
    allele_haplotype = []
    for _var_pos, _var_id in Var_list[gene_name]:
        if allele_name in Links[_var_id]:
            _var = Vars[gene_name][_var_id]
            allele_haplotype.append("%s-%d-%s-%s" % (_var[0], _var[1], _var[2], _var_id))                                
    allele_haplotype = ','.join(allele_haplotype)
    return allele_haplotype


"""
"""
def calculate_allele_coverage(allele_haplotype,
                              N_haplotypes,
                              exons,
                              partial,
                              exonic_only,
                              output):
    _var_count = {}
    for read_haplotypes in N_haplotypes.values():
        for read_haplotype in read_haplotypes:
            if haplotype_cmp(allele_haplotype, read_haplotype) <= 0:
                _, assembled = assemble_two_haplotypes(allele_haplotype.split(','), read_haplotype.split(','))
            else:
                _, assembled = assemble_two_haplotypes(read_haplotype.split(','), allele_haplotype.split(','))
            read_vars = read_haplotype.split(',')
            for read_var in read_vars:
                _type, _left, _data, _id = read_var.split('-')
                if _type in ["left", "right", "unknown"]:
                    continue
                if _id not in _var_count:
                    _var_count[_id] = 1
                else:
                    _var_count[_id] += 1
    total_var, covered_var = 0, 0
    for allele_var in allele_haplotype.split(',')[1:-1]:
        _type, _left, _data, _id = allele_var.split('-')
        _left = int(_left)
        _count = 0

        if partial and \
                exonic_only and \
                not var_in_exon([_type, _left, _data], exons):
            continue
        
        total_var += 1
        if _id in _var_count:
            _count = _var_count[_id]
            covered_var += 1
        if output:
            print "\t %d %s %s (%s - %d)" % (_left, _type, _data, _id, _count)
            
    return covered_var, total_var


"""
HISAT-genotype's mpileup
"""
def get_mpileup(alignview_cmd,
                ref_seq,
                base_locus,
                vars,
                allow_discordant):
    ref_seq_len = len(ref_seq)
    mpileup = []
    for i in range(ref_seq_len):
        mpileup.append([[], {}])
        
    proc = subprocess.Popen(alignview_cmd,
                            stdout=subprocess.PIPE,
                            stderr=open("/dev/null", 'w'))

    prev_pos = -1
    cigar_re = re.compile('\d+\w')
    for line in proc.stdout:
        line = line.strip()
        cols = line.split()
        read_id, flag, _, pos, _, cigar_str = cols[:6]
        read_seq = cols[9]
        flag, pos = int(flag), int(pos)
        # Unalined?
        if flag & 0x4 != 0:
            continue
        pos -= (base_locus + 1)
        if pos < 0:
            continue

        # Concordantly mapped?
        if flag & 0x2 != 0:
            concordant = True
        else:
            concordant = False

        if not allow_discordant and not concordant:
            continue                    

        read_pos, left_pos = 0, pos
        right_pos = left_pos
        cigars = cigar_re.findall(cigar_str)
        cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
        for i in range(len(cigars)):
            cigar_op, length = cigars[i]
            if cigar_op in "MD":
                for j in range(length):
                    if cigar_op == 'M':
                        read_nt = read_seq[read_pos + j]
                    else:
                        read_nt = 'D'
                    if read_nt not in mpileup[right_pos + j][1]:
                        mpileup[right_pos + j][1][read_nt] = 1
                    else:
                        mpileup[right_pos + j][1][read_nt] += 1

            if cigar_op in "MND":
                right_pos += length

            if cigar_op in "MIS":
                read_pos += length

    # Choose representative bases or 'D'
    for i in range(len(mpileup)):
        nt_dic = mpileup[i][1]
        num_nt = sum(nt_dic.values())
        nt_set = []
        if num_nt >= 20:
            for nt, count in nt_dic.items():
                if nt not in "ACGT":
                    continue
                if count >= num_nt * 0.25 or count >= 8:
                    nt_set.append(nt)
        mpileup[i][0] = nt_set

    # Sort variants
    var_list = [[] for i in range(len(mpileup))]
    for var_id, value in vars.items():
        var_type, var_pos, var_data = value
        assert var_pos < len(var_list)
        var_list[var_pos].append([var_id, var_type, var_data])

    # Assign known or unknown variants
    skip_i, prev_del_var_id = -1, ""
    for i in range(len(mpileup)):
        nt_dic = mpileup[i][1]
        ref_nt = ref_seq[i]
        var_set = set()
        for nt, count in nt_dic.items():
            if nt == 'D':
                if i <= skip_i:
                    assert prev_del_var_id != ""
                    var_set.add(prev_del_var_id)
                    continue
                for var_id, var_type, var_data in var_list[i]:
                    if var_type != "deletion":
                        continue
                    del_len = int(var_data)
                    del_exist = True
                    for j in range(i + 1, i + del_len):
                        assert j < len(mpileup)
                        nt_dic2 = mpileup[j][1]
                        if 'D' not in nt_dic2:
                            del_exist = False
                            break
                    if del_exist:
                        var_set.add(var_id)
                        prev_del_var_id = var_id
                        skip_i = i + del_len - 1
                        break                                                
            elif nt != 'N' and nt != ref_nt:
                assert nt in "ACGT"
                id = "unknown"
                for var_id, var_type, var_data in var_list[i]:
                    if var_type != "single":
                        continue
                    if nt == var_data:
                        id = var_id
                        break
                var_set.add(id)
                        
        mpileup[i].append(var_set)

    return mpileup


"""
"""
def error_correct(ref_seq,
                  read_seq,
                  read_pos,
                  mpileup,
                  Vars,
                  Var_list,
                  cmp_list,
                  debug = False):
    if debug:
        print cmp_list
        print read_seq

    i = 0
    while i < len(cmp_list):
        type, left, length = cmp_list[i][:3]
        assert length > 0
        if type == "match":
            middle_cmp_list = []
            last_j = 0
            for j in range(length):
                read_bp, ref_bp = read_seq[read_pos + j], ref_seq[left + j]
                assert left + j < len(mpileup)
                nt_set = mpileup[left + j][0]
                if len(nt_set) > 0 and read_bp not in nt_set:
                    read_bp = 'N' if len(nt_set) > 1 else nt_set[0]                    
                    read_seq = read_seq[:read_pos + j] + read_bp + read_seq[read_pos + j + 1:]
                    assert read_bp != ref_bp
                    new_cmp = ["mismatch", left + j, 1, "unknown"]
                    if read_bp != 'N':
                        var_idx = lower_bound(Var_list, left + j)
                        while var_idx < len(Var_list):
                            var_pos, var_id = Var_list[var_idx]
                            if var_pos > left + j:
                                break
                            if var_pos == left + j:
                                var_type, _, var_data = Vars[var_id]
                                if var_type == "single" and read_bp == var_data:
                                    new_cmp[3] = var_id
                                    break                                                        
                            var_idx += 1
                    if j > last_j:
                        middle_cmp_list.append(["match", left + last_j, j- last_j])
                    middle_cmp_list.append(new_cmp)
                    last_j = j + 1
            if last_j < length:
                middle_cmp_list.append(["match", left + last_j, length - last_j])

            assert len(middle_cmp_list) > 0
            cmp_list = cmp_list[:i] + middle_cmp_list + cmp_list[i+1:]
            i += (len(middle_cmp_list) - 1)
        else:
            assert type == "mismatch"
            read_bp, ref_bp = read_seq[read_pos], ref_seq[left]
            assert left < len(mpileup)
            nt_set = mpileup[left][0]

            if debug:
                print left, read_bp, ref_bp, mpileup[left]

            if len(nt_set) > 0 and read_bp not in nt_set:
                read_bp = 'N' if len(nt_set) > 1 else nt_set[0]
                read_seq = read_seq[:read_pos] + read_bp + read_seq[read_pos+1:]
                if read_bp == 'N':
                    cmp_list[i][3] = "unknown"
                elif read_bp == ref_bp:
                    cmp_list[i] = ["match", left, 1]
                else:
                    cmp_list[i][3] = "unknown"
                    var_idx = lower_bound(Var_list, left)
                    while var_idx < len(Var_list):
                        var_pos, var_id = Var_list[var_idx]
                        if var_pos > left:
                            break
                        if var_pos == left:
                            var_type, _, var_data = Vars[var_id]
                            if var_type == "single" and read_bp == var_data:
                                cmp_list[i][3] = var_id
                                break                                                        
                        var_idx += 1

                if debug:
                    print left, read_bp, ref_bp, mpileup[left]
                    print cmp_list[i]

        read_pos += length
        i += 1

    # Combine matches
    i = 0
    while i < len(cmp_list):
        type, left, length = cmp_list[i][:3]
        if type == "match" and i + 1 < len(cmp_list):
            type2, left2, length2 = cmp_list[i+1][:3]
            if type2 == "match":
                cmp_list[i] = [type, left, length + length2]
                cmp_list = cmp_list[:i+1] + cmp_list[i+2:]
                continue
        i += 1

    if debug:
        print cmp_list
        print read_seq
                            
    return cmp_list, read_seq


"""
"""
def typing(ex_path,
           simulation,
           reference_type,
           hla_list,
           partial,
           partial_alleles,
           refHLAs,
           HLAs,
           HLA_names,
           HLA_lengths,
           refHLA_loci,
           Vars,
           Var_list,
           Links,
           aligners,
           num_editdist,
           assembly,
           error_correction,
           allow_discordant,
           display_alleles,
           fastq,
           read_fname,
           alignment_fname,
           num_frag_list,
           threads,
           best_alleles,
           verbose):
    if simulation:
        test_passed = {}

    for aligner, index_type in aligners:
        if index_type == "graph":
            print >> sys.stderr, "\n\t\t%s %s on %s" % (aligner, index_type, reference_type)
        else:
            print >> sys.stderr, "\n\t\t%s %s" % (aligner, index_type)

        remove_alignment_file = False
        if alignment_fname == "":
            # Align reads, and sort the alignments into a BAM file
            remove_alignment_file = True
            if simulation:
                alignment_fname = "hla_output.bam"
            else:
                alignment_fname = read_fname[0].split('/')[-1]
                alignment_fname = alignment_fname.split('.')[0] + ".bam"
                
            align_reads(ex_path,
                        aligner,
                        simulation,
                        index_type,
                        read_fname,
                        fastq,
                        threads,
                        alignment_fname,
                        verbose)
            
        for test_HLA_names in hla_list:
            if simulation:
                gene = test_HLA_names[0].split('*')[0]
            else:
                gene = test_HLA_names
            ref_allele = refHLAs[gene]
            ref_seq = HLAs[gene][ref_allele]
            ref_exons = refHLA_loci[gene][-1]

            novel_var_count = 0        
            gene_vars, gene_var_list = deepcopy(Vars[gene]), deepcopy(Var_list[gene])
            var_count = {}
            def add_novel_var(gene_vars,
                              gene_var_list,
                              novel_var_count,
                              var_type,
                              var_pos,
                              var_data):
                var_idx = lower_bound(gene_var_list, var_pos)
                while var_idx < len(gene_var_list):
                    pos_, id_ = gene_var_list[var_idx]
                    if pos_ > var_pos:
                        break
                    if pos_ == var_pos:
                        type_, _, data_ = gene_vars[id_]
                        assert type_ != var_type or data_ != var_data
                        if type_ != var_type:
                            if var_type == "insertion":
                                break
                            elif var_type == "single" and type_ == "deletion":
                                break
                        else:
                            if var_data < data_:
                                break
                    var_idx += 1
                var_id = "nv%d" % novel_var_count
                assert var_id not in gene_vars
                gene_vars[var_id] = [var_type, var_pos, var_data]
                gene_var_list.insert(var_idx, [var_pos, var_id])                
                return var_id, novel_var_count + 1

            if not os.path.exists(alignment_fname + ".bai"):
                os.system("samtools index %s" % alignment_fname)
            # Read alignments
            alignview_cmd = ["samtools",
                             "view",
                             alignment_fname]
            base_locus = 0
            if index_type == "graph":
                if reference_type == "gene":
                    alignview_cmd += [ref_allele]
                else:
                    assert reference_type in ["chromosome", "genome"]
                    _, chr, left, right, _ = refHLA_loci[gene]
                    base_locus = left
                    alignview_cmd += ["%s:%d-%d" % (chr, left + 1, right + 1)]

                mpileup = get_mpileup(alignview_cmd,
                                      ref_seq,
                                      base_locus,
                                      gene_vars,
                                      allow_discordant)

                bamview_proc = subprocess.Popen(alignview_cmd,
                                                stdout=subprocess.PIPE,
                                                stderr=open("/dev/null", 'w'))

                sort_read_cmd = ["sort", "-k", "1,1", "-s"] # -s for stable sorting
                alignview_proc = subprocess.Popen(sort_read_cmd,
                                                  stdin=bamview_proc.stdout,
                                                  stdout=subprocess.PIPE,
                                                  stderr=open("/dev/null", 'w'))
            else:
                alignview_proc = subprocess.Popen(alignview_cmd,
                                             stdout=subprocess.PIPE,
                                             stderr=open("/dev/null", 'w'))

            # List of nodes that represent alleles
            allele_vars = {}
            for var_id, allele_list in Links.items():
                for allele_id in allele_list:
                    if allele_id not in HLAs[gene]:
                        continue
                    if allele_id not in allele_vars:
                        allele_vars[allele_id] = [var_id]
                    else:
                        allele_vars[allele_id].append(var_id)

            # Extract variants that are within exons
            exon_vars = get_exonic_vars(gene_vars, ref_exons)

            # Compare two alleles
            allele_nodes = {}
            for allele_id, var_ids in allele_vars.items():
                seq = list(ref_seq)  # sequence that node represents
                var = ["" for i in range(len(ref_seq))]  # how sequence is related to backbone
                for var_id in var_ids:
                    assert var_id in gene_vars
                    var_type, var_pos, var_data = gene_vars[var_id]
                    assert var_pos >= 0 and var_pos < len(ref_seq)
                    if var_type == "single":
                        seq[var_pos] = var_data
                        var[var_pos] = var_id
                    elif var_type == "deletion":
                        assert var_type == "deletion"
                        del_len = int(var_data)
                        assert var_pos + del_len <= len(ref_seq)
                        seq[var_pos:var_pos + del_len] = ['D'] * del_len
                        var[var_pos:var_pos + del_len] = [var_id] * del_len
                    else:
                        # DK - to be implemented for insertions
                        None

                qual = ' ' * len(seq)
                allele_nodes[allele_id] = assembly_graph.Node(allele_id,
                                                              0,
                                                              seq,
                                                              qual,
                                                              var,
                                                              ref_seq,
                                                              gene_vars,
                                                              mpileup,
                                                              simulation)

            true_allele_nodes = {}
            if simulation:
                for allele_name in test_HLA_names:
                    true_allele_nodes[allele_name] = allele_nodes[allele_name]

            display_allele_nodes = {}
            for display_allele in display_alleles:
                if display_allele not in allele_nodes:
                    continue
                display_allele_nodes[display_allele] = allele_nodes[display_allele]

            # Assembly graph
            asm_graph = assembly_graph.Graph(ref_seq,
                                             gene_vars,
                                             ref_exons,
                                             partial_alleles,
                                             true_allele_nodes,
                                             display_allele_nodes,
                                             simulation)

            # Choose allele representives from those that share the same exonic sequences
            allele_reps, allele_rep_groups = get_rep_alleles(Links, exon_vars)
            allele_rep_set = set(allele_reps.values())

            # For checking alternative alignments near the ends of alignments
            Alts_left, Alts_right = get_alternatives(ref_seq, gene_vars, gene_var_list, verbose)

            # Count alleles
            HLA_counts, HLA_cmpt = {}, {}
            HLA_gen_counts, HLA_gen_cmpt = {}, {}
            num_reads, total_read_len = 0, 0

            # For debugging purposes
            debug_allele_names = set(test_HLA_names) if simulation and verbose >= 2 else set()

            # Read information
            prev_read_id = None
            prev_right_pos = 0
            prev_lines = []
            left_read_ids, right_read_ids = set(), set()
            if index_type == "graph":
                # nodes for reads
                read_nodes = []
                read_vars_list = []
                
                # Cigar regular expression
                cigar_re = re.compile('\d+\w')
                for line in alignview_proc.stdout:
                    line = line.strip()
                    cols = line.split()
                    read_id, flag, chr, pos, mapQ, cigar_str = cols[:6]
                    node_read_id = orig_read_id = read_id
                    if simulation:
                        read_id = read_id.split('|')[0]
                    read_seq, read_qual = cols[9], cols[10]
                    total_read_len += len(read_seq)
                    flag, pos = int(flag), int(pos)
                    pos -= (base_locus + 1)
                    if pos < 0:
                        continue

                    # Unalined?
                    if flag & 0x4 != 0:
                        if simulation and verbose >= 2:
                            print "Unaligned"
                            print "\t", line                            
                        continue

                    # Concordantly mapped?
                    if flag & 0x2 != 0:
                        concordant = True
                    else:
                        concordant = False

                    NM, Zs, MD, NH = "", "", "", ""
                    for i in range(11, len(cols)):
                        col = cols[i]
                        if col.startswith("Zs"):
                            Zs = col[5:]
                        elif col.startswith("MD"):
                            MD = col[5:]
                        elif col.startswith("NM"):
                            NM = int(col[5:])
                        elif col.startswith("NH"):
                            NH = int(col[5:])

                    if NM > num_editdist:
                        continue

                    # Only consider unique alignment
                    if NH > 1:
                        continue

                    # Concordantly aligned mate pairs
                    if not allow_discordant and not concordant:
                        continue

                    # Left read?
                    if flag & 0x40 != 0:
                        if read_id in left_read_ids:
                            continue
                        left_read_ids.add(read_id)
                        if not simulation:
                            node_read_id += '|L'
                    else: # Right read?
                        assert flag & 0x80 != 0
                        if read_id in right_read_ids:
                            continue
                        right_read_ids.add(read_id)
                        if not simulation:
                            node_read_id += '|R'

                    if Zs:
                        Zs = Zs.split(',')             

                    assert MD != ""
                    MD_str_pos, MD_len = 0, 0
                    Zs_pos, Zs_i = 0, 0
                    for _i in range(len(Zs)):
                        Zs[_i] = Zs[_i].split('|')
                        Zs[_i][0] = int(Zs[_i][0])
                    if Zs_i < len(Zs):
                        Zs_pos += Zs[Zs_i][0]
                    read_pos, left_pos = 0, pos
                    right_pos = left_pos
                    cigars = cigar_re.findall(cigar_str)
                    cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
                    cmp_list = []

                    likely_misalignment = False

                    # Extract variants w.r.t backbone from CIGAR string
                    softclip = [0, 0]
                    for i in range(len(cigars)):
                        cigar_op, length = cigars[i]
                        if cigar_op == 'M':
                            first = True
                            MD_len_used = 0
                            cmp_list_i = len(cmp_list)
                            while True:
                                if not first or MD_len == 0:
                                    if MD[MD_str_pos].isdigit():
                                        num = int(MD[MD_str_pos])
                                        MD_str_pos += 1
                                        while MD_str_pos < len(MD):
                                            if MD[MD_str_pos].isdigit():
                                                num = num * 10 + int(MD[MD_str_pos])
                                                MD_str_pos += 1
                                            else:
                                                break
                                        MD_len += num
                                # Insertion or full match followed
                                if MD_len >= length:
                                    MD_len -= length
                                    if length > MD_len_used:
                                        cmp_list.append(["match", right_pos + MD_len_used, length - MD_len_used])
                                    break
                                first = False
                                read_base = read_seq[read_pos + MD_len]
                                MD_ref_base = MD[MD_str_pos]
                                MD_str_pos += 1
                                assert MD_ref_base in "ACGT"
                                if MD_len > MD_len_used:
                                    cmp_list.append(["match", right_pos + MD_len_used, MD_len - MD_len_used])

                                _var_id = "unknown"
                                if read_pos + MD_len == Zs_pos and Zs_i < len(Zs):
                                    assert Zs[Zs_i][1] == 'S'
                                    _var_id = Zs[Zs_i][2]
                                    Zs_i += 1
                                    Zs_pos += 1
                                    if Zs_i < len(Zs):
                                        Zs_pos += Zs[Zs_i][0]
                                else:
                                    # Search for a known (yet not indexed) variant or a novel variant
                                    ref_pos = right_pos + MD_len
                                    var_idx = lower_bound(gene_var_list, ref_pos)
                                    while var_idx < len(gene_var_list):
                                        var_pos, var_id = gene_var_list[var_idx]
                                        if var_pos > ref_pos:
                                            break
                                        if var_pos == ref_pos:
                                            var_type, _, var_data = gene_vars[var_id]
                                            if var_type == "single" and var_data == read_base:
                                                _var_id = var_id
                                                break
                                        var_idx += 1

                                cmp_list.append(["mismatch", right_pos + MD_len, 1, _var_id])
                                MD_len_used = MD_len + 1
                                MD_len += 1
                                # Full match
                                if MD_len == length:
                                    MD_len = 0
                                    break

                            # Correction for sequencing errors and update for cmp_list
                            if error_correction:
                                assert cmp_list_i < len(cmp_list)
                                new_cmp_list, read_seq = error_correct(ref_seq,
                                                                       read_seq,
                                                                       read_pos,
                                                                       mpileup,
                                                                       gene_vars,
                                                                       gene_var_list,
                                                                       cmp_list[cmp_list_i:],
                                                                       node_read_id == "#HSQ1008:176:D0UYCACXX:4:1304:19006:96208|R")
                                cmp_list = cmp_list[:cmp_list_i] + new_cmp_list                            

                        elif cigar_op == 'I':
                            _var_id = "unknown"
                            if read_pos == Zs_pos and Zs_i < len(Zs):
                                assert Zs[Zs_i][1] == 'I'
                                _var_id = Zs[Zs_i][2]
                                Zs_i += 1
                                if Zs_i < len(Zs):
                                    Zs_pos += Zs[Zs_i][0]
                            else:
                                # Search for a known (yet not indexed) variant or a novel variant
                                var_idx = lower_bound(gene_var_list, right_pos)
                                while var_idx < len(gene_var_list):
                                    var_pos, var_id = gene_var_list[var_idx]
                                    if var_pos > right_pos:
                                        break
                                    if var_pos == right_pos:
                                        var_type, _, var_data = gene_vars[var_id]
                                        if var_type == "insertion" and len(var_data) == length:
                                            _var_id = var_id
                                            break
                                    var_idx += 1                            
                            cmp_list.append(["insertion", right_pos, length, _var_id])
                            if 'N' in read_seq[read_pos:read_pos+length]:
                                likely_misalignment = True
                                
                        elif cigar_op == 'D':
                            if MD[MD_str_pos] == '0':
                                MD_str_pos += 1
                            assert MD[MD_str_pos] == '^'
                            MD_str_pos += 1
                            while MD_str_pos < len(MD):
                                if not MD[MD_str_pos] in "ACGT":
                                    break
                                MD_str_pos += 1
                            _var_id = "unknown"
                            if read_pos == Zs_pos and Zs_i < len(Zs):
                                assert Zs[Zs_i][1] == 'D'
                                _var_id = Zs[Zs_i][2]
                                Zs_i += 1
                                if Zs_i < len(Zs):
                                    Zs_pos += Zs[Zs_i][0]
                            else:
                                # Search for a known (yet not indexed) variant or a novel variant
                                var_idx = lower_bound(gene_var_list, right_pos)
                                while var_idx < len(gene_var_list):
                                    var_pos, var_id = gene_var_list[var_idx]
                                    if var_pos > right_pos:
                                        break
                                    if var_pos == right_pos:
                                        var_type, _, var_data = gene_vars[var_id]
                                        if var_type == "deletion" and int(var_data) == length:
                                            _var_id = var_id
                                            break
                                    var_idx += 1

                            cmp_list.append(["deletion", right_pos, length, _var_id])

                            # Check if this deletion is artificial alignment
                            assert right_pos < mpileup
                            del_count, nt_count = 0, 0
                            for nt, count in mpileup[right_pos][1].items():
                                if nt == 'D':
                                    del_count += count
                                else:
                                    nt_count += count
                            # DK - debugging purposes
                            if del_count * 6 < nt_count: # and nt_count >= 15:
                                likely_misalignment = True
                            
                        elif cigar_op == 'S':
                            if i == 0:
                                softclip[0] = length
                            else:
                                assert i + 1 == len(cigars)
                                softclip[1] = length
                        else:                    
                            assert cigar_op == 'N'
                            assert False
                            cmp_list.append(["intron", right_pos, length])

                        if cigar_op in "MND":
                            right_pos += length

                        if cigar_op in "MIS":
                            read_pos += length                     

                    # Remove softclip in cigar and modify read_seq and read_qual accordingly
                    if sum(softclip) > 0:
                        if softclip[0] > 0:
                            cigars = cigars[1:]
                            read_seq = read_seq[softclip[0]:]
                            read_qual = read_qual[softclip[0]:]
                        if softclip[1] > 0:
                            cigars = cigars[:-1]
                            read_seq = read_seq[:-softclip[1]]
                            read_qual = read_qual[:-softclip[1]]

                        cigar_str = ""
                        for type, length in cigars:
                            cigar_str += str(length)
                            cigar_str += type
                    
                    if right_pos > len(ref_seq):
                        continue

                    if likely_misalignment:
                        continue

                    # Add novel variants
                    read_pos = 0
                    for cmp_i in range(len(cmp_list)):
                        type_, pos_, length_ = cmp_list[cmp_i][:3]
                        if type_ != "match":
                            var_id_ = cmp_list[cmp_i][3]
                            if var_id_ == "unknown":
                                add = True
                                if type_ == "mismatch":
                                    data_ = read_seq[read_pos]
                                    if data_ == 'N':
                                        add = False
                                elif type_ == "deletion":
                                    data_ = str(length_)
                                else:
                                    assert type_ == "insertion"
                                    data_ = read_seq[read_pos:read_pos + length_]
                                if add:
                                    var_id, novel_var_count = add_novel_var(gene_vars,
                                                                            gene_var_list,
                                                                            novel_var_count,
                                                                            type_ if type_ != "mismatch" else "single",
                                                                            pos_,
                                                                            data_)
                                    cmp_list[cmp_i][3] = var_id
                            if var_id not in var_count:
                                var_count[var_id] = 1
                            else:
                                var_count[var_id] += 1
                                
                        if type_ != "deletion":
                            read_pos += length_

                    # Count the number of reads aligned uniquely with some constraints
                    num_reads += 1

                    def add_stat(HLA_cmpt, HLA_counts, HLA_count_per_read, include_alleles = set()):
                        max_count = max(HLA_count_per_read.values())
                        cur_cmpt = set()
                        for allele, count in HLA_count_per_read.items():
                            if count < max_count:
                                continue

                            if len(include_alleles) > 0 and allele not in include_alleles:
                                continue
                            
                            cur_cmpt.add(allele)                    
                            if allele not in HLA_counts:
                                HLA_counts[allele] = 1
                            else:
                                HLA_counts[allele] += 1

                        if len(cur_cmpt) == 0:
                            return ""

                        # DK - for debugging purposes                            
                        alleles = ["", ""]
                        # alleles = ["A*24:36N", "A*24:359N"]
                        allele1_found, allele2_found = False, False
                        if alleles[0] != "":
                            for allele, count in HLA_count_per_read.items():
                                if count < max_count:
                                    continue
                                if allele == alleles[0]:
                                    allele1_found = True
                                elif allele == alleles[1]:
                                    allele2_found = True
                            if allele1_found != allele2_found:
                                print alleles[0], HLA_count_per_read[alleles[0]]
                                print alleles[1], HLA_count_per_read[alleles[1]]
                                if allele1_found:
                                    print ("%s\tread_id %s - %d vs. %d]" % (alleles[0], prev_read_id, max_count, HLA_count_per_read[alleles[1]]))
                                else:
                                    print ("%s\tread_id %s - %d vs. %d]" % (alleles[1], prev_read_id, max_count, HLA_count_per_read[alleles[0]]))
                                print read_seq

                        cur_cmpt = sorted(list(cur_cmpt))
                        cur_cmpt = '-'.join(cur_cmpt)
                        if not cur_cmpt in HLA_cmpt:
                            HLA_cmpt[cur_cmpt] = 1
                        else:
                            HLA_cmpt[cur_cmpt] += 1

                        return cur_cmpt

                    if read_id != prev_read_id:
                        if prev_read_id != None:
                            cur_cmpt = add_stat(HLA_cmpt, HLA_counts, HLA_count_per_read, allele_rep_set)
                            add_stat(HLA_gen_cmpt, HLA_gen_counts, HLA_gen_count_per_read)
                            for read_id_, read_node in read_nodes:
                                asm_graph.add_node(read_id_,
                                                   read_node,
                                                   simulation)
                            read_nodes, read_var_list = [], []

                            if verbose >= 2:
                                cur_cmpt = cur_cmpt.split('-')
                                if not(set(cur_cmpt) & set(test_HLA_names)):
                                    print "%s are chosen instead of %s" % ('-'.join(cur_cmpt), '-'.join(test_HLA_names))
                                    for prev_line in prev_lines:
                                        print "\t", prev_line

                            prev_lines = []

                        HLA_count_per_read, HLA_gen_count_per_read = {}, {}
                        for HLA_name in HLA_names[gene]:
                            if HLA_name.find("BACKBONE") != -1:
                                continue
                            HLA_count_per_read[HLA_name] = 0
                            HLA_gen_count_per_read[HLA_name] = 0

                    prev_lines.append(line)

                    def add_count(count_per_read, var_id, add):
                        alleles = Links[var_id]
                        if verbose >= 2:
                            if add > 0 and not (set(alleles) & debug_allele_names):
                                print "Add:", add, debug_allele_names, "-", var_id
                                print "\t", line
                                print "\t", alleles
                            if add < 0 and set(alleles) & debug_allele_names:
                                print "Add:", add, debug_allele_names, "-", var_id
                                print "\t", line

                        for allele in alleles:
                            count_per_read[allele] += add

                    # Decide which allele(s) a read most likely came from
                    for var_id, data in gene_vars.items():
                        if var_id == "unknown" or var_id.startswith("nv"):
                            continue
                        var_type, var_pos, var_data = data
                        if var_type != "deletion":
                            continue
                        if left_pos >= var_pos and right_pos <= var_pos + int(var_data):
                            if var_id in exon_vars:
                                add_count(HLA_count_per_read, var_id, -1)
                            add_count(HLA_gen_count_per_read, var_id, -1)

                    # Node
                    read_node_pos, read_node_seq, read_node_qual, read_node_var = -1, [], [], []
                    read_vars = []

                    # Positive and negative evidence
                    positive_vars, negative_vars = set(), set()

                    # Sanity check - read length, cigar string, and MD string
                    ref_pos, read_pos, cmp_cigar_str, cmp_MD = left_pos, 0, "", ""
                    cigar_match_len, MD_match_len = 0, 0

                    cmp_list_left, cmp_list_right = identify_ambigious_diffs(gene_vars,
                                                                             Alts_left,
                                                                             Alts_right,
                                                                             cmp_list,
                                                                             verbose)

                    # Deletions at 5' and 3' ends
                    for var_id, data in gene_vars.items():
                        var_type, var_pos, var_data = data
                        if var_type != "deletion":
                            continue
                        if left_pos >= var_pos and right_pos <= var_pos + int(var_data):
                            negative_vars.add(var_id)
                    
                    cmp_i = 0
                    while cmp_i < len(cmp_list):
                        cmp = cmp_list[cmp_i]
                        type, length = cmp[0], cmp[2]
                        # Disable the following sanity check due to error correction
                        # if num_editdist == 0 and type in ["mismatch", "deletion", "insertion"]:
                        #     assert cmp[3] != "unknown"

                        if type in ["match", "mismatch"]:
                            if read_node_pos < 0:
                                read_node_pos = ref_pos

                        if type == "match":
                            read_node_seq += list(read_seq[read_pos:read_pos+length])
                            read_node_qual += list(read_qual[read_pos:read_pos+length])
                            read_node_var += ([''] * length)
                            
                            var_idx = lower_bound(gene_var_list, ref_pos)
                            while var_idx < len(gene_var_list):
                                var_pos, var_id = gene_var_list[var_idx]
                                if ref_pos + length <= var_pos:
                                    break
                                if ref_pos <= var_pos:
                                    var_type, _, var_data = gene_vars[var_id]
                                    if var_type == "insertion":
                                        if ref_pos < var_pos and ref_pos + length > var_pos + len(var_data):
                                            negative_vars.add(var_id)
                                    elif var_type == "deletion":
                                        del_len = int(var_data)
                                        if ref_pos < var_pos and ref_pos + length > var_pos + del_len:
                                            # Check if this might be one of the two tandem repeats (the same left coordinate)
                                            cmp_left, cmp_right = cmp[1], cmp[1] + cmp[2]
                                            test1_seq1 = ref_seq[cmp_left:cmp_right]
                                            test1_seq2 = ref_seq[cmp_left:var_pos] + ref_seq[var_pos + del_len:cmp_right + del_len]
                                            # Check if this happens due to small repeats (the same right coordinate - e.g. 19 times of TTTC in DQA1*05:05:01:02)
                                            cmp_left -= read_pos
                                            cmp_right += (len(read_seq) - read_pos - cmp[2])
                                            test2_seq1 = ref_seq[cmp_left+int(var_data):cmp_right]
                                            test2_seq2 = ref_seq[cmp_left:var_pos] + ref_seq[var_pos+int(var_data):cmp_right]
                                            if test1_seq1 != test1_seq2 and test2_seq1 != test2_seq2:
                                                negative_vars.add(var_id)
                                    else:
                                        negative_vars.add(var_id)
                                var_idx += 1
                            read_pos += length
                            ref_pos += length
                            cigar_match_len += length
                            MD_match_len += length
                        elif type == "mismatch":
                            var_id = cmp[3]
                            read_base, qual = read_seq[read_pos], read_qual[read_pos]
                            read_node_seq += [read_base]
                            read_node_qual += [qual]
                            read_node_var.append(var_id)
                            if var_id != "unknown":
                                if cmp_i >= cmp_list_left and cmp_i <= cmp_list_right:
                                    positive_vars.add(var_id)
                            
                            cmp_MD += ("%d%s" % (MD_match_len, ref_seq[ref_pos]))
                            MD_match_len = 0
                            cigar_match_len += 1
                            read_pos += 1
                            ref_pos += 1
                        elif type == "insertion":
                            var_id = cmp[3]
                            ins_len = length
                            ins_seq = read_seq[read_pos:read_pos+ins_len]
                            if var_id != "unknown" or not var_id.startswith("nv"):
                                if cmp_i >= cmp_list_left and cmp_i <= cmp_list_right:
                                    # Require at least 5bp match before and after a deletion
                                    if read_pos >= 5 and read_pos + 5 <= len(read_seq):
                                        positive_vars.add(var_id)
                            read_node_seq += ["I%s" % nt for nt in ins_seq]
                            read_node_qual += list(read_qual[read_pos:read_pos+ins_len])
                            read_node_var += ([var_id] * ins_len)                                        
                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            read_pos += length
                            cmp_cigar_str += ("%dI" % length)
                        elif type == "deletion":
                            var_id = cmp[3]
                            alt_match = False
                            del_len = length
                            read_node_seq += (['D'] * del_len)
                            read_node_qual += ([''] * del_len)
                            if var_id != "unknown" or not var_id.statswith("nv"):
                                if cmp_i >= cmp_list_left and cmp_i <= cmp_list_right:
                                    # Require at least 5bp match before and after a deletion
                                    if read_pos >= 5 and read_pos + 5 <= len(read_seq):
                                        positive_vars.add(var_id)

                            if len(read_node_seq) > len(read_node_var):
                                assert len(read_node_seq) == len(read_node_var) + del_len
                                read_node_var += ([var_id] * del_len)

                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            cmp_MD += ("%d" % MD_match_len)
                            MD_match_len = 0
                            cmp_cigar_str += ("%dD" % length)
                            cmp_MD += ("^%s" % ref_seq[ref_pos:ref_pos+length])
                            ref_pos += length
                        else:
                            assert type == "intron"
                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            cmp_cigar_str += ("%dN" % length)
                            ref_pos += length

                        cmp_i += 1
                
                    if cigar_match_len > 0:
                        cmp_cigar_str += ("%dM" % cigar_match_len)
                    cmp_MD += ("%d" % MD_match_len)
                    # Sanity check
                    if read_pos != len(read_seq) or \
                            cmp_cigar_str != cigar_str:
                            # cmp_MD != MD: # Disabled due to error correction
                        print >> sys.stderr, "Error:", cigar_str, MD
                        print >> sys.stderr, "\tcomputed:", cmp_cigar_str, cmp_MD
                        print >> sys.stderr, "\tcmp list:", cmp_list
                        assert False

                    # Node
                    read_nodes.append([node_read_id,
                                       assembly_graph.Node(node_read_id,
                                                           read_node_pos,
                                                           read_node_seq,
                                                           read_node_qual,
                                                           read_node_var,
                                                           ref_seq,
                                                           gene_vars,
                                                           mpileup,
                                                           simulation)])

                    for positive_var in positive_vars:
                        if positive_var == "unknown" or positive_var.startswith("nv"):
                            continue
                        if positive_var in exon_vars:
                            add_count(HLA_count_per_read, positive_var, 1)
                        add_count(HLA_gen_count_per_read, positive_var, 1)
                    for negative_var in negative_vars:
                        if negative_var == "unknown" or negative_var.startswith("nv"):
                            continue
                        if negative_var in exon_vars:
                            add_count(HLA_count_per_read, negative_var, -1)
                        add_count(HLA_gen_count_per_read, negative_var, -1)

                    prev_read_id = read_id
                    prev_right_pos = right_pos

                if num_reads <= 0:
                    continue

                print >> sys.stderr, "\t\t\tNumber of reads aligned: %d" % num_reads

                if prev_read_id != None:
                    add_stat(HLA_cmpt, HLA_counts, HLA_count_per_read, allele_rep_set)
                    add_stat(HLA_gen_cmpt, HLA_gen_counts, HLA_gen_count_per_read)
                    for read_id_, read_node in read_nodes:
                        asm_graph.add_node(read_id_,
                                           read_node,
                                           simulation)
                    read_nodes, read_var_list = [], []
                
            else:
                assert index_type == "linear"
                def add_alleles(alleles):
                    if not allele in HLA_counts:
                        HLA_counts[allele] = 1
                    else:
                        HLA_counts[allele] += 1

                    cur_cmpt = sorted(list(alleles))
                    cur_cmpt = '-'.join(cur_cmpt)
                    if not cur_cmpt in HLA_cmpt:
                        HLA_cmpt[cur_cmpt] = 1
                    else:
                        HLA_cmpt[cur_cmpt] += 1

                prev_read_id, prev_AS = None, None
                alleles = set()
                for line in alignview_proc.stdout:
                    cols = line[:-1].split()
                    read_id, flag, allele = cols[:3]
                    flag = int(flag)
                    if flag & 0x4 != 0:
                        continue
                    if not allele.startswith(gene):
                        continue
                    if allele.find("BACKBONE") != -1:
                        continue

                    AS = None
                    for i in range(11, len(cols)):
                        col = cols[i]
                        if col.startswith("AS"):
                            AS = int(col[5:])
                    assert AS != None
                    if read_id != prev_read_id:
                        if alleles:
                            if aligner == "hisat2" or \
                                    (aligner == "bowtie2" and len(alleles) < 10):
                                add_alleles(alleles)
                            alleles = set()
                        prev_AS = None
                    if prev_AS != None and AS < prev_AS:
                        continue
                    prev_read_id = read_id
                    prev_AS = AS
                    alleles.add(allele)

                if alleles:
                    add_alleles(alleles)

            HLA_counts = [[allele, count] for allele, count in HLA_counts.items()]
            def HLA_count_cmp(a, b):
                if a[1] != b[1]:
                    return b[1] - a[1]
                assert a[0] != b[0]
                if a[0] < b[0]:
                    return -1
                else:
                    return 1
            HLA_counts = sorted(HLA_counts, cmp=HLA_count_cmp)
            for count_i in range(len(HLA_counts)):
                count = HLA_counts[count_i]
                if simulation:
                    found = False
                    for test_HLA_name in test_HLA_names:
                        if count[0] == test_HLA_name:
                            print >> sys.stderr, "\t\t\t*** %d ranked %s (count: %d)" % (count_i + 1, test_HLA_name, count[1])
                            found = True
                            """
                            if count_i > 0 and HLA_counts[0][1] > count[1]:
                                print >> sys.stderr, "Warning: %s ranked first (count: %d)" % (HLA_counts[0][0], HLA_counts[0][1])
                                assert False
                            else:
                                test_passed += 1
                            """
                    if count_i < 5 and not found:
                        print >> sys.stderr, "\t\t\t\t%d %s (count: %d)" % (count_i + 1, count[0], count[1])
                else:
                    print >> sys.stderr, "\t\t\t\t%d %s (count: %d)" % (count_i + 1, count[0], count[1])
                    if count_i >= 9:
                        break
            print >> sys.stderr

            # Calculate the abundance of representative alleles on exonic sequences
            HLA_prob = single_abundance(HLA_cmpt, HLA_lengths[gene])

            # Incorporate non representative alleles (full length alleles)
            gen_alleles = set()
            gen_prob_sum = 0.0
            for prob_i in range(len(HLA_prob)):
                allele, prob = HLA_prob[prob_i][:2]
                if prob_i >= 10 and prob < 0.03:
                    break
                if allele in partial_alleles:
                    continue
                
                gen_prob_sum += prob
                for allele2 in allele_rep_groups[allele]:
                    gen_alleles.add(allele2)
            if len(gen_alleles) > 0:
                HLA_gen_cmpt2 = {}
                for cmpt, value in HLA_gen_cmpt.items():
                    cmpt2 = []
                    for allele in cmpt.split('-'):
                        if allele in gen_alleles:
                            cmpt2.append(allele)
                    if len(cmpt2) == 0:
                        continue
                    cmpt2 = '-'.join(cmpt2)
                    if cmpt2 not in HLA_gen_cmpt2:
                        HLA_gen_cmpt2[cmpt2] = value
                    else:
                        HLA_gen_cmpt2[cmpt2] += value
                HLA_gen_cmpt = HLA_gen_cmpt2
                HLA_gen_prob = single_abundance(HLA_gen_cmpt, HLA_lengths[gene])

                HLA_combined_prob = {}
                for allele, prob in HLA_prob:
                    assert allele not in HLA_combined_prob
                    if allele in gen_alleles:
                        HLA_combined_prob[allele] = 0.0
                    else:
                        HLA_combined_prob[allele] = prob
                for allele, prob in HLA_gen_prob:
                    HLA_combined_prob[allele] = prob * gen_prob_sum
                HLA_prob = [[allele, prob] for allele, prob in HLA_combined_prob.items()]
                HLA_prob = sorted(HLA_prob, cmp=HLA_prob_cmp)

            if index_type == "graph" and assembly:
                allele_node_order = []
                if not simulation:
                    predicted_allele_nodes = {}
                    for allele_name, prob in HLA_prob:
                        assert allele_name in allele_nodes
                        assert allele_name not in predicted_allele_nodes
                        predicted_allele_nodes[allele_name] = allele_nodes[allele_name]
                        allele_node_order.append([allele_name, prob])
                        if len(predicted_allele_nodes) >= 2:
                            break
                    asm_graph.set_allele_nodes(predicted_allele_nodes)
                    asm_graph.allele_node_order = allele_node_order

                # Filter out nodes
                asm_graph.filter_nodes()

                # Generate edges
                asm_graph.generate_edges(0.5)

                # Start drawing assembly graph
                asm_graph.begin_draw("assembly_graph")

                # Draw assembly graph
                begin_y = asm_graph.draw(0, "Initial graph")
                begin_y += 200

                # Reduce graph
                asm_graph.reduce()

                # Draw assembly graph
                begin_y = asm_graph.draw(begin_y, "Unitigs")
                begin_y += 200

                # Further reduce graph with mate pairs
                asm_graph.assemble_with_mates()

                # Draw assembly graph
                begin_y = asm_graph.draw(begin_y, "Graph with mate pairs")
                begin_y += 200

                # DK - debugging purposes
                """

                asm_graph.assemble_with_alleles(allele_nodes)

                # Draw assembly graph
                begin_y = asm_graph.draw(begin_y, "Graph with alleles")


                """

                # End drawing assembly graph
                asm_graph.end_draw()

                # Compare two alleles
                if simulation and len(test_HLA_names) == 2:
                    allele_name1, allele_name2 = test_HLA_names
                    print >> sys.stderr, allele_name1, "vs.", allele_name2
                    asm_graph.print_node_comparison(asm_graph.allele_nodes)

                def compare_alleles(vars1, vars2):
                    skip = True
                    var_i, var_j = 0, 0
                    while var_i < len(vars1) and var_j < len(vars2):
                        cmp_var_id, node_var_id = vars1[var_i], vars2[var_j]
                        if cmp_var_id == node_var_id:
                            skip = False
                            var = gene_vars[cmp_var_id]
                            print >> sys.stderr, cmp_var_id, var, "\t\t\t", mpileup[var[1]]
                            var_i += 1; var_j += 1
                            continue
                        cmp_var, node_var = gene_vars[cmp_var_id], gene_vars[node_var_id]
                        if cmp_var[1] <= node_var[1]:
                            if not skip:
                                if (var_i > 0 and var_i + 1 < len(vars1)) or cmp_var[0] != "deletion":
                                    print >> sys.stderr, "***", cmp_var_id, cmp_var, "==", "\t\t\t", mpileup[cmp_var[1]]
                            var_i += 1
                        else:
                            if not skip:
                                print >> sys.stderr, "*** ==", node_var_id, node_var, "\t\t\t", mpileup[node_var[1]]
                            var_j += 1
                    
                tmp_nodes = asm_graph.nodes
                print >> sys.stderr, "Number of tmp nodes:", len(tmp_nodes)
                count = 0
                for id, node in tmp_nodes.items():
                    count += 1
                    if count > 10:
                        break
                    node_vars = node.get_var_ids(gene_vars)
                    node.print_info(); print >> sys.stderr
                    if node.id in asm_graph.to_node:
                        for id2, at in asm_graph.to_node[node.id]:
                            print >> sys.stderr, "\tat %d ==> %s" % (at, id2)

                    if simulation:
                        cmp_HLA_names = test_HLA_names
                    else:
                        cmp_HLA_names = [allele_name for allele_name, _ in allele_node_order]
                        
                    alleles, cmp_vars, max_common = "", [], -sys.maxint
                    for cmp_HLA_name in cmp_HLA_names:
                        tmp_vars = allele_nodes[cmp_HLA_name].get_var_ids(gene_vars, node.left, node.right)
                        tmp_common = len(set(node_vars) & set(tmp_vars))
                        tmp_common -= len(set(node_vars) | set(tmp_vars))
                        if max_common < tmp_common:
                            max_common = tmp_common
                            alleles = [[cmp_HLA_name, tmp_vars]]
                        elif max_common == tmp_common:
                            alleles.append([cmp_HLA_name, tmp_vars])

                    for allele_name, cmp_vars in alleles:
                        print >> sys.stderr, "vs.", allele_name
                        compare_alleles(cmp_vars, node_vars)

                    print >> sys.stderr
                    print >> sys.stderr

            success = [False for i in range(len(test_HLA_names))]
            found_list = [False for i in range(len(test_HLA_names))]
            for prob_i in range(len(HLA_prob)):
                prob = HLA_prob[prob_i]
                found = False
                _allele_rep = prob[0]
                """
                if partial and exonic_only:
                    _fields = _allele_rep.split(':')
                    if len(_fields) == 4:
                        _allele_rep = ':'.join(_fields[:-1])
                """

                if simulation:
                    for name_i in range(len(test_HLA_names)):
                        test_HLA_name = test_HLA_names[name_i]
                        if prob[0] == test_HLA_name:
                            rank_i = prob_i
                            while rank_i > 0:
                                if prob == HLA_prob[rank_i - 1][1]:
                                    rank_i -= 1
                                else:
                                    break
                            print >> sys.stderr, "\t\t\t*** %d ranked %s (abundance: %.2f%%)" % (rank_i + 1, test_HLA_name, prob[1] * 100.0)
                            if rank_i < len(success):
                                success[rank_i] = True
                            found_list[name_i] = True
                            found = True
                    # DK - for debugging purposes
                    if not False in found_list and prob_i >= 10:
                        break
                if not found:
                    print >> sys.stderr, "\t\t\t\t%d ranked %s (abundance: %.2f%%)" % (prob_i + 1, _allele_rep, prob[1] * 100.0)
                    if best_alleles and prob_i < 2:
                        print >> sys.stdout, "SingleModel %s (abundance: %.2f%%)" % (_allele_rep, prob[1] * 100.0)
                if not simulation and prob_i >= 9:
                    break
                if prob_i >= 19:
                    break
            print >> sys.stderr
            
            # DK - debugging purposes
            # sys.exit(1)

            # DK - for debugging purposes
            if False and (len(test_HLA_names) == 2 or not simulation):
                HLA_prob = joint_abundance(HLA_cmpt, HLA_lengths[gene])
                if len(HLA_prob) <= 0:
                    continue
                success = [False]
                for prob_i in range(len(HLA_prob)):
                    allele_pair, prob = HLA_prob[prob_i]
                    allele1, allele2 = allele_pair.split('-')
                    if best_alleles and prob_i < 1:
                        print >> sys.stdout, "PairModel %s (abundance: %.2f%%)" % (allele_pair, prob * 100.0)
                    if simulation:
                        if allele1 in test_HLA_names and allele2 in test_HLA_names:
                            rank_i = prob_i
                            while rank_i > 0:
                                if HLA_prob[rank_i-1][1] == prob:                                        
                                    rank_i -= 1
                                else:
                                    break
                            print >> sys.stderr, "\t\t\t*** %d ranked %s (abundance: %.2f%%)" % (rank_i + 1, allele_pair, prob * 100.0)
                            if rank_i == 0:
                                success[0] = True
                            break
                    print >> sys.stderr, "\t\t\t\t%d ranked %s (abundance: %.2f%%)" % (prob_i + 1, allele_pair, prob * 100.0)
                    if not simulation and prob_i >= 9:
                        break
                print >> sys.stderr

                # Li's method
                """
                li_hla = os.path.join(ex_path, "li_hla/hla")
                if os.path.exists(li_hla):
                    li_hla_cmd = [li_hla,
                                  "hla",
                                  "hla_input.bam",
                                  "-b", "%s*BACKBONE" % gene]
                    li_hla_proc = subprocess.Popen(li_hla_cmd,
                                                   stdout=subprocess.PIPE,
                                                   stderr=open("/dev/null", 'w'))

                    # read in the result of Li's hla
                    for line in li_hla_proc.stdout:
                        allele1, allele2, score = line.strip().split()
                        score = float(score)
                        if simulation:
                            if allele1 in test_HLA_names and allele2 in test_HLA_names:
                                print >> sys.stderr, "\t\t\t*** 1 ranked %s-%s (score: %.2f)" % (allele1, allele2, score)
                                success[0] = True
                            else:
                                print >> sys.stderr, "\t\t\tLiModel fails"
                        if best_alleles:
                            print >> sys.stdout, "LiModel %s-%s (score: %.2f)" % (allele1, allele2, score)
                    li_hla_proc.communicate()
                """

            if simulation and not False in success:
                aligner_type = "%s %s" % (aligner, index_type)
                if not aligner_type in test_passed:
                    test_passed[aligner_type] = 1
                else:
                    test_passed[aligner_type] += 1

        if remove_alignment_file and not simulation:
            os.system("rm %s*" % (alignment_fname))

    if simulation:
        return test_passed

    
"""
"""
def read_HLA_alleles(fname, HLAs):
    for line in open(fname):
        if line.startswith(">"):
            HLA_name = line.strip().split()[0][1:]
            HLA_gene = HLA_name.split('*')[0]
            if not HLA_gene in HLAs:
                HLAs[HLA_gene] = {}
            if not HLA_name in HLAs[HLA_gene]:
                HLAs[HLA_gene][HLA_name] = ""
        else:
            HLAs[HLA_gene][HLA_name] += line.strip()
    return HLAs


"""
"""
def read_HLA_vars(fname, reference_type):
    Vars, Var_list = {}, {}
    for line in open(fname):
        var_id, var_type, allele, pos, data = line.strip().split('\t')
        pos = int(pos)
        if reference_type != "gene":
            allele, dist = None, 0
            for tmp_gene, values in refHLA_loci.items():
                allele_name, chr, left, right, exons = values
                if allele == None:
                    allele = allele_name
                    dist = abs(pos - left)
                else:
                    if dist > abs(pos - left):
                        allele = allele_name
                        dist = abs(pos - left)
            
        gene = allele.split('*')[0]
        if not gene in Vars:
            Vars[gene] = {}
            assert not gene in Var_list
            Var_list[gene] = []
            
        assert not var_id in Vars[gene]
        left = 0
        if reference_type != "gene":
            _, _, left, _, _ = refHLA_loci[gene]
        Vars[gene][var_id] = [var_type, pos - left, data]
        Var_list[gene].append([pos - left, var_id])
        
    for gene, in_var_list in Var_list.items():
        Var_list[gene] = sorted(in_var_list)

    return Vars, Var_list


"""
"""
def read_HLA_links(fname):
    Links = {}
    for line in open(fname):
        var_id, alleles = line.strip().split('\t')
        alleles = alleles.split()
        assert not var_id in Links
        Links[var_id] = alleles

    return Links


"""
"""
def construct_allele_seq(backbone_seq, var_ids, Vars):
    allele_seq = list(backbone_seq)
    for id in var_ids:
        assert id in Vars
        type, pos, data = Vars[id]
        assert pos < len(allele_seq)
        if type == "single":
            assert allele_seq[pos] != data
            allele_seq[pos] = data
        else:
            assert type == "deletion"
            del_len = int(data)
            assert pos + del_len <= len(allele_seq)
            for i in range(pos, pos + del_len):
                allele_seq[i] = '.'

    allele_seq = ''.join(allele_seq)
    allele_seq = allele_seq.replace('.', '')
    return allele_seq


"""
"""
def test_HLA_genotyping(base_fname,
                        reference_type,
                        hla_list,
                        partial,
                        aligners,
                        read_fname,
                        alignment_fname,
                        threads,
                        simulate_interval,
                        read_len,
                        fragment_len,
                        best_alleles,
                        num_editdist,
                        perbase_errorrate,
                        perbase_snprate,
                        skip_fragment_regions,
                        assembly,
                        error_correction,
                        discordant,
                        display_alleles,
                        verbose,
                        daehwan_debug):
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(test_HLA_genotyping))
    ex_path = os.path.dirname(curr_script)

    # Clone a git repository, IMGTHLA
    if not os.path.exists("IMGTHLA"):
        HLA_typing.clone_IMGTHLA_database()

    simulation = (read_fname == [] and alignment_fname == "")

    def check_files(fnames):
        for fname in fnames:
            if not os.path.exists(fname):
                return False
        return True

    # Download human genome and HISAT2 index
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]
    if not check_files(HISAT2_fnames):
        typing_common.download_genome_and_index(ex_path)

    # Check if the pre-existing files (hla*) are compatible with the current parameter setting
    if os.path.exists("hla.ref"):
        left = 0
        HLA_genes = []
        BACKBONE = False
        for line in open("hla.ref"):
            HLA_name = line.strip().split()[0]
            if HLA_name.find("BACKBONE") != -1:
                BACKBONE = True
            HLA_gene = HLA_name.split('*')[0]
            HLA_genes.append(HLA_gene)
        delete_hla_files = False
        if reference_type == "gene":
            if not BACKBONE:
                delete_hla_files = True
        elif reference_type in ["chromosome", "genome"]:
            if BACKBONE:
                delete_hla_files = True
        else:
            assert False
        if len(hla_list) == 0:
            hla_list = HLA_genes
        if not set(hla_list).issubset(set(HLA_genes)):
            delete_hla_files = True
        if delete_hla_files:
            os.system("rm hla*")

    # Extract HLA variants, backbone sequence, and other sequeces  
    if len(base_fname) > 0:
        base_fname = "_" + base_fname
    base_fname = "hla" + base_fname
    
    HLA_fnames = [base_fname + "_backbone.fa",
                  base_fname + "_sequences.fa",
                  base_fname + ".ref",
                  base_fname + ".snp",
                  base_fname + ".index.snp",
                  base_fname + ".haplotype",
                  base_fname + ".link"]
    
    if verbose >= 1:
        print >> sys.stderr, HLA_fnames
    
    if not check_files(HLA_fnames):
        extract_hla_script = os.path.join(ex_path, "hisatgenotype_extract_vars.py")
        extract_cmd = [extract_hla_script,
                       "--reference-type", reference_type]
        if len(hla_list) > 0:
            extract_cmd += ["--hla-list", ','.join(hla_list)]

        if len(base_fname) > 3:
            extract_cmd += ["--base", base_fname]

        if not partial:
            extract_cmd += ["--no-partial"]
        extract_cmd += ["--inter-gap", "30",
                        "--intra-gap", "50"]

        # DK - debugging purposes
        extract_cmd += ["--min-var-freq", "0.1"]
        
        # extract_cmd += ["--leftshift"]
        # extract_cmd += ["--ext-seq", "300"]
        if verbose >= 1:
            print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
        proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()
        
        if not check_files(HLA_fnames):
            print >> sys.stderr, "Error: hisatgenotype_extract_vars failed!"
            sys.exit(1)

    for aligner, index_type in aligners:
        if aligner == "hisat2":
            # Build HISAT2 graph indexes based on the above information
            if index_type == "graph":
                HLA_hisat2_graph_index_fnames = ["%s.graph.%d.ht2" % (base_fname, i+1) for i in range(8)]
                if not check_files(HLA_hisat2_graph_index_fnames):
                    hisat2_build = os.path.join(ex_path, "hisat2-build")
                    build_cmd = [hisat2_build,
                                 "-p", str(threads),
                                 "--snp", "%s.index.snp" % base_fname,
                                 "--haplotype", "%s.haplotype" % base_fname,
                                 "%s_backbone.fa" % base_fname,
                                 "%s.graph" % base_fname]
                    if verbose >= 1:
                        print >> sys.stderr, "\tRunning:", ' '.join(build_cmd)
                    proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
                    proc.communicate()        
                    if not check_files(HLA_hisat2_graph_index_fnames):
                        print >> sys.stderr, "Error: indexing HLA failed!  Perhaps, you may have forgotten to build hisat2 executables?"
                        sys.exit(1)
            # Build HISAT2 linear indexes based on the above information
            else:
                assert index_type == "linear"
                HLA_hisat2_linear_index_fnames = ["%s.linear.%d.ht2" % (base_fname, i+1) for i in range(8)]
                if reference_type == "gene" and not check_files(HLA_hisat2_linear_index_fnames):
                    hisat2_build = os.path.join(ex_path, "hisat2-build")
                    build_cmd = [hisat2_build,
                                 "%s_backbone.fa,%s_sequences.fa" % (base_fname, base_fname),
                                 "%s.linear" % base_fname]
                    proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
                    proc.communicate()        
                    if not check_files(HLA_hisat2_linear_index_fnames):
                        print >> sys.stderr, "Error: indexing HLA failed!"
                        sys.exit(1)
        else:
            assert aligner == "bowtie2" and index_type == "linear"
            # Build Bowtie2 indexes based on the above information
            HLA_bowtie2_index_fnames = ["%s.%d.bt2" % (base_fname, i+1) for i in range(4)]
            HLA_bowtie2_index_fnames += ["%s.rev.%d.bt2" % (base_fname, i+1) for i in range(2)]
            if reference_type == "gene" and not check_files(HLA_bowtie2_index_fnames):
                build_cmd = ["bowtie2-build",
                             "%s_backbone.fa,%s_sequences.fa" % (base_fname, base_fname),
                             base_fname]
                proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'))
                proc.communicate()        
                if not check_files(HLA_bowtie2_index_fnames):
                    print >> sys.stderr, "Error: indexing HLA failed!"
                    sys.exit(1)

    # Read partial alleles from hla.data (temporary)
    partial_alleles = set()
    for line in open("IMGTHLA/hla.dat"):
        if not line.startswith("DE"):
            continue
        allele_name = line.split()[1][:-1]
        if allele_name.startswith("HLA-"):
            allele_name = allele_name[4:]
        gene = allele_name.split('*')[0]
        if line.find("partial") != -1:
            partial_alleles.add(allele_name)

    # Read HLA alleles (names and sequences)
    refHLAs, refHLA_loci = {}, {}
    for line in open("hla.ref"):
        HLA_name, chr, left, right, length, exon_str = line.strip().split()
        HLA_gene = HLA_name.split('*')[0]
        assert not HLA_gene in refHLAs
        refHLAs[HLA_gene] = HLA_name
        left, right = int(left), int(right)
        exons = []
        for exon in exon_str.split(','):
            exon_left, exon_right = exon.split('-')
            exons.append([int(exon_left), int(exon_right)])
        refHLA_loci[HLA_gene] = [HLA_name, chr, left, right, exons]
    HLAs = {}
    if len(hla_list) == 0:
        hla_list = refHLA_loci.keys()

    if reference_type == "gene":
        read_HLA_alleles(base_fname + "_backbone.fa", HLAs)
    read_HLA_alleles(base_fname + "_sequences.fa", HLAs)

    # HLA gene alleles
    HLA_names = {}
    for HLA_gene, data in HLAs.items():
        HLA_names[HLA_gene] = list(data.keys())

    # HLA gene allele lengths
    HLA_lengths = {}
    for HLA_gene, HLA_alleles in HLAs.items():
        HLA_lengths[HLA_gene] = {}
        for allele_name, seq in HLA_alleles.items():
            HLA_lengths[HLA_gene][allele_name] = len(seq)

    # Read HLA variants, and link information
    Vars, Var_list = read_HLA_vars("%s.snp" % base_fname, reference_type)
    Links = read_HLA_links("%s.link" % base_fname)
    # Test HLA typing
    test_list = []
    if simulation:
        basic_test, pair_test = True, False
        if daehwan_debug:
            if "basic_test" in daehwan_debug:
                basic_test, pair_test = True, False
            else:
                basic_test, pair_test = False, True

        test_passed = {}
        test_list = []
        genes = list(set(hla_list) & set(HLA_names.keys()))
        if basic_test:
            for gene in genes:
                HLA_gene_alleles = HLA_names[gene]
                for HLA_name in HLA_gene_alleles:
                    if HLA_name.find("BACKBONE") != -1:
                        continue
                    test_list.append([[HLA_name]])
        if pair_test:
            test_size = 500
            allele_count = 2
            for test_i in range(test_size):
                test_pairs = []
                for gene in genes:
                    HLA_gene_alleles = []

                    for allele in HLA_names[gene]:
                        if allele.find("BACKBONE") != -1:
                            continue

                        if "full" in debug:
                            if allele in partial_alleles:
                                continue

                        HLA_gene_alleles.append(allele)
                    nums = [i for i in range(len(HLA_gene_alleles))]
                    random.shuffle(nums)
                    test_pairs.append(sorted([HLA_gene_alleles[nums[i]] for i in range(allele_count)]))
                test_list.append(test_pairs)

        # DK - for debugging purposes
        # test_list = [[["A*01:01:01:01"]], [["A*32:29"]]]
        # test_list = [[["A*01:01:01:01", "A*03:01:01:01"]]]
        # test_list = [[["A*24:36N", "A*30:03"]]]
        # test_list = [[["A*26:25N"]]] # for allele that includes an insertion
        # test_list = [[["A*24:36N"]]] # for normal allele?
        # test_list = [[["A*02:01:21"]], [["A*03:01:01:01"]], [["A*03:01:01:04"]], [["A*02:521"]]]
        for test_i in range(len(test_list)):
            if "test_id" in daehwan_debug:
                daehwan_test_ids = daehwan_debug["test_id"].split('-')
                if str(test_i + 1) not in daehwan_test_ids:
                    continue

            print >> sys.stderr, "Test %d" % (test_i + 1), str(datetime.now())
            test_HLA_list = test_list[test_i]
            num_frag_list = typing_common.simulate_reads(HLAs,
                                                         test_HLA_list,
                                                         Vars,
                                                         Links,
                                                         simulate_interval,
                                                         read_len,
                                                         fragment_len,
                                                         perbase_errorrate,
                                                         perbase_snprate,
                                                         skip_fragment_regions)

            assert len(num_frag_list) == len(test_HLA_list)
            for i_ in range(len(test_HLA_list)):
                test_HLA_names = test_HLA_list[i_]
                num_frag_list_i = num_frag_list[i_]
                assert len(num_frag_list_i) == len(test_HLA_names)
                for j_ in range(len(test_HLA_names)):
                    test_HLA_name = test_HLA_names[j_]
                    gene = test_HLA_name.split('*')[0]
                    test_HLA_seq = HLAs[gene][test_HLA_name]
                    seq_type = "partial" if test_HLA_name in partial_alleles else "full"
                    print >> sys.stderr, "\t%s - %d bp (%s sequence, %d pairs)" % (test_HLA_name, len(test_HLA_seq), seq_type, num_frag_list_i[j_])

            if "single-end" in daehwan_debug:
                read_fname = ["hla_input_1.fa"]
            else:
                read_fname = ["hla_input_1.fa", "hla_input_2.fa"]

            fastq = False
            tmp_test_passed = typing(ex_path,
                                     simulation,
                                     reference_type,
                                     test_HLA_list,
                                     partial,
                                     partial_alleles,
                                     refHLAs,
                                     HLAs,                       
                                     HLA_names,
                                     HLA_lengths,
                                     refHLA_loci,
                                     Vars,
                                     Var_list,
                                     Links,
                                     aligners,
                                     num_editdist,
                                     assembly,
                                     error_correction,
                                     discordant,
                                     display_alleles,
                                     fastq,
                                     read_fname,
                                     alignment_fname,
                                     num_frag_list,
                                     threads,
                                     best_alleles,
                                     verbose)

            for aligner_type, passed in tmp_test_passed.items():
                if aligner_type in test_passed:
                    test_passed[aligner_type] += passed
                else:
                    test_passed[aligner_type] = passed

                print >> sys.stderr, "\t\tPassed so far: %d/%d (%.2f%%)" % (test_passed[aligner_type], test_i + 1, (test_passed[aligner_type] * 100.0 / (test_i + 1)))


        for aligner_type, passed in test_passed.items():
            print >> sys.stderr, "%s:\t%d/%d passed (%.2f%%)" % (aligner_type, passed, len(test_list), passed * 100.0 / len(test_list))
    
    else: # With real reads or BAMs
        print >> sys.stderr, "\t", ' '.join(hla_list)
        fastq = True
        typing(ex_path,
               simulation,
               reference_type,
               hla_list,
               partial,
               partial_alleles,
               refHLAs,
               HLAs,                       
               HLA_names,
               HLA_lengths,
               refHLA_loci,
               Vars,
               Var_list,
               Links,
               aligners,
               num_editdist,
               assembly,
               error_correction,
               discordant,
               display_alleles,
               fastq,
               read_fname,
               alignment_fname,
               [],
               threads,
               best_alleles,
               verbose)

        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='test HLA genotyping')
    parser.add_argument("--base",
                        dest="base_fname",
                        type=str,
                        default="",
                        help="base filename for backbone HLA sequence, HLA variants, and HLA linking info")
    parser.add_argument("--reference-type",
                        dest="reference_type",
                        type=str,
                        default="gene",
                        help="Reference type: gene, chromosome, and genome (default: gene)")
    parser.add_argument("--hla-list",
                        dest="hla_list",
                        type=str,
                        default="",
                        help="A comma-separated list of HLA genes (default: empty, all HLA genes in IMGT/HLA database)")
    parser.add_argument('--no-partial',
                        dest='partial',
                        action='store_false',
                        help='Include partial alleles (e.g. A_nuc.fasta)')
    parser.add_argument("--aligner-list",
                        dest="aligners",
                        type=str,
                        default="hisat2.graph,hisat2.linear,bowtie2.linear",
                        help="A comma-separated list of aligners (default: hisat2.graph,hisat2.linear,bowtie2.linear)")
    parser.add_argument("--reads",
                        dest="read_fname",
                        type=str,
                        default="",
                        help="Fastq read file name")
    parser.add_argument("--alignment",
                        dest="alignment_fname",
                        type=str,
                        default="",
                        help="BAM file name")
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    parser.add_argument("--simulate-interval",
                        dest="simulate_interval",
                        type=int,
                        default=10,
                        help="Reads simulated at every these base pairs (default: 10)")
    parser.add_argument("--read-len",
                        dest="read_len",
                        type=int,
                        default=100,
                        help="Length of simulated reads (default: 100)")
    parser.add_argument("--fragment-len",
                        dest="fragment_len",
                        type=int,
                        default=250,
                        help="Length of fragments (default: 250)")
    parser.add_argument("--best-alleles",
                        dest="best_alleles",
                        action='store_true',
                        help="")
    parser.add_argument("--random-seed",
                        dest="random_seed",
                        type=int,
                        default=1,
                        help="A seeding number for randomness (default: 1)")
    parser.add_argument("--num-editdist", "--num-mismatch",
                        dest="num_editdist",
                        type=int,
                        default=0,
                        help="Maximum number of mismatches per read alignment to be considered (default: 0)")
    parser.add_argument("--perbase-errorrate",
                        dest="perbase_errorrate",
                        type=float,
                        default=0.0,
                        help="Per basepair error rate in percentage when simulating reads (default: 0.0)")
    parser.add_argument("--perbase-snprate",
                        dest="perbase_snprate",
                        type=float,
                        default=0.0,
                        help="Per basepair SNP rate in percentage when simulating reads (default: 0.0)")
    parser.add_argument("--skip-fragment-regions",
                        dest="skip_fragment_regions",
                        type=str,
                        default="",
                        help="A comma-separated list of regions from which no reads originate, e.g., 500-600,1200-1400 (default: None).")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')
    parser.add_argument('--verbose-level',
                        dest='verbose_level',
                        type=int,
                        default=0,
                        help='also print some statistics to stderr (default: 0)')
    parser.add_argument("--debug",
                        dest="debug",
                        type=str,
                        default="",
                        help="e.g., test_id:10,read_id:10000,basic_test")
    parser.add_argument("--no-assembly",
                        dest="assembly",
                        action="store_false",
                        help="Perform assembly")
    parser.add_argument("--no-error-correction",
                        dest="error_correction",
                        action="store_false",
                        help="Correct sequencing errors")
    parser.add_argument("--discordant",
                        dest="discordant",
                        action="store_true",
                        help="Allow discordantly mapped pairs or singletons")    
    parser.add_argument("--display-alleles",
                        dest="display_alleles",
                        type=str,
                        default="",
                        help="A comma-separated list of alleles to display in HTML (default: empty)")

    args = parser.parse_args()
    if not args.reference_type in ["gene", "chromosome", "genome"]:
        print >> sys.stderr, "Error: --reference-type (%s) must be one of gene, chromosome, and genome." % (args.reference_type)
        sys.exit(1)
    if args.hla_list == "":
        hla_list = []
    else:
        hla_list = args.hla_list.split(',')
    if args.aligners == "":
        print >> sys.stderr, "Error: --aligners must be non-empty."
        sys.exit(1)    
    args.aligners = args.aligners.split(',')
    for i in range(len(args.aligners)):
        args.aligners[i] = args.aligners[i].split('.')
    if args.read_fname:
        args.read_fname = args.read_fname.split(',')
    else:
        args.read_fname = []
    if args.alignment_fname != "" and \
            not os.path.exists(args.alignment_fname):
        print >> sys.stderr, "Error: %s doesn't exist." % args.alignment_fname
        sys.exit(1)

    if args.verbose and args.verbose_level == 0:
        args.verbose_level = 1
        
    debug = {}
    if args.debug != "":
        for item in args.debug.split(','):
            if ':' in item:
                key, value = item.split(':')
                debug[key] = value
            else:
                debug[item] = 1

    if not args.partial:
        print >> sys.stderr, "Warning: --no-partial will be no longer supported!"

    if args.read_len * 2 > args.fragment_len:
        print >> sys.stderr, "Warning: fragment might be too short (%d)" % (args.fragment_len)

    skip_fragment_regions = []
    if args.skip_fragment_regions != "":
        prev_left, prev_right = -1, -1
        for region in args.skip_fragment_regions.split(','):
            left, right = region.split('-')
            left, right = int(left), int(right)
            assert left < right
            assert prev_right < left
            prev_left, prev_right = left, right
            skip_fragment_regions.append([left, right])

    if args.display_alleles == "":
        display_alleles = []
    else:
        display_alleles = args.display_alleles.split(',')

    random.seed(args.random_seed)
    test_HLA_genotyping(args.base_fname,
                        args.reference_type,
                        hla_list,
                        args.partial,
                        args.aligners,
                        args.read_fname,
                        args.alignment_fname,
                        args.threads,
                        args.simulate_interval,
                        args.read_len,
                        args.fragment_len,
                        args.best_alleles,
                        args.num_editdist,
                        args.perbase_errorrate,
                        args.perbase_snprate,
                        skip_fragment_regions,
                        args.assembly,
                        args.error_correction,
                        args.discordant,
                        display_alleles,
                        args.verbose_level,
                        debug)

