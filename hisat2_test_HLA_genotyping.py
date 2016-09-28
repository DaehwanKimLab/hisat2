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
from argparse import ArgumentParser, FileType


"""
"""
def simulate_reads(HLAs,
                   test_HLA_list,
                   simulate_interval):
    HLA_reads_1, HLA_reads_2 = [], []
    for test_HLA_names in test_HLA_list:
        gene = test_HLA_names[0].split('*')[0]

        # Simulate reads from two HLA alleles
        def simulate_reads_impl(seq, simulate_interval = 1, frag_len = 250, read_len = 100):
            comp_table = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
            reads_1, reads_2 = [], []
            for i in range(0, len(seq) - frag_len + 1, simulate_interval):
                reads_1.append(seq[i:i+read_len])
                tmp_read_2 = reversed(seq[i+frag_len-read_len:i+frag_len])
                read_2 = ""
                for s in tmp_read_2:
                    if s in comp_table:
                        read_2 += comp_table[s]
                    else:
                        read_2 += s
                reads_2.append(read_2)
            return reads_1, reads_2

        for test_HLA_name in test_HLA_names:
            HLA_seq = HLAs[gene][test_HLA_name]
            tmp_reads_1, tmp_reads_2 = simulate_reads_impl(HLA_seq, simulate_interval)
            HLA_reads_1 += tmp_reads_1
            HLA_reads_2 += tmp_reads_2

    # Write reads into a fasta read file
    def write_reads(reads, idx):
        read_file = open('hla_input_%d.fa' % idx, 'w')
        for read_i in range(len(reads)):
            print >> read_file, ">%d" % (read_i + 1)
            print >> read_file, reads[read_i]
        read_file.close()
    write_reads(HLA_reads_1, 1)
    write_reads(HLA_reads_2, 2)


"""
Align reads, and sort the alignments into a BAM file
"""
def align_reads(ex_path,
                aligner,
                index_type,
                read_fname,
                fastq,
                threads,
                out_fname,
                verbose):
    if aligner == "hisat2":
        hisat2 = os.path.join(ex_path, "hisat2")
        aligner_cmd = [hisat2,
                       "--no-unal",
                       "--mm"]
        if index_type == "linear":
            aligner_cmd += ["-k", "10"]
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

    # normalize2(HLA_prob, HLA_length)
    normalize(HLA_prob)
    def next_prob(HLA_cmpt, HLA_prob, HLA_length):
        HLA_prob_next = {}
        for cmpt, count in HLA_cmpt.items():
            alleles = cmpt.split('-')
            alleles_prob = 0.0
            for allele in alleles:
                assert allele in HLA_prob
                alleles_prob += HLA_prob[allele]
            for allele in alleles:
                if allele not in HLA_prob_next:
                    HLA_prob_next[allele] = 0.0
                HLA_prob_next[allele] += (float(count) * HLA_prob[allele] / alleles_prob)
        # normalize2(HLA_prob_next, HLA_length)
        normalize(HLA_prob_next)
        return HLA_prob_next

    diff, iter = 1.0, 0
    while diff > 0.0001 and iter < 1000:
        HLA_prob_next = next_prob(HLA_cmpt, HLA_prob, HLA_length)
        diff = prob_diff(HLA_prob, HLA_prob_next)
        HLA_prob = HLA_prob_next

        # DK - for debugging purposes
        if iter % 10 == 0 and False:
            print "iter", iter
            for allele, prob in HLA_prob.items():
                if prob >= 0.01:
                    print "\t", allele, prob
        
        iter += 1
    for allele, prob in HLA_prob.items():
        allele_len = HLA_length[allele]
        HLA_prob[allele] /= float(allele_len)
    normalize(HLA_prob)
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
def assemble_two_haplotypes(a_vars, b_vars, debug = False):
    const_deletion_relax = 10
    def get_var_info(var):
        var_type, var_left, var_data = var.split('-')[:3]
        var_left = int(var_left)
        if var_type == "unknown":
            var_right = int(var_data)
        elif var_type == "deletion":
            var_right = var_left + int(var_data) - 1
        else:
            var_right = var_left
        assert var_left <= var_right
        if var_type in ["left", "right"]:
            assert var_left == var_right
        return var_type, var_left, var_right

    a_vars, b_vars = a_vars[:], b_vars[:]
    c_vars = []
    a_i, b_i = 0, 0
    while a_i < len(a_vars) or b_i < len(b_vars):
        if a_i == len(a_vars):
            c_vars += b_vars[b_i:]
            break
        elif b_i == len(b_vars):
            c_vars += a_vars[a_i:]
            break
        else:
            a_var = a_vars[a_i]
            b_var = b_vars[b_i]
            if a_var == b_var:
                a_i += 1
                b_i += 1
                c_vars.append(a_var)
                continue

            a_type, a_left, a_right = get_var_info(a_var)
            b_type, b_left, b_right = get_var_info(b_var)

            if debug:
                print >> sys.stderr, "a_var:", a_i, a_var
                print >> sys.stderr, "b_var:", b_i, b_var
                print >> sys.stderr, "c_vars:", c_vars

            if a_type == "left":
                assert a_i == 0
            elif a_type == "right":
                assert a_i + 1 == len(a_vars)
            if b_type == "left":
                assert b_i == 0
            elif b_type == "right":
                assert b_i + 1 == len(b_vars)

            if a_type in ["left", "right"]:
                if a_type == "left":
                    assert a_left <= b_left
                    c_vars.append("left-%d-%d-left" % (a_left, a_left))
                    a_i += 1
                    while a_i < len(a_vars):
                        a_var = a_vars[a_i]
                        a_type, a_left, a_right = get_var_info(a_var)
                        if a_right < b_left:
                            if a_i + 1 == len(a_vars):
                                if b_type in ["left", "right"]:
                                    b_i += 1
                                    break
                                break
                            else:
                                c_vars.append(a_var)
                        else:
                            if b_type == "left":
                                if a_type == "unknown":
                                    if a_left < b_left and a_right > b_left:
                                        a_vars[a_i] = "unknown-%d-%d-unknown" % (b_left, a_right)
                                        assert b_left <= a_right
                                elif a_type == "deletion":
                                    relax = const_deletion_relax if b_type == "deletion" else 0
                                    if a_right < b_left + relax:
                                        c_vars.append(a_var)
                                        a_i += 1
                                b_i += 1
                            break
                        a_i += 1
                else:
                    assert a_type == "right"
                    a_i += 1
                    while b_i < len(b_vars):
                        b_var = b_vars[b_i]
                        b_type, b_left, b_right = get_var_info(b_var)
                        if b_type in ["left", "right"]:
                            if b_type == "left":
                                if a_right < b_left:
                                    return c_vars, False
                            else:
                                c_vars.append("right-%d-%d-right" % (b_right, b_right))
                        elif b_type == "unknown":
                            if b_left < a_right and b_right > a_right:
                                c_vars.append("unknown-%d-%d-unknown" % (a_right + 1, b_right))
                        else:
                            relax = const_deletion_relax if b_type == "deletion" else 0
                            if b_left + relax < a_right:
                                return c_vars, False
                            else:
                                break
                        b_i += 1
            elif b_type in ["left", "right"]:
                if b_type == "left":
                    if a_left <= b_left:
                        a_i += 1
                        c_vars.append(a_var)
                    else:
                        b_i += 1
                else:
                    assert b_type == "right"
                    b_i += 1
                    while a_i < len(a_vars):
                        a_var = a_vars[a_i]
                        a_type, a_left, a_right = get_var_info(a_var)
                        relax = 1
                        if a_left > b_right + relax:
                            break
                        else:
                            if a_type == "unknown":
                                if a_right > b_right:
                                    c_vars.append("unknown-%d-%d-unknown" % (b_right + 1, a_right))
                                    a_i += 1
                            elif a_type == "deletion":
                                relax = const_deletion_relax if b_type == "deletion" else 0
                                if a_right < b_right + relax:
                                    c_vars.append(a_var)
                                    a_i += 1
                            else:
                                return c_vars, False
                            break
                        a_i += 1

            else:
                assert a_type not in ["left", "right"] and b_type not in ["left", "right"]
                if a_type == "unknown" and b_type == "unknown":
                    if (a_left <= b_left and a_right >= b_left) or \
                            (b_left <= a_left and b_right >= a_left):
                        c_vars.append("unknown-%d-%d-unknown" % (max(a_left, b_left), min(a_right, b_right)))
                    a_i += 1
                    b_i += 1
                    if a_right < b_right:
                        while a_i < len(a_vars):
                            a_var = a_vars[a_i]
                            a_type, a_left, a_right = get_var_info(a_var)
                            if (a_type == "deletion" and a_right < b_right + 1) or \
                                    (a_type != "deletion" and a_right < b_right):
                                a_i += 1
                                if a_type == "right":
                                    assert a_i == len(a_vars)
                                    c_vars.append("unknown-%d-%d-unknown" % (a_right + 1, b_right))
                                else:
                                    c_vars.append(a_var)
                            else:
                                break                                            
                    else:
                        while b_i < len(b_vars):
                            b_var = b_vars[b_i]
                            b_type, b_left, b_right = get_var_info(b_var)
                            if (b_type == "deletion" and b_right < a_right + 1) or \
                                    (b_type != "deletion" and b_right < a_right):
                                b_i += 1
                                if b_type == "right":
                                    assert b_i == len(b_vars)
                                    c_vars.append("unknown-%d-%d-unknown" % (b_right + 1, a_right))
                                else:
                                    c_vars.append(b_var)
                            else:
                                break
                elif a_type == "unknown" or b_type == "unknown":
                    if a_right < b_right:
                        if b_type == "unknown":
                            # DK - for debugging purposes
                            """
                            if a_type == "deletion":
                                if a_right + const_deletion_relax < b_left:
                                    return c_vars, False
                            else:
                                if a_right < b_left:
                                    return c_vars, False
                            """
                            if a_right < b_left:
                                return c_vars, False
                            
                            a_i += 1
                            b_vars[b_i] = "unknown-%d-%d-unknown" % (a_right + 1, b_right)
                            c_vars.append(a_var)
                        else:
                            assert a_type == "unknown"
                            a_i += 1
                    else:
                        if a_type == "unknown":
                            if b_i == 0 or b_i + 1 == len(b_vars):
                                return c_vars, False                                            
                            b_i += 1
                            if b_right  < a_right:
                                a_vars[a_i] = "unknown-%d-%d-unknown" % (b_right + 1, a_right)
                            else:
                                a_i += 1
                            c_vars.append(b_var)                                            
                        else:
                            assert b_type == "unknown"
                            b_i += 1
                            if b_right > a_right:
                                c_vars.append("unknown-%d-%d-unknown" % (max(a_right + 1, b_left), b_right))
                else:
                    if a_left == b_left:
                        return c_vars, False
                    if a_left < b_left and b_i == 0:
                        c_vars.append(a_var)
                        a_i += 1
                    else:
                        return c_vars, False

    return c_vars, True


"""
"""
def assemble_alleles(N_haplotype_list, max_allele = 2):
    N_alleles = set()
    alleles = [[N_haplotype_list[0], 1]]
    while len(alleles) > 0:
        a, i = alleles.pop()
        if i == len(N_haplotype_list):
            N_alleles.add(a)
            if len(N_alleles) >= max_allele:
                # print >> sys.stderr, "Warning: too many alleles have been assembled (>= %d)" % len(N_alleles)
                break
        else:
            b = N_haplotype_list[i]
            a_vars = a.split(',')
            assert len(a_vars) > 0
            b_vars = b.split(',')
            assert len(b_vars) > 0

            c_vars, success = assemble_two_haplotypes(a_vars, b_vars)
            c = ','.join(c_vars)

            # DK - for debugging purposes
            if c.find("single-1154-C-hv880,single-1158-T-hv889") != -1 and False:
                print "a:", a
                print "b:", b
                print "c:", ','.join(c_vars), success
                sys.exit(1)
            # if len(c_vars) >= 10:
            #    sys.exit(1)

            # Sanity check
            for c_var in c_vars[1:-1]:
                if c_var.find("left") != -1 or \
                        c_var.find("right") != -1:
                    print >> sys.stderr, "Error containing left or right in intermediate elements!"
                    print >> sys.stderr, "a:", a
                    print >> sys.stderr, "b:", b
                    print >> sys.stderr, "c:", c
                    sys.exit(1)

            if success:
                alleles.append([c, i + 1])
            else:
                # DK - temporary
                if len(alleles) == 0:
                    # DK - for debugging purposes
                    """
                    print "a:", a
                    print "b:", b
                    print "c:", c, success
                    sys.exit(1)
                    """

                    alleles.append([b, i + 1])
                alleles.append([a, i + 1])                                        
    return N_alleles


"""
Example:
   a: left-2925-2925-left,unknown-3015-3161-unknown,single-3166-A-hv1696,single-3240-G-hv1705,right-3251-3251-right
   b: left-1572-1572-left,deletion-1625-1-hv1345,unknown-1663-1754-unknown,right-1844-1844-right
"""
def haplotype_cmp(a, b):
    a_vars = a.split(',')
    assert len(a_vars) > 0
    a_first_var, a_last_var = a_vars[0], a_vars[-1]
    a_left = int(a_first_var.split('-')[1])
    a_right = int(a_last_var.split('-')[1])

    b_vars = b.split(',')
    assert len(b_vars) > 0
    b_first_var, b_last_var = b_vars[0], b_vars[-1]

    b_left = int(b_first_var.split('-')[1])
    b_right = int(b_last_var.split('-')[1])

    if a_left != b_left:
        return a_left - b_left
    else:
        return a_right - b_right                    


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
"""
def HLA_typing(ex_path,
               simulation,
               reference_type,
               hla_list,
               partial,
               refHLAs,
               HLAs,
               HLA_names,
               HLA_lengths,
               refHLA_loci,
               Vars,
               Var_list,
               Links,
               HLAs_default,
               Vars_default,
               Var_list_default,
               Links_default,
               exclude_allele_list,
               aligners,
               num_mismatch,
               assembly,
               concordant_assembly,
               exonic_only,
               fastq,
               read_fname,
               alignment_fname,
               threads,
               enable_coverage,
               best_alleles,
               verbose):
    
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

            allele_reps = {} # allele representatives
            if partial and exonic_only:
                _allele_groups = {}
                _alleles = HLAs[gene].keys()
                for _allele in _alleles:
                    _fields = _allele.split(':')
                    if len(_fields) <= 3:
                        assert _allele not in _allele_groups
                        _allele_groups[_allele] = [_allele]
                    else:
                        assert len(_fields) == 4
                        _allele_rep = ':'.join(_fields[:-1])
                        if _allele_rep not in _allele_groups:
                            _allele_groups[_allele_rep] = []
                        assert _allele not in _allele_groups[_allele_rep]
                        _allele_groups[_allele_rep].append(_allele)
                    
                for _allele_members in _allele_groups.values():
                    assert len(_allele_members) > 0
                    _allele_rep = _allele_members[0]
                    for _allele_member in _allele_members:
                        assert _allele_member not in allele_reps
                        allele_reps[_allele_member] = _allele_rep

            if not os.path.exists(alignment_fname + ".bai"):
                os.system("samtools index %s" % alignment_fname)
            # Read alignments
            alignview_cmd = ["samtools",
                             "view",
                             alignment_fname]
            base_locus = 0
            if index_type == "graph":
                if reference_type == "gene":
                    alignview_cmd += ["%s" % ref_allele]
                else:
                    assert reference_type in ["chromosome", "genome"]
                    _, chr, left, right, _ = refHLA_loci[gene]
                    base_locus = left
                    alignview_cmd += ["%s:%d-%d" % (chr, left + 1, right + 1)]

                bamview_proc = subprocess.Popen(alignview_cmd,
                                                stdout=subprocess.PIPE,
                                                stderr=open("/dev/null", 'w'))

                sort_read_cmd = ["sort", "-k", "1,1", "-s"] # -s for stable sorting
                # sort_read_cmd = ["sort", "-k", "1,1", "-k", "4,4", "-n"]
                alignview_proc = subprocess.Popen(sort_read_cmd,
                                                  stdin=bamview_proc.stdout,
                                                  stdout=subprocess.PIPE,
                                                  stderr=open("/dev/null", 'w'))
            else:
                alignview_proc = subprocess.Popen(alignview_cmd,
                                             stdout=subprocess.PIPE,
                                             stderr=open("/dev/null", 'w'))

            # Check deletions' alternatives
            Del_alts = {}
            for _var_i, _var_id in Var_list[gene]:
                _var_type, _var_pos, _var_data = Vars[gene][_var_id]
                if _var_type != "deletion" or _var_pos == 0:
                    continue
                _var_del_len = int(_var_data)
                Del_alts[_var_id] = [[], []]
                debug = (_var_id == "hv93ee0")
                if debug:
                    print Vars[gene][_var_id]
                _var_j = lower_bound(Var_list[gene], _var_pos)
                _latest_pos = _var_pos - 1
                while _var_j < len(Var_list[gene]):
                    _j_pos, _j_id = Var_list[gene][_var_j]
                    if _j_pos > _var_pos + _var_del_len - 1:
                        break
                    if _j_pos < _var_pos:
                        _var_j += 1
                        continue
                    # Check bases between SNPs
                    while _latest_pos < _j_pos:
                        if ref_seq[_latest_pos + 1] != ref_seq[_var_pos + _var_del_len - 1 - (_latest_pos - _var_pos)]:
                            break
                        _latest_pos += 1
                        Del_alts[_var_id][0].append(["ref", _latest_pos])
                    if _latest_pos + 1 < _j_pos:
                        break
                    if _j_pos == _latest_pos:
                        _var_j += 1
                        continue
                    _j_type, _, _j_data = Vars[gene][_j_id]
                    if _j_type == "single":
                        if debug: print Vars[gene][_j_id]
                        _off = _j_pos - _var_pos
                        if debug: print _var_pos + _off, ref_seq[_var_pos + _off]
                        if debug: print _var_pos + _var_del_len + _off, ref_seq[_var_pos + _var_del_len + _off]
                        if _j_data == ref_seq[_var_pos + _var_del_len + _off]:
                            Del_alts[_var_id][0].append([_j_id, _j_pos])
                            _latest_pos = _j_pos
                        else:
                            # Not having single at the same position?
                            if _var_j + 1 >= len(Var_list[gene]) or Var_list[gene][_var_j+1][0] != _j_pos:
                                break
                    _var_j += 1
                _var_j = lower_bound(Var_list[gene], _var_pos + _var_del_len - 1)
                _latest_pos = _var_pos + _var_del_len
                while _var_j >= 0:
                    _j_pos, _j_id = Var_list[gene][_var_j]
                    if _j_pos <= _var_pos:
                        break
                    if _j_pos >= _var_pos + _var_del_len:
                        _var_j -= 1
                        continue
                    # Check bases between SNPs
                    while _latest_pos > _j_pos:
                        if debug: print _latest_pos - 1, ref_seq[_latest_pos - 1], _latest_pos - 1 - _var_del_len, ref_seq[_latest_pos - 1 - _var_del_len]
                        if ref_seq[_latest_pos - 1] != ref_seq[_latest_pos - 1 - _var_del_len]:
                            break
                        _latest_pos -= 1
                        Del_alts[_var_id][1].append(["ref", _latest_pos])
                    if _latest_pos - 1 > _j_pos:
                        break
                    if _j_pos == _latest_pos:
                        _var_j -= 1
                        continue
                    _j_type, _, _j_data = Vars[gene][_j_id]
                    if _j_type == "single":
                        if debug: print Vars[gene][_j_id]
                        _off = _var_pos + _var_del_len - _j_pos
                        if debug: print _var_pos - _off, ref_seq[_var_pos - _off]
                        if debug: print _j_pos, ref_seq[_j_pos]
                        if _j_data == ref_seq[_var_pos - _off]:
                            Del_alts[_var_id][1].append([_j_id, _j_pos])
                            _latest_pos = _j_pos
                        else:
                            # Not having single at the same position?
                            if _var_j >= 0 or Var_list[gene][_var_j-1][0] != _j_pos:
                                break
                    _var_j -= 1
                

                if debug:
                    print Del_alts
                    print "DK :-)"
                    sys.exit(1)


            # Count alleles
            HLA_counts, HLA_cmpt = {}, {}
            coverage = [0 for i in range(len(ref_seq) + 1)]
            num_reads, total_read_len = 0, 0

            # Novel alleles
            N_vars, N_haplotypes = {}, {}
            # Variants to idenity novel alleles
            N_read_vars = []

            # Novel variants used to construct novel alleles
            Novel_vars, Novel_var_ids = {}, {}
            def add_novel_var(Novel_vars, var):
                var_str = "%s-%d-%s" % (var[0], var[1], var[2])
                if var_str in Novel_var_ids:
                    return Novel_var_ids[var_str]
                    
                var_id = "novel%d" % len(Novel_vars)
                Novel_var_ids[var_str] = var_id
                Novel_vars[var_id] = var
                assert len(Novel_var_ids) == len(Novel_vars)
                return var_id

            # Read information
            prev_read_id = None
            prev_exon = False
            prev_right_pos = 0
            if index_type == "graph":
                # Cigar regular expression
                cigar_re = re.compile('\d+\w')
                for line in alignview_proc.stdout:
                    cols = line.strip().split()
                    read_id, flag, chr, pos, mapQ, cigar_str = cols[:6]
                    read_seq, qual = cols[9], cols[10]
                    num_reads += 1
                    total_read_len += len(read_seq)
                    flag, pos = int(flag), int(pos)
                    pos -= (base_locus + 1)
                    if pos < 0:
                        continue

                    # Unalined?
                    if flag & 0x4 != 0:
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

                    if NM > num_mismatch:
                        continue

                    if NH > 1:
                        concordant = False

                    # DK - for debugging purposes
                    debug = False
                    if read_id in ["2339"] and False:
                        debug = True
                        print "read_id: %s)" % read_id, pos, cigar_str, "NM:", NM, MD, Zs
                        print "            ", read_seq

                    vars = []
                    if Zs:
                        vars = Zs.split(',')

                    assert MD != ""
                    MD_str_pos, MD_len = 0, 0
                    read_pos, left_pos = 0, pos
                    right_pos = left_pos
                    cigars = cigar_re.findall(cigar_str)
                    cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
                    cmp_list = []

                    # Extract variants w.r.t backbone from CIGAR string
                    for i in range(len(cigars)):
                        cigar_op, length = cigars[i]
                        if cigar_op == 'M':
                            # Update coverage
                            if enable_coverage:
                                if right_pos + length < len(coverage):
                                    coverage[right_pos] += 1
                                    coverage[right_pos + length] -= 1
                                elif right_pos < len(coverage):
                                    coverage[right_pos] += 1
                                    coverage[-1] -= 1

                            first = True
                            MD_len_used = 0
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
                                    cmp_list.append(["match", right_pos + MD_len_used, length - MD_len_used])
                                    break
                                first = False
                                read_base = read_seq[read_pos + MD_len]
                                MD_ref_base = MD[MD_str_pos]
                                MD_str_pos += 1
                                # DK - for debugging purposes
                                # assert MD_ref_base in "ACGT"
                                cmp_list.append(["match", right_pos + MD_len_used, MD_len - MD_len_used])
                                cmp_list.append(["mismatch", right_pos + MD_len, 1])
                                MD_len_used = MD_len + 1
                                MD_len += 1
                                # Full match
                                if MD_len == length:
                                    MD_len = 0
                                    break
                        elif cigar_op == 'I':
                            cmp_list.append(["insertion", right_pos, length])
                        elif cigar_op == 'D':
                            if MD[MD_str_pos] == '0':
                                MD_str_pos += 1
                            assert MD[MD_str_pos] == '^'
                            MD_str_pos += 1
                            while MD_str_pos < len(MD):
                                if not MD[MD_str_pos] in "ACGT":
                                    break
                                MD_str_pos += 1
                            cmp_list.append(["deletion", right_pos, length])
                        elif cigar_op == 'S':
                            cmp_list.append(["soft", right_pos, length])
                        else:                    
                            assert cigar_op == 'N'
                            cmp_list.append(["intron", right_pos, length])

                        if cigar_op in "MND":
                            right_pos += length

                        if cigar_op in "MIS":
                            read_pos += length

                    exon = False
                    for exon in ref_exons:
                        exon_left, exon_right = exon
                        if right_pos <= exon_left or pos > exon_right:
                            continue
                        else:
                            exon = True
                            break

                    if right_pos > len(ref_seq):
                        continue

                    def add_stat(HLA_cmpt, HLA_counts, HLA_count_per_read, exon = True):
                        max_count = max(HLA_count_per_read.values())
                        cur_cmpt = set()
                        for allele, count in HLA_count_per_read.items():
                            if count < max_count:
                                continue
                            if allele in exclude_allele_list:
                                continue                                
                            cur_cmpt.add(allele)                    
                            if not allele in HLA_counts:
                                HLA_counts[allele] = 1
                            else:
                                HLA_counts[allele] += 1

                        if len(cur_cmpt) == 0:
                            return

                        # DK - for debugging purposes                            
                        # alleles = ["", ""]
                        alleles = ["A*32:29", "A*03:01:07"]
                        allele1_found, allele2_found = False, False
                        for allele, count in HLA_count_per_read.items():
                            if count < max_count:
                                continue
                            if allele == alleles[0]:
                                allele1_found = True
                            elif allele == alleles[1]:
                                allele2_found = True
                        if allele1_found != allele2_found and allele2_found:
                            print alleles[0], HLA_count_per_read[alleles[0]]
                            print alleles[1], HLA_count_per_read[alleles[1]]
                            if allele1_found:
                                print ("%s\tread_id %s - %d vs. %d]" % (alleles[0], prev_read_id, max_count, HLA_count_per_read[alleles[1]]))
                            else:
                                print ("%s\tread_id %s - %d vs. %d]" % (alleles[1], prev_read_id, max_count, HLA_count_per_read[alleles[0]]))
                            print read_seq
                            print line,

                        cur_cmpt = sorted(list(cur_cmpt))
                        cur_cmpt = '-'.join(cur_cmpt)
                        add = 1
                        if not cur_cmpt in HLA_cmpt:
                            HLA_cmpt[cur_cmpt] = add
                        else:
                            HLA_cmpt[cur_cmpt] += add

                    def add_N_stat(N_vars, N_haplotypes, N_read_vars):
                        haplotype_str, haplotype_canonical_str = "", ""
                        # Sanity check
                        for i in range(len(N_read_vars)):
                            N_read_var = N_read_vars[i]
                            if N_read_var[0] in ["left", "right"]:
                                assert i == 0 or i + 1 == len(N_read_vars)
                            if N_read_var[0] in ["unknown", "left", "right"]:
                                var_str = "%s-%d-%d-%s" % (N_read_var[0], N_read_var[1], N_read_var[2], N_read_var[3])
                            else:
                                var_type, var_pos, var_data, var_id = N_read_var
                                var_str = "%s-%d-%s-%s" % (var_type, var_pos, var_data, var_id)
                                if var_str not in N_vars:
                                    N_vars[var_str] = 1
                                else:
                                    N_vars[var_str] += 1
                            if var_str.find("unknown") == -1 and var_str.find("left") == -1 and var_str.find("right") == -1:
                                if haplotype_canonical_str != "":
                                    haplotype_canonical_str += ","
                                haplotype_canonical_str += var_str                                    
                            if haplotype_str != "":
                                haplotype_str += ","
                            haplotype_str += var_str
                        if len(N_read_vars) >= 1:
                            if haplotype_canonical_str not in N_haplotypes:
                                N_haplotypes[haplotype_canonical_str] = [haplotype_str]
                            else:
                                N_haplotypes[haplotype_canonical_str].append(haplotype_str)

                    if read_id != prev_read_id:
                        if prev_read_id != None:
                            add_stat(HLA_cmpt, HLA_counts, HLA_count_per_read, prev_exon)
                            add_N_stat(N_vars, N_haplotypes, N_read_vars)
                            N_read_vars = []

                        HLA_count_per_read = {}
                        for HLA_name in HLA_names[gene]:
                            if HLA_name.find("BACKBONE") != -1:
                                continue
                            HLA_count_per_read[HLA_name] = 0

                    def add_count(var_id, add):
                        if partial and \
                                exonic_only and \
                                not var_in_exon(Vars[gene][var_id], ref_exons):
                            return
                        assert var_id in Links
                        alleles = Links[var_id]
                        for allele in alleles:
                            if allele.find("BACKBONE") != -1:
                                continue
                            if partial and exonic_only:
                                assert allele in allele_reps
                                if allele != allele_reps[allele]:
                                    continue
                            HLA_count_per_read[allele] += add
                            # DK - for debugging purposes
                            if debug or True:
                                if allele in ["A*32:29", "A*03:01:07"]:
                                    print allele, add, var_id

                    def add_N_var(N_read_vars, N_read_var):
                        if len(N_read_vars) == 0:
                            N_read_vars.append(N_read_var)
                        else:
                            var_type, var_pos, var_data, var_id = N_read_vars[-1]
                            var_type2, var_pos2, var_data, var_id2 = N_read_var
                            if var_pos < var_pos2:
                                N_read_vars.append(N_read_var)

                    # Decide which allele(s) a read most likely came from
                    # also sanity check - read length, cigar string, and MD string
                    for var_id, data in Vars[gene].items():
                        var_type, var_pos, var_data = data
                        if var_type != "deletion":
                            continue
                        if left_pos >= var_pos and right_pos <= var_pos + int(var_data):
                            add_count(var_id, -1)

                    if concordant and concordant_assembly:
                        concordant_second_read = (len(N_read_vars) > 0)
                    else:
                        concordant_second_read = False
                        add_N_stat(N_vars, N_haplotypes, N_read_vars)
                        N_read_vars = []
                         
                    ref_pos, read_pos, cmp_cigar_str, cmp_MD = left_pos, 0, "", ""
                    cigar_match_len, MD_match_len = 0, 0
                    first_diff = True
                    for cmp_i in range(len(cmp_list)):
                        cmp = cmp_list[cmp_i]
                        type = cmp[0]
                        length = cmp[2]
                        unknown_in_pair = []
                        if type in ["match", "mismatch", "deletion", "insertion"]:
                            if concordant_second_read and first_diff:
                                if prev_right_pos < left_pos - 1:
                                    add_N_var(N_read_vars, ["unknown", prev_right_pos, left_pos - 1, "unknown"])
                            first_diff = False

                        if type == "match":
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos + length <= var_pos:
                                    break
                                if ref_pos <= var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "insertion":
                                        if ref_pos < var_pos and ref_pos + length > var_pos + len(var_data):
                                            add_count(var_id, -1)
                                            # DK - for debugging purposes
                                            if debug:
                                                print cmp, var_id, Links[var_id]
                                    elif var_type == "deletion":
                                        del_len = int(var_data)
                                        if ref_pos < var_pos and ref_pos + length > var_pos + del_len:
                                            # DK - for debugging purposes
                                            if debug:
                                                print cmp, var_id, Links[var_id], -1, Vars[gene][var_id]
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
                                                add_count(var_id, -1)
                                    else:
                                        if debug:
                                            print cmp, var_id, Links[var_id], -1
                                        add_count(var_id, -1)
                                var_idx += 1
                            if cmp_i == 0 or (cmp_i == 1 and cmp_list[0][0] == "soft"):
                                if not concordant or not concordant_second_read:
                                    if ref_pos > ref_pos + length - 1:
                                        # DK - should fix the following
                                        None
                                    else:
                                        if length > 20:
                                            add_N_var(N_read_vars, ["left", ref_pos + 10, ref_pos + 10, "left"])                                            
                            if cmp_i + 1 == len(cmp_list) or (cmp_i + 2 == len(cmp_list) and cmp_list[-1][0] == "soft"):
                                if not concordant or concordant_second_read:
                                    if length > 20:
                                        add_N_var(N_read_vars, ["right", ref_pos + length - 1 - 10, ref_pos + length - 1 - 10, "right"])

                            read_pos += length
                            ref_pos += length
                            cigar_match_len += length
                            MD_match_len += length
                        elif type == "mismatch":
                            read_base = read_seq[read_pos]
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            known_var = False
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos < var_pos:
                                    break
                                if ref_pos == var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "single":
                                        if var_data == read_base:
                                            # DK - for debugging purposes
                                            if debug:
                                                print cmp, var_id, 1, var_data, read_base, Links[var_id]

                                            # DK - for debugging purposes
                                            if False:
                                                read_qual = ord(qual[read_pos])
                                                add_count(var_id, (read_qual - 60) / 60.0)
                                            else:
                                                add_count(var_id, 1)
                                                known_var = True
                                                add_N_var(N_read_vars, [var_type, ref_pos, read_base, var_id])
                                        # DK - check out if this routine is appropriate
                                        # else:
                                        #    add_count(var_id, -1)
                                var_idx += 1

                            if not known_var:
                                novel_var_id = add_novel_var(Novel_vars, ["single", ref_pos, read_base])
                                add_N_var(N_read_vars, ["single", ref_pos, read_base, novel_var_id])
                                
                            cmp_MD += ("%d%s" % (MD_match_len, ref_seq[ref_pos]))
                            MD_match_len = 0
                            cigar_match_len += 1
                            read_pos += 1
                            ref_pos += 1
                        elif type == "insertion":
                            ins_seq = read_seq[read_pos:read_pos+length]
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            # DK - for debugging purposes
                            if debug:
                                print left_pos, cigar_str, MD, vars
                                print ref_pos, ins_seq, Var_list[gene][var_idx], Vars[gene][Var_list[gene][var_idx][1]]
                                # sys.exit(1)
                            known_var = False
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos < var_pos:
                                    break
                                if ref_pos == var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "insertion":                                
                                        if var_data == ins_seq:
                                            # DK - for debugging purposes
                                            if debug:
                                                print cmp, var_id, 1, Links[var_id]
                                            add_count(var_id, 1)
                                            known_var = True
                                            add_N_var(N_read_vars, [var_type, ref_pos, ins_seq, var_id])
                                var_idx += 1

                            if not known_var:
                                novel_var_id = add_novel_var(Novel_vars, ["insertion", ref_pos, ins_seq])
                                add_N_var(N_read_vars, ["insertion", ref_pos, ins_seq, novel_var_id])

                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            read_pos += length
                            cmp_cigar_str += ("%dI" % length)
                        elif type == "deletion":
                            alt_match = False
                            del_len = length
                            # Check if this deletion can be removed
                            if (cmp_i == 1 and cmp_list[0][0] == "match") or \
                                    (cmp_i == 2 and cmp_list[0][0] == "soft" and cmp_list[0][1] == "match"):
                                _match_len = cmp_list[cmp_i-1][2]
                                temp_ref_pos = ref_pos
                                while _match_len > 0:
                                    last_bp = ref_seq[temp_ref_pos + del_len - 1]
                                    prev_bp = ref_seq[temp_ref_pos - 1]
                                    if last_bp != prev_bp:
                                        break
                                    temp_ref_pos -= 1
                                    _match_len -= 1
                                if _match_len == 0:
                                    alt_match = True
                            if (cmp_i + 2 == len(cmp_list) and cmp_list[-1][0] == "match") or \
                                    (cmp_i + 3 == len(cmp_list) and cmp_list[-1][0] == "soft" and cmp_list[-2][0] == "match"):
                                _match_len = cmp_list[cmp_i+1][2]
                                temp_ref_pos = ref_pos
                                """
                                while _match_len > 0:
                                    last_bp = ref_seq[temp_ref_pos + del_len - 1]
                                    prev_bp = ref_seq[temp_ref_pos - 1]
                                    if last_bp != prev_bp:
                                        break
                                    temp_ref_pos -= 1
                                    _match_len -= 1
                                if _match_len == 0:
                                    by_chance = True
                                """
                                
                            # Deletions can be shifted bidirectionally by HISAT2
                            temp_ref_pos = ref_pos
                            while temp_ref_pos > 0:
                                last_bp = ref_seq[temp_ref_pos + del_len - 1]
                                prev_bp = ref_seq[temp_ref_pos - 1]
                                if last_bp != prev_bp:
                                    break
                                temp_ref_pos -= 1
                            known_var = False
                            var_idx = lower_bound(Var_list[gene], temp_ref_pos)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if temp_ref_pos < var_pos:
                                    first_bp = ref_seq[temp_ref_pos]
                                    next_bp = ref_seq[temp_ref_pos + del_len]
                                    if first_bp == next_bp:
                                        temp_ref_pos += 1
                                        continue
                                    else:
                                        break
                                if temp_ref_pos == var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "deletion":
                                        var_len = int(var_data)
                                        if var_len == length:
                                            if debug:
                                                print cmp, var_id, 1, Links[var_id]
                                                print ref_seq[var_pos - 10:var_pos], ref_seq[var_pos:var_pos+int(var_data)], ref_seq[var_pos+int(var_data):var_pos+int(var_data)+10]
                                            if not alt_match:
                                                add_count(var_id, 1)
                                            known_var = True
                                            add_N_var(N_read_vars, [var_type, ref_pos, del_len, var_id])
                                var_idx += 1

                            if not known_var:
                                novel_var_id = add_novel_var(Novel_vars, ["deletion", ref_pos, str(del_len)])
                                add_N_var(N_read_vars, ["deleltion", ref_pos, del_len, novel_var_id])

                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            cmp_MD += ("%d" % MD_match_len)
                            MD_match_len = 0
                            cmp_cigar_str += ("%dD" % length)
                            cmp_MD += ("^%s" % ref_seq[ref_pos:ref_pos+length])
                            ref_pos += length
                        elif type == "soft":
                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            read_pos += length
                            cmp_cigar_str += ("%dS" % length)
                        else:
                            assert type == "intron"
                            if cigar_match_len > 0:
                                cmp_cigar_str += ("%dM" % cigar_match_len)
                                cigar_match_len = 0
                            cmp_cigar_str += ("%dN" % length)
                            ref_pos += length                    
                    if cigar_match_len > 0:
                        cmp_cigar_str += ("%dM" % cigar_match_len)
                    cmp_MD += ("%d" % MD_match_len)
                    # Sanity check
                    if read_pos != len(read_seq) or \
                            cmp_cigar_str != cigar_str or \
                            cmp_MD != MD:
                        print >> sys.stderr, "Error:", cigar_str, MD
                        print >> sys.stderr, "\tcomputed:", cmp_cigar_str, cmp_MD
                        print >> sys.stderr, "\tcmp list:", cmp_list
                        assert False

                    prev_read_id = read_id
                    prev_exon = exon
                    prev_right_pos = right_pos

                if num_reads <= 0:
                    continue

                if prev_read_id != None:
                    add_stat(HLA_cmpt, HLA_counts, HLA_count_per_read)
                    add_N_stat(N_vars, N_haplotypes, N_read_vars)
                    N_read_vars = []

                # Delete N vars with low evidence
                var_count_avg = sum(N_vars.values()) / float(len(N_vars))
                N_ignored_vars = set()
                for N_var, count in N_vars.items():
                    _var_type, _, _, _var_id = N_var.split('-')
                    if _var_id.find("novel"):
                        if _var_type in ["deleltion", "insertion"]:
                            if count < var_count_avg / 8.0:
                                N_ignored_vars.add(N_var)
                        else:
                            if count < var_count_avg / 6.0:
                                N_ignored_vars.add(N_var)
                    else:
                        if _var_type in ["deleltion", "insertion"]:
                            if count < var_count_avg / 10.0:
                                N_ignored_vars.add(N_var)
                        else:
                            if count < var_count_avg / 8.0:
                                N_ignored_vars.add(N_var)

                    # DK - for debugging purposes
                    if count <= 3:
                        N_ignored_vars.add(N_var)

                if verbose >= 1:
                    print "Number of N haplotypes: %d (before cleaning)" % len(N_haplotypes)

                cleaned_N_haplotypes = {}
                for _haplotypes in N_haplotypes.values():
                    for _haplotype in _haplotypes:
                        _vars = _haplotype.split(',')
                        cleaned_vars = []
                        cleaned_canonical_vars = []
                        filtered = False
                        for _var in _vars:
                            if _var.find("unknown") != -1 or \
                                    _var.find("left") != -1 or \
                                    _var.find("right") != -1:
                                cleaned_vars.append(_var)
                            elif _var not in N_ignored_vars:
                                cleaned_vars.append(_var)
                                cleaned_canonical_vars.append(_var)
                            else:
                                filtered = True
                                break

                        # DK - for debugging purposes
                        if filtered:
                            continue

                        if len(cleaned_vars) > 0:
                            if cleaned_vars[0].find("unknown") != -1:
                                cleaned_vars = cleaned_vars[1:]
                            elif cleaned_vars[-1].find("unknown") != -1:
                                cleaned_vars = cleaned_vars[:1]

                        if len(cleaned_vars) > 0:
                            assert cleaned_vars[0].find("unknown") == -1 and cleaned_vars[-1].find("unknown") == -1
                            cleaned_haplotype = ','.join(cleaned_vars)
                            cleaned_canonical_haplotype = ','.join(cleaned_canonical_vars)
                            if cleaned_canonical_haplotype == "":
                                cleaned_canonical_haplotype = "nothing-0-0"
                            if cleaned_canonical_haplotype not in cleaned_N_haplotypes:
                                cleaned_N_haplotypes[cleaned_canonical_haplotype] = [cleaned_haplotype]
                            else:
                                cleaned_N_haplotypes[cleaned_canonical_haplotype].append(cleaned_haplotype)

                N_haplotypes = cleaned_N_haplotypes

                def var_cmp(a, b):
                    a_left = int(a.split('-')[1])
                    b_left = int(b.split('-')[1])
                    return a_left - b_left

                N_var_list = list(N_vars.keys())
                N_var_list = sorted(N_var_list, cmp=var_cmp)

                if verbose:
                    print "Number of vars: %d (filtered: %d)" % (len(N_vars), len(N_ignored_vars))
                    for N_var in N_var_list:
                        print "%s: %d (removed: %s)" % (N_var, N_vars[N_var], N_var in N_ignored_vars)
                        
                N_haplotype_list = set()
                for _haplotypes in N_haplotypes.values():
                    for _haplotype in _haplotypes:
                        N_haplotype_list.add(_haplotype)
                N_haplotype_list = list(N_haplotype_list)
                N_haplotype_list = sorted(N_haplotype_list, cmp=haplotype_cmp)

                if verbose >= 2:
                    v_N_haplotype_list = list(N_haplotypes.keys())
                    v_N_haplotype_list = sorted(v_N_haplotype_list, cmp=haplotype_cmp)
                    print "Number of haplotypes: %d" % len(N_haplotypes)
                    for N_haplotype in v_N_haplotype_list:
                        print N_haplotype, len(N_haplotypes[N_haplotype])
                        if verbose >= 3:
                            for read_haplotype in N_haplotypes[N_haplotype]:
                                print "\t\t", read_haplotype

                if assembly:
                    if simulation:
                        num_assemblies = len(test_HLA_names)
                    else:
                        # DK - for debugging purposes
                        num_assemblies = 10
                    N_alleles = assemble_alleles(N_haplotype_list, num_assemblies)
                    if verbose >= 2:
                        print "Number of N alleles: %d" % len(N_alleles)
                        for N_allele in N_alleles:
                            print "\t", N_allele
                else:
                    N_alleles = []

                # Coverage
                # it is not used by the default
                if enable_coverage:
                    assert num_reads > 0
                    read_len = int(total_read_len / float(num_reads))
                    coverage_sum = 0
                    for i in range(len(coverage)):
                        if i > 0:
                            coverage[i] += coverage[i-1]
                        coverage_sum += coverage[i]
                    coverage_avg = coverage_sum / float(len(coverage))
                    assert len(ref_seq) < len(coverage)
                    for i in range(len(ref_seq)):
                        coverage_threshold = 1.0 * coverage_avg
                        if i < read_len:
                            coverage_threshold *= ((i+1) / float(read_len))
                        elif i + read_len > len(ref_seq):
                            coverage_threshold *= ((len(ref_seq) - i) / float(read_len))
                        if coverage[i] >= coverage_threshold:
                            continue
                        pseudo_num_reads = (coverage_threshold - coverage[i]) / read_len
                        var_idx = lower_bound(Var_list[gene], i + 1)
                        if var_idx >= len(Var_list[gene]):
                            var_idx = len(Var_list[gene]) - 1
                        cur_cmpt = set()
                        while var_idx >= 0:
                            var_pos, var_id = Var_list[gene][var_idx]
                            var_type, _, var_data = Vars[gene][var_id]
                            if var_type == "deletion":
                                del_len = int(var_data)
                                if i < var_pos:
                                    break
                                if i + read_len < var_pos + int(var_data):
                                    assert var_id in Links
                                    cur_cmpt = cur_cmpt.union(set(Links[var_id]))
                            var_idx -= 1
                        if cur_cmpt:
                            cur_cmpt = '-'.join(list(cur_cmpt))
                            if not cur_cmpt in HLA_cmpt:
                                HLA_cmpt[cur_cmpt] = 0
                            HLA_cmpt[cur_cmpt] += pseudo_num_reads
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

            HLA_prob = single_abundance(HLA_cmpt, HLA_lengths[gene])

            success = [False for i in range(len(test_HLA_names))]
            found_list = [False for i in range(len(test_HLA_names))]
            for prob_i in range(len(HLA_prob)):
                prob = HLA_prob[prob_i]
                found = False
                allele_haplotype = get_allele(gene, prob[0], Vars_default, Var_list_default, Links_default)
                if allele_haplotype != "":
                    covered, total = calculate_allele_coverage(allele_haplotype, N_haplotypes, ref_exons, partial, exonic_only, verbose >= 3)
                else:
                    covered, total = 0, 0
                _allele_rep = prob[0]
                if partial and exonic_only:
                    _fields = _allele_rep.split(':')
                    if len(_fields) == 4:
                        _allele_rep = ':'.join(_fields[:-1])

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
                            print >> sys.stderr, "\t\t\t*** %d ranked %s (abundance: %.2f%%, vars_covered: %d/%d)" % (rank_i + 1, test_HLA_name, prob[1] * 100.0, covered, total)
                            if rank_i < len(success):
                                success[rank_i] = True
                            found_list[name_i] = True
                            found = True
                    # DK - for debugging purposes
                    if not False in found_list and prob_i >= 10:
                        break
                if not found:
                    print >> sys.stderr, "\t\t\t\t%d ranked %s (abundance: %.2f%%, vars_covered: %d/%d)" % (prob_i + 1, _allele_rep, prob[1] * 100.0, covered, total)
                    if best_alleles and prob_i < 2:
                        print >> sys.stdout, "SingleModel %s (abundance: %.2f%%)" % (_allele_rep, prob[1] * 100.0)
                if not simulation and prob_i >= 9:
                    break
                if prob_i >= 19:
                    break
            print >> sys.stderr

            if verbose >= 2 and index_type == "graph":
                if verbose >= 3:
                    allele_haplotypes = []
                    for HLA_name, abundance in HLA_prob:
                        if abundance >= 0.03:
                            allele_haplotypes.append(get_allele(gene, HLA_name, Vars_default, Var_list_default, Links_default))

                    for read_haplotype in v_N_haplotype_list:
                        for read_haplotype in N_haplotypes[read_haplotype]:
                            assembled = False
                            for allele_haplotype in allele_haplotypes:
                                # DK - for debugging purposes
                                debug = (read_haplotype == "left-199-199-left,single-245-T-hv54,unknown-289-429-unknown,single-448-G-hv99,single-452-G-hv100,single-478-A-hv108,single-492-G-hv112,right-519-519-right") and \
                                    (allele_haplotype.find("single-452-G-hv100,single-478-A-hv108,single-492-G-hv112") != -1)
                                debug = False
                                if haplotype_cmp(allele_haplotype, read_haplotype) <= 0:
                                    _, assembled = assemble_two_haplotypes(allele_haplotype.split(','), read_haplotype.split(','), debug)
                                else:
                                    _, assembled = assemble_two_haplotypes(read_haplotype.split(','), allele_haplotype.split(','), debug)
                                if assembled:
                                    break
                            # DK - for debugging purposes
                            if not assembled:
                                print "Not compatible", read_haplotype

                N_alleles = list(N_alleles)
                for N_allele_i in range(len(N_alleles)):
                    N_allele = N_alleles[N_allele_i]
                    same = True
                    _vars = N_allele.split(',')

                    backbone_seq = HLAs[gene][ref_allele]
                    backbone_seq_default = HLAs_default[gene][ref_allele]
                    assert backbone_seq == backbone_seq_default
                    if backbone_seq == backbone_seq_default:
                        calculated_var_ids = []
                        calculated_vars = set()
                        for _var in _vars:
                            if _var.find("unknown") != -1 or \
                                    _var.find("left") != -1 or \
                                    _var.find("right") != -1:
                                continue
                            _var_type, _var_pos, _var_data, _var_id = _var.split('-')
                            calculated_var_ids.append(_var_id)
                            calculated_vars.add("%s-%s" % (_var_type, _var_pos))

                        cmp_allele_name = ""
                        num_common = 0
                        if simulation:
                            cmp_HLA_names = test_HLA_names
                        else:
                            cmp_HLA_names = []
                            for HLA_name, abundance in HLA_prob:
                                if abundance >= 0.03:
                                    cmp_HLA_names.append(HLA_name)

                        for test_allele_name in cmp_HLA_names:
                            cmp_vars = set()
                            for _var_id in Var_list_default[gene]:
                                _var_pos, _var_id = _var_id
                                if test_allele_name in Links_default[_var_id]:
                                    _var_type, _var_pos, _= Vars_default[gene][_var_id]
                                    cmp_vars.add("%s-%d" % (_var_type, _var_pos))
                            tmp_num_common = len(calculated_vars.intersection(cmp_vars))
                            if num_common < tmp_num_common:
                                num_common = tmp_num_common
                                cmp_allele_name = test_allele_name

                        # DK - for debugging purposes
                        # cmp_allele_name = "A*32:01:01"
                        assert cmp_allele_name != ""

                        print "\t%d) vs. %s" % (N_allele_i + 1, cmp_allele_name)
                        
                        true_var_ids = []
                        for _var_id in Var_list_default[gene]:
                            _var_pos, _var_id = _var_id
                            if cmp_allele_name in Links_default[_var_id]:
                                true_var_ids.append(_var_id)

                        _i, _j, _k = 0, 0, 0
                        while _i < len(calculated_var_ids) or _j < len(true_var_ids):
                            if _i == len(calculated_var_ids):
                                t_var_id = true_var_ids[_j]
                                t_var = Vars_default[gene][t_var_id]
                                print "%4d\t<> %d %s %s (%s)" % (_k, t_var[1], t_var[0], t_var[2], t_var_id)
                                _j += 1
                                _k += 1
                            elif _j == len(true_var_ids):
                                c_var_id = calculated_var_ids[_i]
                                if c_var_id.find("novel") != -1:
                                    c_var = Novel_vars[c_var_id]
                                else:
                                    c_var = Vars[gene][c_var_id]
                                print "%4d\t%d %s %s (%s) <>" % (_k, c_var[1], c_var[0], c_var[2], c_var_id)
                                _i += 1
                                _k += 1
                            else:
                                c_var_id = calculated_var_ids[_i]
                                if c_var_id.find("novel") != -1:
                                    c_var = Novel_vars[c_var_id]
                                else:
                                    c_var = Vars[gene][c_var_id]
                                t_var_id = true_var_ids[_j]
                                t_var = Vars_default[gene][t_var_id]

                                if c_var == t_var:
                                    print "%4d\tsame: %d %s %s (%s)" % (_k, c_var[1], c_var[0], c_var[2], c_var_id)
                                    _i += 1
                                    _j += 1
                                else:
                                    if _i > 0 and _j > 0:
                                        same = False
                                    if c_var[1] <= t_var[1]:
                                        _i += 1
                                        print "%4d\t%d %s %s (%s) <>" % (_k, c_var[1], c_var[0], c_var[2], c_var_id)
                                    else:
                                        _j += 1
                                        print "%4d\t<> %d %s %s (%s)" % (_k, t_var[1], t_var[0], t_var[2], t_var_id)
                                _k += 1
                    else:
                        _var_ids = []
                        for _var in _vars:
                            _var_id = _var.split('-')[2]
                            _var_ids.append(_var_id)
                        calculated_allele_seq = construct_allele_seq(backbone_seq, _var_ids, Vars[gene])
                        true_allele_seq = HLAs_default[gene][test_allele_name]

                        print "calculated allele: %d bp" % len(calculated_allele_seq)
                        print "true allele: %d bp" % len(true_allele_seq)
                        if calculated_allele_seq != true_allele_seq:
                            if calculated_allele_seq.find(true_allele_seq) == -1 and \
                                    true_allele_seq.find(calculated_allele_seq) == -1:
                                same = False

                    if same:
                        print "Same as %s!!" % cmp_allele_name
                    else:
                        print "Different!!"


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

        if remove_alignment_file:
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
                        enable_coverage,
                        best_alleles,
                        exclude_allele_list,
                        default_allele_list,
                        num_mismatch,
                        assembly,
                        concordant_assembly,
                        exonic_only,
                        verbose,
                        daehwan_debug):
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(test_HLA_genotyping))
    ex_path = os.path.dirname(curr_script)

    # Clone a git repository, IMGTHLA
    if not os.path.exists("IMGTHLA"):
        os.system("git clone https://github.com/jrob119/IMGTHLA.git")

    simulation = (read_fname == [] and alignment_fname == "")

    def check_files(fnames):
        for fname in fnames:
            if not os.path.exists(fname):
                return False
        return True

    # Download HISAT2 index
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]
    if not check_files(HISAT2_fnames):
        os.system("wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz; tar xvzf grch38.tar.gz; rm grch38.tar.gz")
        hisat2_inspect = os.path.join(ex_path, "hisat2-inspect")
        os.system("%s grch38/genome > genome.fa" % hisat2_inspect)
        os.system("samtools faidx genome.fa")

    # Check if the pre-existing files (hla*) are compatible with the current parameter setting
    if os.path.exists("hla.ref"):
        left = 0
        HLA_genes = set()
        BACKBONE = False
        for line in open("hla.ref"):
            HLA_name = line.strip().split()[0]
            if HLA_name.find("BACKBONE") != -1:
                BACKBONE = True
            HLA_gene = HLA_name.split('*')[0]
            HLA_genes.add(HLA_gene)
        delete_hla_files = False
        if reference_type == "gene":
            if not BACKBONE:
                delete_hla_files = True
        elif reference_type in ["chromosome", "genome"]:
            if BACKBONE:
                delete_hla_files = True
        else:
            assert False
        if not set(hla_list).issubset(HLA_genes):
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
                  base_fname + ".haplotype",
                  base_fname + ".link",
                  base_fname + "_alleles_excluded.txt"]
    
    # Check if excluded alleles in current files match
    excluded_alleles_match = False
    if(os.path.exists(HLA_fnames[6])):
        afile = open(HLA_fnames[6],'r')
        afile.readline()
        lines = afile.read().split()
        excluded_alleles_match = (set(exclude_allele_list) == set(lines))
        afile.close()
    elif len(exclude_allele_list) == 0:
        excluded_alleles_match = True
        try:
            temp_name = HLA_fnames[6]
            HLA_fnames.remove(HLA_fnames[6])
            os.remove(temp_name)
        except OSError:
            pass
        
    if not excluded_alleles_match:
        print("Creating Allele Exclusion File.\n")
        afile = open(HLA_fnames[6],'w')
        afile.write("Alleles excluded:\n")
        afile.write("\n".join(exclude_allele_list))
        afile.close()

    if verbose >= 1:
        print >> sys.stderr, HLA_fnames
    
    if (not check_files(HLA_fnames)) or (not excluded_alleles_match) :
        extract_hla_script = os.path.join(ex_path, "hisatgenotype_extract_vars.py")
        extract_cmd = [extract_hla_script,
                       "--reference-type", reference_type,
                       "--hla-list", ','.join(hla_list)]

        if len(exclude_allele_list) > 0:
            print exclude_allele_list
            extract_cmd += ["--exclude-allele-list", ",".join(exclude_allele_list)]

        if len(base_fname) > 3:
            extract_cmd += ["--base", base_fname]

        if partial:
            extract_cmd += ["--partial"]
        extract_cmd += ["--inter-gap", "30",
                        "--intra-gap", "50"]
        if verbose >= 1:
            print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
        proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()
        
        if not check_files(HLA_fnames):
            print >> sys.stderr, "Error: extract_HLA_vars failed!"
            sys.exit(1)

    for aligner, index_type in aligners:
        if aligner == "hisat2":
            # Build HISAT2 graph indexes based on the above information
            if index_type == "graph":
                HLA_hisat2_graph_index_fnames = ["%s.graph.%d.ht2" % (base_fname, i+1) for i in range(8)]
                if not check_files(HLA_hisat2_graph_index_fnames) or (not excluded_alleles_match):
                    hisat2_build = os.path.join(ex_path, "hisat2-build")
                    build_cmd = [hisat2_build,
                                 "-p", str(threads),
                                 "--snp", "%s.snp" % base_fname,
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
                if reference_type == "gene" and (not check_files(HLA_hisat2_linear_index_fnames) or (not excluded_alleles_match)):
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
            if reference_type == "gene" and (not check_files(HLA_bowtie2_index_fnames) or (not excluded_alleles_match)):
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
        allele_name = line.split()[1][4:-1]
        gene = allele_name.split('*')[0]
        if line.find("partial") != -1:
            partial_alleles.add(allele_name)

    if len(default_allele_list) > 0:
        if not os.path.exists("Default-HLA/hla_backbone.fa"):
            try:
                os.mkdir("Default-HLA")
            except:
                pass
            #os.chdir(current_path + "/Default-HLA")
            
            extract_hla_script = os.path.join(ex_path, "hisatgenotype_extract_vars.py")
            extract_cmd = [extract_hla_script,
                           "--reference-type", reference_type,
                           "--hla-list", ','.join(hla_list),
                           "--base", "Default-HLA/hla"]

            if partial:
                extract_cmd += ["--partial"]
            extract_cmd += ["--inter-gap", "30",
                            "--intra-gap", "50"]
            if verbose >= 1:
                print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
            proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
            proc.communicate()
            
            if not os.path.exists("Default-HLA/hla_backbone.fa"):
                print >> sys.stderr, "Error: extract_HLA_vars (Default) failed!"
                sys.exit(1)
    
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

    # Construct excluded alleles (Via default backbone data)
    custom_allele_check = False
    if len(default_allele_list) > 0:
        custom_allele_check = True
        HLAs_default = {}
        read_HLA_alleles("Default-HLA/hla_backbone.fa", HLAs_default)
        read_HLA_alleles("Default-HLA/hla_sequences.fa", HLAs_default)
        #HLA_lengths_default = {}

        for HLA_gene, HLA_alleles in HLAs_default.items():
            for allele_name, seq in HLA_alleles.items():
                if allele_name in default_allele_list:
                    HLA_lengths[HLA_gene][allele_name] = len(seq)
        
        #for allele_name, seq in HLAs_default.items():
         #   if allele_name in default_allele_list:
          #      HLA_lengths[allele_name] = len(seq)
            #if (allele_name in default_allele_list):
            #    HLA_lengths_default[allele_name] = len(seq)
    else:
        HLAs_default = HLAs

    # Read HLA variants, and link information
    Vars, Var_list = read_HLA_vars("%s.snp" % base_fname, reference_type)
    Links = read_HLA_links("%s.link" % base_fname)
    Vars_default, Var_list_default, Links_default = {}, {}, {}
    if len(default_allele_list) > 0:
        Vars_default, Var_list_default = read_HLA_vars("Default-HLA/hla.snp", reference_type)
        Links_default = read_HLA_links("Default-HLA/hla.link")
    else:
        Vars_default, Var_list_default = Vars, Var_list
        Links_default = Links

    # DK - for debugging purposes
    """
    print "Links:", len(Links)
    DK_set = set(Links["hv252"])
    for var_id in ["hv252", "hv254", "hv265", "hv269", "hv272", "hv273", "hv274", "hv275", "hv276", "hv277", "hv278", "hv279", "hv280", "hv281", "hv282", "hv283", "hv284", "hv285", "hv286", "hv287", "hv288", "hv289", "hv290", "hv291", "hv292"]:
        print var_id,  set(Links[var_id])
        # print var_id, DK_set
    sys.exit(1)
    """
        
    # Scoring schemes from Sangtae Kim (Illumina)'s implementation
    # Currently not used.
    """
    max_qual_value = 100
    match_score, mismatch_score = [0] * max_qual_value, [0] * max_qual_value
    for qual in range(max_qual_value):
        error_rate = 0.1 ** (qual / 10.0)
        match_score[qual] = math.log(1.000000000001 - error_rate);
        mismatch_score[qual] = math.log(error_rate / 3.0);
    """
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
            if custom_allele_check:
                for allele in default_allele_list:
                    test_list.append([[allele]])
            else:
                for gene in genes:
                    HLA_gene_alleles = HLA_names[gene]
                    for HLA_name in HLA_gene_alleles:
                        if HLA_name.find("BACKBONE") != -1:
                            continue
                        test_list.append([[HLA_name]])
        if pair_test:
            test_size = 500
            allele_count = 2
            if custom_allele_check:
                if (default_allele_list) < allele_count:
                    print >> sys.stderr, "# of default alleles (%d) is at least %d" % (len(defeault_allele_list), allele_count)
                    sys.exit(1)
                    
                for test_i in range(1):
                    random.shuffle(default_allele_list)
                    test_pair = [default_allele_list[:allele_count]]
                    test_list.append(test_pair)
            else:
                for test_i in range(test_size):
                    test_pairs = []
                    for gene in genes:
                        HLA_gene_alleles = []

                        for allele in HLA_names[gene]:
                            if allele.find("BACKBONE") != -1:
                                continue
                            HLA_gene_alleles.append(allele)
                        nums = [i for i in range(len(HLA_gene_alleles))]
                        random.shuffle(nums)
                        test_pairs.append(sorted([HLA_gene_alleles[nums[i]] for i in range(allele_count)]))
                    test_list.append(test_pairs)

        for test_i in range(len(test_list)):
            if "test_id" in daehwan_debug:
                daehwan_test_ids = daehwan_debug["test_id"].split('-')
                if str(test_i + 1) not in daehwan_test_ids:
                    continue

            print >> sys.stderr, "Test %d" % (test_i + 1)
            test_HLA_list = test_list[test_i]
           
            # DK - for debugging purposes
            test_HLA_list = [["A*32:29"]]
            for test_HLA_names in test_HLA_list:
                for test_HLA_name in test_HLA_names:
                    if custom_allele_check:
                        gene = test_HLA_name.split('*')[0]
                        test_HLA_seq = HLAs_default[gene][test_HLA_name]
                        seq_type = "partial" if test_HLA_name in partial_alleles else "full"
                        print >> sys.stderr, "\t%s - %d bp (%s sequence)" % (test_HLA_name, len(test_HLA_seq), seq_type)
                        continue
                    gene = test_HLA_name.split('*')[0]
                    test_HLA_seq = HLAs[gene][test_HLA_name]
                    seq_type = "partial" if test_HLA_name in partial_alleles else "full"
                    print >> sys.stderr, "\t%s - %d bp (%s sequence)" % (test_HLA_name, len(test_HLA_seq), seq_type)
                    
            if custom_allele_check:
                simulate_reads(HLAs_default, test_HLA_list, simulate_interval)
            else:
                simulate_reads(HLAs, test_HLA_list, simulate_interval)

            if "single-end" in daehwan_debug:
                read_fname = ["hla_input_1.fa"]
            else:
                read_fname = ["hla_input_1.fa", "hla_input_2.fa"]

            fastq = False
            tmp_test_passed = HLA_typing(ex_path,
                                         simulation,
                                         reference_type,
                                         test_HLA_list,
                                         partial,
                                         refHLAs,
                                         HLAs,                       
                                         HLA_names,
                                         HLA_lengths,
                                         refHLA_loci,
                                         Vars,
                                         Var_list,
                                         Links,
                                         HLAs_default,
                                         Vars_default,
                                         Var_list_default,
                                         Links_default,
                                         exclude_allele_list,
                                         aligners,
                                         num_mismatch,
                                         assembly,
                                         concordant_assembly,
                                         exonic_only,
                                         fastq,
                                         read_fname,
                                         alignment_fname,
                                         threads,
                                         enable_coverage,
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
        HLA_typing(ex_path,
                   simulation,
                   reference_type,
                   hla_list,
                   partial,
                   refHLAs,
                   HLAs,                       
                   HLA_names,
                   HLA_lengths,
                   refHLA_loci,
                   Vars,
                   Var_list,
                   Links,
                   HLAs_default,
                   Vars_default,
                   Var_list_default,
                   Links_default,
                   exclude_allele_list,
                   aligners,
                   num_mismatch,
                   assembly,
                   concordant_assembly,
                   exonic_only,
                   fastq,
                   read_fname,
                   alignment_fname,
                   threads,
                   enable_coverage,
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
    parser.add_argument("--default-list",
                        dest = "default_allele_list",
                        type=str,
                        default="",
                        help="A comma-separated list of HLA alleles to be tested. Alleles are retrieved from default backbone data (all alleles included in backbone).")
    parser.add_argument("--reference-type",
                        dest="reference_type",
                        type=str,
                        default="gene",
                        help="Reference type: gene, chromosome, and genome (default: gene)")
    parser.add_argument("--hla-list",
                        dest="hla_list",
                        type=str,
                        default="A,B,C,DQA1,DQB1,DRB1",
                        help="A comma-separated list of HLA genes (default: A,B,C,DQA1,DQB1,DRB1)")
    parser.add_argument('--partial',
                        dest='partial',
                        action='store_true',
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
                        default=1,
                        help="Reads simulated at every these base pairs (default: 1)")
    parser.add_argument("--coverage",
                        dest="coverage",
                        action='store_true',
                        help="Experimental purpose (assign reads based on coverage)")
    parser.add_argument("--best-alleles",
                        dest="best_alleles",
                        action='store_true',
                        help="")
    parser.add_argument("--exclude-allele-list",
                        dest="exclude_allele_list",
                        type=str,
                        default="",
                        help="A comma-separated list of alleles to be excluded. Enter a number N to randomly select N alleles for exclusion and N non-excluded alleles for testing (2N tested in total).")
    parser.add_argument("--random-seed",
                        dest="random_seed",
                        type=int,
                        default=0,
                        help="A seeding number for randomness (default: 0)")
    parser.add_argument("--num-mismatch",
                        dest="num_mismatch",
                        type=int,
                        default=0,
                        help="Maximum number of mismatches per read alignment to be considered (default: 0)")
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
    parser.add_argument("--assembly",
                        dest="assembly",
                        action="store_true",
                        help="Perform assembly")
    parser.add_argument("--no-concordant-assembly",
                        dest="concordant_assembly",
                        action="store_false",
                        help="")
    parser.add_argument("--exonic-only",
                        dest="exonic_only",
                        action="store_true",
                        help="Consider exonic regions only")
    parser.add_argument("--novel_allele_detection",
                        dest="novel_allele_detection",
                        action='store_true',
                        help="Change test to detection of new alleles. Report sensitivity and specificity rate at the end.")


    args = parser.parse_args()
    if not args.reference_type in ["gene", "chromosome", "genome"]:
        print >> sys.stderr, "Error: --reference-type (%s) must be one of gene, chromosome, and genome." % (args.reference_type)
        sys.exit(1)
    args.hla_list = args.hla_list.split(',')
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
    
    if len(args.default_allele_list) > 0:
        args.default_allele_list = args.default_allele_list.split(',')
        
    if len(args.exclude_allele_list) > 0:
        if args.exclude_allele_list.strip().isdigit():
            num_alleles = int(args.exclude_allele_list)            
            
            if not os.path.exists("Default-HLA/hla_backbone.fa"):
                curr_script = os.path.realpath(inspect.getsourcefile(test_HLA_genotyping))
                ex_path = os.path.dirname(curr_script)
                extract_hla_script = os.path.join(ex_path, "hisatgenotype_extract_vars.py")
                extract_cmd = [extract_hla_script,
                               "--reference-type", args.reference_type,
                               "--hla-list", ','.join(args.hla_list),
                               "--base", "Default-HLA/hla"]
                if args.partial:
                    extract_cmd += ["--partial"]
                extract_cmd += ["--inter-gap", "30",
                                "--intra-gap", "50"]
                if args.verbose_level >= 1:
                    print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
                proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
                proc.communicate()
                if not os.path.exists("Default-HLA/hla_backbone.fa"):
                    print >> sys.stderr, "Error: extract_HLA_vars (Default) failed!"
                    sys.exit(1)
       
            HLAs_default = {}
            #read_HLA_alleles("Default-HLA/hla_backbone.fa", HLAs_default)
            read_HLA_alleles("Default-HLA/hla_sequences.fa", HLAs_default)
            
            allele_names = list(HLAs_default['A'].keys())
            random.seed(args.random_seed)
            random.shuffle(allele_names)
            args.exclude_allele_list = allele_names[0:num_alleles]
            args.default_allele_list = allele_names[num_alleles:2*num_alleles]
            
            args.default_allele_list = args.default_allele_list + args.exclude_allele_list
            
            # DK - for debugging purposes
            args.default_allele_list = args.exclude_allele_list
        else:
            args.exclude_allele_list = args.exclude_allele_list.split(',')

        if args.num_mismatch == 0:
            args.num_mismatch = 3
        
    debug = {}
    if args.debug != "":
        for item in args.debug.split(','):
            if ':' in item:
                key, value = item.split(':')
                debug[key] = value
            else:
                debug[item] = 1

    random.seed(1)
    test_HLA_genotyping(args.base_fname,
                        args.reference_type,
                        args.hla_list,
                        args.partial,
                        args.aligners,
                        args.read_fname,
                        args.alignment_fname,
                        args.threads,
                        args.simulate_interval,
                        args.coverage,
                        args.best_alleles,
                        args.exclude_allele_list,
                        args.default_allele_list,
                        args.num_mismatch,
                        args.assembly,
                        args.concordant_assembly,
                        args.exonic_only,
                        args.verbose_level,
                        debug)
