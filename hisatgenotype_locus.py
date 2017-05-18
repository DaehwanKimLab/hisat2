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
from hisatgenotype_modules import typing_common, Gene_typing, assembly_graph


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
def Gene_prob_cmp(a, b):
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
def single_abundance(Gene_cmpt,
                     Gene_length):
    def normalize2(prob, length):
        total = 0
        for allele, mass in prob.items():
            assert allele in length
            total += (mass / length[allele])
        for allele, mass in prob.items():
            assert allele in length
            prob[allele] = mass / length[allele] / total

    Gene_prob, Gene_prob_next = {}, {}
    for cmpt, count in Gene_cmpt.items():
        alleles = cmpt.split('-')
        for allele in alleles:
            if allele not in Gene_prob:
                Gene_prob[allele] = 0.0
            Gene_prob[allele] += (float(count) / len(alleles))

    normalize(Gene_prob)
    def next_prob(Gene_cmpt, Gene_prob, Gene_length):
        Gene_prob_next = {}
        for cmpt, count in Gene_cmpt.items():
            alleles = cmpt.split('-')
            alleles_prob = 0.0
            for allele in alleles:
                if allele not in Gene_prob:
                    continue
                alleles_prob += Gene_prob[allele]
            if alleles_prob <= 0.0:
                continue
            for allele in alleles:
                if allele not in Gene_prob:
                    continue
                if allele not in Gene_prob_next:
                    Gene_prob_next[allele] = 0.0
                Gene_prob_next[allele] += (float(count) * Gene_prob[allele] / alleles_prob)
        normalize(Gene_prob_next)
        return Gene_prob_next


    fast_EM = True
    diff, iter = 1.0, 0
    while diff > 0.0001 and iter < 1000:
        Gene_prob_next = next_prob(Gene_cmpt, Gene_prob, Gene_length)
        if fast_EM:
            # Accelerated version of EM - SQUAREM iteration
            #    Varadhan, R. & Roland, C. Scand. J. Stat. 35, 335-353 (2008)
            #    Also, this algorithm is used in Sailfish - http://www.nature.com/nbt/journal/v32/n5/full/nbt.2862.html
            Gene_prob_next2 = next_prob(Gene_cmpt, Gene_prob_next, Gene_length)
            sum_squared_r, sum_squared_v = 0.0, 0.0
            p_r, p_v = {}, {}
            for a in Gene_prob.keys():
                p_r[a] = Gene_prob_next[a] - Gene_prob[a]
                sum_squared_r += (p_r[a] * p_r[a])
                p_v[a] = Gene_prob_next2[a] - Gene_prob_next[a] - p_r[a]
                sum_squared_v += (p_v[a] * p_v[a])
            if sum_squared_v > 0.0:
                gamma = -math.sqrt(sum_squared_r / sum_squared_v)
                for a in Gene_prob.keys():
                    Gene_prob_next2[a] = max(0.0, Gene_prob[a] - 2 * gamma * p_r[a] + gamma * gamma * p_v[a]);
                Gene_prob_next = next_prob(Gene_cmpt, Gene_prob_next2, Gene_length)

        diff = prob_diff(Gene_prob, Gene_prob_next)
        Gene_prob = Gene_prob_next

        # Accelerate convergence
        if iter >= 10:
            Gene_prob2 = {}
            avg_prob = sum(Gene_prob.values()) / len(Gene_prob)
            for allele, prob in Gene_prob.items():
                if prob >= 0.005 or prob > avg_prob:
                    Gene_prob2[allele] = prob
                Gene_prob = Gene_prob2

        # DK - debugging purposes
        if iter % 10 == 0 and False:
            print "iter", iter
            for allele, prob in Gene_prob.items():
                if prob >= 0.01:
                    print >> sys.stderr, "\t", iter, allele, prob, str(datetime.now())
        
        iter += 1
        
    """
    for allele, prob in Gene_prob.items():
        allele_len = Gene_length[allele]
        Gene_prob[allele] /= float(allele_len)
    """
    
    # normalize(Gene_prob)
    normalize2(Gene_prob, Gene_length)
    Gene_prob = [[allele, prob] for allele, prob in Gene_prob.items()]
    Gene_prob = sorted(Gene_prob, cmp=Gene_prob_cmp)
    return Gene_prob


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
                if count >= num_nt * 0.2 or count >= 7:
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
        new_nt_dic = {}
        for nt, count in nt_dic.items():
            var_id = ""
            if nt == 'D':
                if i <= skip_i:
                    assert prev_del_var_id != ""
                    var_id = prev_del_var_id
                else:
                    for var_id_, var_type, var_data in var_list[i]:
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
                            var_id = var_id_
                            prev_del_var_id = var_id
                            skip_i = i + del_len - 1
                            break                                                
            elif nt != 'N' and nt != ref_nt:
                assert nt in "ACGT"
                id = "unknown"
                for var_id_, var_type, var_data in var_list[i]:
                    if var_type != "single":
                        continue
                    if nt == var_data:
                        var_id = var_id_
                        break
            new_nt_dic[nt] = [count, var_id]
                        
        mpileup[i][1] = new_nt_dic

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
                        var_idx = typing_common.lower_bound(Var_list, left + j)
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
                    var_idx = typing_common.lower_bound(Var_list, left)
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
           base_fname,
           locus_list,
           partial,
           partial_alleles,
           refGenes,
           Genes,
           Gene_names,
           Gene_lengths,
           refGene_loci,
           Vars,
           Var_list,
           Links,
           aligners,
           num_editdist,
           assembly,
           output_base,
           error_correction,
           allow_discordant,
           display_alleles,
           fastq,
           read_fname,
           alignment_fname,
           num_frag_list,
           read_len,
           fragment_len,
           threads,
           best_alleles,
           verbose):
    if simulation:
        test_passed = {}
    report_file = open(output_base + ".report", 'w')
    for aligner, index_type in aligners:
        for f_ in [sys.stderr, report_file]:
            if index_type == "graph":
                print >> f_, "\n\t\t%s %s" % (aligner, index_type)
            else:
                print >> f_, "\n\t\t%s %s" % (aligner, index_type)

        remove_alignment_file = False
        if alignment_fname == "":
            # Align reads, and sort the alignments into a BAM file
            remove_alignment_file = True
            if simulation:
                alignment_fname = "%s_output.bam" % base_fname
            else:
                alignment_fname = read_fname[0].split('/')[-1]
                alignment_fname = alignment_fname.split('.')[0] + ".bam"
                
            typing_common.align_reads(aligner,
                                      simulation,
                                      base_fname,
                                      index_type,
                                      read_fname,
                                      fastq,
                                      threads,
                                      alignment_fname,
                                      verbose)
            
        for test_Gene_names in locus_list:
            if simulation:
                gene = test_Gene_names[0].split('*')[0]
            else:
                gene = test_Gene_names
            ref_allele = refGenes[gene]
            ref_seq = Genes[gene][ref_allele]
            ref_exons = refGene_loci[gene][-1]
            
            novel_var_count = 0        
            gene_vars, gene_var_list = deepcopy(Vars[gene]), deepcopy(Var_list[gene])
            gene_del_var_list = []
            for var_id, var_data in gene_vars.items():
                var_type, var_pos, var_data = var_data
                if var_type != "deletion":
                    continue
                gene_del_var_list.append([var_pos, var_pos + int(var_data) - 1, var_id])
            var_count = {}
            def add_novel_var(gene_vars,
                              gene_var_list,
                              novel_var_count,
                              var_type,
                              var_pos,
                              var_data):
                var_idx = typing_common.lower_bound(gene_var_list, var_pos)
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
                alignview_cmd += [ref_allele]
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
            for _, var_id in gene_var_list:
                allele_list = Links[var_id]
                for allele_id in allele_list:
                    if allele_id not in Genes[gene]:
                        continue
                    if allele_id not in allele_vars:
                        allele_vars[allele_id] = [var_id]
                    else:
                        allele_vars[allele_id].append(var_id)

            # Extract variants that are within exons
            exon_vars = get_exonic_vars(gene_vars, ref_exons)

            # Store nodes that represent alleles
            allele_nodes = {}
            def create_allele_node(allele_name):
                if allele_name in allele_nodes:
                    return allele_nodes[allele_name]
                if allele_name in allele_vars:
                    var_ids = allele_vars[allele_name]
                else:
                    var_ids = []
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
                        del_len = int(var_data)
                        assert var_pos + del_len <= len(ref_seq)
                        seq[var_pos:var_pos + del_len] = ['D'] * del_len
                        var[var_pos:var_pos + del_len] = [var_id] * del_len
                    else:
                        # DK - to be implemented for insertions
                        assert var_type == "insertion"

                qual = ' ' * len(seq)
                allele_node = assembly_graph.Node(allele_name,
                                                  0,
                                                  seq,
                                                  qual,
                                                  var,
                                                  ref_seq,
                                                  gene_vars,
                                                  mpileup,
                                                  simulation)
                allele_nodes[allele_name] = allele_node
                return allele_node

            true_allele_nodes = {}
            if simulation:
                for allele_name in test_Gene_names:
                    true_allele_nodes[allele_name] = create_allele_node(allele_name)

            display_allele_nodes = {}
            for display_allele in display_alleles:
                display_allele_nodes[display_allele] = create_allele_node(display_allele)

            # Assembly graph
            asm_graph = assembly_graph.Graph(ref_seq,
                                             gene_vars,
                                             ref_exons,
                                             partial_alleles,
                                             true_allele_nodes,
                                             {}, # predicted_allele_nodes, which is empty for now
                                             display_allele_nodes,
                                             simulation)

            # Choose allele representives from those that share the same exonic sequences
            allele_reps, allele_rep_groups = get_rep_alleles(Links, exon_vars)
            allele_rep_set = set(allele_reps.values())

            # For checking alternative alignments near the ends of alignments
            Alts_left, Alts_right = typing_common.get_alternatives(ref_seq,
                                                                   allele_vars,
                                                                   gene_vars,
                                                                   gene_var_list,
                                                                   verbose)

            def haplotype_alts_list(haplotype_alts, left = True):
                haplotype_list = []
                for haplotype in haplotype_alts.keys():
                    if left:
                        pos = int(haplotype.split('-')[-1])
                    else:
                        pos = int(haplotype.split('-')[0])
                    haplotype_list.append([pos, haplotype])
                return sorted(haplotype_list, cmp = lambda a, b: a[0] - b[0])

            Alts_left_list, Alts_right_list = haplotype_alts_list(Alts_left, True), haplotype_alts_list(Alts_right, False)

            # Count alleles
            Gene_counts, Gene_cmpt = {}, {}
            Gene_gen_counts, Gene_gen_cmpt = {}, {}
            num_reads, total_read_len = 0, 0

            # For debugging purposes
            debug_allele_names = set(test_Gene_names) if simulation and verbose >= 2 else set()

            # Read information
            prev_read_id = None
            prev_right_pos = 0
            prev_lines = []
            left_read_ids, right_read_ids = set(), set()
            if index_type == "graph":
                # nodes for reads
                read_nodes = []
                read_vars_list = []

                # 
                def add_count(count_per_read, ht, add):
                    orig_ht = ht
                    ht = ht.split('-')

                    assert len(ht) >= 2
                    left, right = int(ht[0]), int(ht[-1])
                    assert left <= right

                    ht = ht[1:-1]
                    alleles = set(Genes[gene].keys()) - set([ref_allele])
                    for i in range(len(ht)):
                        var_id = ht[i]
                        if var_id.startswith("nv"):
                            continue
                        alleles &= set(Links[var_id])
                    ht = set(ht)

                    tmp_alleles = set()
                    var_idx = typing_common.lower_bound(gene_var_list, left)
                    for var_idx in range(var_idx, len(gene_var_list)):
                        var_pos, var_id = gene_var_list[var_idx]
                        if var_pos > right:
                            break
                        if var_pos < left or var_id in ht:
                            continue
                        if var_id in ht:
                            continue
                        if var_id.startswith("nv"):
                            continue
                        tmp_alleles |= set(Links[var_id])

                    for var_idx in range(var_idx, len(gene_del_var_list)):
                        var_left, var_right, var_id = gene_del_var_list[var_idx]
                        if var_id in ht:
                            continue
                        if var_right < left or var_left > right:
                            continue
                        tmp_alleles |= set(Links[var_id])
                    alleles -= tmp_alleles
                    
                    for allele in alleles:
                        count_per_read[allele] += add

                    return len(alleles)

                # Identify best pairs
                def choose_pairs(left_positive_hts, right_positive_hts):
                    if len(left_positive_hts) > 0 and \
                       len(right_positive_hts) > 0 and \
                       max(len(left_positive_hts), len(right_positive_hts)) >= 2:
                        expected_inter_dist = fragment_len - read_len * 2
                        best_diff = sys.maxint
                        picked = []                                
                        for left_ht_str in left_positive_hts:
                            left_ht = left_ht_str.split('-')
                            l_left, l_right = int(left_ht[0]), int(left_ht[-1])
                            for right_ht_str in right_positive_hts:
                                right_ht = right_ht_str.split('-')
                                r_left, r_right = int(right_ht[0]), int(right_ht[-1])
                                if l_right < r_right:
                                    inter_dist = r_left - l_right - 1
                                else:
                                    inter_dist = l_left - r_right - 1

                                cur_diff = abs(expected_inter_dist - inter_dist)
                                if best_diff > cur_diff:
                                    best_diff = cur_diff
                                    picked = [[left_ht_str, right_ht_str]]
                                elif best_diff == cur_diff:
                                    picked.append([left_ht_str, right_ht_str])

                        assert len(picked) > 0

                        left_positive_hts, right_positive_hts = set(), set()
                        for left_ht_str, right_ht_str in picked:
                            left_positive_hts.add(left_ht_str)
                            right_positive_hts.add(right_ht_str)

                    return left_positive_hts, right_positive_hts                
                
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
                    is_left_read = flag & 0x40 != 0
                    if is_left_read:
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
                                    var_idx = typing_common.lower_bound(gene_var_list, ref_pos)
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
                                var_idx = typing_common.lower_bound(gene_var_list, right_pos)
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
                            if read_pos == Zs_pos and \
                               Zs_i < len(Zs) and \
                               Zs[Zs_i][1] == 'D':
                                _var_id = Zs[Zs_i][2]
                                Zs_i += 1
                                if Zs_i < len(Zs):
                                    Zs_pos += Zs[Zs_i][0]
                            else:
                                # Search for a known (yet not indexed) variant or a novel variant
                                var_idx = typing_common.lower_bound(gene_var_list, right_pos)
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
                            for nt, value in mpileup[right_pos][1].items():
                                count = value[0]
                                if nt == 'D':
                                    del_count += count
                                else:
                                    nt_count += count
                                    
                            # DK - debugging purposes
                            if base_fname == "hla":
                                if del_count * 6 < nt_count: # and nt_count >= 15:
                                    likely_misalignment = True
                            
                        elif cigar_op == 'S':
                            if i == 0:
                                softclip[0] = length
                                Zs_pos += length
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

                    def add_stat(Gene_cmpt, Gene_counts, Gene_count_per_read, include_alleles = set()):
                        max_count = max(Gene_count_per_read.values())
                        cur_cmpt = set()
                        for allele, count in Gene_count_per_read.items():
                            if count < max_count:
                                continue

                            if len(include_alleles) > 0 and allele not in include_alleles:
                                continue
                            
                            cur_cmpt.add(allele)                    
                            if allele not in Gene_counts:
                                Gene_counts[allele] = 1
                            else:
                                Gene_counts[allele] += 1

                        if len(cur_cmpt) == 0:
                            return ""

                        # DK - for debugging purposes                            
                        alleles = ["", ""]
                        # alleles = ["A*24:36N", "A*24:359N"]
                        allele1_found, allele2_found = False, False
                        if alleles[0] != "":
                            for allele, count in Gene_count_per_read.items():
                                if count < max_count:
                                    continue
                                if allele == alleles[0]:
                                    allele1_found = True
                                elif allele == alleles[1]:
                                    allele2_found = True
                            if allele1_found != allele2_found:
                                print alleles[0], Gene_count_per_read[alleles[0]]
                                print alleles[1], Gene_count_per_read[alleles[1]]
                                if allele1_found:
                                    print ("%s\tread_id %s - %d vs. %d]" % (alleles[0], prev_read_id, max_count, Gene_count_per_read[alleles[1]]))
                                else:
                                    print ("%s\tread_id %s - %d vs. %d]" % (alleles[1], prev_read_id, max_count, Gene_count_per_read[alleles[0]]))
                                print read_seq

                        cur_cmpt = sorted(list(cur_cmpt))
                        cur_cmpt = '-'.join(cur_cmpt)
                        if not cur_cmpt in Gene_cmpt:
                            Gene_cmpt[cur_cmpt] = 1
                        else:
                            Gene_cmpt[cur_cmpt] += 1

                        return cur_cmpt

                    if read_id != prev_read_id:
                        if prev_read_id != None:
                            if base_fname == "codis" and gene == "D18S51":
                                left_positive_hts, right_positive_hts = choose_pairs(left_positive_hts, right_positive_hts)

                            for positive_ht in left_positive_hts | right_positive_hts:
                                if positive_ht == "unknown" or positive_ht.startswith("nv"):
                                    continue
                                if set(positive_ht.split('-')) & exon_vars:
                                    add_count(Gene_count_per_read, positive_ht, 1)
                                add_count(Gene_gen_count_per_read, positive_ht, 1)

                            if base_fname == "hla":
                                cur_cmpt = add_stat(Gene_cmpt, Gene_counts, Gene_count_per_read, allele_rep_set)
                                cur_cmpt_gen = add_stat(Gene_gen_cmpt, Gene_gen_counts, Gene_gen_count_per_read)
                            else:
                                cur_cmpt = add_stat(Gene_gen_cmpt, Gene_gen_counts, Gene_gen_count_per_read)
                            for read_id_, read_node in read_nodes:
                                asm_graph.add_node(read_id_,
                                                   read_node,
                                                   simulation)
                            read_nodes, read_var_list = [], []
                            if verbose >= 2 and base_fname in ["hla", "codis"]:
                                cur_cmpt = cur_cmpt.split('-') if cur_cmpt != "" else set()
                                cur_cmpt_gen = cur_cmpt_gen.split('-') if cur_cmpt_gen != "" else set()
                                show_debug = (partial and cur_cmpt != "" and not set(cur_cmpt) & set(test_Gene_names)) or \
                                              (not partial and cur_cmpt_gen != "" and not set(cur_cmpt_gen) & set(test_Gene_names))
                                              
                                if show_debug:
                                    print "%s are chosen instead of %s" % (cur_cmpt if partial else cur_cmpt_gen, '-'.join(test_Gene_names))

                                    # DK - debugging purpose
                                    print "DK:", cur_cmpt
                                    print "DK:", cur_cmpt_gen
                                    
                                    for prev_line in prev_lines:
                                        print "\t", prev_line

                            prev_lines = []

                        Gene_count_per_read, Gene_gen_count_per_read = {}, {}
                        for Gene_name in Gene_names[gene]:
                            if Gene_name.find("BACKBONE") != -1:
                                continue
                            Gene_count_per_read[Gene_name] = 0
                            Gene_gen_count_per_read[Gene_name] = 0

                    prev_lines.append(line)

                    # Node
                    read_node_pos, read_node_seq, read_node_qual, read_node_var = -1, [], [], []
                    read_vars = []

                    # Positive evidence
                    if prev_read_id == None:
                        left_positive_hts, right_positive_hts = set(), set()
                    elif is_left_read:
                        left_positive_hts = set()
                    else:
                        right_positive_hts = set()

                    # Sanity check - read length, cigar string, and MD string
                    ref_pos, read_pos, cmp_cigar_str, cmp_MD = left_pos, 0, "", ""
                    cigar_match_len, MD_match_len = 0, 0

                    cmp_list_left, cmp_list_right, cmp_left_alts, cmp_right_alts = \
                    typing_common.identify_ambigious_diffs(ref_seq,
                                                           gene_vars,
                                                           Alts_left,
                                                           Alts_right,
                                                           Alts_left_list,
                                                           Alts_right_list,
                                                           cmp_list,
                                                           verbose,
                                                           orig_read_id.startswith("a59|L"))  # debug?

                    mid_ht = []
                    for cmp in cmp_list[cmp_list_left:cmp_list_right+1]:
                        type = cmp[0]
                        if type not in ["mismatch", "deletion", "insertion"]:
                            continue                            
                        var_id = cmp[3]
                        if var_id == "unknown" or var_id.startswith("nv"):
                            continue
                        mid_ht.append(var_id)

                    for l in range(len(cmp_left_alts)):
                        left_ht = cmp_left_alts[l].split('-')
                        left_ht += mid_ht
                        for r in range(len(cmp_right_alts)):
                            right_ht = cmp_right_alts[r].split('-')
                            ht = left_ht + right_ht
                            if len(ht) <= 0:
                                continue
                            ht_str = '-'.join(ht)
                            if is_left_read:
                                left_positive_hts.add(ht_str)
                            else:
                                right_positive_hts.add(ht_str)

                    # DK - debugging purposes
                    DK_debug = False
                    if orig_read_id.startswith("a59|L"):
                        DK_debug = True
                        print line
                        print cmp_list
                        print "positive hts:", left_positive_hts, right_positive_hts
                        print "cmp_list[%d, %d]" % (cmp_list_left, cmp_list_right)

                    cmp_i = 0
                    while cmp_i < len(cmp_list):
                        cmp = cmp_list[cmp_i]
                        type, length = cmp[0], cmp[2]
                        if type in ["match", "mismatch"]:
                            if read_node_pos < 0:
                                read_node_pos = ref_pos

                        if type == "match":
                            read_node_seq += list(read_seq[read_pos:read_pos+length])
                            read_node_qual += list(read_qual[read_pos:read_pos+length])
                            read_node_var += ([''] * length)
                            read_pos += length
                        elif type == "mismatch":
                            var_id = cmp[3]
                            read_base, qual = read_seq[read_pos], read_qual[read_pos]
                            read_node_seq += [read_base]
                            read_node_qual += [qual]
                            read_node_var.append(var_id)
                            read_pos += 1
                        elif type == "insertion":
                            var_id = cmp[3]
                            ins_len = length
                            ins_seq = read_seq[read_pos:read_pos+ins_len]
                            read_node_seq += ["I%s" % nt for nt in ins_seq]
                            read_node_qual += list(read_qual[read_pos:read_pos+ins_len])
                            read_node_var += ([var_id] * ins_len)                                        
                            read_pos += length
                        elif type == "deletion":
                            var_id = cmp[3]
                            alt_match = False
                            del_len = length
                            read_node_seq += (['D'] * del_len)
                            read_node_qual += ([''] * del_len)
                            if len(read_node_seq) > len(read_node_var):
                                assert len(read_node_seq) == len(read_node_var) + del_len
                                read_node_var += ([var_id] * del_len)
                        else:
                            assert type == "intron"

                        cmp_i += 1

                    # Node
                    if assembly:
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

                    # DK - debugging purposes
                    if DK_debug:
                        print "count_per_read:", Gene_gen_count_per_read

                    prev_read_id = read_id
                    prev_right_pos = right_pos

                if num_reads <= 0:
                    continue

                for f_ in [sys.stderr, report_file]:
                    print >> f_, "\t\t\tNumber of reads aligned: %d" % num_reads

                if prev_read_id != None:
                    if base_fname == "codis" and gene == "D18S51":
                        left_positive_hts, right_positive_hts = choose_pairs(left_positive_hts, right_positive_hts)                            
                    for positive_ht in left_positive_hts | right_positive_hts:
                        if positive_ht == "unknown" or positive_ht.startswith("nv"):
                            continue
                        if set(positive_ht.split('-')) & exon_vars:
                            add_count(Gene_count_per_read, positive_ht, 1)
                        add_count(Gene_gen_count_per_read, positive_ht, 1)

                    if base_fname == "hla":
                        add_stat(Gene_cmpt, Gene_counts, Gene_count_per_read, allele_rep_set)
                    add_stat(Gene_gen_cmpt, Gene_gen_counts, Gene_gen_count_per_read)
                    for read_id_, read_node in read_nodes:
                        asm_graph.add_node(read_id_,
                                           read_node,
                                           simulation)
                    read_nodes, read_var_list = [], []
                
            else:
                assert index_type == "linear"
                def add_alleles(alleles):
                    if not allele in Gene_counts:
                        Gene_counts[allele] = 1
                    else:
                        Gene_counts[allele] += 1

                    cur_cmpt = sorted(list(alleles))
                    cur_cmpt = '-'.join(cur_cmpt)
                    if not cur_cmpt in Gene_cmpt:
                        Gene_cmpt[cur_cmpt] = 1
                    else:
                        Gene_cmpt[cur_cmpt] += 1

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

            # DK - debugging purposes
            # if base_fname != "hla":
            if base_fname != "hla" or True:
                Gene_counts = Gene_gen_counts
            Gene_counts = [[allele, count] for allele, count in Gene_counts.items()]
            def Gene_count_cmp(a, b):
                if a[1] != b[1]:
                    return b[1] - a[1]
                assert a[0] != b[0]
                if a[0] < b[0]:
                    return -1
                else:
                    return 1
            Gene_counts = sorted(Gene_counts, cmp=Gene_count_cmp)
            for count_i in range(len(Gene_counts)):
                count = Gene_counts[count_i]
                if simulation:
                    found = False
                    for test_Gene_name in test_Gene_names:
                        if count[0] == test_Gene_name:
                            for f_ in [sys.stderr, report_file]:
                                print >> f_, "\t\t\t*** %d ranked %s (count: %d)" % (count_i + 1, test_Gene_name, count[1])
                            found = True
                            """
                            if count_i > 0 and Gene_counts[0][1] > count[1]:
                                print >> sys.stderr, "Warning: %s ranked first (count: %d)" % (Gene_counts[0][0], Gene_counts[0][1])
                                assert False
                            else:
                                test_passed += 1
                            """
                    if count_i < 5 and not found:
                        for f_ in [sys.stderr, report_file]:
                            print >> f_, "\t\t\t\t%d %s (count: %d)" % (count_i + 1, count[0], count[1])
                else:
                    for f_ in [sys.stderr, report_file]:
                        print >> f_, "\t\t\t\t%d %s (count: %d)" % (count_i + 1, count[0], count[1])
                    if count_i >= 9:
                        break
            for f_ in [sys.stderr, report_file]:
                print >> f_

            # Calculate the abundance of representative alleles on exonic sequences
            if base_fname == "hla":
                # Incorporate non representative alleles (full length alleles)
                Gene_prob = single_abundance(Gene_cmpt, Gene_lengths[gene])
                gen_alleles = set()
                gen_prob_sum = 0.0
                for prob_i in range(len(Gene_prob)):
                    allele, prob = Gene_prob[prob_i][:2]
                    if prob_i >= 10 and prob < 0.03:
                        break
                    if allele in partial_alleles:
                        continue

                    gen_prob_sum += prob
                    for allele2 in allele_rep_groups[allele]:
                        gen_alleles.add(allele2)

                if len(gen_alleles) > 0:
                    Gene_gen_cmpt2 = {}
                    for cmpt, value in Gene_gen_cmpt.items():
                        cmpt2 = []
                        for allele in cmpt.split('-'):
                            if allele in gen_alleles:
                                cmpt2.append(allele)
                        if len(cmpt2) == 0:
                            continue
                        cmpt2 = '-'.join(cmpt2)
                        if cmpt2 not in Gene_gen_cmpt2:
                            Gene_gen_cmpt2[cmpt2] = value
                        else:
                            Gene_gen_cmpt2[cmpt2] += value
                    Gene_gen_cmpt = Gene_gen_cmpt2
                    Gene_gen_prob = single_abundance(Gene_gen_cmpt, Gene_lengths[gene])

                    Gene_combined_prob = {}
                    for allele, prob in Gene_prob:
                        assert allele not in Gene_combined_prob
                        if allele in gen_alleles:
                            Gene_combined_prob[allele] = 0.0
                        else:
                            Gene_combined_prob[allele] = prob
                    for allele, prob in Gene_gen_prob:
                        Gene_combined_prob[allele] = prob * gen_prob_sum
                    Gene_prob = [[allele, prob] for allele, prob in Gene_combined_prob.items()]
                    Gene_prob = sorted(Gene_prob, cmp=Gene_prob_cmp)
            else:
                Gene_prob = single_abundance(Gene_gen_cmpt, Gene_lengths[gene])

            if index_type == "graph" and assembly:
                allele_node_order = []
                predicted_allele_nodes = {}
                for allele_name, prob in Gene_prob:
                    if prob < 0.1: # abundance of 10%
                        break
                    predicted_allele_nodes[allele_name] = create_allele_node(allele_name)
                    allele_node_order.append([allele_name, prob])
                    if len(predicted_allele_nodes) >= 2:
                        break
                asm_graph.predicted_allele_nodes = predicted_allele_nodes
                asm_graph.allele_node_order = allele_node_order

                # Start drawing assembly graph
                asm_graph.begin_draw(output_base)

                # Draw assembly graph
                begin_y = asm_graph.draw(0, "Initial graph")
                begin_y += 200
                
                # Apply De Bruijn graph
                asm_graph.guided_DeBruijn()

                # Draw assembly graph
                begin_y = asm_graph.draw(begin_y, "Asssembly")
                begin_y += 200

                # Draw assembly graph
                asm_graph.nodes = asm_graph.nodes2
                asm_graph.to_node, asm_graph.from_node = {}, {}
                begin_y = asm_graph.draw(begin_y, "Assembly with known alleles")

                # End drawing assembly graph
                asm_graph.end_draw()

                # Compare two alleles
                if simulation and len(test_Gene_names) == 2:
                    allele_name1, allele_name2 = test_Gene_names
                    print >> sys.stderr, allele_name1, "vs.", allele_name2
                    asm_graph.print_node_comparison(asm_graph.true_allele_nodes)

                def compare_alleles(vars1, vars2, print_output = True):
                    skip = True
                    var_i, var_j = 0, 0
                    exon_i = 0
                    mismatches = 0
                    while var_i < len(vars1) and var_j < len(vars2):
                        cmp_var_id, node_var_id = vars1[var_i], vars2[var_j]
                        cmp_var, node_var = gene_vars[cmp_var_id], gene_vars[node_var_id]

                        min_pos = min(cmp_var[1], node_var[1])
                        cmp_var_in_exon, node_var_in_exon = False, False
                        while exon_i < len(ref_exons):
                            exon_left, exon_right = ref_exons[exon_i]
                            if min_pos <= exon_right:
                                if cmp_var[1] >= exon_left and cmp_var[1] <= exon_right:
                                    cmp_var_in_exon = True
                                else:
                                    cmp_var_in_exon = False
                                if node_var[1] >= exon_left and node_var[1] <= exon_right:
                                    node_var_in_exon = True
                                else:
                                    node_var_in_exon = False                                
                                break
                            exon_i += 1
                        
                        if cmp_var_id == node_var_id:
                            skip = False
                            if print_output:
                                if cmp_var_in_exon:
                                    print >> sys.stderr, "\033[94mexon%d\033[00m" % (exon_i + 1),
                                print >> sys.stderr, cmp_var_id, cmp_var, "\t\t\t", mpileup[cmp_var[1]]
                            var_i += 1; var_j += 1
                            continue
                        if cmp_var[1] <= node_var[1]:
                            if not skip:
                                if (var_i > 0 and var_i + 1 < len(vars1)) or cmp_var[0] != "deletion":
                                    if print_output:
                                        if cmp_var_in_exon:
                                            for f_ in [sys.stderr, report_file]:
                                                print >> f_, "\033[94mexon%d\033[00m" % (exon_i + 1),
                                        for f_ in [sys.stderr, report_file]:
                                            print >> f_, "***", cmp_var_id, cmp_var, "==", "\t\t\t", mpileup[cmp_var[1]]
                                    mismatches += 1
                            var_i += 1
                        else:
                            if print_output:
                                if node_var_in_exon:
                                    for f_ in [sys.stderr, report_file]:
                                        print >> f_, "\033[94mexon%d\033[00m" % (exon_i + 1),
                                for f_ in [sys.stderr, report_file]:
                                    print >> f_, "*** ==", node_var_id, node_var, "\t\t\t", mpileup[node_var[1]]
                            mismatches += 1
                            var_j += 1
                            
                    return mismatches
                    
                tmp_nodes = asm_graph.nodes
                print >> sys.stderr, "Number of tmp nodes:", len(tmp_nodes)
                count = 0
                for id, node in tmp_nodes.items():
                    count += 1
                    if count > 10:
                        break
                    node_vars = node.get_var_ids()
                    node.print_info(); print >> sys.stderr
                    if node.id in asm_graph.to_node:
                        for id2, at in asm_graph.to_node[node.id]:
                            print >> sys.stderr, "\tat %d ==> %s" % (at, id2)

                    if simulation:
                        cmp_Gene_names = test_Gene_names
                    else:
                        cmp_Gene_names = [allele_name for allele_name, _ in allele_node_order]
                        
                    alleles, cmp_vars, max_common = [], [], -sys.maxint
                    for cmp_Gene_name in cmp_Gene_names:
                        tmp_vars = allele_nodes[cmp_Gene_name].get_var_ids(node.left, node.right)
                        tmp_common = len(set(node_vars) & set(tmp_vars))
                        tmp_common -= len(set(node_vars) | set(tmp_vars))
                        if max_common < tmp_common:
                            max_common = tmp_common
                            alleles = [[cmp_Gene_name, tmp_vars]]
                        elif max_common == tmp_common:
                            alleles.append([cmp_Gene_name, tmp_vars])

                    for allele_name, cmp_vars in alleles:
                        for f_ in [sys.stderr, report_file]:
                            print >> f_, "vs.", allele_name
                        compare_alleles(cmp_vars, node_vars)

                    print >> sys.stderr
                    print >> sys.stderr


            # Identify alleles that perfectly or closesly match assembled alleles
            for node_name, node in asm_graph.nodes.items():
                vars = set(node.get_var_ids())

                max_allele_names, max_common = [], -sys.maxint
                for allele_name, vars2 in allele_vars.items():
                    vars2 = set(vars2)
                    tmp_common = len(vars & vars2) - len(vars | vars2)
                    if tmp_common > max_common:
                        max_common = tmp_common
                        max_allele_names = [allele_name]                        
                    elif tmp_common == max_common:
                        max_allele_names.append(allele_name)

                for f_ in [sys.stderr, report_file]:
                    print >> f_, "Genomic:", node_name
                    node_vars = node.get_var_ids()
                    min_mismatches = sys.maxint
                    for max_allele_name in max_allele_names:
                        cmp_vars = allele_vars[max_allele_name]
                        cmp_vars = sorted(cmp_vars, cmp=lambda a, b: int(a[2:]) - int(b[2:]))
                        print_output = False
                        tmp_mismatches = compare_alleles(cmp_vars, node_vars, print_output)
                        print >> f_, "\t\t%s:" % max_allele_name, max_common, tmp_mismatches
                        if tmp_mismatches < min_mismatches:
                            min_mismatches = tmp_mismatches
                    if min_mismatches > 0:
                        print >> f_, "Novel allele"
                    else:
                        print >> f_, "Known allele"

            """
            allele_exon_vars = {}
            for allele_name, vars in allele_vars.items():
                allele_exon_vars[allele_name] = set(vars) & exon_vars

            for node_name, node in asm_graph.nodes.items():
                vars = []
                for left, right in ref_exons:
                    vars += node.get_var_ids(left, right)
                vars = set(vars) & exon_vars

                max_allele_names, max_common = [], -sys.maxint
                for allele_name, vars2 in allele_exon_vars.items():
                    tmp_common = len(vars & vars2) - len(vars | vars2)
                    if tmp_common > max_common:
                        max_common = tmp_common
                        max_allele_names = [allele_name]                        
                    elif tmp_common == max_common:
                        max_allele_names.append(allele_name)

                for f_ in [sys.stderr, report_file]:
                    print >> f_, "Exonic:", node_name
                    for max_allele_name in max_allele_names:
                        print >> f_, "\t\t%s:" % max_allele_name, max_common
            """
            
            success = [False for i in range(len(test_Gene_names))]
            found_list = [False for i in range(len(test_Gene_names))]
            for prob_i in range(len(Gene_prob)):
                prob = Gene_prob[prob_i]
                found = False
                _allele_rep = prob[0]
                """
                if partial and exonic_only:
                    _fields = _allele_rep.split(':')
                    if len(_fields) == 4:
                        _allele_rep = ':'.join(_fields[:-1])
                """

                if simulation:
                    for name_i in range(len(test_Gene_names)):
                        test_Gene_name = test_Gene_names[name_i]
                        if prob[0] == test_Gene_name:
                            rank_i = prob_i
                            while rank_i > 0:
                                if prob == Gene_prob[rank_i - 1][1]:
                                    rank_i -= 1
                                else:
                                    break
                            for f_ in [sys.stderr, report_file]:
                                print >> f_, "\t\t\t*** %d ranked %s (abundance: %.2f%%)" % (rank_i + 1, test_Gene_name, prob[1] * 100.0)
                            if rank_i < len(success):
                                success[rank_i] = True
                            found_list[name_i] = True
                            found = True
                    # DK - for debugging purposes
                    if not False in found_list and prob_i >= 10:
                        break
                if not found:
                    for f_ in [sys.stderr, report_file]:
                        print >> f_, "\t\t\t\t%d ranked %s (abundance: %.2f%%)" % (prob_i + 1, _allele_rep, prob[1] * 100.0)
                    if best_alleles and prob_i < 2:
                        for f_ in [sys.stderr, report_file]:
                            print >> f_, "SingleModel %s (abundance: %.2f%%)" % (_allele_rep, prob[1] * 100.0)
                if not simulation and prob_i >= 9:
                    break
                if prob_i >= 19:
                    break
            print >> sys.stderr

            if simulation and not False in success:
                aligner_type = "%s %s" % (aligner, index_type)
                if not aligner_type in test_passed:
                    test_passed[aligner_type] = 1
                else:
                    test_passed[aligner_type] += 1

        if remove_alignment_file and not simulation:
            os.system("rm %s*" % (alignment_fname))

    report_file.close()
    if simulation:
        return test_passed

    
"""
"""
def read_Gene_alleles(fname, Genes):
    for line in open(fname):
        if line.startswith(">"):
            Gene_name = line.strip().split()[0][1:]
            Gene_gene = Gene_name.split('*')[0]
            if not Gene_gene in Genes:
                Genes[Gene_gene] = {}
            if not Gene_name in Genes[Gene_gene]:
                Genes[Gene_gene][Gene_name] = ""
        else:
            Genes[Gene_gene][Gene_name] += line.strip()
    return Genes


"""
"""
def read_Gene_vars(fname):
    Vars, Var_list = {}, {}
    for line in open(fname):
        var_id, var_type, allele, pos, data = line.strip().split('\t')
        pos = int(pos)
        gene = allele.split('*')[0]
        if not gene in Vars:
            Vars[gene] = {}
            assert not gene in Var_list
            Var_list[gene] = []
            
        assert not var_id in Vars[gene]
        left = 0
        Vars[gene][var_id] = [var_type, pos - left, data]
        Var_list[gene].append([pos - left, var_id])
        
    for gene, in_var_list in Var_list.items():
        Var_list[gene] = sorted(in_var_list)

    return Vars, Var_list


"""
"""
def read_Gene_links(fname):
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
def test_Gene_genotyping(base_fname,
                         locus_list,
                         only_locus_list,
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
                         output_base,
                         error_correction,
                         discordant,
                         display_alleles,
                         verbose,
                         debug_instr):
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(test_Gene_genotyping))
    ex_path = os.path.dirname(curr_script)

    if not os.path.exists("hisatgenotype_db"):
        typing_common.clone_hisatgenotype_database()

    simulation = (read_fname == [] and alignment_fname == "")

    # Download human genome and HISAT2 index
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]

    if not typing_common.check_files(HISAT2_fnames):
        typing_common.download_genome_and_index(ex_path)

    # Check if the pre-existing files (hla*) are compatible with the current parameter setting
    if os.path.exists("%s.ref" % base_fname):
        left = 0
        Gene_genes = []
        BACKBONE = False
        for line in open("%s.ref" % base_fname):
            Gene_name = line.strip().split()[0]
            if Gene_name.find("BACKBONE") != -1:
                BACKBONE = True
            Gene_gene = Gene_name.split('*')[0]
            Gene_genes.append(Gene_gene)
        delete_hla_files = False
        if not BACKBONE:
            delete_hla_files = True
        if len(locus_list) == 0:
            locus_list = Gene_genes
        if not set(locus_list).issubset(set(Gene_genes)):
            delete_hla_files = True
        if delete_hla_files:
            os.system("rm %s*" % base_fname)

    # Extract HLA variants, backbone sequence, and other sequeces  
    Gene_fnames = [base_fname + "_backbone.fa",
                   base_fname + "_sequences.fa",
                   base_fname + ".ref",
                   base_fname + ".snp",
                   base_fname + ".index.snp",
                   base_fname + ".haplotype",
                   base_fname + ".link"]
    
    if verbose >= 1:
        print >> sys.stderr, Gene_fnames
    
    if not typing_common.check_files(Gene_fnames):
        extract_cmd = ["hisatgenotype_extract_vars.py"]
        if len(only_locus_list) > 0:
            extract_cmd += ["--locus-list", ','.join(only_locus_list)]

        extract_cmd += ["--base", base_fname]
        if not partial:
            extract_cmd += ["--no-partial"]

        if base_fname == "codis":
            extract_cmd += ["--whole-haplotype"]
        else:
            extract_cmd += ["--inter-gap", "30",
                            "--intra-gap", "50"]

        # DK - debugging purposes
        extract_cmd += ["--min-var-freq", "0.1"]
        
        if base_fname == "codis":
            extract_cmd += ["--leftshift"]
        
        # DK - debugging purposes
        # extract_cmd += ["--ext-seq", "300"]
        if verbose >= 1:
            print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
        proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()
        
        if not typing_common.check_files(Gene_fnames):
            print >> sys.stderr, "Error: hisatgenotype_extract_vars failed!"
            sys.exit(1)

    for aligner, index_type in aligners:
        if aligner == "hisat2":
            # Build HISAT2 graph indexes based on the above information
            if index_type == "graph":
                Gene_hisat2_graph_index_fnames = ["%s.graph.%d.ht2" % (base_fname, i+1) for i in range(8)]
                if not typing_common.check_files(Gene_hisat2_graph_index_fnames):
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
                    if not typing_common.check_files(Gene_hisat2_graph_index_fnames):
                        print >> sys.stderr, "Error: indexing HLA failed!  Perhaps, you may have forgotten to build hisat2 executables?"
                        sys.exit(1)
            # Build HISAT2 linear indexes based on the above information
            else:
                assert index_type == "linear"
                Gene_hisat2_linear_index_fnames = ["%s.linear.%d.ht2" % (base_fname, i+1) for i in range(8)]
                if not typing_common.check_files(Gene_hisat2_linear_index_fnames):
                    hisat2_build = os.path.join(ex_path, "hisat2-build")
                    build_cmd = [hisat2_build,
                                 "%s_backbone.fa,%s_sequences.fa" % (base_fname, base_fname),
                                 "%s.linear" % base_fname]
                    proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
                    proc.communicate()        
                    if not typing_common.check_files(Gene_hisat2_linear_index_fnames):
                        print >> sys.stderr, "Error: indexing HLA failed!"
                        sys.exit(1)
        else:
            assert aligner == "bowtie2" and index_type == "linear"
            # Build Bowtie2 indexes based on the above information
            Gene_bowtie2_index_fnames = ["%s.%d.bt2" % (base_fname, i+1) for i in range(4)]
            Gene_bowtie2_index_fnames += ["%s.rev.%d.bt2" % (base_fname, i+1) for i in range(2)]
            if not typing_common.check_files(Gene_bowtie2_index_fnames):
                build_cmd = ["bowtie2-build",
                             "%s_backbone.fa,%s_sequences.fa" % (base_fname, base_fname),
                             base_fname]
                proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'))
                proc.communicate()        
                if not typing_common.check_files(Gene_bowtie2_index_fnames):
                    print >> sys.stderr, "Error: indexing HLA failed!"
                    sys.exit(1)

    # Read partial alleles from hla.data (temporary)
    partial_alleles = set()
    for line in open("hisatgenotype_db/HLA/hla.dat"):
        if not line.startswith("DE"):
            continue
        allele_name = line.split()[1][:-1]
        if allele_name.startswith("HLA-"):
            allele_name = allele_name[4:]
        gene = allele_name.split('*')[0]
        if line.find("partial") != -1:
            partial_alleles.add(allele_name)

    # Read HLA alleles (names and sequences)
    refGenes, refGene_loci = {}, {}
    for line in open("%s.ref" % base_fname):
        Gene_name, chr, left, right, length, exon_str, strand = line.strip().split()
        Gene_gene = Gene_name.split('*')[0]
        assert not Gene_gene in refGenes
        refGenes[Gene_gene] = Gene_name
        left, right = int(left), int(right)
        exons = []
        for exon in exon_str.split(','):
            exon_left, exon_right = exon.split('-')
            exons.append([int(exon_left), int(exon_right)])
        refGene_loci[Gene_gene] = [Gene_name, chr, left, right, exons]
    Genes = {}
    if len(locus_list) == 0:
        locus_list = refGene_loci.keys()

    read_Gene_alleles(base_fname + "_backbone.fa", Genes)
    read_Gene_alleles(base_fname + "_sequences.fa", Genes)

    # HLA gene alleles
    Gene_names = {}
    for Gene_gene, data in Genes.items():
        Gene_names[Gene_gene] = list(data.keys())

    # HLA gene allele lengths
    Gene_lengths = {}
    for Gene_gene, Gene_alleles in Genes.items():
        Gene_lengths[Gene_gene] = {}
        for allele_name, seq in Gene_alleles.items():
            Gene_lengths[Gene_gene][allele_name] = len(seq)

    # Read HLA variants, and link information
    Vars, Var_list = read_Gene_vars("%s.snp" % base_fname)
    Links = read_Gene_links("%s.link" % base_fname)
    # Test HLA typing
    test_list = []
    if simulation:
        basic_test, pair_test = True, False
        if debug_instr:
            if "basic_test" in debug_instr:
                basic_test, pair_test = True, False
            else:
                basic_test, pair_test = False, True

        test_passed = {}
        test_list = []
        genes = list(set(locus_list) & set(Gene_names.keys()))
        if basic_test:
            for gene in genes:
                Gene_gene_alleles = Gene_names[gene]
                for Gene_name in Gene_gene_alleles:
                    if Gene_name.find("BACKBONE") != -1:
                        continue
                    test_list.append([[Gene_name]])
        if pair_test:
            test_size = 500
            allele_count = 2
            for test_i in range(test_size):
                test_pairs = []
                for gene in genes:
                    Gene_gene_alleles = []

                    for allele in Gene_names[gene]:
                        if allele.find("BACKBONE") != -1:
                            continue

                        if "full" in debug:
                            if allele in partial_alleles:
                                continue

                        Gene_gene_alleles.append(allele)
                    nums = [i for i in range(len(Gene_gene_alleles))]
                    random.shuffle(nums)
                    test_pairs.append(sorted([Gene_gene_alleles[nums[i]] for i in range(allele_count)]))
                test_list.append(test_pairs)

        if "test_list" in debug_instr:
            test_list = [[debug_instr["test_list"].split('-')]]
            
        for test_i in range(len(test_list)):
            if "test_id" in debug_instr:
                test_ids = debug_instr["test_id"].split('-')
                if str(test_i + 1) not in test_ids:
                    continue

            print >> sys.stderr, "Test %d" % (test_i + 1), str(datetime.now())
            test_locus_list = test_list[test_i]
            num_frag_list = typing_common.simulate_reads(Genes,
                                                         base_fname,
                                                         test_locus_list,
                                                         Vars,
                                                         Links,
                                                         simulate_interval,
                                                         read_len,
                                                         fragment_len,
                                                         perbase_errorrate,
                                                         perbase_snprate,
                                                         skip_fragment_regions)

            assert len(num_frag_list) == len(test_locus_list)
            for i_ in range(len(test_locus_list)):
                test_Gene_names = test_locus_list[i_]
                num_frag_list_i = num_frag_list[i_]
                assert len(num_frag_list_i) == len(test_Gene_names)
                for j_ in range(len(test_Gene_names)):
                    test_Gene_name = test_Gene_names[j_]
                    gene = test_Gene_name.split('*')[0]
                    test_Gene_seq = Genes[gene][test_Gene_name]
                    seq_type = "partial" if test_Gene_name in partial_alleles else "full"
                    print >> sys.stderr, "\t%s - %d bp (%s sequence, %d pairs)" % (test_Gene_name, len(test_Gene_seq), seq_type, num_frag_list_i[j_])

            if "single-end" in debug_instr:
                read_fname = ["%s_input_1.fa" % base_fname]
            else:
                read_fname = ["%s_input_1.fa" % base_fname, "%s_input_2.fa" % base_fname]

            fastq = False
            tmp_test_passed = typing(ex_path,
                                     simulation,
                                     base_fname,
                                     test_locus_list,
                                     partial,
                                     partial_alleles,
                                     refGenes,
                                     Genes,                       
                                     Gene_names,
                                     Gene_lengths,
                                     refGene_loci,
                                     Vars,
                                     Var_list,
                                     Links,
                                     aligners,
                                     num_editdist,
                                     assembly,
                                     output_base,
                                     error_correction,
                                     discordant,
                                     display_alleles,
                                     fastq,
                                     read_fname,
                                     alignment_fname,
                                     num_frag_list,
                                     read_len,
                                     fragment_len,
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
        print >> sys.stderr, "\t", ' '.join(locus_list)
        fastq = True
        typing(ex_path,
               simulation,
               base_fname,
               locus_list,
               partial,
               partial_alleles,
               refGenes,
               Genes,                       
               Gene_names,
               Gene_lengths,
               refGene_loci,
               Vars,
               Var_list,
               Links,
               aligners,
               num_editdist,
               assembly,
               output_base,
               error_correction,
               discordant,
               display_alleles,
               fastq,
               read_fname,
               alignment_fname,
               [],
               read_len,
               fragment_len,
               threads,
               best_alleles,
               verbose)


"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='test HLA genotyping')
    parser.add_argument("--base", "--base-fname",
                        dest="base_fname",
                        type=str,
                        default="hla",
                        help="base filename for backbone sequence, variants, and linking info (default: hla)")
    parser.add_argument("--locus-list",
                        dest="locus_list",
                        type=str,
                        default="",
                        help="A comma-separated list of genes (default: empty, all genes)")
    parser.add_argument('--no-partial',
                        dest='partial',
                        action='store_false',
                        help='Include partial alleles (e.g. A_nuc.fasta)')
    parser.add_argument("--aligner-list",
                        dest="aligners",
                        type=str,
                        default="hisat2.graph",
                        help="A comma-separated list of aligners such as hisat2.graph,hisat2.linear,bowtie2.linear (default: hisat2.graph)")
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
                        default=350,
                        help="Length of fragments (default: 350)")
    parser.add_argument("--best-alleles",
                        dest="best_alleles",
                        action='store_true',
                        help="")
    parser.add_argument("--random-seed",
                        dest="random_seed",
                        type=int,
                        default=1,
                        help="A seeding number for randomness (default: 1)")
    parser.add_argument("--num-editdist",
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
    parser.add_argument("--output-base", "--assembly-base",
                        dest="output_base",
                        type=str,
                        default="assembly_graph",
                        help="base file name (default: assembly_graph)")
    parser.add_argument("--no-assembly",
                        dest="assembly",
                        action="store_false",
                        help="Perform assembly")
    parser.add_argument("--no-error-correction",
                        dest="error_correction",
                        action="store_false",
                        help="Correct sequencing errors")
    parser.add_argument("--only-locus-list",
                        dest="only_locus_list",
                        type=str,
                        default="",
                        help="A comma-separated list of genes (default: empty, all genes)")
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
    if args.locus_list == "":
        locus_list = []
    else:
        locus_list = args.locus_list.split(',')
    if args.only_locus_list == "":
        only_locus_list = []
    else:
        locus_list = only_locus_list = args.only_locus_list.split(',')    
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
                fields = item.split(':')
                assert len(fields) >= 2
                key, value = fields[0], ':'.join(fields[1:])
                debug[key] = value
            else:
                debug[item] = 1

    if not args.partial:
        print >> sys.stderr, "Warning: --no-partial should be used for debugging purpose only."

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
    test_Gene_genotyping(args.base_fname,
                         locus_list,
                         only_locus_list,
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
                         args.output_base,
                         args.error_correction,
                         args.discordant,
                         display_alleles,
                         args.verbose_level,
                         debug)

