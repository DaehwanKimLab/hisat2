#!/usr/bin/env python
#
# Copyright 2017, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT-genotype.
#
# HISAT-genotype is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT-genotype is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT-genotype.  If not, see <http://www.gnu.org/licenses/>.
#


import sys, os, subprocess, re
import random
import math, multiprocessing
from datetime import datetime, date, time
from copy import deepcopy
import hisatgenotype_typing_common as typing_common, hisatgenotype_assembly_graph as assembly_graph


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
def get_rep_alleles(Links, exon_vars, in_alleles = None):
    allele_vars = {}
    for var, alleles in Links.items():
        if var not in exon_vars:
            continue
        for allele in alleles:
            if in_alleles != None and allele not in in_alleles:
                continue
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
        print >> sys.stderr, cmp_list
        print >> sys.stderr, read_seq

    num_correction = 0
    i = 0
    while i < len(cmp_list):
        type, left, length = cmp_list[i][:3]
        assert length > 0
        if left >= len(ref_seq):
            break
        if type == "match":
            middle_cmp_list = []
            last_j = 0
            for j in range(length):
                if read_pos + j >= len(read_seq) or \
                   left + j >= len(ref_seq):
                    continue
                
                read_bp, ref_bp = read_seq[read_pos + j], ref_seq[left + j]
                assert left + j < len(mpileup)
                nt_set = mpileup[left + j][0]
                if len(nt_set) > 0 and read_bp not in nt_set:
                    read_bp = 'N' if len(nt_set) > 1 else nt_set[0]                    
                    read_seq = read_seq[:read_pos + j] + read_bp + read_seq[read_pos + j + 1:]
                    assert read_bp != ref_bp
                    new_cmp = ["mismatch", left + j, 1, "unknown"]
                    num_correction += 1
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
                print >> sys.stderr, left, read_bp, ref_bp, mpileup[left]

            if len(nt_set) > 0 and read_bp not in nt_set:
                read_bp = 'N' if len(nt_set) > 1 else nt_set[0]
                read_seq = read_seq[:read_pos] + read_bp + read_seq[read_pos+1:]
                if read_bp == 'N':
                    cmp_list[i][3] = "unknown"
                elif read_bp == ref_bp:
                    cmp_list[i] = ["match", left, 1]
                    num_correction += 1
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
                    print >> sys.stderr, left, read_bp, ref_bp, mpileup[left]
                    print >> sys.stderr, cmp_list[i]

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
        print >> sys.stderr, cmp_list
        print >> sys.stderr, read_seq
                            
    return cmp_list, read_seq, num_correction


"""
"""
def typing(simulation,
           base_fname,
           locus_list,
           genotype_genome,
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
           keep_alignment,
           allow_discordant,
           type_primary_exons,
           remove_low_abundance_alleles,
           display_alleles,
           fastq,
           read_fname,
           alignment_fname,
           num_frag_list,
           read_len,
           fragment_len,
           threads,
           best_alleles,
           verbose,
           assembly_verbose,
           out_dir,
           test_i = 0):

    core_fid = "" # May add to bottom of options
    if simulation:
        test_passed = {}
        report_file = open('%s/%s-%s_test-%s.report' % (out_dir, output_base, base_fname, test_i + 1), 'w')
    else:
        core_fid = read_fname[0].split('/')[-1].split('.')[0]
        report_file = open('%s/%s-%s-%s.report' % (out_dir, output_base, base_fname, core_fid), 'w')
    
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
                # alignment_fname = read_fname[0].split('/')[-1] #Temporary removal for testing
                alignment_fname = "%s.bam" % core_fid # '.'.join(alignment_fname.split('.')[:2])
                    
            typing_common.align_reads(aligner,
                                      simulation,
                                      genotype_genome if genotype_genome != "" else (base_fname + "." + index_type),
                                      index_type,
                                      base_fname,
                                      read_fname,
                                      fastq,
                                      threads,
                                      alignment_fname,
                                      verbose)
            
        for test_Gene_names in locus_list:
            if base_fname == "genome":
                if simulation:
                    region_chr, region_left, region_right = test_Gene_names[0]
                else:
                    region_chr, region_left, region_right = test_Gene_names
                gene = "%s:%d-%d" % (region_chr, region_left, region_right)
            else:
                if simulation:
                    gene = test_Gene_names[0].split('*')[0]
                else:
                    gene = test_Gene_names
                
            ref_allele = refGenes[gene]
            ref_seq = Genes[gene][ref_allele]
            ref_locus = refGene_loci[gene]
            ref_exons, ref_primary_exons = ref_locus[-2], ref_locus[-1]
            novel_var_count = 0        
            gene_vars, gene_var_list = deepcopy(Vars[gene]), deepcopy(Var_list[gene])
            cur_maxright = -1
            gene_var_maxrights = {}
            for var_pos, var_id in gene_var_list:
                var_type, var_pos, var_data = gene_vars[var_id]
                if var_type == "deletion":
                    var_pos = var_pos + int(var_data) - 1
                cur_maxright = max(cur_maxright, var_pos)
                gene_var_maxrights[var_id] = cur_maxright
                    
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
            if genotype_genome != "":
                _, chr, left, right = ref_locus[:4]
                alignview_cmd += ["%s:%d-%d" % (chr, left+1, right+1)]
                base_locus = left

            if index_type == "graph":
                alignview_cmd += [ref_allele]
                mpileup = typing_common.get_mpileup(alignview_cmd,
                                                    ref_seq,
                                                    base_locus,
                                                    gene_vars,
                                                    allow_discordant)

                if base_fname == "codis":
                    pair_interdist = typing_common.get_pair_interdist(alignview_cmd,
                                                                      simulation,
                                                                      verbose)
                else:
                    pair_interdist = None

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
                if var_id not in Links:
                    continue
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
            primary_exon_vars = get_exonic_vars(gene_vars, ref_primary_exons)

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
                                             ref_primary_exons,
                                             partial_alleles,
                                             true_allele_nodes,
                                             {}, # predicted_allele_nodes, which is empty for now
                                             display_allele_nodes,
                                             simulation)

            # Choose allele representives from those that share the same exonic sequences
            allele_reps, allele_rep_groups = get_rep_alleles(Links, exon_vars)
            allele_rep_set = set(allele_reps.values())

            # Choose allele representives from those that share the primary exonic sequences
            primary_exon_allele_reps, primary_exon_allele_rep_groups = get_rep_alleles(Links, primary_exon_vars, allele_rep_set)
            primary_exon_allele_rep_set = set(primary_exon_allele_reps.values())

            # Sanity check
            for exon_allele in primary_exon_allele_reps.keys():
                # DK - debugging purposes
                if exon_allele not in allele_rep_set:
                    print exon_allele, allele_reps[exon_allele], exon_allele in primary_exon_allele_reps.keys()
                    
                assert exon_allele in allele_rep_set
                                    
            # For checking alternative alignments near the ends of alignments
            Alts_left, Alts_right = typing_common.get_alternatives(ref_seq,
                                                                   allele_vars,
                                                                   gene_vars,
                                                                   gene_var_list,
                                                                   verbose >= 2)

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
            Gene_primary_exons_counts, Gene_primary_exons_cmpt = {}, {}
            Gene_exons_counts, Gene_exons_cmpt = {}, {}
            Gene_counts, Gene_cmpt = {}, {}
            num_reads, num_pairs = 0, 0

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
                    if base_fname == "genome" and len(count_per_read) == 1:
                        for allele in count_per_read.keys():
                            count_per_read[allele] = add
                        return
                    
                    orig_ht = ht
                    ht = ht.split('-')

                    assert len(ht) >= 2
                    left, right = int(ht[0]), int(ht[-1])
                    assert left <= right

                    ht = ht[1:-1]
                    alleles = set(Genes[gene].keys()) - set([ref_allele])
                    for i in range(len(ht)):
                        var_id = ht[i]
                        if var_id.startswith("nv") or \
                           var_id not in Links:
                            continue
                        alleles &= set(Links[var_id])
                    ht = set(ht)

                    tmp_alleles = set()
                    var_idx = typing_common.lower_bound(gene_var_list, right + 1)
                    var_idx = min(var_idx, len(gene_var_list) - 1)
                    while var_idx >= 0:
                        _, var_id = gene_var_list[var_idx]
                        if var_id.startswith("nv") or \
                           var_id in ht or \
                           var_id not in Links:
                            var_idx -= 1
                            continue
                        if var_id in gene_var_maxrights and gene_var_maxrights[var_id] < left:
                            break
                        var_type, var_left, var_data = gene_vars[var_id]
                        var_right = var_left
                        if var_type == "deletion":
                            var_right = var_left + int(var_data) - 1
                        if (var_left >= left and var_left <= right) or \
                           (var_right >= left and var_right <= right):
                            tmp_alleles |= set(Links[var_id])
                        var_idx -= 1                        
                    alleles -= tmp_alleles
                    alleles &= set(count_per_read.keys())
                    
                    for allele in alleles:
                        count_per_read[allele] += add

                    return len(alleles)

                # Identify best pairs
                def choose_pairs(left_positive_hts, right_positive_hts):
                    if len(left_positive_hts) > 0 and \
                       len(right_positive_hts) > 0 and \
                       max(len(left_positive_hts), len(right_positive_hts)) >= 2:
                        expected_inter_dist = pair_interdist
                        """
                        if simulation:
                            expected_inter_dist = fragment_len - read_len * 2
                        """
                            
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

                def get_exon_haplotypes(ht, exons):
                    if len(exons) <= 0:
                        return []
                    
                    debug_ht = deepcopy(ht)
                    ht = ht.split('-')
                    assert len(ht) >= 2
                    ht[0], ht[-1] = int(ht[0]), int(ht[-1])
                    exon_hts = []
                    for e_left, e_right in exons:
                        assert len(ht) >= 2
                        ht_left, ht_right = ht[0], ht[-1]
                        if e_left > ht_right or e_right < ht_left:
                            continue

                        new_ht = deepcopy(ht)
                        if ht_left < e_left:
                            split = False
                            for i in range(1, len(new_ht) - 1):
                                var_id = new_ht[i]
                                type, left, data = gene_vars[var_id]
                                if (type != "deletion" and left >= e_left) or \
                                   (type == "deletion" and left - 1 >= e_left):
                                    ht_left = e_left
                                    new_ht = [ht_left] + new_ht[i:]
                                    split = True
                                    break
                                if type == "deletion":
                                    right = left + int(data)
                                    if right >= e_left:
                                        ht_left = right
                                        new_ht = [right] + new_ht[i+1:]
                                        split = True
                                        break
                            if not split:
                                ht_left = e_left
                                new_ht = [ht_left, ht_right]
                        assert ht_left >= e_left
                        if ht_right > e_right:
                            split = False
                            for i in reversed(range(1, len(new_ht) - 1)):
                                var_id = new_ht[i]
                                type, right, data = gene_vars[var_id]
                                if type == "deletion":
                                    right = right + int(data) - 1
                                if (type != "deletion" and right <= e_right) or \
                                   (type == "deletion" and right + 1 <= e_right):
                                    ht_right = e_right
                                    new_ht = new_ht[:i+1] + [ht_right]
                                    split = True
                                    break
                                if type == "deletion":
                                    left = right - int(data)
                                    if left <= e_right:
                                        ht_right = left
                                        new_ht = new_ht[:i] + [ht_right]
                                        split = True
                                        break
                            if not split:
                                ht_right = e_right
                                new_ht = [ht_left, ht_right]

                        if len(new_ht) == 2:
                            new_ht = "%d-%d" % (new_ht[0], new_ht[-1])
                        else:
                            assert len(new_ht) > 2
                            new_ht = "%d-%s-%d" % (new_ht[0], '-'.join(new_ht[1:-1]), new_ht[-1])
                        assert ht_left <= ht_right
                        exon_hts.append(new_ht)

                    return exon_hts

                # Positive evidence for left and right reads
                left_positive_hts, right_positive_hts = set(), set()
                
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
                        Zs_str = Zs
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
                    num_error_correction = 0
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
                                new_cmp_list, read_seq, _num_error_correction = error_correct(ref_seq,
                                                                                              read_seq,
                                                                                              read_pos,
                                                                                              mpileup,
                                                                                              gene_vars,
                                                                                              gene_var_list,
                                                                                              cmp_list[cmp_list_i:],
                                                                                              node_read_id == "aHSQ1008:175:C0JVFACXX:5:1109:17665:21583|L")
                                cmp_list = cmp_list[:cmp_list_i] + new_cmp_list
                                num_error_correction += _num_error_correction

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
                            if right_pos < len(mpileup):
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

                    if sum(softclip) > 0:
                        continue

                    if right_pos > len(ref_seq):
                        continue

                    if num_error_correction > max(1, num_editdist):
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
                                    var_id_, novel_var_count = add_novel_var(gene_vars,
                                                                             gene_var_list,
                                                                             novel_var_count,
                                                                             type_ if type_ != "mismatch" else "single",
                                                                             pos_,
                                                                             data_)
                                    cmp_list[cmp_i][3] = var_id_
                            if var_id_ != "unknown":
                                if var_id_ not in var_count:
                                    var_count[var_id_] = 1
                                else:
                                    var_count[var_id_] += 1
                                
                        if type_ != "deletion":
                            read_pos += length_

                    # Count the number of reads aligned uniquely with some constraints
                    num_reads += 1

                    def add_stat(Gene_cmpt, Gene_counts, Gene_count_per_read, include_alleles = set()):
                        if len(Gene_count_per_read) <= 0:
                            return ""
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

                        if verbose >= 2:
                            alleles = ["", ""]
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
                                    print >> sys.stderr, alleles[0], Gene_count_per_read[alleles[0]]
                                    print >> sys.stderr, alleles[1], Gene_count_per_read[alleles[1]]
                                    if allele1_found:
                                        print >> sys.stderr, ("%s\tread_id %s - %d vs. %d]" % (alleles[0], prev_read_id, max_count, Gene_count_per_read[alleles[1]]))
                                    else:
                                        print >> sys.stderr, ("%s\tread_id %s - %d vs. %d]" % (alleles[1], prev_read_id, max_count, Gene_count_per_read[alleles[0]]))

                        cur_cmpt = sorted(list(cur_cmpt))
                        cur_cmpt = '-'.join(cur_cmpt)
                        if not cur_cmpt in Gene_cmpt:
                            Gene_cmpt[cur_cmpt] = 1
                        else:
                            Gene_cmpt[cur_cmpt] += 1

                        return cur_cmpt

                    if read_id != prev_read_id:
                        if prev_read_id != None:
                            num_pairs += 1
                            # DK - needs more test
                            #      Several alleles go over 100 bps
                            """
                            if base_fname == "codis" and gene == "D18S51":
                                left_positive_hts, right_positive_hts = choose_pairs(left_positive_hts, right_positive_hts)
                            """

                            for positive_ht in left_positive_hts | right_positive_hts:
                                primary_exon_hts = get_exon_haplotypes(positive_ht, ref_primary_exons)
                                for exon_ht in primary_exon_hts:
                                    add_count(Gene_primary_exons_count_per_read, exon_ht, 1)
                                exon_hts = get_exon_haplotypes(positive_ht, ref_exons)
                                for exon_ht in exon_hts:
                                    add_count(Gene_exons_count_per_read, exon_ht, 1)
                                add_count(Gene_count_per_read, positive_ht, 1)

                            # DK - debugging purposes
                            if prev_read_id.startswith("NS500497:33:HY32TBGXX:3:13511:0:56517876") and False:
                                print prev_read_id, left_positive_hts, right_positive_hts
                                max_count = max(Gene_primary_exons_count_per_read.values())
                                for allele, count in Gene_primary_exons_count_per_read.items():
                                    if allele not in primary_exon_allele_rep_set:
                                        continue
                                    if count < max_count:
                                        continue
                                    print allele, count

                            # DK - debugging purposes
                            """
                            debug_allele_id = "TH01*10"
                            assert debug_allele_id in Gene_gen_count_per_read
                            debug_max_read_count = max(Gene_gen_count_per_read.values())
                            debug_read_count = Gene_gen_count_per_read[debug_allele_id]
                            if debug_read_count < debug_max_read_count:
                                print prev_read_id, debug_read_count, debug_max_read_count, Gene_gen_count_per_read
                                print "\t", left_positive_hts, right_positive_hts
                                None
                            if prev_read_id == "HSQ1008:175:C0JVFACXX:5:1109:17665:21583":
                                for line in prev_lines:
                                    print line
                                print "left_positive_hts :", left_positive_hts
                                print "right_positive_hts:", right_positive_hts
                                print "exon:", debug_read_count, "max:", debug_max_read_count
                                print "gen:", Gene_gen_count_per_read[debug_allele_id], "max:", max(Gene_gen_count_per_read.values())

                                for allele_id, count in Gene_count_per_read.items():
                                    if count == debug_max_read_count:
                                        None
                                        # print "allele max:", allele_id, count
                                # sys.exit(1)
                                None
                            """                                

                            cur_cmpt, cur_cmpt_gen = "", ""
                            if base_fname == "hla":
                                cur_primary_exons_cmpt = add_stat(Gene_primary_exons_cmpt, Gene_primary_exons_counts, Gene_primary_exons_count_per_read, primary_exon_allele_rep_set)

                                # DK - debugging purposes
                                # for cmpt, count in Gene_primary_exons_count_per_read.items():
                                if cur_primary_exons_cmpt.find("A*24:145") != -1 and cur_primary_exons_cmpt.find("A*24:02:01") == -1:
                                    print prev_read_id
                                    print cur_primary_exons_cmpt
        

            
                                cur_exons_cmpt = add_stat(Gene_exons_cmpt, Gene_exons_counts, Gene_exons_count_per_read, allele_rep_set)
                                cur_cmpt = add_stat(Gene_cmpt, Gene_counts, Gene_count_per_read)
                            else:
                                cur_cmpt = add_stat(Gene_cmpt, Gene_counts, Gene_count_per_read)
                            for read_id_, read_id_i, read_node in read_nodes:
                                asm_graph.add_node(read_id_,
                                                   read_id_i,
                                                   read_node,
                                                   simulation)
                            read_nodes, read_var_list = [], []
                            if simulation and \
                               verbose >= 2 and \
                               base_fname in ["hla", "codis"]:
                                cur_cmpt = cur_cmpt.split('-') if cur_cmpt != "" else set()
                                cur_cmpt_gen = cur_cmpt_gen.split('-') if cur_cmpt_gen != "" else set()
                                show_debug = (partial and cur_cmpt != "" and not set(cur_cmpt) & set(test_Gene_names)) or \
                                              (not partial and cur_cmpt_gen != "" and not set(cur_cmpt_gen) & set(test_Gene_names))
                                              
                                if show_debug:
                                    print "%s are chosen instead of %s" % (cur_cmpt if partial else cur_cmpt_gen, '-'.join(test_Gene_names))
                                    for prev_line in prev_lines:
                                        print "\t", prev_line

                            prev_lines = []

                        left_positive_hts, right_positive_hts = set(), set()                        
                        Gene_primary_exons_count_per_read, Gene_exons_count_per_read, Gene_count_per_read = {}, {}, {}
                        for allele in Gene_names[gene]:
                            if allele.find("BACKBONE") != -1:
                                continue
                            if base_fname == "genome" and allele.find("GRCh38") != -1:
                                continue
                            if allele in primary_exon_allele_rep_set:
                                Gene_primary_exons_count_per_read[allele] = 0
                            if allele in allele_rep_set:
                                Gene_exons_count_per_read[allele] = 0
                            Gene_count_per_read[allele] = 0

                    prev_lines.append(line)

                    # Remove mismatches due to unknown or novel variants
                    cmp_list2 = []
                    for cmp in cmp_list:
                        cmp = deepcopy(cmp)
                        type, pos, length = cmp[:3]
                        if type == "match":
                            if len(cmp_list2) > 0 and cmp_list2[-1][0] == "match":
                                cmp_list2[-1][2] += length
                            else:
                                cmp_list2.append(cmp)
                        elif type == "mismatch" and \
                             (cmp[3] == "unknown" or cmp[3].startswith("nv")):
                            if len(cmp_list2) > 0 and cmp_list2[-1][0] == "match":
                                cmp_list2[-1][2] += 1
                            else:
                                cmp_list2.append(["match", pos, 1])
                        else:
                            cmp_list2.append(cmp)
                            
                    cmp_list_left, cmp_list_right, cmp_left_alts, cmp_right_alts = \
                    typing_common.identify_ambigious_diffs(ref_seq,
                                                           gene_vars,
                                                           Alts_left,
                                                           Alts_right,
                                                           Alts_left_list,
                                                           Alts_right_list,
                                                           cmp_list2,
                                                           verbose,
                                                           orig_read_id.startswith("HSQ1009:126:D0UUYACXX:4:2212:9787:80992#"))  # debug?

                    mid_ht = []
                    for cmp in cmp_list2[cmp_list_left:cmp_list_right+1]:
                        type = cmp[0]
                        if type not in ["mismatch", "deletion", "insertion"]:
                            continue                            
                        var_id = cmp[3]
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
                    if orig_read_id.startswith("30|R!"):
                        DK_debug = True
                        print line
                        print cmp_list
                        print "positive hts:", left_positive_hts, right_positive_hts
                        print "cmp_list [%d, %d]" % (cmp_list_left, cmp_list_right)

                    if assembly:
                        # Construct multiple candidate realignments for CODIS
                        cmp_llist = []
                        hts = left_positive_hts if is_left_read else right_positive_hts
                        assert len(hts) > 0
                        for ht in hts:
                            cmp_list = []
                            read_pos = 0
                            vars_ = ht.split('-')
                            left_ = int(vars_[0])
                            vars_ = vars_[1:]
                            for var_i in range(len(vars_)):
                                var_id = vars_[var_i]
                                # ref_seq, read_seq
                                if var_i == len(vars_) - 1:
                                    right_ = int(var_id)
                                else:
                                    var_type, var_pos, var_data = gene_vars[var_id]
                                    right_ = var_pos - 1
                                    
                                for pos in range(left_, right_ + 1):
                                    if read_seq[read_pos] != ref_seq[pos]:
                                        if left_ < pos:
                                            cmp_list.append(["match", left_, pos - left_])
                                        cmp_list.append(["mismatch", pos, 1, "unknown"])
                                        left_ = pos + 1
                                    read_pos += 1                                    
                                if left_ <= right_:
                                    cmp_list.append(["match", left_, right_ - left_ + 1])
                                    
                                if var_i == len(vars_) - 1:
                                    left_ = right_ + 1
                                    break

                                if var_type == "single":
                                    cmp_list.append(["mismatch", var_pos, 1, var_id])
                                    left_ = var_pos + 1
                                    read_pos += 1
                                elif var_type == "deletion":
                                    del_len = int(var_data)
                                    cmp_list.append(["deletion", var_pos, del_len, var_id])
                                    left_ = var_pos + del_len                                    
                                else:
                                    assert var_type == "insertion"
                                    cmp_list.append(["insertion", var_pos, len(var_data), var_id])
                                    left_ = var_pos
                                    read_pos += len(var_data)
                                    
                            assert len(cmp_list) > 0
                            cmp_llist.append(cmp_list)

                        for cmp_list_i in range(len(cmp_llist)):
                            # Node
                            cmp_list = cmp_llist[cmp_list_i]
                            read_node_pos, read_node_seq, read_node_qual, read_node_var = -1, [], [], []
                            read_vars = []
                            ref_pos, read_pos = cmp_list[0][1], 0
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
                                elif type == "deletion":
                                    var_id = cmp[3]
                                    del_len = length
                                    read_node_seq += (['D'] * del_len)
                                    read_node_qual += ([''] * del_len)
                                    if len(read_node_seq) > len(read_node_var):
                                        assert len(read_node_seq) == len(read_node_var) + del_len
                                        read_node_var += ([var_id] * del_len)
                                elif type == "insertion":
                                    var_id = cmp[3]
                                    ins_len = length
                                    ins_seq = read_seq[read_pos:read_pos+ins_len]
                                    read_node_seq += ["I%s" % nt for nt in ins_seq]
                                    read_node_qual += list(read_qual[read_pos:read_pos+ins_len])
                                    read_node_var += ([var_id] * ins_len)                                        
                                    read_pos += length
                                else:
                                    assert type == "intron"
                                cmp_i += 1

                            read_nodes.append([node_read_id,
                                               cmp_list_i,
                                               assembly_graph.Node(node_read_id,
                                                                   read_node_pos,
                                                                   read_node_seq,
                                                                   read_node_qual,
                                                                   read_node_var,
                                                                   ref_seq,
                                                                   gene_vars,
                                                                   mpileup,
                                                                   simulation)])

                    prev_read_id = read_id
                    prev_right_pos = right_pos

                if prev_read_id != None:
                    num_pairs += 1
                    if base_fname == "codis" and gene == "D18S51":
                        left_positive_hts, right_positive_hts = choose_pairs(left_positive_hts, right_positive_hts)                            
                    for positive_ht in left_positive_hts | right_positive_hts:
                        primary_exon_hts = get_exon_haplotypes(positive_ht, ref_primary_exons)
                        for exon_ht in primary_exon_hts:
                            add_count(Gene_primary_exons_count_per_read, exon_ht, 1)
                        exon_hts = get_exon_haplotypes(positive_ht, ref_exons)
                        for exon_ht in exon_hts:
                            add_count(Gene_exons_count_per_read, exon_ht, 1)
                        add_count(Gene_count_per_read, positive_ht, 1)

                    if base_fname == "hla":
                        add_stat(Gene_primary_exons_cmpt, Gene_primary_exons_counts, Gene_primary_exons_count_per_read, primary_exon_allele_rep_set)
                        add_stat(Gene_exons_cmpt, Gene_exons_counts, Gene_exons_count_per_read, allele_rep_set)
                    add_stat(Gene_cmpt, Gene_counts, Gene_count_per_read)
                    for read_id_, read_id_i, read_node in read_nodes:
                        asm_graph.add_node(read_id_,
                                           read_id_i,
                                           read_node,
                                           simulation)
                    read_nodes, read_var_list = [], []

                if num_reads <= 0:
                    continue

                for f_ in [sys.stderr, report_file]:
                    print >> f_, "\t\t\t%d reads and %d pairs are aligned" % (num_reads, num_pairs)
                
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
                perform_typing_primary_exon = False
                # Incorporate representive alleles for primary exons (experimental feature)
                if perform_typing_primary_exon:
                    Gene_prob = primary_exon_prob = typing_common.single_abundance(Gene_primary_exons_cmpt)
                    primary_exon_alleles = set()
                    primary_exon_prob_sum = 0.0
                    for prob_i in range(len(primary_exon_prob)):
                        allele, prob = primary_exon_prob[prob_i][:2]
                        if len(primary_exon_allele_rep_groups[allele]) <= 1:
                            continue
                        primary_exon_prob_sum += prob
                        primary_exon_alleles |= set(primary_exon_allele_rep_groups[allele])

                    # Incorporate representative alleles for exons
                    if len(primary_exon_alleles) > 0:
                        Gene_exons_cmpt2 = {}
                        for cmpt, value in Gene_exons_cmpt.items():
                            cmpt2 = []
                            for allele in cmpt.split('-'):
                                if allele in primary_exon_alleles:
                                    cmpt2.append(allele)
                            if len(cmpt2) == 0:
                                continue
                            cmpt2 = '-'.join(cmpt2)
                            if cmpt2 not in Gene_exons_cmpt2:
                                Gene_exons_cmpt2[cmpt2] = value
                            else:
                                Gene_exons_cmpt2[cmpt2] += value
                        exon_prob = typing_common.single_abundance(Gene_exons_cmpt2,
                                                                   remove_low_abundance_alleles)
                        exon_prob2 = {}
                        for allele, prob in primary_exon_prob:
                            if allele not in primary_exon_alleles:
                                exon_prob2[allele] = prob
                        for allele, prob in exon_prob:
                            exon_prob2[allele] = prob * primary_exon_prob_sum
                        exon_prob = [[allele, prob] for allele, prob in exon_prob2.items()]
                        Gene_prob = exon_prob = sorted(exon_prob, cmp=typing_common.Gene_prob_cmp)
                else:
                    # Incorporate representative alleles for exons
                    Gene_prob =  exon_prob = typing_common.single_abundance(Gene_exons_cmpt,
                                                                            remove_low_abundance_alleles)

                exon_alleles = set()
                exon_prob_sum = 0.0
                for prob_i in range(len(exon_prob)):
                    allele, prob = exon_prob[prob_i][:2]
                    if prob_i >= 10 and prob < 0.03:
                        break
                    if len(allele_rep_groups[allele]) <= 1:
                        continue

                    exon_prob_sum += prob
                    exon_alleles |= set(allele_rep_groups[allele])

                # Incorporate full-length alleles, non-representative alleles
                if len(exon_alleles) > 0:
                    Gene_cmpt2 = {}
                    for cmpt, value in Gene_cmpt.items():
                        cmpt2 = []
                        for allele in cmpt.split('-'):
                            if allele in exon_alleles:
                                cmpt2.append(allele)
                        if len(cmpt2) == 0:
                            continue
                        cmpt2 = '-'.join(cmpt2)
                        if cmpt2 not in Gene_cmpt2:
                            Gene_cmpt2[cmpt2] = value
                        else:
                            Gene_cmpt2[cmpt2] += value
                    Gene_cmpt = Gene_cmpt2
                    Gene_prob = typing_common.single_abundance(Gene_cmpt,
                                                               True,
                                                               Gene_lengths[gene])

                    Gene_combined_prob = {}
                    for allele, prob in exon_prob:
                        if allele not in exon_alleles:
                            Gene_combined_prob[allele] = prob

                    for allele, prob in Gene_prob:
                        Gene_combined_prob[allele] = prob * exon_prob_sum
                                            
                    Gene_prob = [[allele, prob] for allele, prob in Gene_combined_prob.items()]
                    Gene_prob = sorted(Gene_prob, cmp=typing_common.Gene_prob_cmp)
            else:
                if len(Gene_cmpt.keys()) <= 1:
                    Gene_prob = []
                    if len(Gene_cmpt.keys()) == 1:
                        Gene_prob = [[Gene_cmpt.keys()[0], 1.0]]
                else:
                    Gene_prob = typing_common.single_abundance(Gene_cmpt)

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
                asm_graph.calculate_coverage()
                
                # Start drawing assembly graph
                if simulation:
                    asm_graph.begin_draw("%s/%s-%d.%s.%s" % (out_dir, output_base, test_i+1, base_fname, gene))
                else:
                    asm_graph.begin_draw("%s/%s.%s.%s" % (out_dir, output_base, base_fname, gene))

                # Draw assembly graph
                begin_y = asm_graph.draw(0, "a. Read alignment")
                begin_y += 200
                
                # Apply De Bruijn graph
                asm_graph.guided_DeBruijn(assembly_verbose)

                # Draw assembly graph
                begin_y = asm_graph.draw(begin_y, "b. Assembly")
                begin_y += 200

                # Draw assembly graph
                asm_graph.nodes = asm_graph.nodes2
                asm_graph.to_node, asm_graph.from_node = {}, {}
                begin_y = asm_graph.draw(begin_y, "c. Assembly with known alleles")

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
                    allele_seq, mismatches = list(ref_seq), 0
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

                            var_type, var_pos, var_data = cmp_var
                            if var_type == "single":
                                allele_seq[var_pos] = var_data
                            elif var_type == "deletion":
                                allele_seq[var_pos:var_pos+int(var_data)] = '.' * int(var_data)
                            else:
                                assert var_type == "insertion"
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

                    allele_exons = ref_exons[:]
                    allele_seq = ''.join(allele_seq)
                    del_counts = []
                    for del_i in range(len(allele_seq)):
                        del_count = 0 if del_i == 0 else del_counts[-1]
                        if allele_seq[del_i] == '.':
                            del_count += 1
                        del_counts.append(del_count)
                    for exon_i in range(len(allele_exons)):
                        exon_left, exon_right = allele_exons[exon_i]
                        exon_left -= del_counts[exon_left]
                        exon_right -= del_counts[exon_right]
                        allele_exons[exon_i] = [exon_left, exon_right]
                        
                    allele_seq = allele_seq.replace('.', '')
                    return allele_seq, allele_exons, mismatches
                    
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
                            allele_seq, allele_exons, allele_mm = compare_alleles(cmp_vars, node_vars)
                            print >> f_, "\t\tallele sequence (%d bps):" % len(allele_seq), allele_seq
                            print >> f_, "\t\texons (zero-based offset):", allele_exons

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
                        _, _, tmp_mismatches = compare_alleles(cmp_vars, node_vars, print_output)
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

            if simulation:
                success = [False for i in range(len(test_Gene_names))]
                found_list = [False for i in range(len(test_Gene_names))]
            for prob_i in range(len(Gene_prob)):
                prob = Gene_prob[prob_i]
                if prob[1] < 0.01:
                    break
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

                # DK - debugging purposes
                """
                # ref_allele_node_ = create_allele_node("A*03:01:01:01")
                ref_allele_node_ = create_allele_node("DQA1*01:02:01:01")
                cmp_node_ = create_allele_node(_allele_rep)
                count_ = 0
                for i_ in range(len(ref_allele_node_.seq)):
                    if assembly_graph.get_major_nt(ref_allele_node_.seq[i_]) != assembly_graph.get_major_nt(cmp_node_.seq[i_]):
                        count_ += 1
                print "\t\t\t\t\tDK:", count_, len(ref_allele_node_.seq)
                vars1, vars2 = allele_vars["DQA1*01:02:01:01"], allele_vars[_allele_rep]
                print "\t\t\t\t\tDK:", set(vars1) - set(vars2), set(vars2) - set(vars1)
                """

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

        if not keep_alignment and remove_alignment_file:
            os.system("rm %s*" % (alignment_fname))

    report_file.close()
    if simulation:
        return test_passed

    
"""
"""
def read_backbone_alleles(genotype_genome, refGene_loci, Genes):
    for gene_name in refGene_loci:
        allele_name, chr, left, right = refGene_loci[gene_name][:4]
        seq_extract_cmd = ["samtools",
                           "faidx",
                           "%s.fa" % genotype_genome,
                           "%s:%d-%d" % (chr, left+1, right+1)]

        length = right - left + 1
        proc = subprocess.Popen(seq_extract_cmd, stdout=subprocess.PIPE, stderr=open("/dev/null", 'w'))
        seq = ""
        for line in proc.stdout:
            line = line.strip()
            if line.startswith('>'):
                continue
            seq += line
        assert len(seq) == length
        assert gene_name not in Genes
        Genes[gene_name] = {}
        Genes[gene_name][allele_name] = seq

        
"""
"""
def read_Gene_alleles_from_vars(Vars, Var_list, Links, Genes):
    for gene_name in Genes:
        # Assert there is only one allele per gene, which is a backbone allele
        assert len(Genes[gene_name]) == 1
        backbone_allele_name, backbone_seq = Genes[gene_name].items()[0]
        gene_vars, gene_var_list = Vars[gene_name], Var_list[gene_name]
        allele_vars = {}
        for _, var_id in gene_var_list:
            if var_id not in Links:
                continue
            for allele_name in Links[var_id]:
                if allele_name not in allele_vars:
                    allele_vars[allele_name] = []
                allele_vars[allele_name].append(var_id)

        for allele_name, vars in allele_vars.items():
            seq = ""
            prev_pos = 0
            for var_id in vars:
                type, pos, data = gene_vars[var_id]
                assert prev_pos <= pos
                if pos > prev_pos:
                    seq += backbone_seq[prev_pos:pos]
                if type == "single":
                    prev_pos = pos + 1
                    seq += data
                elif type == "deletion":
                    prev_pos = pos + int(data)
                else:
                    assert type == "insertion"
                    seq += data
                    prev_pos = pos
            if prev_pos < len(backbone_seq):
                seq += backbone_seq[prev_pos:]
            Genes[gene_name][allele_name] = seq

        if len(Genes[gene_name]) <= 1:
            Genes[gene_name]["%s*GRCh38" % gene_name] = backbone_seq
            
    
"""
"""
def read_Gene_alleles(fname, Genes):
    for line in open(fname):
        if line.startswith(">"):
            allele_name = line.strip().split()[0][1:]
            gene_name = allele_name.split('*')[0]
            if not gene_name in Genes:
                Genes[gene_name] = {}
            if not allele_name in Genes[gene_name]:
                Genes[gene_name][allele_name] = ""
        else:
            Genes[gene_name][allele_name] += line.strip()
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
        Vars[gene][var_id] = [var_type, pos, data]
        Var_list[gene].append([pos, var_id])
        
    for gene, in_var_list in Var_list.items():
        Var_list[gene] = sorted(in_var_list)

    return Vars, Var_list


"""
"""
def read_Gene_vars_genotype_genome(fname, refGene_loci):
    loci = {}
    for gene, values in refGene_loci.items():
        allele_name, chr, left, right = values[:4]
        if chr not in loci:
            loci[chr] = []
        loci[chr].append([allele_name, left, right])
        
    Vars, Var_list = {}, {}
    for line in open(fname):
        var_id, var_type, var_chr, pos, data = line.strip().split('\t')
        if var_chr not in loci:
            continue
        pos = int(pos)
        found = False
        for allele_name, left, right in loci[var_chr]:
            if pos >= left and pos <= right:
                found = True
                break
        if not found:
            continue
        
        gene = allele_name.split('*')[0]
        if not gene in Vars:
            Vars[gene] = {}
            assert not gene in Var_list
            Var_list[gene] = []
            
        assert not var_id in Vars[gene]
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
def genotyping_locus(base_fname,
                     locus_list,
                     genotype_genome,
                     only_locus_list,
                     partial,
                     aligners,
                     read_fname,
                     fastq,
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
                     keep_alignment,
                     discordant,
                     type_primary_exons,
                     remove_low_abundance_alleles,
                     display_alleles,
                     verbose,
                     assembly_verbose,
                     out_dir,
                     debug_instr):
    assert isinstance(base_fname, basestring)
    if ',' in base_fname:
        print >> sys.stderr, "Don't use list in hisatgenotype_locus"
        exit(1)
        
    simulation = (read_fname == [] and alignment_fname == "")
    if genotype_genome == "":
        if not os.path.exists("hisatgenotype_db"):
            typing_common.clone_hisatgenotype_database()

        # Download human genome and HISAT2 index
        HISAT2_fnames = ["grch38",
                         "genome.fa",
                         "genome.fa.fai"]
        
        try: lock.acquire()
        except: pass
        if not typing_common.check_files(HISAT2_fnames):
            typing_common.download_genome_and_index()
        try: lock.release()
        except: pass

    # Check if the pre-existing files (hla*) are compatible with the current parameter setting
    if genotype_genome != "":
        if os.path.exists("%s.locus" % base_fname):
            left = 0
            Gene_genes = []
            BACKBONE = False
            for line in open("%s.locus" % base_fname):
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

    # Extract variants, backbone sequence, and other sequeces  
    if genotype_genome != "":
        genome_fnames = [genotype_genome + ".fa",
                         genotype_genome + ".fa.fai",
                         genotype_genome + ".locus",
                         genotype_genome + ".snp",
                         genotype_genome + ".index.snp",
                         genotype_genome + ".haplotype",
                         genotype_genome + ".link",
                         genotype_genome + ".clnsig",
                         genotype_genome + ".coord",
                         genotype_genome + ".allele",
                         genotype_genome + ".partial"]
        for i in range(8):
            genome_fnames.append(genotype_genome + ".%d.ht2" % (i+1))

        if not typing_common.check_files(genome_fnames):
            print >> sys.stderr, "Error: some of the following files are not available:", ' '.join(genome_fnames)
            sys.exit(1)
    else:
        try: lock.acquire()
        except: pass
        typing_common.extract_database_if_not_exists(base_fname,
                                                     only_locus_list,
                                                     30,              # inter_gap
                                                     50,              # intra_gap
                                                     partial,
                                                     verbose >= 1)        
        for aligner, index_type in aligners:
            typing_common.build_index_if_not_exists(base_fname,
                                                    aligner,
                                                    index_type,
                                                    threads,
                                                    verbose >= 1)
        try: lock.release()
        except: pass

    # Read alleles
    alleles = set()
    if genotype_genome != "":
        for line in open("%s.allele" % genotype_genome):
            family, allele_name = line.strip().split('\t')
            if family == base_fname:
                alleles.add(allele_name)
    else:
        for line in open("%s.allele" % base_fname):
            alleles.add(line.strip())

    # Read partial alleles
    partial_alleles = set()
    if genotype_genome != "":
        for line in open("%s.partial" % genotype_genome):
            family, allele_name = line.strip().split('\t')
            if family == base_fname:
                partial_alleles.add(allele_name)

    else:
        for line in open("%s.partial" % base_fname):
            partial_alleles.add(line.strip())

    # Read alleles (names and sequences)
    refGenes, refGene_loci = {}, {}
    if base_fname == "genome":
        for chr, left, right in locus_list:
            region_name = "%s:%d-%d" % (chr, left, right)
            refGenes[region_name] = region_name
            refGene_loci[region_name] = [region_name, chr, left, right, []]
    else:
        for line in open("%s.locus" % (genotype_genome if genotype_genome != "" else base_fname)):
            fields = line.strip().split()
            if genotype_genome != "" and base_fname != fields[0].lower():
                continue
            if genotype_genome != "":
                _, Gene_name, chr, left, right, exon_str, strand = fields
            else:
                Gene_name, chr, left, right, _, exon_str, strand = fields
            Gene_gene = Gene_name.split('*')[0]
            assert not Gene_gene in refGenes
            refGenes[Gene_gene] = Gene_name
            left, right = int(left), int(right)
            exons, primary_exons = [], []
            for exon in exon_str.split(','):
                primary = exon.endswith('p')
                if primary:
                    exon = exon[:-1]
                exon_left, exon_right = exon.split('-')
                exon_left, exon_right = int(exon_left), int(exon_right)
                exons.append([exon_left, exon_right])
                if primary:
                    primary_exons.append([exon_left, exon_right])
            refGene_loci[Gene_gene] = [Gene_name, chr, left, right, exons, primary_exons]
    Genes = {}
    if len(locus_list) == 0:
        locus_list = refGene_loci.keys()

    # Read variants, and link information
    if genotype_genome:
        Vars, Var_list = read_Gene_vars_genotype_genome("%s.snp" % genotype_genome, refGene_loci)
        Links = read_Gene_links("%s.link" % genotype_genome)
    else:
        Vars, Var_list = read_Gene_vars("%s.snp" % base_fname)
        Links = read_Gene_links("%s.link" % base_fname)

    # Some loci may have only one allele such as AMELX and AMELY
    for gene_name in refGene_loci.keys():
        if gene_name in Vars:
            continue
        Vars[gene_name], Var_list[gene_name], Links[gene_name] = {}, [], {}        

    # Read allele sequences
    if genotype_genome != "":
        read_backbone_alleles(genotype_genome, refGene_loci, Genes)
        read_Gene_alleles_from_vars(Vars, Var_list, Links, Genes)        
    else:
        read_Gene_alleles(base_fname + "_backbone.fa", Genes)
        read_Gene_alleles_from_vars(Vars, Var_list, Links, Genes)

    # alleles corresponding to backbones
    for allele in alleles:
        locus = allele.split('*')[0]
        assert locus in Genes
        if allele not in Genes[locus]:
            Genes[locus][allele] = Genes[locus]["%s*BACKBONE" % locus]

    # Sanity Check
    if os.path.exists(base_fname + "_backbone.fa") and \
       os.path.exists(base_fname + "_sequences.fa"):
        Genes2 = {}
        read_Gene_alleles(base_fname + "_backbone.fa", Genes2)
        read_Gene_alleles(base_fname + "_sequences.fa", Genes2)
        for gene_name, alleles in Genes.items():
            assert gene_name in Genes2
            for allele_name, allele_seq in alleles.items():
                assert allele_name in Genes2[gene_name]
                allele_seq2 = Genes2[gene_name][allele_name]
                assert allele_seq == allele_seq2, 'Problems with: %s - length %d vs %d' % (allele_name, len(allele_seq), len(allele_seq2))

    # alleles names
    Gene_names = {}
    for Gene_gene, data in Genes.items():
        Gene_names[Gene_gene] = list(data.keys())

    # allele lengths
    Gene_lengths = {}
    for Gene_gene, Gene_alleles in Genes.items():
        Gene_lengths[Gene_gene] = {}
        for allele_name, seq in Gene_alleles.items():
            Gene_lengths[Gene_gene][allele_name] = len(seq)

    # Test typing
    test_list = []
    if simulation:
        basic_test, pair_test = True, False
        if debug_instr and "pair" in debug_instr:
            basic_test, pair_test = False, True

        test_passed = {}
        test_list = []
        genes = list(set(locus_list) & set(Gene_names.keys()))
        if basic_test:
            for gene in genes:
                Gene_gene_alleles = Gene_names[gene]
                for allele in Gene_gene_alleles:
                    if allele.find("BACKBONE") != -1:
                        continue
                    test_list.append([[allele]])
                random.shuffle(test_list)
        if pair_test:
            test_size = 200
            allele_count = 2
            for test_i in range(test_size):
                test_pairs = []
                for gene in genes:
                    Gene_gene_alleles = []

                    for allele in Gene_names[gene]:
                        if allele.find("BACKBONE") != -1:
                            continue

                        if "full" in debug_instr:
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
                                                         skip_fragment_regions,
                                                         out_dir)

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
            tmp_test_passed = typing(simulation,
                                     base_fname,
                                     test_locus_list,
                                     genotype_genome,
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
                                     keep_alignment,
                                     discordant,
                                     type_primary_exons,
                                     remove_low_abundance_alleles,
                                     display_alleles,
                                     fastq,
                                     read_fname,
                                     alignment_fname,
                                     num_frag_list,
                                     read_len,
                                     fragment_len,
                                     threads,
                                     best_alleles,
                                     verbose,
                                     assembly_verbose,
                                     out_dir,
                                     test_i)

            for aligner_type, passed in tmp_test_passed.items():
                if aligner_type in test_passed:
                    test_passed[aligner_type] += passed
                else:
                    test_passed[aligner_type] = passed

                print >> sys.stderr, "\t\tPassed so far: %d/%d (%.2f%%)" % (test_passed[aligner_type], test_i + 1, (test_passed[aligner_type] * 100.0 / (test_i + 1)))


        for aligner_type, passed in test_passed.items():
            print >> sys.stderr, "%s:\t%d/%d passed (%.2f%%)" % (aligner_type, passed, len(test_list), passed * 100.0 / len(test_list))
    
    else: # With real reads or BAMs
        if base_fname == "genome":
            print >> sys.stderr, "\t", locus_list
        else:
            print >> sys.stderr, "\t", ' '.join(locus_list)
        typing(simulation,
               base_fname,
               locus_list,
               genotype_genome,
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
               keep_alignment,
               discordant,
               type_primary_exons,
               remove_low_abundance_alleles,
               display_alleles,
               fastq,
               read_fname,
               alignment_fname,
               [],
               read_len,
               fragment_len,
               threads,
               best_alleles,
               verbose,
               assembly_verbose,
               out_dir)

