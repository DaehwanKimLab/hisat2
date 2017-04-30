#!/usr/bin/env python

import sys, os
import math
import random
from copy import deepcopy
from datetime import datetime


"""
"""
def reverse_complement(seq):
    comp_table = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    rc_seq = ""
    for s in reversed(seq):
        if s in comp_table:
            rc_seq += comp_table[s]
        else:
            rc_seq += s
    return rc_seq


"""
"""
def read_genome(genome_file):
    chr_dic, chr_names, chr_full_names = {}, [], []
    chr_name, chr_full_name, sequence = "", "", ""
    for line in genome_file:
        if line.startswith(">"):
            if chr_name and sequence:
                chr_dic[chr_name] = sequence
                chr_names.append(chr_name)
            chr_full_name = line.strip()[1:]
            chr_name = line.strip().split()[0][1:]
            chr_full_names.append(chr_full_name)
            sequence = ""
        else:
            sequence += line.strip()
    if chr_name and sequence:
        chr_dic[chr_name] = sequence
        chr_names.append(chr_name)
        chr_full_names.append(chr_full_name)
    return chr_dic, chr_names, chr_full_names


"""
"""
def read_allele_sequences(fname):
    allele_seqs = {}
    allele_name, sequence = "", ""
    for line in open(fname):
        if line.startswith(">"):
            if allele_name != "" and allele_name not in allele_seqs:
                allele_seqs[allele_name] = sequence
            allele_name = line.strip()[1:]
            sequence = ""
        else:
            sequence += line.strip()
    if allele_name != "" and allele_name not in allele_seqs:
        allele_seqs[allele_name] = sequence
    return allele_seqs


"""
"""
def read_variants(fname):
    allele_vars = {}
    for line in open(fname):
        var_id, type, allele_name, left, data = line.strip().split()
        left = int(left)
        if type == "deletion":
            data = int(data)
        if allele_name not in allele_vars:
            allele_vars[allele_name] = []
        allele_vars[allele_name].append([left, type, data, var_id])
    return allele_vars


"""
"""
def read_haplotypes(fname):
    allele_haplotypes = {}
    for line in open(fname):
        haplotype_id, allele_name, left, right, vars = line.strip().split()
        vars = vars.split(',')
        left, right = int(left), int(right)
        if allele_name not in allele_haplotypes:
            allele_haplotypes[allele_name] = []
        allele_haplotypes[allele_name].append([left, right, vars])
    return allele_haplotypes


"""
"""
def read_links(fname):
    links = []
    for line in open(fname):
        var_id, allele_names = line.strip().split('\t')
        links.append([var_id, allele_names])
    return links


"""
Compare two variants
"""
def compare_vars(a, b):
    a_pos, a_type, a_data = a[:3]
    b_pos, b_type, b_data = b[:3]

    if a_pos != b_pos:
        return a_pos - b_pos
    if a_type != b_type:
         if a_type == 'I':
             return -1
         elif b_type == 'I':
             return 1
         if a_type == 'S':
             return -1
         else:
             return 1
    if a_data < b_data:
        return -1
    elif a_data > b_data:
        return 1
    else:
        return 0


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
"""
def check_files(fnames):
    for fname in fnames:
        if not os.path.exists(fname):
            return False
    return True
   
    
"""
Download GRCh38 human reference and HISAT2 indexes
"""
def download_genome_and_index(ex_path):
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]
    if not check_files(HISAT2_fnames):
        os.system("wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz; tar xvzf grch38.tar.gz; rm grch38.tar.gz")
        hisat2_inspect = os.path.join(ex_path, "hisat2-inspect")
        os.system("%s grch38/genome > genome.fa" % hisat2_inspect)
        os.system("samtools faidx genome.fa")


"""
"""
def clone_hisatgenotype_database():
    os.system("git clone https://github.com/infphilo/hisatgenotype_db.git")
    
    # Check out one particular revision just to have the same data across multiple computers        
    # revision = "d3b559b34b96ff9e7f0d97476222d8e4cdee63ad" # Revision on November 16, 2016
    # os.system("cd IMGTHLA; git checkout %s; cd .." % revision)


"""
Simulate reads from alleles with headers (>) filled with mapping information.
  For an example, see hisat2_test_HLA_genotyping.py.
"""
def simulate_reads(seq_dic,                       # seq_dic["A"]["A*24:36N"] = "ACGTCCG ..."
                   base_fname,                    # hla, codis, cyp, or so on
                   allele_list,                   # ["A*32:29", "B*07:02:01"]
                   Vars,                          # Vars["A"]["hv326"] = ["single", 604, "C"]
                   Links,
                   simulate_interval = 1,
                   read_len = 100,
                   frag_len = 250,
                   perbase_errorrate = 0.0,
                   perbase_snprate = 0.0,
                   skip_fragment_regions = []):
    reads_1, reads_2 = [], []
    num_pairs = []
    for allele_names in allele_list:
        gene = allele_names[0].split('*')[0]
        num_pairs.append([])

        # Introduce SNPs into allele sequences
        def introduce_snps(seq):
            seq = list(seq)
            for i in range(len(seq)):
                if random.random() * 100 < perbase_snprate:
                    if seq[i] == 'A':
                        alt_bases = ['C', 'G', 'T']
                    elif seq[i] == 'C':
                        alt_bases = ['A', 'G', 'T']
                    elif seq[i] == 'G':
                        alt_bases = ['A', 'C', 'T']
                    else:
                        assert seq[i] == 'T'
                        alt_bases = ['A', 'C', 'G']
                    random.shuffle(alt_bases)
                    alt_base = alt_bases[0]
                    seq[i] = alt_base
            seq = ''.join(seq)
            return seq

        # Simulate reads from two alleles
        def simulate_reads_impl(seq,
                                seq_map,
                                ex_seq_map,
                                ex_seq,
                                ex_desc,
                                simulate_interval = 1,
                                read_len = 100,
                                frag_len = 250,
                                perbase_errorrate = 0.0,
                                skip_fragment_regions = []):
            # Introduce sequencing errors
            def introduce_seq_err(read_seq, pos):
                read_seq = list(read_seq)
                for i in range(read_len):
                    map_pos = seq_map[pos + i]
                    if ex_desc[map_pos] != "":
                        continue
                    if random.random() * 100 < perbase_errorrate:
                        if read_seq[i] == 'A':
                            alt_bases = ['C', 'G', 'T']
                        elif read_seq[i] == 'C':
                            alt_bases = ['A', 'G', 'T']
                        elif read_seq[i] == 'G':
                            alt_bases = ['A', 'C', 'T']
                        else:
                            assert read_seq[i] == 'T'
                            alt_bases = ['A', 'C', 'G']
                        random.shuffle(alt_bases)
                        alt_base = alt_bases[0]
                        read_seq[i] = alt_base
                read_seq = ''.join(read_seq)
                return read_seq                            
                            
            # Get read alignment, e.g., 260|R_483_61M5D38M23D1M_46|S|hv154,3|S|hv162,10|D|hv185,38|D|hv266
            def get_info(read_seq, pos):
                info = "%d_" % (seq_map[pos] + 1)
                total_match, match, sub_match = 0, 0, 0
                var_str = ""
                ins_len, ins_var = 0, ""
                for i in range(pos, pos + read_len):
                    map_i = ex_seq_map[i]
                    assert ex_seq[map_i] != 'D'
                    total_match += 1
                    match += 1
                    if ex_seq[map_i] == 'I':
                        if ins_var != "":
                            assert ins_var == ex_desc[map_i]
                        ins_var = ex_desc[map_i]
                        ins_len += 1
                    elif ins_var != "":
                        if var_str != "":
                            var_str += ','
                        var_str += ("%s|I|%s" % (sub_match, ins_var))
                        ins_len, ins_var = 0, ""
                        sub_match = 0
                    if ex_seq[map_i] != 'I':
                        if ex_desc[map_i] != "" or read_seq[i-pos] != ex_seq[map_i]:
                            if var_str != "":
                                var_str += ','
                            var_str += ("%d|S|%s" % (sub_match, ex_desc[map_i] if ex_desc[map_i] != "" else "unknown"))
                            sub_match = 0
                        else:
                            sub_match += 1
                    if i + 1 < pos + read_len and ex_seq[map_i+1] == 'D':
                        assert match > 0
                        info += ("%dM" % match)
                        match = 0
                        del_len = 1
                        while map_i + 1 + del_len < len(ex_seq):
                            if ex_seq[map_i + 1 + del_len] != 'D':
                                break
                            del_len += 1
                        info += ("%dD" % del_len)
                        if var_str != "":
                            var_str += ','
                        var_str += ("%s|D|%s" % (sub_match, ex_desc[map_i + 1]))
                        sub_match = 0
                assert match > 0
                info += ("%dM" % match)
                assert total_match == read_len
                if var_str:
                    info += "_"
                    info += var_str                
                return info
                
            comp_table = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
            reads_1, reads_2 = [], []
            for i in range(0, len(seq) - frag_len + 1, simulate_interval):
                if len(skip_fragment_regions) > 0:
                    skip = False
                    for skip_left, skip_right in skip_fragment_regions:
                        if i <= skip_right and i + frag_len > skip_left:
                            skip = True
                            break
                    if skip:
                        continue
                        
                pos1 = i
                seq1 = seq[pos1:pos1+read_len]
                if perbase_errorrate > 0.0:
                    seq1 = introduce_seq_err(seq1, pos1)
                info1 = get_info(seq1, pos1)
                reads_1.append([seq1, info1])
                
                pos2 = i + frag_len - read_len
                seq2 = seq[pos2:pos2+read_len]
                if perbase_errorrate > 0.0:
                    seq2 = introduce_seq_err(seq2, pos2)                
                info2 = get_info(seq2, pos2)
                tmp_read_2 = reversed(seq2)
                read_2 = ""
                for s in tmp_read_2:
                    if s in comp_table:
                        read_2 += comp_table[s]
                    else:
                        read_2 += s
                reads_2.append([read_2, info2])
            return reads_1, reads_2

        # for each allele in a list of alleles such as ['A*32:29', 'B*07:02:01']
        for allele_name in allele_names:
            allele_seq = seq_dic[gene][allele_name]
            backbone_seq = seq_dic[gene]["%s*BACKBONE" % gene]
            allele_ex_seq = list(backbone_seq)
            allele_ex_desc = [''] * len(allele_ex_seq)
            allele_seq_map = [i for i in range(len(allele_seq))]
            allele_ex_seq_map = [i for i in range(len(allele_seq))]

            if perbase_snprate > 0:
                HLA_seq = introduce_snps(allele_seq)

            # Extract variants included in each allele
            var_ids = []
            for var_id, allele_list in Links.items():
                if allele_name in allele_list:
                    var_ids.append(var_id)

            def var_cmp(a, b):
                assert a.startswith("hv") and b.startswith("hv")
                return int(a[2:]) - int(b[2:])
            var_ids = sorted(var_ids, cmp=var_cmp)

            # Build annotated sequence for the allele w.r.t backbone sequence
            add_pos = 0
            for var_id in var_ids:
                var_type, var_pos, var_data = Vars[gene][var_id]
                var_pos += add_pos
                if var_type == "single":
                    allele_ex_seq[var_pos] = var_data
                    allele_ex_desc[var_pos] = var_id
                elif var_type == "deletion":
                    del_len = int(var_data)
                    assert var_pos + del_len <= len(allele_ex_seq)
                    allele_ex_seq[var_pos:var_pos+del_len] = ['D'] * del_len
                    allele_ex_desc[var_pos:var_pos+del_len] = [var_id] * del_len
                else:
                    assert var_type == "insertion"
                    ins_len = len(var_data)
                    allele_ex_seq = allele_ex_seq[:var_pos] + (['I'] * ins_len) + allele_ex_seq[var_pos:]
                    allele_ex_desc = allele_ex_desc[:var_pos] + ([var_id] * ins_len) + allele_ex_desc[var_pos:]
                    add_pos += ins_len
            allele_ex_seq = ''.join(allele_ex_seq)
            assert len(backbone_seq) + add_pos == len(allele_ex_seq)            

            # Build mapping from the allele to the annotated sequence
            prev_j, minus_pos = 0, 0
            for i in range(len(allele_seq)):
                for j in range(prev_j, len(allele_ex_seq)):
                    if allele_ex_seq[j] != 'D':
                        if allele_ex_seq[j] == 'I':
                            minus_pos += 1
                        break
                allele_seq_map[i] = j - minus_pos
                allele_ex_seq_map[i] = j
                prev_j = j + 1

            # DK - debugging purposes
            """
            for t in range(0, len(allele_ex_seq), 100):
                print t, allele_ex_seq[t:t+100]
                print t, '-'.join(allele_ex_desc[t:t+100])
                print t, allele_seq_map[t:t+100]
            print "allele_seq length:", len(allele_seq)
            print len(allele_ex_seq), "vs.", len(seq_dic[gene]["A*BACKBONE"]), "vs.", len(allele_seq_map)
            print allele_ex_seq[1943:1946]
            print allele_ex_desc[1943:1946]
            sys.exit(1)
            """
            
            tmp_reads_1, tmp_reads_2 = simulate_reads_impl(allele_seq,
                                                           allele_seq_map,
                                                           allele_ex_seq_map,
                                                           allele_ex_seq,
                                                           allele_ex_desc,
                                                           simulate_interval,
                                                           read_len,
                                                           frag_len,
                                                           perbase_errorrate,
                                                           skip_fragment_regions)
            reads_1 += tmp_reads_1
            reads_2 += tmp_reads_2
            num_pairs[-1].append(len(tmp_reads_1))

    # Write reads into a FASTA file
    def write_reads(reads, idx):
        read_file = open('%s_input_%d.fa' % (base_fname, idx), 'w')
        for read_i in range(len(reads)):
            print >> read_file, ">%d|%s_%s" % (read_i + 1, "LR"[idx-1], reads[read_i][1])
            print >> read_file, reads[read_i][0]
        read_file.close()
    write_reads(reads_1, 1)
    write_reads(reads_2, 2)

    return num_pairs


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
                               del_len,
                               other_del_len,
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

        if del_len == other_del_len:
            return
                
        var_type, var_pos, var_data = Vars[var_id]
        if left: # Look in left direction
            if var_j < 0:
                return
            j_pos, j_id = Var_list[var_j]
            alt_del = []
            if var_id != j_id and j_pos < var_pos + del_len:
                prev_latest_pos = latest_pos
                # Check bases between SNPs
                while latest_pos - max(0, del_len - other_del_len) > 0:
                    if ref_seq[latest_pos - 1] != ref_seq[latest_pos - 1 - del_len + other_del_len]:
                        break
                    latest_pos -= 1
                    add_alt(Alts, alt_list, var_id, str(latest_pos))
                if latest_pos - 1 > j_pos:
                    return
                j_type, _, j_data = Vars[j_id]
                if j_type == "deletion":
                    j_del_len = int(j_data)
                if j_type == "single" and j_pos == latest_pos - 1:
                    j_cmp_pos = j_pos - del_len + other_del_len
                    if debug:
                        print Vars[j_id]
                        print j_pos, ref_seq[j_pos]
                        print j_cmp_pos, ref_seq[j_cmp_pos]
                    if j_data == ref_seq[j_cmp_pos]:
                        add_alt(Alts, alt_list, var_id, j_id)
                        latest_pos = j_pos
                elif j_type == "deletion" and j_pos + j_del_len - 1 == prev_latest_pos - 1:
                    alt_list2 = alt_list[:] + [j_id]
                    latest_pos2 = j_pos
                    alt_del = [alt_list2, latest_pos2]
                
            get_alternatives_recur(ref_seq,
                                   Vars,
                                   Var_list,
                                   Alts,
                                   var_id,
                                   del_len,
                                   other_del_len,
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
                                       del_len,
                                       other_del_len + j_del_len,
                                       left,
                                       alt_list2,
                                       var_j - 1,
                                       latest_pos2,
                                       debug)
                # Remove this Deletion if not supported by additional bases?
                assert alt_idx < len(Alts[var_id])
                if Alts[var_id][alt_idx][-1] == j_id:
                    Alts[var_id] = Alts[var_id][:alt_idx] + Alts[var_id][alt_idx+1:]
              
        else: # Look in right direction
            if var_j >= len(Var_list):
                return
            j_pos, j_id = Var_list[var_j]
            alt_del = []
            if var_id != j_id and j_pos >= var_pos:
                # Check bases between SNPs
                prev_latest_pos = latest_pos
                while latest_pos + 1 + max(0, del_len - other_del_len) < len(ref_seq):
                    if ref_seq[latest_pos + 1] != ref_seq[latest_pos + 1 + del_len - other_del_len]:
                        break

                    # DK - debugging purposes
                    if debug:
                        pos2_ = latest_pos + 1 + del_len - other_del_len
                        print "DK: latest_pos:", latest_pos + 1, pos2_
                        print "DK: var_pos:", var_pos, "del_len:", del_len, "other_del_len:", other_del_len
                        print "DK:", ref_seq[latest_pos + 1], ref_seq[pos2_]
                    
                    latest_pos += 1
                    add_alt(Alts, alt_list, var_id, str(latest_pos))

                if latest_pos + 1 < j_pos:
                    return               
                
                j_type, _, j_data = Vars[j_id]
                if j_type == "single" and j_pos == latest_pos + 1:
                    j_cmp_pos = j_pos + del_len - other_del_len
                    if debug:
                        print Vars[j_id]
                        print j_pos, ref_seq[j_pos]
                        print j_cmp_pos, ref_seq[j_cmp_pos]

                    if j_data == ref_seq[j_cmp_pos]:
                        add_alt(Alts, alt_list, var_id, j_id)
                        latest_pos = j_pos
                elif j_type == "deletion" and j_pos == prev_latest_pos + 1:                        
                    j_del_len = int(j_data)
                    alt_list2 = alt_list[:] + [j_id]
                    latest_pos2 = j_pos + j_del_len - 1
                    alt_del = [alt_list2, latest_pos2]

            get_alternatives_recur(ref_seq,
                                   Vars,
                                   Var_list,
                                   Alts,
                                   var_id,
                                   del_len,
                                   other_del_len,
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
                                       del_len,
                                       other_del_len + j_del_len,
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
        debug = (var_id == "hv454a")
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
                                   del_len,
                                   0,
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
                               del_len,
                               0,
                               False, # right
                               alt_list,
                               var_j,
                               latest_pos,
                               debug)

        if debug:
            print "DK :-)"
            sys.exit(1)

    def assert_print_alts(Alts, dir):
        def get_seq_pos(alt_list):
            seq = ""
            seq_left, seq_right = -1, -1
            for i in range(len(alt_list)):
                alt = alt_list[i]
                if alt.isdigit():
                    assert i + 1 == len(alt_list)
                    if dir == "left":
                        if i == 0:
                            seq = alt
                            break

                        alt = int(alt)
                        seq = ref_seq[seq_left-alt+1:seq_left+1] + seq
                        seq_left -= alt
                    else:
                        alt = int(alt)
                        seq += ref_seq[seq_right:seq_right+alt]
                        seq_right += alt
                    break

                var_type, var_pos, var_data = Vars[alt]
                if dir == "left" and var_type == "deletion":
                    var_pos = var_pos + int(var_data) - 1

                if i == 0:
                    if dir == "left":
                        seq_left, seq_right = var_pos, var_pos
                    else:
                        seq_left, seq_right = var_pos, var_pos
                       
                if dir == "left":
                    assert seq_left >= var_pos
                    if i > 0:
                        seq = ref_seq[var_pos+1:seq_left+1] + seq
                    if var_type == "single":
                        seq = var_data + seq
                        seq_left = var_pos - 1
                    elif var_type == "deletion":
                        seq_left = var_pos - int(var_data)
                    else:
                        assert var_type == "insertion"
                        seq = var_data + seq
                else:
                    assert seq_right <= var_pos
                    if i > 0:
                        seq += ref_seq[seq_right:var_pos]
                    if var_type == "single":
                        seq += var_data
                        seq_right = var_pos + 1
                    elif var_type == "deletion":
                        seq_right = var_pos + int(var_data)
                    else:
                        assert var_type == "insertion"
                        seq += var_data
                        
            return seq, seq_left, seq_right
        
        for alt_list1, alt_list2 in Alts.items():
            if verbose >= 2: print >> sys.stderr, "\t", dir, ":", alt_list1, alt_list2
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
            if verbose >= 2: print >> sys.stderr, out_str

            for i in range(len(alt_list2)):
                alt_list3 = alt_list2[i]
                seq1, seq1_left, seq1_right = get_seq_pos(alt_list1)
                seq2, seq2_left, seq2_right = get_seq_pos(alt_list3)
                if seq1.isdigit():
                    assert not seq2.isdigit()
                    seq1_left, seq1_right = seq2_right - int(seq1), seq2_right
                    seq1 = ref_seq[seq1_left+1:seq1_right+1]
                elif seq2.isdigit():
                    seq2_left, seq2_right = seq1_right - int(seq2), seq1_right
                    seq2 = ref_seq[seq2_left+1:seq2_right+1]
                    
                if dir == "left":
                    if seq1_right < seq2_right:
                        seq1 += ref_seq[seq1_right+1:seq2_right+1]
                    elif seq2_right < seq1_right:
                        seq2 += ref_seq[seq2_right+1:seq1_right+1]
                else:
                    if seq1_left < seq2_left:
                        seq2 = ref_seq[seq1_left:seq2_left] + seq2
                    elif seq2_left < seq1_left:
                        seq1 = ref_seq[seq2_left:seq1_left] + seq1
                seq1_len, seq2_len = len(seq1), len(seq2)
                if seq1_len != seq2_len:
                    len_diff = abs(seq1_len - seq2_len)
                    if dir == "left":
                        if seq1_len < seq2_len:
                            seq1 = ref_seq[seq1_left-len_diff+1:seq1_left+1] + seq1
                        else:
                            seq2 = ref_seq[seq2_left-len_diff+1:seq2_left+1] + seq2
                    else:
                        if seq1_len < seq2_len:
                            seq1 += ref_seq[seq1_right:seq1_right+len_diff]
                        else:
                            seq2 += ref_seq[seq2_right:seq2_right+len_diff]
                if verbose >= 3:
                    print >> sys.stderr, "\t\t", alt_list1, alt_list3
                    print >> sys.stderr, "\t\t\t", seq1, seq1_left, seq1_right
                    print >> sys.stderr, "\t\t\t", seq2, seq2_left, seq2_right
                assert seq1 == seq2            
            
    assert_print_alts(Alts_left, "left")
    assert_print_alts(Alts_right, "right")

    return Alts_left, Alts_right


"""
Identify ambigious differences that may account for other alleles,
  given a list of differences (cmp_list) between a read and a potential allele   
"""
def identify_ambigious_diffs(Vars,
                             Alts_left,
                             Alts_right,
                             cmp_list,
                             verbose,
                             debug = False):
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

            # DK - debugging purposes
            if debug:
                print "DK: var_id:", var_id
                print "DK: cmp_list:", cmp_list
                print "DK: cmp_right:", cmp_right
                # sys.exit(1)
    
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

                # DK - debugging purposes
                if debug:
                    print "DK: id_str:", id_str
                
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
Identify alternative haplotypes
   insertions are not considered...

"""
def get_alternatives2(ref_seq,     # GATAACTAGATACATGAGATAGATTTGATAGATAGATAGATACATACATACATACATACATACAGGATAGATAACTAGG...
                      allele_vars, # {'VWA*20(22)': ['hv231', 'hv245'], "VWA*16(18')": ['hv235', 'hv250', 'hv256'], ...}
                      Vars,        # {'hv241': ['deletion', 529, '52'], 'hv240': ['deletion', 529, '48'], ... }
                      Var_list,    # [[529, 'hv230'], [529, 'hv231'], [529, 'hv232'], [529, 'hv233'], ...]
                      verbose):
    haplotype_alts_left, haplotype_alts_right = {}, {}
    second_order_haplotypes = set()
    for allele_name, vars in allele_vars.items():
        for v in range(len(vars) - 1):
            ht = vars[v] + "-" + vars[v+1]
            second_order_haplotypes.add(ht)

    rev_Var_list = []
    for _, var_id in Var_list:
        var_type, var_pos, var_data = Vars[var_id]
        if var_type == "deletion":
            var_pos = var_pos + int(var_data) - 1
        elif var_type == "insertion":
            var_pos += 1
        rev_Var_list.append([var_pos, var_id])
    rev_Var_list = sorted(rev_Var_list, cmp=lambda a, b: a[0] - b[0])

    def nextbases(haplotype,
                  left = True,
                  exclude_list = []):
        if left:
            pos = int(haplotype[0]) - 1
        else:
            pos = haplotype[-1] + 1
        assert pos >= 0 and pos < len(ref_seq)

        if left:
            bases = [[[pos] + haplotype[1:], ref_seq[pos]]]
            prev_id = None
            if len(haplotype) > 2:
                prev_id = haplotype[1]        

            var_i = lower_bound(rev_Var_list, pos + 1)
            for var_j in reversed(range(0, var_i)):
                _, var_id = rev_Var_list[var_j]
                var_type, var_pos, var_data = Vars[var_id]
                if var_type == "deletion":
                    if var_pos == 0:
                        continue
                    var_pos = var_pos + int(var_data) - 1
                if var_pos > pos:
                    continue
                if var_pos < pos:
                    break
                if var_id in exclude_list:
                    continue
                if prev_id:
                    second_ht = var_id + "-" + prev_id
                    if second_ht not in second_order_haplotypes:
                        continue

                if var_type == "single":
                    bases.append([[var_pos, var_id] + haplotype[1:], var_data])
                elif var_type == "deletion":
                    bases2 = nextbases([var_pos - int(var_data) + 1, var_id] + haplotype[1:],
                                       left,
                                       exclude_list)
                    bases += bases2
                else:
                    assert var_type == "insertion"
        else:
            bases = [[haplotype[:-1] + [pos], ref_seq[pos]]]
            prev_id = None
            if len(haplotype) > 2:
                prev_id = haplotype[-2]       

            var_i = lower_bound(Var_list, pos)
            for var_j in range(var_i, len(Var_list)):
                _, var_id = Var_list[var_j]
                var_type, var_pos, var_data = Vars[var_id]
                if var_pos < pos:
                    continue
                if var_pos > pos:
                    break
                if var_id in exclude_list:
                    continue
                if prev_id:
                    second_ht = prev_id + "-" + var_id
                    if second_ht not in second_order_haplotypes:
                        continue

                if var_type == "single":
                    bases.append([haplotype[:-1] + [var_id, var_pos], var_data])
                elif var_type == "deletion":
                    bases2 = nextbases(haplotype[:-1] + [var_id, var_pos + int(var_data) - 1],
                                       left,
                                       exclude_list)
                    bases += bases2
                else:
                    assert var_type == "insertion"

        return bases

    def get_haplotype_seq(haplotype):
        seq = ""
        pos = int(haplotype[0])
        for i in range(1, len(haplotype) - 1):
            var_id = haplotype[i]
            var_type, var_pos, var_data = Vars[var_id]
            if pos < var_pos:
                seq += ref_seq[pos:var_pos]
            if var_type == "single":
                seq += var_data
                pos = var_pos + 1
            elif var_type == "deletion":
                pos = var_pos + int(var_data)
            else:
                assert var_type == "insertion"
                seq += var_data
                pos = var_pos
            
        last_pos = int(haplotype[-1]) + 1
        assert pos <= last_pos
        if pos < last_pos:
            seq += ref_seq[pos:last_pos]                
        return seq

    def get_alternative_recur(var_orig_id,
                              haplotype,
                              haplotype_alt,
                              left = True,
                              dep = 0):
        bases1 = nextbases(haplotype,
                           left)
        bases2 = nextbases(haplotype_alt,
                           left,
                           [var_orig_id]) # exclude

        found = False
        for base1 in bases1:
            next_haplotype, bp = base1
            for base2 in bases2:
                next_haplotype_alt, bp2 = base2
                if left:
                    if next_haplotype[0] == next_haplotype_alt[0]:
                        continue
                else:
                    if next_haplotype[-1] == next_haplotype_alt[-1]:
                        continue
                if bp != bp2:
                    continue

                found = True
                get_alternative_recur(var_orig_id,
                                      next_haplotype,
                                      next_haplotype_alt,
                                      left,
                                      dep + 1)            
    
        if dep > 0:
            if not found:
                def to_haplotype_str(haplotype):
                    if len(haplotype) <= 2:
                        haplotype = "%d-%d" % (haplotype[0], haplotype[1])
                    else:
                        haplotype = "%d-%s-%d" % (haplotype[0], '-'.join(haplotype[1:-1]), haplotype[-1])
                    return haplotype

                haplotype, haplotype_alt = to_haplotype_str(haplotype), to_haplotype_str(haplotype_alt)
                haplotype_alts = haplotype_alts_left if left else haplotype_alts_right
                if haplotype not in haplotype_alts:
                    haplotype_alts[haplotype] = set()
                assert haplotype_alt not in haplotype_alts[haplotype]
                haplotype_alts[haplotype].add(haplotype_alt)



    # Search alternative haplotypes in both left and right directions
    for var_i in range(len(Var_list)):
        _, var_id = Var_list[var_i]
        var_type, var_pos, var_data = Vars[var_id]
        if var_pos == 0:
            continue
        if var_type != "deletion":
            continue
        del_len = int(var_data)
        if var_pos + del_len >= len(ref_seq):
            continue

        # Left direction
        get_alternative_recur(var_id,
                              [var_pos, var_id, var_pos + del_len - 1],
                              [var_pos + del_len, var_pos + del_len - 1])

        # Right direction    
        get_alternative_recur(var_id,
                              [var_pos, var_id, var_pos + del_len - 1],
                              [var_pos, var_pos - 1],
                              False)


    # Print alternative haplotypes / Sanity check
    def print_haplotype_alts(haplotype_alts):
        for haplotype, haplotype_set in haplotype_alts.items():
            if verbose: print "%s:" % haplotype, haplotype_set
            haplotype_seq = get_haplotype_seq(haplotype.split('-'))
            for haplotype_alt in haplotype_set:
                haplotype_alt_seq = get_haplotype_seq(haplotype_alt.split('-'))
                assert haplotype_seq == haplotype_alt_seq            

    if verbose: print "number of left haplotypes:", len(haplotype_alts_left)
    print_haplotype_alts(haplotype_alts_left)
    if verbose: print "number of right haplotypes:", len(haplotype_alts_right)
    print_haplotype_alts(haplotype_alts_right)
    
    return haplotype_alts_left, haplotype_alts_right



"""
Identify alternative haplotypes
   insertions are not considered...

"""
def get_alternatives3(ref_seq,     # GATAACTAGATACATGAGATAGATTTGATAGATAGATAGATACATACATACATACATACATACAGGATAGATAACTAGG...
                      allele_vars, # {'VWA*20(22)': ['hv231', 'hv245'], "VWA*16(18')": ['hv235', 'hv250', 'hv256'], ...}
                      Vars,        # {'hv241': ['deletion', 529, '52'], 'hv240': ['deletion', 529, '48'], ... }
                      Var_list,    # [[529, 'hv230'], [529, 'hv231'], [529, 'hv232'], [529, 'hv233'], ...]
                      verbose):
    haplotype_alts_left, haplotype_alts_right = {}, {}

    allele_vars_list = []
    for allele_name, vars in allele_vars.items():
        allele_vars_list.append([allele_name, vars])
    allele_vars = allele_vars_list

    def print_allele(vars, var_i, left = True):
        if left:
            pos = -1
            for i in range(0, var_i + 1):
                var_id = vars[i]
                var_type, var_pos, var_data = Vars[var_id]
                if i == 0:
                    assert var_pos >= 0
                    print ref_seq[var_pos-50:var_pos],
                else:
                    assert pos > 0
                    if pos < var_pos:
                        print ref_seq[var_pos+1:pos+1],
                if var_type in ["single", "insertion"]:
                    pos = var_pos + 1
                else:
                    assert var_type == "deletion"
                    pos = var_pos + int(var_data)

                print Vars[var_id],
                if i == var_i:
                    print
        else:
            pos = -1
            for i in range(var_i, len(vars)):
                var_id = vars[i]
                var_type, var_pos, var_data = Vars[var_id]
                if i > var_i:
                    assert pos > 0
                    if pos < var_pos:
                        print ref_seq[pos:var_pos],
                if var_type in ["single", "insertion"]:
                    pos = var_pos + 1
                else:
                    assert var_type == "deletion"
                    pos = var_pos + int(var_data)

                print Vars[var_id],

                if i + 1 == len(vars):
                    assert pos >= 0
                    print ref_seq[pos:pos+50]

    def next_base(vars, var_i, pos, left = True):
        base = ""
        if left:
            assert pos >= 0
            while var_i >= 0:
                var_type, var_pos, var_data = Vars[vars[var_i]]
                if var_type == "deletion":
                    var_pos = var_pos + int(var_data) - 1
                if pos > var_pos:
                    break
                assert pos == var_pos
                var_i -= 1
                if var_type == "single":
                    base = var_data
                    pos -= 1
                    break
                elif var_type == "insertion":
                    continue
                else:
                    assert var_type == "deletion"
                    pos -= int(var_data)
                    continue
        else:
            assert pos < len(ref_seq)
            while var_i < len(vars):
                var_type, var_pos, var_data = Vars[vars[var_i]]
                if pos < var_pos:
                    break
                assert pos == var_pos
                var_i += 1
                if var_type == "single":
                    base = var_data
                    pos += 1
                    break
                elif var_type == "insertion":
                    continue
                else:
                    assert var_type == "deletion"
                    pos += int(var_data)
                    continue

        if base == "" and \
           pos >= 0 and \
           pos < len(ref_seq):
            base = ref_seq[pos]
            pos = pos - 1 if left else pos + 1

        return var_i, pos, base

    def get_seq(vars, var_i, pos, pos_end, left = True):
        haplotype, seq = str(pos), ""
        if left:
            for i in reversed(range(0, var_i + 1)):
                var_id = vars[i]
                var_type, var_pos, var_data = Vars[var_id]
                if var_type == "single":
                    pos_next = var_pos - 1
                elif var_type == "deletion":
                    pos_next = var_pos - 1
                    var_pos = var_pos + int(var_data) - 1
                else:
                    assert var_type == "insertion"
                    continue
                if pos_next < pos_end:
                    break

                if pos > var_pos:
                    seq = ref_seq[var_pos+1:pos+1] + seq
                haplotype = var_id + "-" + haplotype
                if var_type == "single":
                    seq = var_data + seq
                elif var_type == "insertion":
                    continue
                else:
                    assert var_type == "deletion"
                pos = pos_next

            seq = ref_seq[pos_end+1:pos+1] + seq
            haplotype = ("%d-" % (pos_end + 1)) + haplotype
        else:
            for i in range(var_i, len(vars)):
                var_id = vars[i]
                var_type, var_pos, var_data = Vars[var_id]
                if var_type == "single":
                    pos_next = var_pos + 1
                elif var_type == "deletion":
                    pos_next = var_pos + int(var_data)
                else:
                    assert var_type == "insertion"
                    continue
                if pos_next > pos_end:
                    break
                
                if pos < var_pos:
                    seq += ref_seq[pos:var_pos]
                haplotype = haplotype + "-" + var_id
                if var_type == "single":
                    seq += var_data
                elif var_type == "insertion":
                    continue
                else:
                    assert var_type == "deletion"
                pos = pos_next

            seq += ref_seq[pos:pos_end]
            haplotype += ("-%d" % (pos_end - 1))

        return haplotype, seq

    def is_searched(haplotype_alts, vars, var_i, left = True):
        assert var_i < len(vars)
        var_id = vars[var_i]
        var_type, var_pos, var_data = Vars[var_id]
        assert var_type == "deletion"
        del_len = int(var_data)
        if var_pos not in haplotype_alts:
            return False
        haplotypes = haplotype_alts[var_pos]

        for haplotype in haplotypes:
            cmp_vars, cmp_var_i, cmp_var_i2 = haplotype
            num_vars = cmp_var_i2 - cmp_var_i + 1
            if left:
                None
            else:
                if cmp_vars[cmp_var_i:cmp_var_i+num_vars] == \
                   vars[var_i:var_i+num_vars]:
                    return True
        return False

    _haplotype_alts = {}

    # DK - initial - left direction
    for allele_i in range(len(allele_vars)):
        allele_name, vars = allele_vars[allele_i]
        var_j_list = [sys.maxint for _ in range(len(allele_vars))]
        for var_i in range(len(vars)):
            var_id = vars[var_i]
            var_type, var_pos, var_data = Vars[var_id]
            if var_type != "deletion":
                continue
            del_len = int(var_data)
            var_pos = var_pos + del_len - 1
            if is_searched(_haplotype_alts, vars, var_i):
                continue

            for allele_j in range(len(allele_vars)):
                if allele_j == allele_i:
                    continue
                allele_name2, vars2 = allele_vars[allele_j]
                if var_j_list[allele_j] > len(vars2):
                    var_j_list[allele_j] = len(vars2) - 1

                skip = False
                var_j = var_j_list[allele_j]
                while var_j >= 0:
                    var_id2 = vars2[var_j]
                    if var_id2 == var_id:
                        skip = True
                        break
                    var_type2, var_pos2, var_data2 = Vars[var_id2]
                    if var_type2 == "deletion":
                        var_pos2 = var_pos2 + int(var_data2) - 1
                    if var_pos2 <= var_pos:
                        break
                    var_j_list[allele_j] = var_j
                    var_j -= 1
                    if var_type2 in ["single", "insertion"]:
                        continue
                    assert var_type2 == "deletion"
                    del_len2 = int(var_data2)
                    if var_pos2 - del_len2 >= var_pos:
                        continue
                    else:
                        skip = True
                        break
                if skip:
                    continue

                # DK - debugging purposes
                # print allele_name,  ":",; print_allele(vars, var_i)
                # print allele_name2, ":",; print_allele(vars2, var_j)

                seq_len = 0
                init_pos1 = var_pos
                pos1, pos2 = var_pos, var_pos
                var_i2, var_j2 = var_i, var_j
                while pos1 >= 0 and pos2 >= 0:
                    var_i2_, pos1_, base1 = next_base(vars,  var_i2, pos1)
                    var_j2_, pos2_, base2 = next_base(vars2, var_j2, pos2)
                    if base1 != base2:
                        break

                    var_i2, pos1 = var_i2_, pos1_
                    var_j2, pos2 = var_j2_, pos2_

                if pos1 < init_pos1:
                    haplotype1, seq1 = get_seq(vars,  var_i, init_pos1, pos1)
                    haplotype2, seq2 = get_seq(vars2, var_j, init_pos1, pos2)
                    assert seq1 == seq2

    _haplotype_alts = {}

    # DK - initial - right direction
    max_time = None
    for allele_i in range(len(allele_vars)):
        # DK - debugging purposes
        start_time = datetime.now()
        allele_name, vars = allele_vars[allele_i]
        var_j_list = [0 for _ in range(len(allele_vars))]
        for var_i in range(len(vars)):
            var_id = vars[var_i]
            var_type, var_pos, var_data = Vars[var_id]
            if var_type != "deletion":
                continue
            if var_pos <= 0:
                continue
            del_len = int(var_data)

            # DK - debugging purposes
            if allele_name == "A*24:359N" and var_pos == 1600:
                continue
            
            if is_searched(_haplotype_alts, vars, var_i, False):
                continue

            max_var_i2 = var_i
            for allele_j in range(len(allele_vars)):
                if allele_j == allele_i:
                    continue
                allele_name2, vars2 = allele_vars[allele_j]

                skip = False
                var_j = var_j_list[allele_j]
                var_j = 0
                while var_j < len(vars2):
                    var_id2 = vars2[var_j]
                    if var_id2 == var_id:
                        skip = True
                        break
                    var_type2, var_pos2, var_data2 = Vars[var_id2]
                    if var_pos2 >= var_pos:
                        break
                    var_j_list[allele_j] = var_j
                    var_j += 1
                    if var_type2 in ["single", "insertion"]:
                        continue
                    assert var_type2 == "deletion"
                    del_len2 = int(var_data2)
                    if var_pos2 + del_len2 > var_pos:
                        skip = True
                        break
                if skip:
                    continue

                # DK - debugging purposes
                # print allele_name,  ":",; print_allele(vars, var_i, False)
                # print allele_name2, ":",; print_allele(vars2, var_j, False)
                            
                seq_len = 0
                init_pos1 = var_pos
                pos1, pos2 = var_pos, var_pos
                var_i2, var_j2 = var_i, var_j
                while pos1 < len(ref_seq) and pos2 < len(ref_seq):
                    var_i2_, pos1_, base1 = next_base(vars,  var_i2, pos1, False)
                    var_j2_, pos2_, base2 = next_base(vars2, var_j2, pos2, False)
                    if base1 != base2:
                        break

                    var_i2, pos1 = var_i2_, pos1_
                    var_j2, pos2 = var_j2_, pos2_

                if init_pos1 < pos1:
                    haplotype1, seq1 = get_seq(vars,  var_i, init_pos1, pos1, False)
                    haplotype2, seq2 = get_seq(vars2, var_j, init_pos1, pos2, False)

                    if allele_name == "A*24:359N":
                        print "R ", allele_name, "vs.", allele_name2
                        print "\t", allele_name, haplotype1
                        print "\t", allele_name2, haplotype2

                    # DK - debugging purposes
                    if seq1 != seq2:
                        print "\t", seq1, allele_name, pos1, ref_seq[pos1-4:pos1+1]
                        print "\t", seq2, allele_name2, pos2, ref_seq[pos2-4:pos2+1]
                    
                    assert seq1 == seq2

                    if max_var_i2 < var_i2:
                        max_var_i2 = var_i2

                    if haplotype1 not in haplotype_alts_right:
                        haplotype_alts_right[haplotype1] = set()
                    haplotype_alts_right[haplotype1].add(haplotype2)

            if max_var_i2 > var_i:
                if var_pos not in _haplotype_alts:
                    _haplotype_alts[var_pos] = []
                _haplotype_alts[var_pos].append([vars, var_i, max_var_i2])


        tmp_time = datetime.now() - start_time
        if not max_time or tmp_time > max_time or allele_name == "A*24:359N":
            max_time = tmp_time
            print >> sys.stderr, max_time, allele_name

            
            # DK - debugging purposes
            # sys.exit(1)

    # DK - debugging purposes
    print >> sys.stderr, "num of left haplotypes:", len(haplotype_alts_left)
    for haplotype, haplotype_set in haplotype_alts_left.items():
        print >> sys.stderr,  haplotype, haplotype_set
    print >> sys.stderr, "num of right haplotypes:", len(haplotype_alts_right)
    for haplotype, haplotype_set in haplotype_alts_right.items():
        print >> sys.stderr,  haplotype, haplotype_set

    return haplotype_alts_left, haplotype_alts_right

