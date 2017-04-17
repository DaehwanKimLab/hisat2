#!/usr/bin/env python

import sys, os
import math
import random
from copy import deepcopy


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

