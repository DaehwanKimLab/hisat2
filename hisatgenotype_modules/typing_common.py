#!/usr/bin/env python

import sys, os, subprocess
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
            query_name = "%d|%s_%s" % (read_i + 1, "LR"[idx-1], reads[read_i][1])
            if len(query_name) > 254:
                query_name = query_name[:254]
            print >> read_file, ">%s" % query_name
            print >> read_file, reads[read_i][0]
        read_file.close()
    write_reads(reads_1, 1)
    write_reads(reads_2, 2)

    return num_pairs


"""
Align reads, and sort the alignments into a BAM file
"""
def align_reads(aligner,
                simulation,
                base_fname,
                index_type,
                read_fname,
                fastq,
                threads,
                out_fname,
                verbose):
    if aligner == "hisat2":
        aligner_cmd = [aligner, "--mm"]
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
            aligner_cmd += ["--haplotype"]
            if base_fname == "codis":
                aligner_cmd += ["--enable-codis"]
        aligner_cmd += ["-x", "%s.%s" % (base_fname, index_type)]
    elif aligner == "bowtie2":
        aligner_cmd = [aligner,
                       "--no-unal",
                       "-k", "10",
                       "-x", base_fname]
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

   INPUT: see the function's parameters below
   OUPUT: 529-hv8-hv22-606: set(['529-hv13-570', '529-hv4-hv18-590', '529-hv2-hv16-582'])
          529-hv3-hv17-598: set(['529-hv6-hv21-hv26-610'])
"""
def get_alternatives(ref_seq,     # GATAACTAGATACATGAGATAGATTTGATAGATAGATAGATACATACATACATACATACATACAGGATAGATAACTAGG...
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
        if pos < 0 or pos >= len(ref_seq):
            return []

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
                if bp != bp2:
                    continue

                # Todo: implement a routine to handle haplotypes ending with the same coordinate
                if left:
                    left1, left2 = int(next_haplotype[0]), int(next_haplotype_alt[0])
                    if left1 == left2:
                        continue
                else:
                    right1, right2 = int(next_haplotype[-1]), int(next_haplotype_alt[-1])
                    if right1 == right2:
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
                haplotype_alts[haplotype].add(haplotype_alt)

                if haplotype_alt not in haplotype_alts:
                    haplotype_alts[haplotype_alt] = set()
                haplotype_alts[haplotype_alt].add(haplotype)

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
            if verbose: print "\t%s:" % haplotype, haplotype_set
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
Identify ambigious differences that may account for other alleles,
  given a list of differences (cmp_list) between a read and a potential allele   
"""
def identify_ambigious_diffs(ref_seq,
                             Vars,
                             Alts_left,
                             Alts_right,
                             Alts_left_list,
                             Alts_right_list,
                             cmp_list,
                             verbose,
                             debug = False):
    cmp_left, cmp_right = 0, len(cmp_list) - 1
    left, right = cmp_list[0][1], cmp_list[-1][1] + cmp_list[-1][2] - 1
    left_alt_set, right_alt_set = set(), set()

    def get_haplotype_and_seq(cmp_list):
        ht, seq = [], ""
        for i in range(len(cmp_list)):
            cmp_i = cmp_list[i]
            type, pos, length = cmp_i[:3]
            if len(cmp_i) <= 3:
                var_id = ""
            else:
                var_id = cmp_i[3]
            if type == "match":
                seq += ref_seq[pos:pos+length]
            elif type == "mismatch":
                seq += ref_seq[pos]
            elif type == "insertion":
                None
                # seq += data
            else:
                assert type == "deletion"

            if var_id != "" and var_id != "unknown":
                ht.append(var_id)
        return ht, seq

    # Left direction
    found = False
    for i in reversed(range(len(cmp_list))):
        cmp_i = cmp_list[i]
        type, cur_left, length = cmp_i[:3]
        var_id = cmp_i[3] if type in ["mismatch", "deletion"] else ""
        if type in ["match", "deletion"]:
            cur_right = cur_left + length - 1
        else:
            cur_right = cur_left

        cur_ht, cur_seq = get_haplotype_and_seq(cmp_list[:i+1])
        if len(cur_ht) == 0:
            cur_ht_str = str(left)
        else:
            cur_ht_str = "%d-%s" % (left, '-'.join(cur_ht))

        ht_i = lower_bound(Alts_left_list, cur_right + 1)
        for ht_j in reversed(range(0, min(ht_i + 1, len(Alts_left_list)))):
            ht_pos, ht = Alts_left_list[ht_j]
            if ht_pos < cur_left:
                break
            if ht_pos > cur_right:
                continue

            if len(cur_ht) > 0:
                if ht.find('-'.join(cur_ht)) == -1:
                    continue

            ht = ht.split('-')[:-1]
            if len(cur_ht) + 1 == len(ht):
                ht_pos = int(ht[0])
            else:
                var_id2 = ht[len(ht) - len(cur_ht) - 1]
                ht_type, ht_pos, ht_data = Vars[var_id2]
                if ht_type == "deletion":
                    ht_pos = ht_pos + int(ht_data) - 1
                    
            if left <= ht_pos:
                continue

            found = True

            if debug:
                print cmp_list[:i+1]
                print "\t", cur_ht, "vs", Alts_left_list[ht_j], ht_pos

            _, rep_ht = Alts_left_list[ht_j]

            if debug:
                print "\t\t", rep_ht, Alts_left[rep_ht], cur_seq

            for alt_ht_str in Alts_left[rep_ht]:
                alt_ht = alt_ht_str.split('-')
                alt_ht_left, alt_ht_right = int(alt_ht[0]), int(alt_ht[-1])
                assert alt_ht_right <= cur_right
                seq_pos = cur_right - alt_ht_right
                cur_pos = alt_ht_right
                part_alt_ht = []
                alt_ht = alt_ht[1:-1]
                for var_id_ in reversed(alt_ht):
                    var_type_, var_pos_, var_data_ = Vars[var_id_]
                    if var_type_ == "deletion":
                        del_len = int(var_data_)
                        var_pos_ = var_pos_ + del_len - 1
                    assert var_pos_ <= cur_pos
                    next_seq_pos = seq_pos + (cur_pos - var_pos_)
                    if next_seq_pos >= len(cur_seq):
                        break
                    if var_type_ == "single":
                        next_seq_pos += 1
                        cur_pos = var_pos_ - 1
                    elif var_type_ == "deletion":
                        cur_pos = var_pos_ - del_len
                    else:
                        assert var_type_ == "insertion"
                        assert False

                    if next_seq_pos >= len(cur_seq):
                        break
                    part_alt_ht.insert(0, var_id_)
                    seq_pos = next_seq_pos

                if len(part_alt_ht) > 0:
                    seq_left = len(cur_seq) - seq_pos - 1
                    part_alt_ht_str = "%d-%s" % (cur_pos - seq_left, '-'.join(part_alt_ht))
                    left_alt_set.add(part_alt_ht_str)
                        
                if debug:
                    print "\t\t", cur_left, alt_ht_str

        if found:
            cmp_left = i + 1
            left_alt_set.add(cur_ht_str)
            break

    if not found:
        left_alt_set.add(str(left))

    # Right direction
    found = False
    for i in range(0, len(cmp_list)):
        cmp_i = cmp_list[i]
        type, cur_left, length = cmp_i[:3]
        var_id = cmp_i[3] if type in ["mismatch", "deletion"] else ""
        if type in ["match", "deletion"]:
            cur_right = cur_left + length - 1
        else:
            cur_right = cur_left

        cur_ht, cur_seq = get_haplotype_and_seq(cmp_list[i:])
        if len(cur_ht) == 0:
            cur_ht_str = str(right)
        else:
            cur_ht_str = "%s-%d" % ('-'.join(cur_ht), right)

        ht_i = lower_bound(Alts_right_list, cur_left)
        for ht_j in range(ht_i, len(Alts_right_list)):
            ht_pos, ht = Alts_right_list[ht_j]
            if ht_pos > cur_right:
                break
            if ht_pos < cur_left:
                continue

            if len(cur_ht) > 0:
                if ht.find('-'.join(cur_ht)) == -1:
                    continue

            ht = ht.split('-')[1:]
            if len(cur_ht) + 1 == len(ht):
                ht_pos = int(ht[-1])
            else:
                var_id2 = ht[len(cur_ht)]
                _, ht_pos, _ = Vars[var_id2]

            if right >= ht_pos:
                continue

            found = True
            _, rep_ht = Alts_right_list[ht_j]

            # DK - debugging purposes
            if debug:
                print "DK1:", cmp_i, cmp_list
                print "DK2:", rep_ht, Alts_right[rep_ht]
                print "DK3:", left, right, ht_pos

            for alt_ht_str in Alts_right[rep_ht]:
                alt_ht = alt_ht_str.split('-')
                alt_ht_left, alt_ht_right = int(alt_ht[0]), int(alt_ht[-1])
                assert cur_left <= alt_ht_left
                seq_pos = alt_ht_left - cur_left
                cur_pos = alt_ht_left
                part_alt_ht = []
                alt_ht = alt_ht[1:-1]
                for var_id_ in alt_ht:
                    var_type_, var_pos_, var_data_ = Vars[var_id_]
                    assert var_pos_ >= cur_pos
                    next_seq_pos = seq_pos + (var_pos_ - cur_pos)
                    if next_seq_pos >= len(cur_seq):
                        break
                    
                    if var_type_ == "single":
                        next_seq_pos += 1
                        cur_pos = var_pos_ + 1
                    elif var_type_ == "deletion":
                        cur_pos = var_pos_ + int(var_data_)
                    else:
                        assert var_type_ == "insertion"
                        assert False

                    if next_seq_pos >= len(cur_seq):
                        break
                    part_alt_ht.append(var_id_)
                    seq_pos = next_seq_pos

                if len(part_alt_ht) > 0:
                    seq_left = len(cur_seq) - seq_pos - 1
                    part_alt_ht_str = "%s-%d" % ('-'.join(part_alt_ht), cur_pos + seq_left)
                    right_alt_set.add(part_alt_ht_str)
                        
        if found:
            cmp_right = i - 1
            right_alt_set.add(cur_ht_str)
            break

    if not found:
        right_alt_set.add(str(right))

    if cmp_right < cmp_left:
        cmp_left = 0
        left_alt_set = set([str(left)])

    # Sanity check
    ht_set_ = set()
    for ht in left_alt_set:
        ht = '-'.join(ht.split('-')[1:])
        if ht == "":
            continue
        if ht in ht_set_:
            print >> sys.stderr, "Error) %s should not be in" % ht, ht_set_
            assert False
        ht_set_.add(ht)
    for ht in right_alt_set:
        ht = '-'.join(ht.split('-')[:-1])
        if ht == "":
            continue
        if ht in ht_set_:
            print >> sys.stderr, "Error) %s should not be in" % ht, ht_set_
            assert False
        ht_set_.add(ht)

    # DK - debugging purposes
    if debug:
        print "cmp_list_range: [%d, %d]" % (cmp_left, cmp_right)
        print "left  alt set:", left_alt_set
        print "right alt set:", right_alt_set
    
    return cmp_left, cmp_right, list(left_alt_set), list(right_alt_set)

