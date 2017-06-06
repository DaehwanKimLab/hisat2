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
import math
import random
from copy import deepcopy
from datetime import datetime


##################################################
#   Sequence processing routines
##################################################


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


##################################################
#   Alleles, variants, haplotypes, etc.
##################################################


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


##################################################
#   Database releated routines
##################################################

    
"""
Download GRCh38 human reference and HISAT2 indexes
"""
def download_genome_and_index():
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]
    if not check_files(HISAT2_fnames):
        os.system("wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz; tar xvzf grch38.tar.gz; rm grch38.tar.gz")
        os.system("hisat2-inspect grch38/genome > genome.fa")
        os.system("samtools faidx genome.fa")


"""
"""
def clone_hisatgenotype_database():
    os.system("git clone https://github.com/infphilo/hisatgenotype_db.git")
    
    # Check out one particular revision just to have the same data across multiple computers        
    # revision = "d3b559b34b96ff9e7f0d97476222d8e4cdee63ad" # Revision on November 16, 2016
    # os.system("cd IMGTHLA; git checkout %s; cd .." % revision)


"""
"""
def extract_database_if_not_exists(base,
                                   locus_list,
                                   inter_gap = 30,
                                   intra_gap = 50,
                                   partial = True,
                                   verbose = False):
    fnames = [base + "_backbone.fa",
              base + "_sequences.fa",
              base + ".locus",
              base + ".snp",
              base + ".index.snp",
              base + ".haplotype",
              base + ".link",
              base + ".partial"]
    if check_files(fnames):
        return

    extract_cmd = ["hisatgenotype_extract_vars.py"]
    extract_cmd += ["--base", base]
    if len(locus_list) > 0:
        extract_cmd += ["--locus-list", ','.join(locus_list)]    
    if not partial:
        extract_cmd += ["--no-partial"]
    if base == "codis":
        extract_cmd += ["--whole-haplotype"]
    else:
        extract_cmd += ["--inter-gap", str(inter_gap),
                        "--intra-gap", str(intra_gap)]
    if base == "hla":
        extract_cmd += ["--min-var-freq", "0.1"]

    if base == "codis":
        extract_cmd += ["--leftshift"]

    # DK - debugging purposes
    # extract_cmd += ["--ext-seq", "300"]
    if verbose:
        print >> sys.stderr, "\tRunning:", ' '.join(extract_cmd)
    proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
    proc.communicate()

    if not check_files(fnames):
        print >> sys.stderr, "Error: hisatgenotype_extract_vars failed!"
        sys.exit(1)

        
"""
"""
def build_index_if_not_exists(base,
                              aligner,
                              index_type,
                              threads = 1,
                              verbose = False):
    if aligner == "hisat2":
        # Build HISAT2 graph indexes based on the above information
        if index_type == "graph":
            hisat2_graph_index_fnames = ["%s.graph.%d.ht2" % (base, i+1) for i in range(8)]
            if not check_files(hisat2_graph_index_fnames):
                build_cmd = ["hisat2-build",
                             "-p", str(threads),
                             "--snp", "%s.index.snp" % base,
                             "--haplotype", "%s.haplotype" % base,
                             "%s_backbone.fa" % base,
                             "%s.graph" % base]
                if verbose:
                    print >> sys.stderr, "\tRunning:", ' '.join(build_cmd)
                proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
                proc.communicate()        
                if not check_files(hisat2_graph_index_fnames):
                    print >> sys.stderr, "Error: indexing HLA failed!  Perhaps, you may have forgotten to build hisat2 executables?"
                    sys.exit(1)
        # Build HISAT2 linear indexes based on the above information
        else:
            assert index_type == "linear"
            hisat2_linear_index_fnames = ["%s.linear.%d.ht2" % (base, i+1) for i in range(8)]
            if not check_files(hisat2_linear_index_fnames):
                build_cmd = ["hisat2-build",
                             "%s_backbone.fa,%s_sequences.fa" % (base, base),
                             "%s.linear" % base]
                proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
                proc.communicate()        
                if not check_files(hisat2_linear_index_fnames):
                    print >> sys.stderr, "Error: indexing HLA failed!"
                    sys.exit(1)                    
    else:
        # Build Bowtie2 indexes based on the above information
        assert aligner == "bowtie2" and index_type == "linear"        
        bowtie2_index_fnames = ["%s.%d.bt2" % (base, i+1) for i in range(4)]
        bowtie2_index_fnames += ["%s.rev.%d.bt2" % (base, i+1) for i in range(2)]
        if not tcheck_files(bowtie2_index_fnames):
            build_cmd = ["bowtie2-build",
                         "%s_backbone.fa,%s_sequences.fa" % (base, base),
                         base]
            proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'))
            proc.communicate()        
            if not check_files(bowtie2_index_fnames):
                print >> sys.stderr, "Error: indexing HLA failed!"
                sys.exit(1)

                    

##################################################
#   Read simulation and alignment
##################################################


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
                index_name,
                index_type,
                base_fname,
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
            aligner_cmd += ["--haplotype"]
            if base_fname == "codis":
                aligner_cmd += ["--enable-codis"]        
    elif aligner == "bowtie2":
        aligner_cmd = [aligner,
                       "--no-unal",
                       "-k", "10"]
    else:
        assert False
    aligner_cmd += ["-x", index_name]
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
                    if right_pos + j < len(mpileup):
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
def get_pair_interdist(alignview_cmd,
                       simulation,
                       verbose):
    bamview_proc = subprocess.Popen(alignview_cmd,
                                    stdout=subprocess.PIPE,
                                    stderr=open("/dev/null", 'w'))
    sort_read_cmd = ["sort", "-k", "1,1", "-s"] # -s for stable sorting
    alignview_proc = subprocess.Popen(sort_read_cmd,
                                      stdin=bamview_proc.stdout,
                                      stdout=subprocess.PIPE,
                                      stderr=open("/dev/null", 'w'))

    dist_list = []
    prev_read_id = None
    cigar_re = re.compile('\d+\w')
    reads = []
    for line in alignview_proc.stdout:
        line = line.strip()
        cols = line.split()
        read_id, flag, _, pos, _, cigar_str = cols[:6]
        read_seq = cols[9]
        flag, pos = int(flag), int(pos)
        # Unalined?
        if flag & 0x4 != 0:
            continue

        if simulation:
            read_id = read_id.split('|')[0]

        # Concordantly mapped?
        if flag & 0x2 != 0:
            concordant = True
        else:
            concordant = False

        NH, YT = sys.maxint, ""
        for i in range(11, len(cols)):
             col = cols[i]
             if col.startswith("NH"):
                 NH = int(col[5:])
             elif col.startswith("YT"):
                 YT = col[5:]
        if NH > 1 or YT != "CP":
            continue

        if prev_read_id != None and read_id != prev_read_id:
            if len(reads) == 2:
                left1, right1 = reads[0]
                left2, right2 = reads[1]
                if left1 <= left2:
                    dist = left2 - right1 - 1
                else:
                    dist = left1 - right2 - 1
                dist_list.append(dist)
            reads = []

        left_pos = right_pos =  pos
        cigars = cigar_re.findall(cigar_str)
        cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
        for i in range(len(cigars)):
            cigar_op, length = cigars[i]
            if cigar_op in "MND":
                right_pos += length

        reads.append([left_pos, right_pos - 1])
        
        prev_read_id = read_id

    dist_list = sorted(dist_list)
    dist_avg = sum(dist_list) / max(1, len(dist_list))
    if len(dist_list) > 0:
        dist_median = dist_list[len(dist_list)/2]
    else:
        dist_median = -1

    return dist_median


##################################################
#   Statistical routines
##################################################


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
def single_abundance(Gene_cmpt, Gene_length, exonic = False):
    def normalize(prob):
        total = sum(prob.values())
        for allele, mass in prob.items():
            prob[allele] = mass / total        

    def normalize_len(prob, length):
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
    if exonic:
        normalize(Gene_prob)
    else:
        normalize_len(Gene_prob, Gene_length)

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
        if exonic:
            normalize(Gene_prob_next)
        else:
            normalize_len(Gene_prob_next, Gene_length)
        return Gene_prob_next

    def select_alleles(Gene_prob):
        if len(Gene_prob) == 0:
            return Gene_prob
        Gene_prob2 = {}
        max_prob = max(Gene_prob.values())
        for allele, prob in Gene_prob.items():
            if prob >= max_prob / 10.0:
                Gene_prob2[allele] = prob
        return Gene_prob2

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
            Gene_prob = select_alleles(Gene_prob)

        # DK - debugging purposes
        if iter % 10 == 0 and False:
            print >> sys.stderr, "iter", iter
            for allele, prob in Gene_prob.items():
                if prob >= 0.01:
                    print >> sys.stderr, "\t", iter, allele, prob
        
        iter += 1
        
    Gene_prob = select_alleles(Gene_prob)
    if exonic:
        normalize(Gene_prob)
    else:
        normalize_len(Gene_prob, Gene_length)
    Gene_prob = [[allele, prob] for allele, prob in Gene_prob.items()]
    Gene_prob = sorted(Gene_prob, cmp=Gene_prob_cmp)
    return Gene_prob


##################################################
#   Realignment, alternative alignments
##################################################


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
        i_found = False
        cmp_i = cmp_list[i]
        type, cur_left, length = cmp_i[:3]
        var_id = cmp_i[3] if type in ["mismatch", "deletion"] else ""

        # DK - debugging purpose
        if type in ["mismatch", "deletion", "insertion"]:
            if not var_id.startswith("hv"):
                continue
        
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
                    
            if left < ht_pos:
                continue

            i_found = True

            if debug:
                print cmp_list[:i+1]
                print "\t", cur_ht, "vs", Alts_left_list[ht_j], ht_pos

            _, rep_ht = Alts_left_list[ht_j]

            if debug:
                print "DK1:", cmp_i, cmp_list
                print "DK2:", rep_ht, Alts_left[rep_ht]
                print "DK3:", left, right, ht_pos

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
                        next_cur_pos = var_pos_ - 1
                    elif var_type_ == "deletion":
                        next_cur_pos = var_pos_ - del_len
                    else:
                        assert var_type_ == "insertion"
                        assert False

                    part_alt_ht.insert(0, var_id_)
                    if next_seq_pos >= len(cur_seq):
                        break
                    seq_pos, cur_pos = next_seq_pos, next_cur_pos

                if len(part_alt_ht) > 0:
                    seq_left = len(cur_seq) - seq_pos - 1
                    part_alt_ht_str = ""
                    if found:
                        var_id_list = []
                        for j in range(i + 1, cmp_left):
                            cmp_j = cmp_list[j]
                            if cmp_j[0] in ["mismatch", "deletion", "insertion"]:
                                var_id_ = cmp_j[3]
                                if var_id_.startswith("hv"):
                                    var_id_list.append(var_id_)
                        if len(var_id_list) > 0:
                            part_alt_ht_str = '-' + '-'.join(var_id_list)
                    part_alt_ht_str = ("%d-%s" % (cur_pos - seq_left, '-'.join(part_alt_ht))) + part_alt_ht_str
                    left_alt_set.add(part_alt_ht_str)
                        
                if debug:
                    print "\t\t", cur_left, alt_ht_str

        if i_found:
            if not found:
                cmp_left = i + 1
                left_alt_set.add(cur_ht_str)
            found = True

    if not found:
        left_alt_set.add(str(left))

    # Right direction
    found = False
    for i in range(0, len(cmp_list)):
        i_found = False
        cmp_i = cmp_list[i]
        type, cur_left, length = cmp_i[:3]
        var_id = cmp_i[3] if type in ["mismatch", "deletion"] else ""

        # DK - debugging purpose
        if type in ["mismatch", "deletion", "insertion"]:
            if not var_id.startswith("hv"):
                continue

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

            if right > ht_pos:
                continue

            i_found = True
            _, rep_ht = Alts_right_list[ht_j]

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
                        next_cur_pos = var_pos_ + 1
                    elif var_type_ == "deletion":
                        next_cur_pos = var_pos_ + int(var_data_)
                    else:
                        assert var_type_ == "insertion"
                        assert False

                    part_alt_ht.append(var_id_)
                    if next_seq_pos >= len(cur_seq):
                        break
                    seq_pos, cur_pos = next_seq_pos, next_cur_pos

                if len(part_alt_ht) > 0:
                    seq_left = len(cur_seq) - seq_pos - 1
                    assert seq_left >= 0
                    part_alt_ht_str = ""
                    if found:
                        var_id_list = []
                        for j in range(cmp_right + 1, i):
                            cmp_j = cmp_list[j]
                            if cmp_j[0] in ["mismatch", "deletion", "insertion"]:
                                var_id_ = cmp_j[3]
                                if var_id_.startswith("hv"):
                                    var_id_list.append(var_id_)
                        if len(var_id_list) > 0:
                            part_alt_ht_str = '-'.join(var_id_list) + '-'
                    part_alt_ht_str += ("%s-%d" % ('-'.join(part_alt_ht), cur_pos + seq_left))
                    right_alt_set.add(part_alt_ht_str)
                        
        if i_found:            
            if not found:
                cmp_right = i - 1
                right_alt_set.add(cur_ht_str)
            found = True

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

            # DK - debugging purposes
            print "DK: cmp_list_range: [%d, %d]" % (cmp_left, cmp_right)
            print "DK: cmp_list:", cmp_list
            print "DK: left_alt_set:", left_alt_set, "right_alt_set:", right_alt_set
            
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

    if debug:
        print "cmp_list_range: [%d, %d]" % (cmp_left, cmp_right)
        print "left  alt set:", left_alt_set
        print "right alt set:", right_alt_set
    
    return cmp_left, cmp_right, list(left_alt_set), list(right_alt_set)

